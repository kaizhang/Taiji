{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Pipeline (pipeline) where

import           Bio.Data.Bed
import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Utils
import           Bio.GO.GREAT
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Instances
import           Bio.Pipeline.NGS
import           Bio.Pipeline.ScanMotifs
import           Bio.Utils.Misc            (readInt)
import           Control.Arrow             (first, second, (&&&), (***))
import           Control.Lens
import           Data.Binary               (decodeFile, encodeFile)
import qualified Data.ByteString.Char8     as B
import           Data.Dynamic
import           Data.Function             (on)
import           Data.List
import qualified Data.Map                  as M
import           Data.Maybe
import           Data.Ord
import qualified Data.Text                 as T

import           Scientific.Workflow

(!) :: Typeable a => M.Map T.Text Dynamic -> T.Text -> a
(!) m x = fromJust $ fromDynamic $ fromJust $ m ^.at x

config :: M.Map T.Text Dynamic
config = [ ("input", toDyn ("data/input.yaml" :: FilePath))
         , ("outputDir", toDyn ("output/" :: FilePath))
         , ("motifFile", toDyn ("data/cisbp.meme" :: FilePath))
         , ("genome", toDyn ("/home/kai/data/genome/mm10.index" :: FilePath))
         , ("gencode", toDyn ("/home/kai/project/TRM/data/gencode.vM9.annotation.gff3" :: FilePath))
         ]

mm10Index = "/home/kai/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa"

gencodeActiveGenes :: FilePath   -- ^ gencode file
                   -> [BED3]     -- ^ feature that is used to determine the activity
                                 -- of promoters, e.g., H3K4me3 peaks or ATAC-seq peaks
                   -> IO [((B.ByteString, Int, Bool), B.ByteString)]
gencodeActiveGenes input peaks = do
    c <- B.readFile input
    let promoters = map ( \((chr, i, isForward), geneName) ->
            if isForward
                then BED chr (i-5000) (i+1000) (Just geneName)
                     (Just $ fromIntegral i) (Just True)
                else BED chr (i-1000) (i+5000) (Just geneName)
                     (Just $ fromIntegral i) (Just False)
            ) $ map f $ filter g $ map (B.split '\t') $ B.lines c
    return $ map (\x -> ( (chrom x, truncate $ fromJust $ bedScore x,
        fromJust $ bedStrand x), fromJust $ bedName x )) $ fst $ unzip $
        filter (not . null . snd) $ intersectBedWith id promoters peaks
  where
    f xs = let name = last $ B.split '=' $ head $
                      filter (B.isPrefixOf "gene_name") $ B.split ';' $ last xs
               start = readInt $ xs !! 3
               end = readInt $ xs !! 4
               strand = if xs !! 6 == "+" then True else False
          in ((xs !! 0, if strand then start else end, strand), name)
    g xs = not $ B.isPrefixOf "#" (head xs) || (xs !! 2 /= "transcript")


getDomains :: FilePath -> [Experiment] -> IO [Experiment]
getDomains dir = mapM $ \e -> do
    let [fl] = filter ((==NarrowPeakFile) . (^.format)) $ e^.files
        output = dir ++ "/" ++ T.unpack (e^.eid) ++ "_gene_reg_domains.bed"
        newFile = format .~ BedFile $
                  info .~ [] $
                  keywords .~ ["gene regulatory domain"] $
                  location .~ output $ fl
    peaks <- readBed' $ fl^.location
    activePromoters <- gencodeActiveGenes (config!"gencode") peaks
    writeBed' output $ map (\(b, x) ->
        BED (chrom b) (chromStart b) (chromEnd b) (Just x) Nothing Nothing) $
        getRegulatoryDomains (BasalPlusExtension 5000 1000 1000000) activePromoters
    return $ files %~ (newFile:) $ e

linkGeneToTFs :: FilePath -> (Experiment, FilePath) -> IO Experiment
linkGeneToTFs dir (e, tfbs) = do
    sites <- readBed' tfbs
    let [peakFl] = filter ((==NarrowPeakFile) . (^.format)) $ e^.files
        [domainFl] = filter ((==["gene regulatory domain"]) . (^.keywords)) $ e^.files
    genes <- readBed' $ domainFl^.location
    peaks <- readBed' $ peakFl^.location
    let result :: [(B.ByteString, [(B.ByteString, [BED])])]
        result = map (second (map ((head *** id) . unzip) . groupBy ((==) `on` fst) .
            sortBy (comparing fst) . map (getTFName &&& id))) $ M.toList $
            M.fromListWith (++) $ map (first (fromJust . bedName)) $
            findRegulators genes peaks sites
        output = dir ++ "/" ++ T.unpack (e^.eid) ++ ".assign"
        newFile = format .~ Other $
                  info .~ [] $
                  keywords .~ ["gene-TF assignment"] $
                  location .~ output $ peakFl
    encodeFile output result
    return $ files .~ [newFile] $ e
  where
    getTFName = head . B.split '+' . fromJust . bedName

findRegulators :: [BED]          -- ^ gene domains
               -> [NarrowPeak]   -- ^ open chromatin
               -> [BED]          -- ^ TFBS
               -> [(BED, [BED])]
findRegulators genes peaks sites = filter (not . null . snd) $
    intersectBedWith id genes $ intersectBed sites peaks'
  where
    peaks' = flip map peaks $ \p -> let c = chromStart p + (fromJust . _npPeak) p
                                    in BED3 (chrom p) (c - 50) (c+ 50)

printEdgeList :: FilePath -> Experiment -> IO ()
printEdgeList dir e = do
    let [fl] = filter (\x -> x^.keywords == ["gene-TF assignment"]) $ e^.files
    result <- decodeFile $ fl^.location :: IO [(B.ByteString, [(B.ByteString, [BED])])]
    let output = dir ++ "/" ++ T.unpack (e^.celltype) ++ "_" ++
            T.unpack (e^.target) ++ ".tsv"
    B.writeFile output $ B.unlines $ flip map result $ \(a,b) ->
        B.intercalate "\t" [a, B.intercalate "," $ fst $ unzip b]

pipeline :: Builder ()
pipeline = do
    node "ex00" [| \() -> fmap fromJust $ readExp $ config!"input" |] $
        submitToRemote .= Just False
    node "ex01" [| bwaAlign (config!"outputDir") mm10Index (return ())|] $
        batch .= 1
    node "ex02" [| filterBam (config!"outputDir") |] $
        batch .= 1
    node "ex03" [| removeDuplicates "/home/kai/software/picard-tools-1.140/picard.jar" (config!"outputDir") |] $
        batch .= 1
    node "ex04" [| bamToBed (config!"outputDir") |] $
        batch .= 1
    node "ex05" [| mapM $ \x -> return (x, Nothing) |] $
        batch .= 1
    node "ex06" [| callPeaks (config!"outputDir") (return ()) |] $
        batch .= 1
    node "ex07" [| \exps -> do
        let f x = map (^.location) $ filter ((==NarrowPeakFile) . (^.format)) $ x^.files
        scanMotifs (config!"genome") (config!"motifFile") 1e-5
            ((config!"outputDir") ++ "/TFBS.bed") (concatMap f (exps :: [Experiment]))
        |] $ return ()
    path ["ex00", "ex01", "ex02", "ex03", "ex04", "ex05", "ex06", "ex07"]

    node "ex08" [| getDomains (config!"outputDir") |] $ batch .= 1
    ["ex06"] ~> "ex08"

    node "ex09" [| \(x,y) -> return $ zip x $ repeat y |] $ do
        label .= "prepare input"
        submitToRemote .= Just False
    ["ex08", "ex07"] ~> "ex09"
    node "ex10" [| sequence . map (linkGeneToTFs $ config!"outputDir") |] $ batch .= 1
    node "ex11" [| sequence . map (printEdgeList $ config!"outputDir") |] $ batch .= 1
    path ["ex09", "ex10", "ex11"]
