--------------------------------------------------------------------------------
-- Network construction
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Component.Network (builder) where

import           Bio.Data.Bed
import           Bio.Data.Experiment.Types
import           Bio.GO.GREAT
import           Bio.Pipeline.Instances            ()
import           Bio.Utils.Misc                    (readInt)
import           Conduit
import           Control.Arrow                     (first, second, (&&&), (***))
import           Control.Lens
import           Data.Binary                       (decodeFile, encodeFile)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (mk, original)
import           Data.Function                     (on)
import qualified Data.HashMap.Strict               as M
import           Data.List
import           Data.Maybe
import           Data.Ord
import qualified Data.Text                         as T
import           Scientific.Workflow

import           Constants
import           Types

builder :: Builder ()
builder = do
    node "ass00" [| \x -> getDomains <$> netOutput <*>
        getConfig' "annotation" <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    ["peak01"] ~> "ass00"
    node "ass01" [| \(x,y) -> return $ zip x $ repeat y |] $ do
        label .= "prepare input"
        submitToRemote .= Just False
    ["ass00", "peak02"] ~> "ass01"
    node "ass02" [| \xs -> do
        dir <- netOutput
        liftIO $ mapM (linkGeneToTFs dir) xs
        |] $ batch .= 1 >> stateful .= True
    node "ass03" [| \xs -> do
        dir <- netOutput
        liftIO $ mapM (printEdgeList dir) xs
        |] $ batch .= 1 >> stateful .= True
    path ["ass01", "ass02", "ass03"]

getDomains :: FilePath   -- output
           -> FilePath   -- annotation
           -> [Experiment ATAC_Seq] -> IO [Experiment ATAC_Seq]
getDomains dir anno = mapM $ \e -> do
    let [fl] = filter ((==NarrowPeakFile) . (^.format)) $ e^.files
        output = dir ++ "/" ++ T.unpack (e^.eid) ++ "_gene_reg_domains.bed"
        newFile = format .~ BedFile $
                  info .~ [] $
                  keywords .~ ["gene regulatory domain"] $
                  location .~ output $ fl
    peaks <- readBed' $ fl^.location
    activePromoters <- gencodeActiveGenes anno peaks
    writeBed' output $ getGeneDomains GREAT activePromoters
    return $ files %~ (newFile:) $ e

linkGeneToTFs :: FilePath -> (Experiment ATAC_Seq, FilePath) -> IO (Experiment ATAC_Seq)
linkGeneToTFs dir (e, tfbs) = do
    let [peakFl] = filter ((==NarrowPeakFile) . (^.format)) $ e^.files
        [domainFl] = filter ((==["gene regulatory domain"]) . (^.keywords)) $ e^.files
    genes <- readBed' $ domainFl^.location
    peaks <- readBed' $ peakFl^.location

    regulators <- findRegulators genes peaks tfbs
    let result :: [Linkage]
        result = map (second (map ((head *** id) . unzip) . groupBy ((==) `on` fst) .
            sortBy (comparing fst) . map (getTFName &&& id))) $ M.toList $
            M.fromListWith (++) $ map (first (mk . fromJust . bedName)) $
            regulators
        output = dir ++ "/" ++ T.unpack (e^.eid) ++ ".assign"
        newFile = format .~ Other $
                  info .~ [] $
                  keywords .~ ["gene-TF assignment"] $
                  location .~ output $ peakFl
    encodeFile output result
    return $ files .~ [newFile] $ e
  where
    getTFName = mk . head . B.split '+' . fromJust . bedName

findRegulators :: [BED]          -- ^ gene domains
               -> [NarrowPeak]   -- ^ open chromatin
               -> FilePath          -- ^ TFBS
               -> IO [(BED, [BED])]
findRegulators genes peaks tfbs = do
    tfbs' <- readBed tfbs =$= intersectBed peaks' $$ sinkList
    yieldMany genes =$= intersectBedWith id tfbs' =$=
        filterC (not . null .snd) $$ sinkList
  where
    peaks' = flip map peaks $ \p -> let c = chromStart p + (fromJust . _npPeak) p
                                    in BED3 (chrom p) (c - 50) (c+ 50)

data Method = GREAT  -- ^ using GREAT strategy

getGeneDomains :: Method -> [((B.ByteString, Int, Bool), B.ByteString)]
               -> [BED]   -- ^ genes' regulatory domains. A single gene can
                          -- potentially have multiple domains.
getGeneDomains GREAT = map ( \(b, x) ->
    BED (chrom b) (chromStart b) (chromEnd b) (Just x) Nothing Nothing ) .
    getRegulatoryDomains (BasalPlusExtension 5000 1000 1000000)

-- | Identify active genes by overlapping their promoters with activity indicators.
gencodeActiveGenes :: FilePath   -- ^ gencode file in GTF format
                   -> [BED3]     -- ^ feature that is used to determine the activity
                                 -- of promoters, e.g., H3K4me3 peaks or ATAC-seq peaks
                   -> IO [((B.ByteString, Int, Bool), B.ByteString)]  -- ^ chr, tss, strand and gene name.
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
        runIdentity $ yieldMany promoters =$= intersectBedWith id peaks =$=
        filterC (not . null . snd) $$ sinkList
  where
    f xs = let name = B.filter (/='\"') $ last $ B.split ' ' $ head $
                      filter (B.isInfixOf "gene_name") $ B.split ';' $ last xs
               start = readInt $ xs !! 3
               end = readInt $ xs !! 4
               strand = if xs !! 6 == "+" then True else False
          in ((xs !! 0, if strand then start else end, strand), name)
    g xs = not $ B.isPrefixOf "#" (head xs) || (xs !! 2 /= "transcript")

printEdgeList :: FilePath -> Experiment a -> IO ()
printEdgeList dir e = do
    let [fl] = filter (\x -> x^.keywords == ["gene-TF assignment"]) $ e^.files
    result <- decodeFile $ fl^.location :: IO [Linkage]
    let output = dir ++ "/" ++ T.unpack (e^.eid) ++ "_network.tsv"
    B.writeFile output $ B.unlines $ flip map result $ \(a,b) ->
        B.intercalate "\t" [original a, B.intercalate "," $ map original $ fst $ unzip b]
