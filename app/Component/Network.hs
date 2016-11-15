--------------------------------------------------------------------------------
-- Network construction
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE GADTs #-}

module Component.Network (builder) where

import           Bio.Data.Bed
import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Utils (mergeExps)
import           Bio.GO.GREAT
import           Bio.Pipeline.Instances            ()
import           Bio.Pipeline.Utils (mapOfFiles)
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
    node "Get_reg_regions" [| \x -> getDomains <$> netOutput <*>
        getConfig' "annotation" <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    ["ATAC_callpeaks"] ~> "Get_reg_regions"
    node "Link_TF_gene_prepare" [| \(peaks, reg_domains, y) ->
        return $ zip (mergeExps $ peaks ++ reg_domains) $ repeat y
        |] $ label .= "prepare input" >> submitToRemote .= Just False
    ["ATAC_callpeaks", "Get_reg_regions", "Find_TF_sites"] ~> "Link_TF_gene_prepare"
    node "Link_TF_gene" [| \xs -> do
        dir <- netOutput
        liftIO $ mapM (linkGeneToTFs dir) xs
        |] $ batch .= 1 >> stateful .= True
    node "Output_network" [| \xs -> do
        dir <- netOutput
        liftIO $ mapM (printEdgeList dir) xs
        |] $ batch .= 1 >> stateful .= True
    path ["Link_TF_gene_prepare", "Link_TF_gene", "Output_network"]

getDomains :: FilePath   -- output
           -> FilePath   -- annotation
           -> [ATACSeq] -> IO [ATACSeq]
getDomains dir anno = mapOfFiles fn
  where
    fn e r (Single fl) = if r^.number /= 0 || fl^.format /= NarrowPeakFile
        then return []
        else do
            let output = dir ++ "/" ++ T.unpack (e^.eid) ++ "_gene_reg_domains.bed"
                newFile = Single $ format .~ BedFile $
                    keywords .~ ["gene regulatory domain"] $
                    location .~ output $ emptyFile
            peaks <- readBed' $ fl^.location
            activePromoters <- gencodeActiveGenes anno peaks
            writeBed' output $ getGeneDomains GREAT activePromoters
            return [newFile]
    fn _ _ _ = return []

linkGeneToTFs :: FilePath -> (ATACSeq, FilePath) -> IO ATACSeq
linkGeneToTFs dir (e, tfbs) = do
    let fls = e^..replicates.folded.filtered ((==0) . (^.number)).files.folded._Single
        [peakFl] = filter ((==NarrowPeakFile) . (^.format)) fls
        [domainFl] = filter ((==["gene regulatory domain"]) . (^.keywords)) fls
    genes <- readBed' $ domainFl^.location
    peaks <- readBed' $ peakFl^.location

    regulators <- findRegulators genes peaks tfbs
    let result :: [Linkage]
        result = map (second (map ((head *** id) . unzip) . groupBy ((==) `on` fst) .
            sortBy (comparing fst) . map (getTFName &&& id))) $ M.toList $
            M.fromListWith (++) $ map (first (mk . fromJust . bedName)) $
            regulators
        output = dir ++ "/" ++ T.unpack (e^.eid) ++ ".assign"
        newFile = Single $ format .~ Other $
                  info .~ [] $
                  keywords .~ ["gene-TF assignment"] $
                  location .~ output $ peakFl
    encodeFile output result
    return $ replicates .~ [files .~ [newFile] $ emptyReplicate] $ e
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

printEdgeList :: Experiment e => FilePath -> e -> IO ()
printEdgeList dir e = do
    let [fl] = e^..replicates.folded.filtered ((==0) . (^.number)).files.folded.
            _Single.filtered ((==["gene-TF assignment"]) . (^.keywords))
    result <- decodeFile $ fl^.location :: IO [Linkage]
    let output = dir ++ "/" ++ T.unpack (e^.eid) ++ "_network.tsv"
    B.writeFile output $ B.unlines $ flip map result $ \(a,b) ->
        B.intercalate "\t" [original a, B.intercalate "," $ map original $ fst $ unzip b]
