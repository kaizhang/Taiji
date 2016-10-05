-- Assign TFs to their targets
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
module Assign
    ( getDomains
    , linkGeneToTFs
    , GeneName
    , Linkage
    ) where

import           Bio.Data.Bed
import           Bio.Data.Experiment.Types
import           Bio.GO.GREAT
import           Bio.Pipeline.Instances    ()
import           Bio.Utils.Misc            (readInt)
import           Control.Arrow             (first, second, (&&&), (***))
import           Control.Lens
import           Data.Binary               (Binary (..), decodeFile, encodeFile)
import qualified Data.ByteString.Char8     as B
import           Data.CaseInsensitive      (mk)
import           Data.Function             (on)
import           Data.List
import qualified Data.Map                  as M
import           Data.Maybe
import           Data.Ord
import qualified Data.Text                 as T

import           Types

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
    sites <- readBed' tfbs
    let [peakFl] = filter ((==NarrowPeakFile) . (^.format)) $ e^.files
        [domainFl] = filter ((==["gene regulatory domain"]) . (^.keywords)) $ e^.files
    genes <- readBed' $ domainFl^.location
    peaks <- readBed' $ peakFl^.location
    let result :: [Linkage]
        result = map (second (map ((head *** id) . unzip) . groupBy ((==) `on` fst) .
            sortBy (comparing fst) . map (getTFName &&& id))) $ M.toList $
            M.fromListWith (++) $ map (first (mk . fromJust . bedName)) $
            findRegulators genes peaks sites
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
               -> [BED]          -- ^ TFBS
               -> [(BED, [BED])]
findRegulators genes peaks sites = filter (not . null . snd) $
    intersectBedWith id genes $ intersectBed sites peaks'
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
        filter (not . null . snd) $ intersectBedWith id promoters peaks
  where
    f xs = let name = B.filter (/='\"') $ last $ B.split ' ' $ head $
                      filter (B.isInfixOf "gene_name") $ B.split ';' $ last xs
               start = readInt $ xs !! 3
               end = readInt $ xs !! 4
               strand = if xs !! 6 == "+" then True else False
          in ((xs !! 0, if strand then start else end, strand), name)
    g xs = not $ B.isPrefixOf "#" (head xs) || (xs !! 2 /= "transcript")
