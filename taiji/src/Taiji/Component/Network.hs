--------------------------------------------------------------------------------
-- Network construction
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE GADTs             #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Taiji.Component.Network where

import           Bio.Data.Bed
import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Utils
import           Bio.GO.GREAT
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Instances    ()
import           Bio.Utils.Misc            (readInt)
import           Conduit
import           Control.Arrow             (first, second, (&&&), (***))
import           Control.Lens
import           Data.Binary               (decodeFile, encodeFile)
import qualified Data.ByteString.Char8     as B
import           Data.CaseInsensitive      (mk, original)
import           Data.Function             (on)
import qualified Data.HashMap.Strict       as M
import qualified Data.IntervalMap.Strict   as IM
import           Data.List
import           Data.List.Ordered         (nubSort)
import           Data.Maybe
import           Data.Ord
import qualified Data.Set                  as S
import qualified Data.Text                 as T
import           Scientific.Workflow
import           System.IO.Temp            (withTempFile)

import           Taiji.Constants
import           Taiji.Types

builder :: Builder ()
builder = do
    node "ATAC_find_active_promoter" [| \e -> do
        anno <- getConfig' "annotation"
        liftIO $ withTempFile "./" "tmp_macs2_file." $ \tmp _ -> do
            let [peakFl] = e^..replicates.folded
                    .filtered ((==0) . (^.number)).files.folded
                    .filtered (\fl -> formatIs BedGZip fl || formatIs BedFile fl)

            _ <- callPeaks tmp peakFl Nothing $
                pair .= pairedEnd e >> cutoff .= QValue 0.1
            peaks <- readBed' tmp :: IO [BED3]
            genes <- gencodeActiveGenes anno peaks
            return (e^.eid, genes)
        |] $ do
            batch .= 1
            stateful .= True
            note .= "Identify active promoters by overlapping annotated promoter regions with active chromatin."
    path ["ATAC_combine_reps", "ATAC_find_active_promoter"]

    node "Link_TF_gene_prepare" [| \(promoters, oriInput, peaks, tfbs) -> do
        let hic = M.fromList $ map (\x -> (x^.groupName, x)) $ oriInput^._4
        return $ ContextData tfbs $ flip map peaks $ \p ->
            (p, M.lookup (p^.groupName) hic, fromJust $ lookup (p^.eid) promoters)
        |] $ do
            submitToRemote .= Just False
            note .= "Prepare for linking TFs and genes."
    [ "ATAC_find_active_promoter", "Initialization", "ATAC_callpeaks"
        , "Find_TF_sites_merge" ] ~> "Link_TF_gene_prepare"

    node "Link_TF_gene" [| \(ContextData tfbs (peak, hic, promoters)) ->
        linkGeneToTFs <$> netOutput <*> return tfbs <*> return peak <*>
            return hic <*> return promoters >>= liftIO
        |] $ do
            batch .= 1
            stateful .= True
            note .= "First link TFs to enhancers and promoters, and then link to genes according to enhancer-gene and promoter-gene assignments."
    node "Output_network" [| \x -> printEdgeList <$> netOutput <*>
        return x >>= liftIO
        |] $ do
            batch .= 1
            stateful .= True
            note .= "Save network structures to files."
    path ["Link_TF_gene_prepare", "Link_TF_gene", "Output_network"]

linkGeneToTFs :: FilePath   -- ^ Output dir
              -> FilePath   -- ^ TFBS
              -> ATACSeq    -- ^ peaks
              -> Maybe HiC        -- ^ 3D interaction
              -> [((B.ByteString, Int, Bool), B.ByteString)]  -- ^ Genes
              -> IO ATACSeq
linkGeneToTFs dir tfbs e hic activeGenes = do
    peaks <- readBed' $ peakFl^.location
    regulators <- runResourceT $ readBed tfbs $$ findRegulators activeGenes loops peaks
    let result :: [Linkage]
        result = map (second (
            map ((head *** id) . unzip) .
            groupBy ((==) `on` fst) .
            sortBy (comparing fst) .
            map (getTFName &&& id)
            )) $ map (first mk) regulators
        output = dir ++ "/" ++ T.unpack (e^.eid) ++ ".assign"
        newFile = Single $ format .~ Other $
                  info .~ [] $
                  tags .~ ["gene-TF assignment"] $
                  location .~ output $ peakFl
    encodeFile output result
    return $ replicates .~ [files .~ [newFile] $ emptyReplicate] $ e
  where
    fls = e^..replicates.folded.filtered ((==0) . (^.number)).files.folded._Single
    [peakFl] = filter ((==NarrowPeakFile) . (^.format)) fls
    getTFName = mk . head . B.split '+' . fromJust . bedName
    loops = case hic of
        Nothing -> Nothing
        Just x -> let [fl] = x^..replicates.folded.files.folded._Single.
                        filtered (elem "loops" . (^.tags))
                  in Just $ read3DContact $ fl^.location

findRegulators :: Monad m
               => [((B.ByteString, Int, Bool), B.ByteString)]  -- ^ Genes
               -> Maybe (Source m (BED3, BED3))  -- ^ 3D contacts
               -> [NarrowPeak]                   -- ^ Open chromatin
               -> Sink BED m [(B.ByteString, [BED])]  -- ^ Take a stream of TFBS as the input
findRegulators genes contacts peaks = do
    -- Overlap TFBS with open chromatin regions
    tfbs <- intersectBed peaks' =$= sinkList

    (assign3D, rest) <- case contacts of
        Just contacts' -> do
            let tfbsTree = bedToTree (++) $ map (\x -> (x, [x])) tfbs
            assignments <- lift $ contacts' =$=
                get3DRegulatoryDomains genes 5000 1000 =$=
                mapC ( \(bed, gene) ->
                    (gene, concat $ IM.elems $ intersecting tfbsTree bed) ) $$
                foldlC (\m (gene, tfs) ->
                    M.insertWith S.union gene (S.fromList tfs) m) M.empty
            let assigned = foldl1' S.union $ M.elems assignments
                unassigned = filter (not . (`S.member` assigned)) tfbs
            return (assignments, unassigned)
        Nothing -> return (M.empty, tfbs)

    assign2D <- fmap (M.fromListWith S.union) $ yieldMany regDomains =$=
        intersectBedWith S.fromList rest =$= filterC (not . null . snd) =$=
        mapC (first (fromJust . bedName)) $$ sinkList

    return $ M.toList $ fmap S.toList $ M.unionWith S.union assign2D assign3D
  where
    regDomains = map ( \(b, x) ->
        BED (chrom b) (chromStart b) (chromEnd b) (Just x) Nothing Nothing ) $
        getRegulatoryDomains (BasalPlusExtension 5000 1000 1000000) genes
    peaks' = flip map peaks $ \p -> let c = chromStart p + (fromJust . _npPeak) p
                                    in BED3 (chrom p) (c - 50) (c + 50)


-- | Identify active genes by overlapping their promoters with activity indicators.
gencodeActiveGenes :: BEDLike b
                   => FilePath   -- ^ gencode file in GTF format
                   -> [b]        -- ^ feature that is used to determine the activity
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
    return $ nubSort $ map (\x -> ( (chrom x, truncate $ fromJust $ bedScore x,
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
            _Single.filtered ((==["gene-TF assignment"]) . (^.tags))
    result <- decodeFile $ fl^.location :: IO [Linkage]
    let output = dir ++ "/" ++ T.unpack (e^.eid) ++ "_network.tsv"
    B.writeFile output $ B.unlines $ flip map result $ \(a,b) ->
        B.intercalate "\t" [original a, B.intercalate "," $ map original $ fst $ unzip b]
