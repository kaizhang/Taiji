--------------------------------------------------------------------------------
-- ATAC-seq data processing
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Taiji.Component.ATACSeq (builder) where

import           Bio.Data.Bed              (chrom, chromEnd, chromStart)
import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Utils
import           Bio.Pipeline.NGS
import           Control.Lens
import           Control.Monad.IO.Class    (liftIO)
import           Data.Maybe                (catMaybes)
import           Scientific.Workflow

import           Taiji.Constants
import           Taiji.Tools

builder :: Builder ()
builder = do
    node "Get_ATAC_data" [| \input -> do
        dir <- downloadOutput
        liftIO $ mapM (downloadData dir) $ input^._1
        |] $ do
            submitToRemote .= Just False
            stateful .= True
            note .= "Extract (or download) ATAC-seq data."
    path ["Initialization", "Get_ATAC_data"]

    node "ATAC_alignment_prepare" [| \input ->
        return $ concatMap splitExpByFile $ filterExpByFile
            (\x -> formatIs FastqFile x || formatIs FastqGZip x) input
        |] $ do
            submitToRemote .= Just False
            note .= "Grab FASTQ files; Prepare for ATAC-seq alignment."
    node "ATAC_alignment" [| \x -> bwaAlign <$> atacOutput <*>
        bwaIndex <*> return (bwaCores .= 4) <*> return x >>= liftIO
        |] $ do
            batch .= 1
            stateful .= True
            -- remoteParam .= "-pe smp 4"
            remoteParam .= "--ntasks-per-node=4"
            note .= "Input: FASTQ.\nOutput: BAM.\nUse BWA to align ATAC-seq data. Default options: \"-l 32 -k 2 -q 5\"."
    path [ "Get_ATAC_data", "ATAC_alignment_prepare", "ATAC_alignment"]

    node "ATAC_bam_filt_prepare" [| \(oriInput, input) -> do
        let filtInput = concatMap splitExpByFile $ filterExpByFile
                (\x -> formatIs BamFile x && not ("filtered" `elem` x^._Single.tags))
                oriInput
        return $ filtInput ++ input
        |] $ do
            submitToRemote .= Just False
            note .= "Prepare for BAM file filtering."
    node "ATAC_bam_filt" [| \x -> filterBam <$> atacOutput <*>
        return x >>= liftIO
        |] $ do
            batch .= 1
            stateful .= True
            note .= "Remove low quality reads using samtools. Default options: \"-F 0x70c -q 30\"."
    node "ATAC_remove_dups" [| \x -> removeDuplicates <$>
        getConfig' "picard" <*> atacOutput <*> return x >>= liftIO
        |] $ do
            batch .= 1
            stateful .= True
            note .= "Remove duplicated reads using picard."
    [ "Get_ATAC_data", "ATAC_alignment"] ~> "ATAC_bam_filt_prepare"
    path ["ATAC_bam_filt_prepare", "ATAC_bam_filt", "ATAC_remove_dups"]

    node "ATAC_makeBED_prepare" [| \(oriInput, input) -> do
        let filtInput = concatMap splitExpByFile $ filterExpByFile
                (\x -> formatIs BamFile x && "filtered" `elem` x^._Single.tags)
                oriInput
        return $ filtInput ++ input
        |] $ do
            submitToRemote .= Just False
            note .= "Prepare for BED file conversion."
    node "ATAC_makeBED" [| \e -> do
        prefix <- atacOutput
        let processWith f = replicates.traverse %%~
                id files (fmap catMaybes . mapM (f prefix)) $ e
        liftIO $ if pairedEnd e
            then processWith ( \x -> sortedBam2BedPE x
                (\(a,b) -> abs (chromStart a - chromEnd b) < 5000 &&
                chrom a /= "chrM" && chrom b /= "chrM") )
            else processWith (\x -> bam2Bed x ((/="chrM") . chrom))
        |] $ do
            batch .= 1
            stateful .= True
            note .= "Convert BAM files to BED files."
    ["Get_ATAC_data", "ATAC_remove_dups"] ~> "ATAC_makeBED_prepare"
    path [ "ATAC_makeBED_prepare", "ATAC_makeBED"]

    node "ATAC_combine_reps_prepare" [| \(oriInput, es) -> do
        let filtInput = filterExpByFile
                (\x -> formatIs BedGZip x || formatIs BedFile x) oriInput
        return $ mergeExps $ es ++ filtInput
        |] $ do
            submitToRemote .= Just False
            note .= "Prepare for combining BED files."
    node "ATAC_combine_reps" [| \x -> mergeReplicatesBed <$>
        atacOutput <*> return x >>= liftIO
        |] $ do
            batch .= 1
            stateful .= True
            note .= "Concatenate BED files from multiple replicates into a single file."
    ["Get_ATAC_data", "ATAC_makeBED"] ~> "ATAC_combine_reps_prepare"
    path ["ATAC_combine_reps_prepare", "ATAC_combine_reps"]

    node "ATAC_callpeaks_prepare" [| \(input1, input2) -> return $
        concatMap splitExpByFile $ filterExpByFile (formatIs BedGZip) $
        mergeExps $ input1 ++ input2
        |] $ do
            submitToRemote .= Just False
            note .= "Prepare for peak calling."
    ["ATAC_combine_reps_prepare", "ATAC_combine_reps"] ~> "ATAC_callpeaks_prepare"
