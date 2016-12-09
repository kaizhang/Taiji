--------------------------------------------------------------------------------
-- ATAC-seq data processing
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Component.ATACSeq (builder) where

import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Utils
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.NGS
import           Bio.Pipeline.Utils (mapOfFiles)
import           Control.Lens
import           Control.Monad.IO.Class    (liftIO)
import           Data.Maybe                (catMaybes)
import qualified Data.Text                 as T
import           Scientific.Workflow

import           Constants

builder :: Builder ()
builder = do
    node "Get_ATAC_data" [| return . (^._1) |] $ do
        submitToRemote .= Just False
        label .= "Get ATAC-seq data"
    path ["Initialization", "Get_ATAC_data"]

    node "ATAC_alignment_prepare" [| \input ->
        return $ concatMap splitExpByFile $ filterExpByFile
            (\x -> formatIs FastqFile x || formatIs FastqGZip x) input
        |] $ submitToRemote .= Just False
    node "ATAC_alignment" [| \x -> bwaAlign <$> atacOutput <*>
        bwaIndex <*> return (bwaCores .= 4) <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True >> remoteParam .= "-pe smp 4"
    path [ "Get_ATAC_data", "ATAC_alignment_prepare", "ATAC_alignment"]

    node "ATAC_bam_filt_prepare" [| \(oriInput, input) -> do
        let filtInput = concatMap splitExpByFile $ filterExpByFile
                (\x -> formatIs BamFile x && not ("filtered" `elem` x^._Single.tags))
                oriInput
        return $ filtInput ++ input
        |] $ submitToRemote .= Just False
    node "ATAC_bam_filt" [| \x -> filterBam <$> atacOutput <*>
        return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    node "ATAC_remove_dups" [| \x -> removeDuplicates <$>
        getConfig' "picard" <*> atacOutput <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    [ "Get_ATAC_data", "ATAC_alignment"] ~> "ATAC_bam_filt_prepare"
    path ["ATAC_bam_filt_prepare", "ATAC_bam_filt", "ATAC_remove_dups"]

    node "ATAC_makeBED_prepare" [| \(oriInput, input) -> do
        let filtInput = concatMap splitExpByFile $ filterExpByFile
                (\x -> formatIs BamFile x && "filtered" `elem` x^._Single.tags)
                oriInput
        return $ filtInput ++ input
        |] $ submitToRemote .= Just False
    node "ATAC_makeBED" [| \e -> do
        prefix <- atacOutput
        let processWith f = replicates.traverse %%~
                id files (fmap catMaybes . mapM (f prefix)) $ e
        liftIO $ if pairedEnd e
            then processWith sortedBam2BedPE
            else processWith bam2Bed
        |] $ batch .= 1 >> stateful .= True
    node "ATAC_combine_reps_prepare" [| return . mergeExps |] $
        submitToRemote .= Just False
    node "ATAC_combine_reps" [| \x -> mergeReplicatesBed <$>
        atacOutput <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    [ "Get_ATAC_data", "ATAC_remove_dups"] ~> "ATAC_makeBED_prepare"
    path [ "ATAC_makeBED_prepare", "ATAC_makeBED", "ATAC_combine_reps_prepare"
         , "ATAC_combine_reps" ]

    node "ATAC_callpeaks_prepare" [| \(input1, input2) -> return $
        concatMap splitExpByFile $ filterExpByFile (formatIs BedGZip) $
        mergeExps $ input1 ++ input2
        |] $ submitToRemote .= Just False
    node "ATAC_callpeaks" [| mapOfFiles $ \e r fl -> do
        dir <- atacOutput
        if formatIs BedGZip fl || formatIs BedFile fl
            then do
                let output = printf "%s/%s_rep%d_p_0.01.narrowPeak" dir
                        (T.unpack $ e^.eid) (r^.number)
                result <- liftIO $ callPeaks output fl Nothing $
                    pair .= pairedEnd e >> cutoff .= PValue 0.01
                return [result]
            else return []
        |] $ batch .= 1 >> stateful .= True
    node "ATAC_IDR_prepare" [| return . mergeExps |] $ submitToRemote .= Just False
    node "ATAC_IDR" [| \e -> do
        dir <- atacOutput
        let [merged] = e^..replicates.folded.filtered (\r -> r^.number == 0)
                .files.folded.filtered (formatIs NarrowPeakFile)
            peakFiles = e^..replicates.folded.filtered (\r -> r^.number /= 0)
                .files.folded.filtered (formatIs NarrowPeakFile)
            output = printf "%s/%s_idr_0.05.narrowPeak" dir (T.unpack $ e^.eid)
        r <- liftIO $ idrMultiple peakFiles merged 0.05 output
        return $ replicates .~ [files .~ [r] $ emptyReplicate] $ e
        |] $ batch .= 1 >> stateful .= True
    ["ATAC_combine_reps_prepare", "ATAC_combine_reps"] ~> "ATAC_callpeaks_prepare"
    path ["ATAC_callpeaks_prepare", "ATAC_callpeaks", "ATAC_IDR_prepare", "ATAC_IDR"]
