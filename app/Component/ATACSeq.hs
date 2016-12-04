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
import           Control.Lens
import           Control.Monad.IO.Class    (liftIO)
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
    [ "Get_ATAC_data", "ATAC_remove_dups"] ~> "ATAC_makeBED_prepare"

    node "ATAC_makeBED" [| \x -> bam2Bed <$> atacOutput <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True

    node "ATAC_combine_reps_prepare" [| return . mergeExps |] $
        submitToRemote .= Just False
    node "ATAC_combine_reps" [| \x -> mergeReplicatesBed <$>
        atacOutput <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True

    node "ATAC_callpeaks_prepare" [| mapM $ \x -> return (x, Nothing) |] $
        submitToRemote .= Just False
    node "ATAC_callpeaks" [| \x -> callPeaks <$> atacOutput <*>
        return (return ()) <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True

    path [ "ATAC_makeBED_prepare", "ATAC_makeBED", "ATAC_combine_reps_prepare"
         , "ATAC_combine_reps", "ATAC_callpeaks_prepare", "ATAC_callpeaks" ]
