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

    node "ATAC_alignment" [| \x -> do
        x' <- bwaAlign <$> atacOutput <*> bwaIndex <*> return (bwaCores .= 4) <*>
            return x >>= liftIO
        let fn (Single fl) = fl^.format == BamFile
            fn _ = False
        return $ filterExpByFile fn $ mergeExps $ x ++ x'
        |] $ batch .= 1 >> stateful .= True >> remoteParam .= "-pe smp 4"
    node "ATAC_bam_filtering" [| \x -> filterBam <$> atacOutput <*> return x
        >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    node "ATAC_remove_dups" [| \x -> removeDuplicates <$>
        getConfig' "picard" <*> atacOutput <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    node "ATAC_make_bed" [| \x -> bam2Bed <$> atacOutput <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    node "ATAC_combine_reps" [| \x -> mergeReplicatesBed <$> atacOutput <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    node "ATAC_callpeaks_prepare" [| mapM $ \x -> return (x, Nothing) |] $
        submitToRemote .= Just False
    node "ATAC_callpeaks" [| \x -> callPeaks <$> atacOutput <*>
        return (return ()) <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    path [ "Get_ATAC_data", "ATAC_alignment", "ATAC_bam_filtering"
         , "ATAC_remove_dups", "ATAC_make_bed", "ATAC_combine_reps"
         , "ATAC_callpeaks_prepare", "ATAC_callpeaks" ]
