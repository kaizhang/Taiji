--------------------------------------------------------------------------------
-- Identification of transcription factor binding sites
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Taiji.Component.TFBS (builder) where

import           Bio.Data.Bed
import           Bio.Data.Experiment.Types
import           Bio.Motif
import           Bio.Seq.IO                (withGenome)
import           Conduit
import           Control.Lens
import           Control.Monad
import           Control.Monad.IO.Class    (liftIO)
import           Data.Default
import           Data.List.Split           (chunksOf)
import qualified Data.Text                 as T
import           Scientific.Workflow
import           Shelly                    (fromText, rm, shelly)

import           Taiji.Constants

builder :: Builder ()
builder = do
    node "Find_TF_sites_prepare" [| \es -> do
        let f x = x^..replicates.folded.filtered ((==0) . (^.number)).
                files.folded._Single.filtered ((==NarrowPeakFile) . (^.format)).
                location
        peaks <- liftIO $ mapM readBed' $ concatMap f es
        openChromatin <- (++ "/openChromatin.bed") <$> tfbsOutput
        liftIO $ mergeBed (concat (peaks :: [[BED3]])) $$ writeBed openChromatin

        motifFile <- getConfig' "motifFile"
        motifs <- liftIO $ readMEME motifFile
        return $ ContextData openChromatin $ zip [1..] $ chunksOf 200 motifs
        |] $ stateful .= True >> submitToRemote .= Just False

    node "Find_TF_sites" [| \(ContextData openChromatin (idx, motifs)) -> do
        let p = 1e-5
        output <- (++ ("/tmp." ++ show (idx::Int))) <$> tfbsOutput
        genome <- getConfig' "seqIndex"
        liftIO $ withGenome genome $ \g ->
            (readBed openChromatin :: Source IO BED3) =$=
            motifScan g motifs def p =$= getMotifScore g motifs def =$=
            getMotifPValue (Just (1 - p * 10)) motifs def $$ writeBed output
        return output
        |] $ stateful .= True >> batch .= 1
    node "Find_TF_sites_merge" [| \xs -> do
        output <- (++ "/TFBS_open_chromatin_union.bed") <$> tfbsOutput
        liftIO $ do
            runResourceT $ mapM_ sourceFileBS xs $$ sinkFile output
            shelly $ mapM_ (rm . fromText . T.pack) xs
        return output
        |] $ stateful .= True >> submitToRemote .= Just False
    path [ "ATAC_callpeaks", "Find_TF_sites_prepare", "Find_TF_sites"
         , "Find_TF_sites_merge" ]

    node "Output_TF_sites" [| \(es, tfbs) -> do
        dir <- tfbsOutput
        liftIO $ forM_ es $ \e -> do
            let [fl] = e^..replicates.folded.filtered ((==0) . (^.number)).files.
                    folded._Single.filtered ((==NarrowPeakFile) . (^.format)).location
                output = dir ++ "/" ++ T.unpack (e^.eid) ++ "_tfbs.bed"
            peaks <- readBed' fl :: IO [BED3]
            (readBed tfbs :: Source IO BED) =$= intersectBed peaks $$
                writeBed output
        |] $ stateful .= True
    ["ATAC_callpeaks", "Find_TF_sites_merge"] ~> "Output_TF_sites"
