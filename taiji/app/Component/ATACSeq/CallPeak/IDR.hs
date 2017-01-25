-- Call Peaks using IDR framework
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Component.ATACSeq.CallPeak.IDR (builder) where

import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Utils
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.NGS
import           Bio.Pipeline.Utils        (mapOfFiles)
import           Control.Lens
import           Control.Monad.IO.Class    (liftIO)
import           Data.Maybe                (catMaybes)
import qualified Data.Text                 as T
import           Scientific.Workflow
import           Shelly                    (fromText, mkdir_p, shelly)

import           Constants


builder :: Builder ()
builder = do
    node "ATAC_IDR_call_relaxed_peaks" [| mapOfFiles $ \e r fl -> do
        dir <- fmap (++"/Peaks_Relax/") atacOutput
        shelly $ mkdir_p $ fromText $ T.pack dir
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
    node "ATAC_callpeaks" [| \e -> do
        dir <- atacOutput
        let [merged] = e^..replicates.folded.filtered (\r -> r^.number == 0)
                .files.folded.filtered (formatIs NarrowPeakFile)
            peakFiles = e^..replicates.folded.filtered (\r -> r^.number /= 0)
                .files.folded.filtered (formatIs NarrowPeakFile)
            output = printf "%s/%s_idr_0.05.narrowPeak" dir (T.unpack $ e^.eid)
        r <- liftIO $ idrMultiple peakFiles merged 0.05 output
        return $ replicates .~ [files .~ [r] $ emptyReplicate] $ e
        |] $ batch .= 1 >> stateful .= True
    path [ "ATAC_callpeaks_prepare", "ATAC_IDR_call_relaxed_peaks"
        , "ATAC_IDR_prepare", "ATAC_callpeaks" ]
