--------------------------------------------------------------------------------
-- Identification of transcription factor binding sites
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Component.TFBS (builder) where

import           Bio.Data.Bed
import           Bio.Data.Experiment.Types
import           Bio.Pipeline.ScanMotifs   (scanMotifs)
import           Conduit
import           Control.Lens
import           Control.Monad
import           Control.Monad.IO.Class    (liftIO)
import qualified Data.Text                 as T
import           Scientific.Workflow

import           Constants

builder :: Builder ()
builder = do
    node "peak02" [| \xs -> do
        let f x = map (^.location) $ filter ((==NarrowPeakFile) . (^.format)) $ x^.files
        scanMotifs <$> (getConfig' "seqIndex") <*> (getConfig' "motifFile") <*>
            return 1e-5 <*> ((++ "/TFBS_open_chromatin_union.bed") <$> tfbsDir ) <*>
            return (concatMap f xs) >>= liftIO
        |] $ stateful .= True
    ["peak01"] ~> "peak02"

    node "peak03" [| \(es, tfbs) -> do
        dir <- tfbsDir
        liftIO $ forM_ es $ \e -> do
            let [fl] = map (^.location) $
                    filter ((==NarrowPeakFile) . (^.format)) $ e^.files
                output = dir ++ "/" ++ T.unpack (e^.eid) ++ "_tfbs.bed"
            peaks <- readBed' fl :: IO [BED3]
            (readBed tfbs :: Source IO BED) =$= intersectBed peaks $$
                writeBed output
        |] $ stateful .= True
    ["peak01", "peak02"] ~> "peak03"
