--------------------------------------------------------------------------------
-- ATAC-seq data processing
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Component.ATACSeq (builder) where

import           Bio.Pipeline.NGS
import           Bio.Pipeline.CallPeaks
import           Control.Monad.IO.Class    (liftIO)
import           Control.Lens
import           Scientific.Workflow

import           Constants

builder :: Builder ()
builder = do
    node "atac00" [| return . (^._1) |] $ do
        submitToRemote .= Just False
        label .= "Get ATAC-seq data"
    path ["init00", "atac00"]

    node "align00" [| \x -> bwaAlign <$> atacSeqDir <*>
        (getConfig' "bwaIndex") <*> return (bwaCores .= 4) <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True >> remoteParam .= "-pe smp 4"
    node "align01" [| \x -> filterBam <$> atacSeqDir <*> return x
        >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    node "align02" [| \x -> removeDuplicates <$>
        getConfig' "picard" <*> atacSeqDir  <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    node "align03" [| \x -> bamToBed <$> atacSeqDir <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    node "align04" [| \x -> mergeReplicatesBed <$> atacSeqDir <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    node "peak00" [| mapM $ \x -> return (x, Nothing) |] $
        batch .= 1
    node "peak01" [| \x -> callPeaks <$> atacSeqDir <*>
        return (return ()) <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    path ["atac00", "align00", "align01", "align02", "align03", "align04", "peak00", "peak01"]
