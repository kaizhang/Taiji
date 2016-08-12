{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE OverloadedStrings #-}
module Builder (graph) where

import           Bio.Data.Bed
import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Utils
import           Bio.GO.GREAT
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Instances
import           Bio.Pipeline.NGS
import           Bio.Pipeline.ScanMotifs
import           Bio.Utils.Misc            (readInt)
import           Control.Arrow             (first, second, (&&&), (***))
import           Control.Lens
import           Data.Binary               (decodeFile, encodeFile)
import qualified Data.ByteString.Char8     as B
import           Data.Dynamic
import           Data.Function             (on)
import           Data.List
import           Data.Maybe
import           Data.Ord
import qualified Data.Text                 as T
import           Scientific.Workflow
import qualified Data.HashMap.Strict as M

import Config
import Assign
import Network

readData :: () -> IO [Experiment ATAC_Seq]
readData _ = fmap fromJust $ readExp $ config!"input"

graph :: Builder ()
graph = do
    node "init00" 'readData $ do
        submitToRemote .= Just False
        label .= "Parse metadata"

    node "align00" [| bwaAlign (config!"outputDir") (config!"genome") (return ()) |] $
        batch .= 1
    node "align01" [| filterBam (config!"outputDir") |] $ batch .= 1
    node "align02" [| removeDuplicates "/home/kai/software/picard-tools-1.140/picard.jar" (config!"outputDir") |] $
        batch .= 1
    node "align03" [| bamToBed (config!"outputDir") |] $
        batch .= 1
    path ["init00", "align00", "align01", "align02", "align03"]

    node "peak00" [| mapM $ \x -> return (x, Nothing) |] $
        batch .= 1
    node "peak01" [| callPeaks (config!"outputDir") (return ()) |] $
        batch .= 1
    node "peak02" [| \exps -> do
        let f x = map (^.location) $ filter ((==NarrowPeakFile) . (^.format)) $ x^.files
        scanMotifs (config!"genomeIndex") (config!"motifFile") 1e-5
            ((config!"outputDir") ++ "/TFBS.bed") (concatMap f (exps :: [Experiment ATAC_Seq]))
        |] $ return ()
    path ["align03", "peak00", "peak01", "peak02"]

    node "ass00" [| getDomains (config!"outputDir") |] $ batch .= 1
    ["peak01"] ~> "ass00"
    node "ass01" [| \(x,y) -> return $ zip x $ repeat y |] $ do
        label .= "prepare input"
        submitToRemote .= Just False
    ["ass00", "peak02"] ~> "ass01"
    node "ass02" [| mapM (linkGeneToTFs $ config!"outputDir") |] $ batch .= 1
    node "ass03" [| mapM (printEdgeList $ config!"outputDir") |] $ batch .= 1
    path ["ass01", "ass02", "ass03"]

    node "net00" 'pageRank $ return ()
    node "net01" [| writeTSV "ranks.tsv" |] $ submitToRemote .= Just False
    path ["ass02", "net00", "net01"]
