{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE OverloadedStrings #-}
module Builder (graph) where

import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Utils
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Instances
import           Bio.Pipeline.NGS
import           Bio.Pipeline.ScanMotifs
import           Control.Arrow             (first, second, (&&&), (***))
import           Control.Lens
import qualified Data.ByteString.Char8     as B
import           Data.Function             (on)
import           Data.List
import           Data.Maybe
import Data.Aeson.Types (Result(..), fromJSON)
import Data.Yaml (Object, decodeFile)
import           Data.Ord
import qualified Data.Text                 as T
import           Scientific.Workflow hiding (Success)
import qualified Data.HashMap.Strict as M

import Config
import Assign
import Network

readData :: () -> IO ( [Experiment ATAC_Seq]
                     , [Experiment ChIP_Seq]
                     , [Experiment RNA_Seq] )
readData _ = do
    Just dat <- decodeFile $ config!"input" :: IO (Maybe Object)
    return ( parse $ M.lookup "atac-seq" dat
           , parse $ M.lookup "chip-seq" dat
           , parse $ M.lookup "rna-seq" dat
           )
  where
    parse x = case x of
        Nothing -> []
        Just x' -> case fromJSON x' of
            Error msg -> error msg
            Success r -> r

graph :: Builder ()
graph = do
    node "init00" 'readData $ do
        submitToRemote .= Just False
        label .= "Parse metadata"

    node "atac00" [| \x -> return $ x^._1 |] $ do
        submitToRemote .= Just False
        label .= "Get ATAC-seq data"
    path ["init00", "atac00"]

    node "align00" [| bwaAlign (config!"outputDir") (config!"genome") (return ()) |] $
        batch .= 1
    node "align01" [| filterBam (config!"outputDir") |] $ batch .= 1
    node "align02" [| removeDuplicates "/home/kai/software/picard-tools-1.140/picard.jar" (config!"outputDir") |] $
        batch .= 1
    node "align03" [| bamToBed (config!"outputDir") |] $
        batch .= 1
    path ["atac00", "align00", "align01", "align02", "align03"]

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
