-- This module contains network-related analysis
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE OverloadedStrings #-}

module Network where

import           Bio.Data.Bed
import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Utils
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Instances
import           Bio.Pipeline.NGS
import           Bio.Pipeline.ScanMotifs
import Control.Monad
import           Control.Lens
import           Data.Binary               (decodeFile, encodeFile)
import qualified Data.ByteString.Char8     as B
import           Data.Function             (on)
import           Data.List
import qualified Data.HashMap.Strict                  as M
import           Data.Maybe
import           Data.Ord
import Data.List.Ordered (nubSort)
import qualified Data.Text                 as T
import           IGraph
import           IGraph.Structure
import Data.Double.Conversion.ByteString (toShortest)

import           Scientific.Workflow

pageRank :: [Experiment] -> IO ([T.Text], [B.ByteString], [[Double]])
pageRank es = do
    results <- forM es $ \e -> do
        gr <- buildNet e
        let labs = map (nodeLab gr) $ nodes gr
        return $ zip labs $ pagerank gr Nothing 0.85
    let genes = nubSort $ concatMap (fst . unzip) results
        expNames = map (^.eid) es
        ranks = flip map results $ \xs ->
            let geneRanks = M.fromList xs
            in flip map genes $ \g -> M.lookupDefault 0 g geneRanks
    return (expNames, genes, ranks)

buildNet :: Experiment -> IO (LGraph D B.ByteString ())
buildNet e = do
    let [fl] = filter (\x -> x^.keywords == ["gene-TF assignment"]) $ e^.files
    result <- decodeFile $ fl^.location :: IO [(B.ByteString, [(B.ByteString, [BED])])]
    return $ fromLabeledEdges $ flip concatMap result $ \(a, b) ->
        zip (zip (repeat a) $ fst $ unzip b) $ repeat ()

writeTSV :: FilePath -> ([T.Text], [B.ByteString], [[Double]]) -> IO ()
writeTSV output (colNames, rowNames, dat) = B.writeFile output $ B.unlines $
    B.intercalate "\t" ("Gene" : map (B.pack . T.unpack) colNames) :
    map (\(a,b) -> B.intercalate "\t" $ a : map toShortest b)
        (zip rowNames $ transpose dat)


graph :: Builder ()
graph = do
    node "net00" 'pageRank $ return ()
    node "net01" [| writeTSV "ranks.tsv" |] $ submitToRemote .= Just False
    path ["ex10", "net00", "net01"]
