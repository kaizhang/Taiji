-- This module contains network-related analysis
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Network where

import           Bio.Data.Bed
import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Utils
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Instances
import           Bio.Pipeline.NGS
import           Bio.Pipeline.ScanMotifs
import           Control.Lens
import           Control.Monad
import           Data.Binary                       (decodeFile, encodeFile)
import qualified Data.ByteString.Char8             as B
import           Data.Double.Conversion.ByteString (toShortest)
import           Data.Function                     (on)
import qualified Data.HashMap.Strict               as M
import           Data.List
import           Data.List.Ordered                 (nubSort)
import           Data.Maybe
import           Data.Ord
import qualified Data.Text                         as T
import           IGraph
import           IGraph.Structure                  (pagerank,
                                                    personalizedPagerank)
import           Scientific.Workflow

import           Expression

pageRank :: [Experiment a] -> IO ([T.Text], [B.ByteString], [[Double]])
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
    return (expNames, genes, transpose ranks)

personalizedPageRank :: (FilePath, [Experiment a])
                     -> IO ([T.Text], [B.ByteString], [[Double]])
personalizedPageRank (rnaseq, es) = do
    rnaseqData <- readExpression rnaseq

    results <- forM es $ \e -> do
        gr <- buildNet e
        let labs = map (nodeLab gr) $ nodes gr
            geneExpression = fromJust $ lookup (B.pack $ T.unpack $ e^.celltype) rnaseqData
            weights = map (\x -> exp (M.lookupDefault (-10) x geneExpression)) labs
        return $ zip labs $ personalizedPagerank gr weights Nothing 0.85
    let genes = nubSort $ concatMap (fst . unzip) results
        expNames = map (^.eid) es
        ranks = flip map results $ \xs ->
            let geneRanks = M.fromList xs
            in flip map genes $ \g -> M.lookupDefault 0 g geneRanks
    return (expNames, genes, transpose ranks)

buildNet :: Experiment a -> IO (LGraph D B.ByteString ())
buildNet e = do
    let [fl] = filter (\x -> x^.keywords == ["gene-TF assignment"]) $ e^.files
    result <- decodeFile $ fl^.location :: IO [(B.ByteString, [(B.ByteString, [BED])])]
    return $ fromLabeledEdges $ flip concatMap result $ \(a, b) ->
        zip (zip (repeat a) $ fst $ unzip b) $ repeat ()

printEdgeList :: FilePath -> Experiment a -> IO ()
printEdgeList dir e = do
    let [fl] = filter (\x -> x^.keywords == ["gene-TF assignment"]) $ e^.files
    result <- decodeFile $ fl^.location :: IO [(B.ByteString, [(B.ByteString, [BED])])]
    let output = dir ++ "/" ++ T.unpack (e^.celltype) ++ "_" ++
            T.unpack (e^.target) ++ "_associations.tsv"
    B.writeFile output $ B.unlines $ flip map result $ \(a,b) ->
        B.intercalate "\t" [a, B.intercalate "," $ fst $ unzip b]

writeTSV :: FilePath -> ([T.Text], [B.ByteString], [[Double]]) -> IO ()
writeTSV output (colNames, rowNames, dat) = B.writeFile output $ B.unlines $
    B.intercalate "\t" ("Gene" : map (B.pack . T.unpack) colNames) :
    map (\(a,b) -> B.intercalate "\t" $ a : map toShortest b)
        (zip rowNames dat)
