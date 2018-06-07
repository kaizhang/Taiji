--------------------------------------------------------------------------------
-- PageRank analysis
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE DataKinds #-}

module Taiji.Component.Rank where

import           Bio.Data.Experiment.Types
import           Bio.Pipeline.Instances            ()
import           Bio.Utils.Functions               (scale)
import           Bio.Utils.Misc                    (readDouble)
import           Conduit
import Control.Monad (replicateM)
import Control.Monad.ST (runST)
import           Control.Lens                      hiding (pre)
import           Data.Binary                       (decodeFile)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (mk, original)
import           Data.Double.Conversion.ByteString (toShortest)
import qualified Data.HashMap.Strict               as M
import qualified Data.HashSet                      as S
import           Data.List
import           Data.List.Ordered                 (nubSort)
import           Data.Maybe
import qualified Data.Text                         as T
import qualified Data.Vector.Unboxed               as U
import           IGraph
import           IGraph.Algorithms (pagerank, rewireEdges)
import           Scientific.Workflow
import Data.Vector.Algorithms.Search

import           Taiji.Constants
import           Taiji.Types

builder :: Builder ()
builder = do
    node "PageRank_prepare" [| \(exps, gene_expr) -> case gene_expr of
        Nothing -> return $ map (\x -> (Nothing, x)) exps
        Just e -> do
            rnaseqData <- readExpression e 1
            return $ flip map exps $ \x ->
                (fmap M.toList $ lookup (B.pack $ T.unpack $ fromJust $ x^.groupName) rnaseqData, x)
        |] $ do
            submitToRemote .= Just False
            note .= "Prepare to run PageRank algorithm."
    node "PageRank" [| \(gene_expr, e) -> do
        r <- pageRank (fmap M.fromList gene_expr) e
        return (fromJust $ e^.groupName, r)
        |] $ do
            batch .= 1
            note .= "Running personalized PageRank algorithm."
    node "Output_ranks" [| \results -> do
        let genes = nubSort $ concatMap (fst . unzip) $ snd $ unzip results
            (groupNames, ranks) = unzip $ flip map results $ \(name, xs) ->
                let geneRanks = M.fromList xs
                in (name, flip map genes $ \g -> M.lookupDefault (0,1) g geneRanks)
        dir <- rankOutput
        let filename1 = dir ++ "/TFRanks.tsv"
        liftIO $ outputData filename1 (groupNames, map original genes
            , (map.map) fst $ transpose ranks)

        let filename2 = dir ++ "/TFRanks_P_values.tsv"
        liftIO $ outputData filename2 (groupNames, map original genes
            , (map.map) snd $ transpose ranks)

        return filename1
        |] $ do
            submitToRemote .= Just False
            stateful .= True
            note .= "Save TF ranks to files."
    ["Link_TF_gene", "Output_expression"] ~> "PageRank_prepare"
    path ["PageRank_prepare", "PageRank", "Output_ranks"]


-- | Read RNA expression data
readExpression :: FilePath
               -> Double    -- ^ Threshold to call a gene as non-expressed
               -> IO [( B.ByteString    -- ^ cell type
                     , M.HashMap GeneName (Double, Double)  -- ^ absolute value and z-score
                     )]
readExpression fl cutoff = do
    c <- B.readFile fl
    let ((_:header):dat) = map (B.split '\t') $ B.lines c
        rowNames = map (mk . head) dat
        dataTable = map (map readDouble . tail) dat
    return $ zipWith (\a b -> (a, M.fromList $ zip rowNames b)) header $
        transpose $ zipWith zip dataTable $ map computeZscore dataTable
  where
    computeZscore xs
        | all (<cutoff) xs || all (==head xs) xs = replicate (length xs) (-10)
        | otherwise = U.toList $ scale $ U.fromList xs

pageRank :: Experiment e
         => Maybe (M.HashMap GeneName (Double, Double))   -- ^ Expression data
         -> e
         -> IO [(GeneName, (Double, Double))]
pageRank expr e = do
    gr <- buildNet e
    let labs = map (nodeLab gr) $ nodes gr
        tfs = S.fromList $ filter (not . null . pre gr) $ nodes gr
    (ranks, randRanks) <- case expr of
        Just expr' -> do
            let lookupExpr x = M.lookupDefault (0.01,-10) x expr'
                nodeWeights = exp . snd . lookupExpr
                gr' = emap (\((_, to), _) -> sqrt $ fst $ lookupExpr $ nodeLab gr to) gr
                pr = pagerank gr' 0.85 (Just nodeWeights) (Just id)
            mg <- thaw gr'
            rpr <- replicateM 5 $ do
                rewireEdges mg 1 False False
                g <- unsafeFreeze mg
                return $ pagerank g 0.85 (Just nodeWeights) (Just id)
            return (pr, U.fromList $ sort $ concat rpr)
        Nothing -> do
            let pr = pagerank gr 0.85 Nothing Nothing
            mg <- thaw gr
            rpr <- replicateM 10 $ do
                rewireEdges mg 1 False False
                g <- unsafeFreeze mg
                return $ pagerank g 0.85 Nothing Nothing
            return (pr, U.fromList $ sort $ concat rpr)
    return $ flip mapMaybe (zip [0..] ranks) $ \(i, rank) ->
        if i `S.member` tfs
            then Just ( nodeLab gr i
                      , ( rank, 1 - (fromIntegral $ bisect randRanks rank) /
                            fromIntegral (U.length randRanks) ) )
            else Nothing

bisect :: U.Vector Double -> Double -> Int
bisect v x = runST $ do
    v' <- U.unsafeThaw v
    binarySearch v' x

buildNet :: Experiment e => e -> IO (Graph 'D GeneName ())
buildNet e = do
    let [fl] = e^..replicates.folded.filtered ((==0) . (^.number)).files.folded.
            _Single.filtered ((==["gene-TF assignment"]) . (^.tags))
    result <- decodeFile $ fl^.location :: IO [Linkage]
    return $ fromLabeledEdges $ flip concatMap result $ \(a, b) ->
        zip (zip (repeat a) $ fst $ unzip b) $ repeat ()

outputData :: FilePath
           -> ([T.Text], [B.ByteString], [[Double]])
           -> IO ()
outputData filename (cts, geneNames, dat) = do
    let allData = zipWith toBS geneNames dat
        header = B.pack $ T.unpack $ T.intercalate "\t" $ "Gene" : cts
    B.writeFile filename $ B.unlines $ header : allData
  where
    toBS nm xs = B.intercalate "\t" $ nm : map toShortest xs
