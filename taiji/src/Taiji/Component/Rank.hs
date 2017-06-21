--------------------------------------------------------------------------------
-- PageRank analysis
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Taiji.Component.Rank where

import           Bio.Data.Experiment.Types
import           Bio.Pipeline.Instances            ()
import           Bio.Utils.Functions               (scale)
import           Bio.Utils.Misc                    (readDouble)
import           Conduit
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
import           IGraph.Structure                  (pagerank,
                                                    personalizedPagerank)
import           Scientific.Workflow
import           Statistics.Sample                 (meanVarianceUnb)

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
                in (name, flip map genes $ \g -> M.lookupDefault 0 g geneRanks)
            dataTable = (groupNames, map original genes, transpose ranks)

        dir <- rankOutput
        liftIO $ outputData dir dataTable $ getMetrics dataTable
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
         -> IO [(GeneName, Double)]
pageRank expr e = do
    gr <- buildNet e
    let labs = map (nodeLab gr) $ nodes gr
        tfs = S.fromList $ filter (not . null . pre gr) $ nodes gr
        ranks = case expr of
            Just expr' ->
                let lookupExpr x = M.lookupDefault (0.01,-10) x expr'
                    nodeWeights = map (exp . snd . lookupExpr) labs
                    edgeWeights = map (sqrt . fst . lookupExpr . nodeLab gr . snd) $
                        edges gr
                in personalizedPagerank gr nodeWeights (Just edgeWeights) 0.85
            Nothing -> pagerank gr Nothing 0.85
    return $ flip mapMaybe (zip [0..] ranks) $ \(i, rank) ->
        if i `S.member` tfs
            then Just (nodeLab gr i, rank)
            else Nothing

buildNet :: Experiment e => e -> IO (LGraph D GeneName ())
buildNet e = do
    let [fl] = e^..replicates.folded.filtered ((==0) . (^.number)).files.folded.
            _Single.filtered ((==["gene-TF assignment"]) . (^.tags))
    result <- decodeFile $ fl^.location :: IO [Linkage]
    return $ fromLabeledEdges $ flip concatMap result $ \(a, b) ->
        zip (zip (repeat a) $ fst $ unzip b) $ repeat ()

data Metrics = Metrics
    { averageRank   :: Double
    , variability   :: Double
    , logFoldChange :: Double
    }

getMetrics :: ([T.Text], [B.ByteString], [[Double]]) -> [(B.ByteString, Metrics)]
getMetrics (cts, geneNames, dat) = zip geneNames $ map (f . U.fromList) dat
  where
    f xs = let (m, v) = meanVarianceUnb xs
               fold = log $ max 1 $ U.maximum xs / U.minimum xs
           in Metrics m (sqrt v / m) fold

outputData :: FilePath
           -> ([T.Text], [B.ByteString], [[Double]])
           -> [(B.ByteString, Metrics)] -> IO ()
outputData dir (cts, geneNames, dat) metrics = do
    let filteredData = map toBS $ filter f $ zip metrics dat
        allData = map toBS $ zip metrics dat
        header = B.pack $ T.unpack $ T.intercalate "\t" $ "Gene" : cts
    B.writeFile (dir++"GeneRank_filtered.tsv") $ B.unlines $ header : filteredData
    B.writeFile (dir++"GeneRank_all.tsv") $ B.unlines $ header : allData
  where
    f ((_, Metrics m v fold), _) = m > 0.0001 && fold > 1.5
    toBS ((nm, _), xs) = B.intercalate "\t" $ nm : map toShortest xs
