--------------------------------------------------------------------------------
-- PageRank analysis
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Component.Rank (builder) where

import           Bio.Data.Experiment.Types
import           Bio.Pipeline.Instances            ()
import           Bio.Utils.Functions               (scale)
import           Bio.Utils.Misc                    (readDouble)
import           Conduit
import           Control.Lens                      hiding (pre)
import           Control.Monad
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
import           System.IO                         (hPutStrLn, stderr)

import           Constants
import           Types

builder :: Builder ()
builder = do
    node "PageRank" [| \(x, gene_expr) -> liftIO $ case gene_expr of
            Nothing -> do
                hPutStrLn stderr "Running PageRank..."
                pageRank x
            Just e -> do
                hPutStrLn stderr "Running personalized PageRank..."
                personalizedPageRank (e, x)
        |] $ return ()
    ["Link_TF_gene", "Output_expression"] ~> "PageRank"

    node "Output_ranks" [| \x -> do
        dir <- rankOutput
        liftIO $ outputData dir x (getMetrics x)
        |] $ submitToRemote .= Just False >> stateful .= True
    ["PageRank"] ~> "Output_ranks"


-- | Read RNA expression data
readExpression :: FilePath
               -> IO [( B.ByteString    -- ^ cell type
                     , M.HashMap GeneName (Double, Double)  -- ^ absolute value and z-score
                     )]
readExpression fl = do
    c <- B.readFile fl
    let ((_:header):dat) = map (B.split '\t') $ B.lines c
        rowNames = map (mk . head) dat
        dataTable = map (map readDouble . tail) dat
    return $ zipWith (\a b -> (a, M.fromList $ zip rowNames b)) header $
        transpose $ zipWith zip dataTable $ map computeZscore dataTable
  where
    computeZscore xs
        | all (<1) xs || all (==head xs) xs = replicate (length xs) (-10)
        | otherwise = U.toList $ scale $ U.fromList xs

pageRank :: Experiment e => [e] -> IO ([T.Text], [B.ByteString], [[Double]])
pageRank es = do
    results <- forM es $ \e -> do
        gr <- buildNet e
        let tfs = S.fromList $ filter (not . null . pre gr) $ nodes gr
        return $ flip mapMaybe (zip [0..] $ pagerank gr Nothing 0.85) $ \(i, rank) ->
            if i `S.member` tfs
                then Just (nodeLab gr i, rank)
                else Nothing
    let genes = nubSort $ concatMap (fst . unzip) results
        expNames = map (fromJust . (^.groupName)) es
        ranks = flip map results $ \xs ->
            let geneRanks = M.fromList xs
            in flip map genes $ \g -> M.lookupDefault 0 g geneRanks
    return (expNames, map original genes, transpose ranks)

personalizedPageRank :: Experiment e
                     => (FilePath, [e])
                     -> IO ([T.Text], [B.ByteString], [[Double]])
personalizedPageRank (rnaseq, es) = do
    rnaseqData <- readExpression rnaseq

    results <- forM es $ \e -> do
        gr <- buildNet e
        let lookupExpr x = M.lookupDefault (0.01,-10) x $ fromJust $ lookup
                (B.pack $ T.unpack $ fromJust $ e^.groupName) rnaseqData
            nodeWeights = map (exp . snd . lookupExpr) labs
            edgeWeights = map (sqrt . fst . lookupExpr . nodeLab gr . snd) $ edges gr
            labs = map (nodeLab gr) $ nodes gr
            tfs = S.fromList $ filter (not . null . pre gr) $ nodes gr
            ranks = personalizedPagerank gr nodeWeights (Just edgeWeights) 0.85
        return $ flip mapMaybe (zip [0..] ranks) $ \(i, rank) ->
            if i `S.member` tfs
                then Just (nodeLab gr i, rank)
                else Nothing

    let genes = nubSort $ concatMap (fst . unzip) results
        expNames = map (fromJust . (^.groupName)) es
        ranks = flip map results $ \xs ->
            let geneRanks = M.fromList xs
            in flip map genes $ \g -> M.lookupDefault 0 g geneRanks
    return (expNames, map original genes, transpose ranks)

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
