{-# LANGUAGE FlexibleContexts #-}
module Main where

import           AI.Clustering.Hierarchical
import           Bio.Utils.Functions        (ihs', scale)
import           Control.Arrow              (first)
import Control.Monad
import qualified Data.ByteString.Char8      as B
import qualified Data.CaseInsensitive       as CI
import           Data.Function
import qualified Data.HashMap.Strict        as HM
import           Data.List                  (groupBy, sort, sortBy)
import           Data.List.Split            (splitOn)
import qualified Data.Matrix                as M
import           Data.Maybe
import           Data.Ord                   (compare, comparing)
import           Data.Semigroup             ((<>))
import qualified Data.Vector                as V
import qualified Data.Vector.Generic        as G
import qualified Data.Vector.Unboxed        as U
import           Diagrams.Backend.Cairo     (renderCairo)
import           Diagrams.Prelude           (dims2D, height, width)
import           Options.Applicative
import           Statistics.Sample          (meanVarianceUnb)
import           System.Environment

import           Taiji.Visualize
import           Taiji.Visualize.Data

data Options = Options
    { input          :: FilePath
    , output         :: FilePath
    , rowNamesFilter :: Maybe FilePath
    , colNamesFilter :: Maybe FilePath
    , expression     :: FilePath
    , pValue         :: Double
    }

parser :: Parser Options
parser = Options
    <$> strArgument
      ( metavar "INPUT"
     <> help "PageRank result" )
    <*> strArgument (metavar "OUTPUT")
    <*> (optional . strOption)
      ( long "rowNamesFilter" )
    <*> (optional . strOption)
      ( long "colNamesFilter" )
    <*> strOption
      ( long "expression" )
    <*> option auto
      ( long "p-value"
     <> value 1e-5
     <> help "P-value for calling cell-type-specific TFs. (default: 1e-5)"
      )

defaultMain :: Options -> IO ()
defaultMain opts = do
    Table r c oriData <- fmap colReorder $ readData (input opts) $ expression opts

    let pvalues = M.fromRows $ map (G.convert . pooledPValue grps . G.convert) $
           M.toRows $ fst $ M.unzip oriData
        grps = [4,6,6,6,6,4,4,4,6,4,4]
    --let pvalues = M.fromRows $ map (G.convert . pValueGaussian . G.convert) $
     --       M.toRows $ fst $ M.unzip oriData
        normalizedData = uncurry M.zip $
            first (M.fromRows . map scale . M.toRows) $ M.unzip oriData
        dataTable = Table r c $ M.zip normalizedData pvalues

    dataTable' <- case rowNamesFilter opts of
        Nothing -> return dataTable
        Just fl -> do
            names <- lines <$> readFile fl
            return $ filterRows (filterByName names) dataTable

    forM_ (zip3 (scanl1 (+) $ 0 : grps) grps ["neural-tube", "forebrain", "midbrain", "hindbrain", "heart", "intestine", "kidney", "limb", "liver", "lung", "stomach"]) $ \(a,b,o) -> do

        let table' = rowReorder $
               filterRows (\(_,x) -> V.any ((<=0.01) . snd) $ V.take b $ V.drop a x) $ dataTable'
        -- let table' = reorderColumns (sortBy (comparing (last . splitOn "_" . fst))) $
        -- let table' = reorderColumns (orderByName ["neural-tube", "forebrain", "midbrain", "hindbrain", "heart", "intestine", "kidney", "limb", "liver", "lung", "stomach"]) $
        --       reorderRows (\x -> flatten $ hclust Ward (V.fromList x) distFn) dataTable'
        -- let table' = reorderColumns (sortBy (comparing (last . splitOn "_" . fst))) $
        --        reorderRows (\x -> flatten $ hclust Ward (V.fromList x) distFn) dataTable'
            w = width dia
            h = height dia
            n = fromIntegral $ M.cols (matrix table') * 50
            dia = spotPlot (pValue opts) table'
            -- dia = heatmap table'
        renderCairo (o ++ ".svg") (dims2D n (n*(h/w))) dia
  where
    colReorder = reorderColumns (orderByName ["neural-tube", "forebrain", "midbrain", "hindbrain", "heart", "intestine", "kidney", "limb", "liver", "lung", "stomach"])
    rowReorder = reorderRows (\x -> flatten $ hclust Ward (V.fromList x) distFn)
      where
        distFn x y = euclidean (fst $ V.unzip $ fst $ V.unzip $ snd x)
            (fst $ V.unzip $ fst $ V.unzip $ snd y)

main :: IO ()
main = defaultMain =<< execParser opts
  where
    opts = info (parser <**> helper)
      ( fullDesc
     <> header "Taiji-viz - a visualization tool for Taiji" )
