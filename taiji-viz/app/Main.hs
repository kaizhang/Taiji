{-# LANGUAGE FlexibleContexts #-}
module Main where

import           AI.Clustering.Hierarchical
import           Bio.Utils.Functions        (scale)
import           Control.Arrow              (first)
import qualified Data.Matrix                as M
import           Data.Semigroup             ((<>))
import qualified Data.Vector                as V
import qualified Data.Vector.Generic        as G
import           Diagrams.Backend.Cairo     (renderCairo)
import           Diagrams.Prelude           (dims2D, height, width)
import           Options.Applicative

import           Taiji.Visualize
import           Taiji.Visualize.Data

import qualified Data.ByteString.Char8 as B
import Data.List (intercalate)
import Control.Monad

data Options = Options
    { input          :: FilePath
    , output         :: FilePath
    , rowNamesFilter :: Maybe FilePath
    , colNamesFilter :: Maybe FilePath
    , expression     :: FilePath
    , pValue         :: Double
    , cv             :: Double
    , minRank        :: Double
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
    <*> option auto
      ( long "cv"
     <> value 1
     <> help "TFs with coefficient of variance less than the specified value will be removed. (default: 1)"
      )
    <*> option auto
      ( long "min"
     <> value 1e-4
     <> help "lowerBound of TF rank."
      )

defaultMain :: Options -> IO ()
defaultMain opts = do
    Table r c oriData <- fmap colReorder $ readData (input opts)
        (expression opts) $ DataFiltOpts (>=(cv opts)) (minRank opts)

        {-
    let pvalues = M.fromRows $ map (G.convert . pooledPValue grps . G.convert) $
           M.toRows $ fst $ M.unzip oriData
        grps = [6, 8, 8, 8, 8, 4, 4, 6, 7, 4, 5, 4]
        -}

    let pvalues = M.fromRows $ map (G.convert . pValueGaussian . G.convert) $
            M.toRows $ fst $ M.unzip oriData
        normalizedData = uncurry M.zip $
            first (M.fromRows . map scale . M.toRows) $ M.unzip oriData
        dataTable = Table r c $ M.zip normalizedData pvalues

    dataTable' <- case rowNamesFilter opts of
        Nothing -> return dataTable
        Just fl -> do
            names <- lines <$> readFile fl
            return $ filterRows (filterByName names) dataTable

    let table' = rowReorder dataTable'
        -- let table' = reorderColumns (sortBy (comparing (last . splitOn "_" . fst))) $
        w = width dia
        h = height dia
        n = fromIntegral $ M.cols (matrix table') * 50
        dia = spotPlot (pValue opts) table'
    renderCairo (output opts) (dims2D n (n*(h/w))) dia

    -- Output processed data table
    B.putStrLn $ B.pack $ intercalate "\t" $ "TF" : colNames table'
    forM_ (zip (rowNames table') $ M.toLists $ matrix table') $ \(x,y) -> do
        B.putStrLn $ B.pack $ intercalate "\t" $ x : map (show . fst . fst) y
    -- output pvalue
    -- B.putStrLn $ B.pack $ intercalate "\t" $ "TF" : colNames table'
    -- forM_ (zip (rowNames table') $ M.toLists $ matrix table') $ \(x,y) -> do
        -- B.putStrLn $ B.pack $ intercalate "\t" $ x : map (show . snd) y
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
