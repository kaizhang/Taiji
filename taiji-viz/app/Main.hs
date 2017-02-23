{-# LANGUAGE BangPatterns #-}
module Main where

import           AI.Clustering.Hierarchical
import           Bio.Utils.Functions         (ihs')
import qualified Data.ByteString.Char8       as B
import qualified Data.CaseInsensitive        as CI
import           Data.Function
import qualified Data.HashMap.Strict         as HM
import           Data.List
import Data.List.Split (splitOn)
import qualified Data.Matrix                 as M
import           Data.Maybe
import           Data.Ord                    (compare, comparing)
import           Data.Semigroup              ((<>))
import qualified Data.Vector                 as V
import qualified Data.Vector.Unboxed         as U
import           Diagrams.Backend.Rasterific (renderRasterific)
import           Diagrams.Prelude (dims2D, width, height)
import           Options.Applicative
import           Statistics.Sample           (meanVarianceUnb)
import           System.Environment
import           Taiji.Visualize

data Options = Options
    { input          :: FilePath
    , output         :: FilePath
    , rowNamesFilter :: Maybe FilePath
    , colNamesFilter :: Maybe FilePath
    , expression     :: FilePath
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

defaultMain :: Options -> IO ()
defaultMain opts = do
    rank <- (fmap ihs' . readTSV) <$> B.readFile (input opts)
    expr <- (fmap ihs' . readTSV) <$> B.readFile (expression opts)

    let (labels, xs) = unzip $ map unzip $ groupBy ((==) `on` (fst.fst)) $ sort $
            HM.toList $ HM.intersectionWith (,) rank expr
        rowlab = map (B.unpack . CI.original) $ fst $ unzip $ map head labels
        collab = map (B.unpack . CI.original) $ snd $ unzip $ head $ labels
        dataTable = filterRows f $ Table rowlab collab $ M.fromLists xs

    dataTable' <- case rowNamesFilter opts of
        Nothing -> return dataTable
        Just fl -> do
            names <- lines <$> readFile fl
            return $ filterRows (filterByName names) dataTable

    --let table' = reorderColumns (orderByName ["neural-tube", "forebrain", "midbrain", "hindbrain", "heart", "intestine", "kidney", "limb", "liver", "lung", "stomach"]) $
     --       reorderRows (\x -> flatten $ hclust Ward (V.fromList x) distFn) dataTable'
    let table' = reorderColumns (sortBy (comparing (last . splitOn "_" . fst))) $
            reorderRows (\x -> flatten $ hclust Ward (V.fromList x) distFn) dataTable'
        w = width dia
        h = height dia
        n = fromIntegral $ length (head xs) * 50
        dia = spotPlot table'
    renderRasterific (output opts) (dims2D n (n*(h/w))) dia
  where
    distFn x y = euclidean (fst $ V.unzip $ snd x) (fst $ V.unzip $ snd y)
    f (_, xs) = V.any (>1e-4) xs' && case () of
        _ | n >= 5 -> let (m, v) = meanVarianceUnb $ fst $ V.unzip xs
                      in sqrt v / m > 1
          | otherwise -> V.maximum xs' / V.minimum xs' >= 2.5
      where
        n = V.length xs
        xs' = fst $ V.unzip xs


main :: IO ()
main = defaultMain =<< execParser opts
  where
    opts = info (parser <**> helper)
      ( fullDesc
     <> header "Taiji-viz - a visualization tool for Taiji" )

orderByName :: [String] -> ReodrderFn a
orderByName prefix = sortBy $ \(a,_) (b,_) ->
    let idx1 = findIdx a
        idx2 = findIdx b
    in case () of
        _ | isJust idx1 && isJust idx2 -> case compare idx1 idx2 of
                LT -> LT
                GT -> GT
                EQ -> compare a b
          | otherwise -> compare a b
  where
    findIdx x = go prefix 0
      where
        go (y:ys) !i | isPrefixOf y x = Just i
                     | otherwise = go ys (i+1)
        go _ _ = Nothing

filterByName :: [String] -> FilterFn a
filterByName xs = \(x, _) -> x `elem` xs
