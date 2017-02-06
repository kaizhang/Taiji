module Main where

import           AI.Clustering.Hierarchical
import           Bio.Utils.Functions         (ihs')
import qualified Data.ByteString.Char8       as B
import           Data.Function
import qualified Data.HashMap.Strict         as M
import           Data.List
import qualified Data.Vector                 as V
import qualified Data.Vector.Unboxed         as U
import           Diagrams.Backend.Rasterific
import           Diagrams.Prelude
import           Statistics.Sample           (meanVarianceUnb)
import           System.Environment
import           Taiji.Visualize

main :: IO ()
main = do
    [f1,f2] <- getArgs
    c1 <- B.readFile f1
    c2 <- B.readFile f2

    let m1 = readTSV c1
        m2 = readTSV c2
        (labels, xs) = unzip $ map unzip $ groupBy ((==) `on` (fst.fst)) $ sort $
            M.toList $ M.intersectionWith (,) m1 m2
        rowlab = map B.unpack $ fst $ unzip $ map head labels
        collab = map B.unpack $ snd $ unzip $ head $ labels
        distFn x y = euclidean (V.fromList $ fst $ unzip $ snd x)
            (V.fromList $ fst $ unzip $ snd y)
        (rowlab', xs') = unzip $ flatten $ hclust Ward
            (V.fromList $ filter f $ zip rowlab $ map (map (\(x,y) -> (ihs' x, ihs' y))) xs)
            distFn
        (collab', xs'') = unzip $ flatten $ hclust Ward
            (V.fromList $ zip collab $ transpose xs') distFn
        w = width dia
        h = height dia
        n = fromIntegral $ length (head xs) * 50
        dia = spotPlot rowlab' collab' $ transpose xs''
    renderRasterific "out.png" (dims2D n (n*(h/w))) dia
  where
    f (_, xs) = let (m, v) = meanVarianceUnb $ V.fromList $ fst $ unzip xs
                in sqrt v / m > 1
