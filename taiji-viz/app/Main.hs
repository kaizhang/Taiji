module Main where

import           AI.Clustering.Hierarchical
import           Bio.Utils.Functions         (ihs')
import qualified Data.ByteString.Char8       as B
import           Data.Function
import qualified Data.HashMap.Strict         as HM
import           Data.List
import qualified Data.Matrix                 as M
import           Data.Ord                    (comparing)
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

    let m1 = fmap ihs' $ readTSV c1
        m2 = fmap ihs' $ readTSV c2
        (labels, xs) = unzip $ map unzip $ groupBy ((==) `on` (fst.fst)) $ sort $
            HM.toList $ HM.intersectionWith (,) m1 m2
        rowlab = map B.unpack $ fst $ unzip $ map head labels
        collab = map B.unpack $ snd $ unzip $ head $ labels
        table = filterRows f $ Table rowlab collab $ M.fromLists xs
        distFn x y = euclidean (fst $ V.unzip $ snd x) (fst $ V.unzip $ snd y)
        table' = reorderColumns (sortBy (comparing fst)) $
            reorderRows (\x -> flatten $ hclust Ward (V.fromList x) distFn) table
        w = width dia
        h = height dia
        n = fromIntegral $ length (head xs) * 50
        dia = spotPlot table'
    renderRasterific "out.png" (dims2D n (n*(h/w))) dia
  where
    f (_, xs) = let (m, v) = meanVarianceUnb $ fst $ V.unzip xs
                in sqrt v / m > 1
