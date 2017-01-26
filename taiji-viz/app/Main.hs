module Main where

import qualified Data.ByteString.Char8 as B
import Taiji.Visualize
import System.Environment
import qualified Data.HashMap.Strict as M
import qualified Data.Vector.Unboxed as U
import Diagrams.Prelude
import Diagrams.Backend.Rasterific
import Data.List
import Data.Function

main :: IO ()
main = do
    [f1,f2] <- getArgs
    c1 <- B.readFile f1
    c2 <- B.readFile f2

    let m1 = readTSV c1
        m2 = readTSV c2
        (labels, xs) = unzip $ map unzip $ groupBy ((==) `on` (fst.fst)) $ sort $
            M.toList $ M.intersectionWith (,) m1 m2
        w = width dia
        h = height dia
        dia = spotPlot (map B.unpack $ fst $ unzip $ map head labels)
            (map B.unpack $ snd $ unzip $ head $ labels) xs
    renderRasterific "out.png" (dims2D 800 (800*(h/w))) dia
