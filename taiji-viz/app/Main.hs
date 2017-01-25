module Main where

import qualified Data.ByteString as B
import Taiji.Visualize
import System.Environment
import qualified Data.Vector.Unboxed as U
import Diagrams.Prelude
import Diagrams.Backend.Rasterific

main :: IO ()
main = do
    [fl] <- getArgs
    c <- B.readFile fl

    let xs' = map logFoldChange $ readTSV c
        w = width dia
        h = height dia
        dia = spotPlot xs' xs'
    renderRasterific "out.png" (dims2D 200 (200*(h/w))) dia
