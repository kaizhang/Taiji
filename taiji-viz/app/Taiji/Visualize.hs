{-# LANGUAGE GADTs #-}
{-# LANGUAGE FlexibleContexts #-}

module Taiji.Visualize where

import qualified Data.ByteString.Char8 as B
import Diagrams.Prelude
import Diagrams.Backend.Rasterific
import qualified Data.Vector as V
import Data.ByteString.Lex.Fractional
import Data.Maybe
import Bio.Utils.Functions (ihs')

spotPlot :: [[Double]]  -- ^ size
         -> [[Double]]  -- ^ color
         -> Diagram B
spotPlot ss cs = vcat $ map (hcat . map mkCircle) ss
  where
    mkCircle x = circle x # lw 0 # fc red

readTSV :: B.ByteString -> [[Double]]
readTSV input = map f $ tail $ B.lines input
  where
    f = map (fst . fromJust . readSigned readDecimal) . tail . B.split '\t'

logFoldChange :: [Double] -> [Double]
logFoldChange xs = map (ihs' . (/ m)) xs
  where
    m = maximum xs
