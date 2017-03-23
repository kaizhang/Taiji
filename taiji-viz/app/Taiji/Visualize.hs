{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs            #-}

module Taiji.Visualize where

import           Bio.Utils.Functions            (ihs')
import           Bio.Utils.Misc                 (readDouble)
import qualified Data.ByteString.Char8          as B
import           Data.ByteString.Lex.Fractional
import qualified Data.CaseInsensitive           as CI
import           Data.Colour                    (blend)
import qualified Data.HashMap.Strict            as M
import           Data.List.Split                (chunksOf)
import qualified Data.Matrix                    as M
import           Data.Maybe
import qualified Data.Vector                    as V
import           Diagrams.Backend.Rasterific
import           Diagrams.Prelude

import Taiji.Visualize.Data

spotPlot :: Table (Double, Double) -> Diagram B
spotPlot (Table rowlab collab xs) = vsep 1 $ (header:) $ zipWith g rowlab $
    map draw $ scaleData $ M.toLists xs
  where
    scaleData dat = chunksOf n $ zip r $ linearMap (4, 16) e
      where
        (r, e) = unzip $ concat dat
        n = length $ head dat
    draw dat = hsep 1 $ zipWith mkCircle ranks' expr
      where
        ranks' = linearMap (0, 1) ranks
        (ranks, expr) = unzip dat
    header = alignR $ hsep 1 $ map (\x ->
        (alignB $ center $ scale 10 $ texterific x # rotate (90 @@ deg)) <> box) collab
    mkCircle x y = withEnvelope box $ circle y # lw 0 # fc (blend x red white)
    g lab x = alignR $ center (scale 10 $ texterific lab) ||| strutX 5 ||| x
    box = circle 15 # lw 0 :: Diagram B

logFoldChange :: [Double] -> [Double]
logFoldChange xs = map (ihs' . (/ m)) xs
  where
    m = maximum xs

linearMap :: (Double, Double) -> [Double] -> [Double]
linearMap (lo, hi) xs = map f xs
  where
    f x = lo + (x - min') / (max' - min') * (hi - lo)
    min' = minimum xs
    max' = maximum xs
