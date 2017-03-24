{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs            #-}

module Taiji.Visualize where

import           Bio.Utils.Functions            (ihs')
import           Bio.Utils.Misc                 (readDouble)
import qualified Data.ByteString.Char8          as B
import qualified Data.Vector.Generic         as G
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

spotPlot :: Double -> Table (Double, Double) -> Diagram B
spotPlot cutoff (Table rowlab collab xs) = vsep 1 $ (header:) $ zipWith g rowlab $
    map draw $ M.toRows $ M.zip3 ranks expr' pvalues
  where
    draw dat = hsep 1 $ V.toList $ V.zipWith3 mkCircle r' e p
      where
        r' = linearMap (0, 1) r
        (r, e, p) = V.unzip3 dat
    header = alignR $ hsep 1 $ map (\x ->
        (alignB $ center $ scale 10 $ texterific x # rotate (90 @@ deg)) <> box) collab
    mkCircle x y p
        | p <= cutoff = withEnvelope box $ circle y # lw 0 # fc (blend x red white) # showOrigin
        | otherwise = withEnvelope box $ circle y # lw 0 # fc (blend x red white)
    g lab x = alignR $ center (scale 10 $ texterific lab) ||| strutX 5 ||| x
    box = circle 15 # lw 0 :: Diagram B
    expr' = M.fromVector (M.dim expr) $ linearMap (4, 16) $ M.flatten expr
    pvalues = M.fromRows $ map (G.convert . pValueGaussian . G.convert) $ M.toRows ranks
    (ranks, expr) = M.unzip xs

logFoldChange :: [Double] -> [Double]
logFoldChange xs = map (ihs' . (/ m)) xs
  where
    m = maximum xs

linearMap :: (Double, Double) -> V.Vector Double -> V.Vector Double
linearMap (lo, hi) xs = V.map f xs
  where
    f x = lo + (x - min') / (max' - min') * (hi - lo)
    min' = V.minimum xs
    max' = V.maximum xs
