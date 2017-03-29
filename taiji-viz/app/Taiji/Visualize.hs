{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs            #-}

module Taiji.Visualize where

import           Bio.Utils.Functions            (ihs', scale)
import           Bio.Utils.Misc                 (readDouble)
import           Control.Arrow                  (first)
import qualified Data.ByteString.Char8          as B
import           Data.ByteString.Lex.Fractional
import qualified Data.CaseInsensitive           as CI
import           Data.Colour                    (blend)
import qualified Data.HashMap.Strict            as M
import           Data.List.Split                (chunksOf)
import qualified Data.Matrix                    as M
import           Data.Maybe
import qualified Data.Vector                    as V
import qualified Data.Vector.Generic            as G
import           Diagrams.Backend.Cairo         (B)
import           Diagrams.Prelude               hiding (scale)
import           Diagrams.TwoD.Text
import           Graphics.SVGFonts              (textSVG)

import           Taiji.Visualize.Data

heatmap :: Table (Double, Double) -> Diagram B
heatmap (Table rowlab collab xs) = vcat $ (header:) $ map (hcat . map mkRect) ranks'
  where
    mkRect x = rect 15 2 # lw 0 # fc (blend x red white)
    header = hcat $ map (\x ->
        (alignB $ textBounded x # rotate (90 @@ deg)) <> box) collab
    ranks' = M.toLists (M.fromVector (M.dim xs) $ linearMap (0, 1) $
        M.flatten $ fst $ M.unzip xs :: M.Matrix Double)
    textBounded x = stroke (textSVG x 15) # lw 0 # fc black
    box = rect 15 3 # lw 0 :: Diagram B

spotPlot :: Double -> Table ((Double, Double), Double) -> Diagram B
spotPlot cutoff (Table rowlab collab xs)
    | null rowlab = error "Nothing to plot"
    | otherwise = vsep 1 $ (header:) $ zipWith g rowlab $
        map (hsep 1 . map mkCircle) $ M.toLists $ M.zip3 ranks' expr' pvalues
  where
    header = alignR $ hsep 1 $ map (\x ->
        (alignB $ textBounded x # rotate (90 @@ deg)) <> box) collab
    mkCircle (x, y, p)
        | p <= cutoff = withEnvelope box $ circle y # lw 0 # fc (blend x red white) # showOrigin
        | otherwise = withEnvelope box $ circle y # lw 0 # fc (blend x red white)
    g lab x = alignR $ textBounded lab ||| strutX 5 ||| x
    box = circle 15 # lw 0 :: Diagram B
    expr' = M.fromVector (M.dim expr) $ linearMap (4, 16) $ M.flatten expr
    ranks' = M.fromVector (M.dim ranks) $ linearMap (0, 1) $ M.flatten ranks
    ((ranks, expr), pvalues) = first M.unzip $ M.unzip xs
    textBounded x = stroke (textSVG x 15) # lw 0 # fc black

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
