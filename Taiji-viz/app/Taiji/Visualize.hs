{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs            #-}

module Taiji.Visualize where

import           Bio.Utils.Functions    (ihs')
import           Control.Arrow          (first)
import           Data.Colour            (blend)
import qualified Data.Matrix            as M
import qualified Data.Vector            as V
import           Diagrams.Backend.Cairo (B)
import           Diagrams.Prelude
import           Graphics.SVGFonts      (textSVG)
import Text.Printf

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
    | otherwise = center plot === strutY 30 === center legend
  where
    plot = vsep 1 $ (header:) $ zipWith g rowlab $
        map (hsep 1 . map mkCircle) $ M.toLists $ M.zip3 ranks' expr' pvalues
    header = alignR $ hsep 1 $ map (\x ->
        (alignB $ textBounded x # rotate (90 @@ deg)) <> box) collab
    mkCircle (x, y, p) = withEnvelope box $ circle y # lw 0 # fc (blend x red white)
    g lab x = alignR $ textBounded lab ||| strutX 5 ||| x
    box = circle 15 # lw 0 :: Diagram B
    expr' = M.fromVector (M.dim expr) $ linearMap (4, 16) $ M.flatten expr
    ranks' = M.fromVector (M.dim ranks) $ linearMap (0, 1) $ M.flatten ranks
    ((ranks, expr), pvalues) = first M.unzip $ M.unzip xs
    textBounded x = stroke (textSVG x 15) # lw 0 # fc black

    legend = rank_legend ||| strutX 50 ||| expr_legend
      where
        expr_legend = center (hsep 1 $ map (\(x,s) -> withEnvelope box
            (circle s # lw 0 # fc black) === strutY 2 === textBounded (printf "%.2f" x)) $
            zip [min_expr, step_expr .. max_expr] $ V.toList $ linearMap (4, 16) $
            V.fromList [min_expr, step_expr .. max_expr] ) === strutY 2 ===
            textBounded "log expression level"
        rank_legend = vsep 2 $
            [ rect 100 25 # lw 0 # fillTexture gradient
            , center $ position $ zip [0^&0, 50^&0, 100^&0] $ map (textBounded . printf "%.2f")
                [min_rank, min_rank + step_rank, max_rank]
            , textBounded "normalized rank score" ]
        gradient = mkLinearGradient stops ((-50) ^& 0) (50 ^& 0) GradPad
        stops = mkStops [(white, 0, 1), (red, 1, 1)]
        min_expr = V.minimum $ M.flatten expr
        max_expr = V.maximum $ M.flatten expr
        step_expr = (max_expr - min_expr) / 3
        min_rank = V.minimum $ M.flatten ranks
        max_rank = V.maximum $ M.flatten ranks
        step_rank = (max_rank - min_rank) / 2

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
