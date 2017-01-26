{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs            #-}

module Taiji.Visualize where

import           Bio.Utils.Functions            (ihs')
import           Bio.Utils.Misc                 (readDouble)
import qualified Data.ByteString.Char8          as B
import           Data.ByteString.Lex.Fractional
import           Data.Colour                    (blend)
import qualified Data.HashMap.Strict            as M
import           Data.Maybe
import qualified Data.Vector                    as V
import           Diagrams.Backend.Rasterific
import           Diagrams.Prelude

spotPlot :: [String]  -- ^ row labels
         -> [String]  -- ^ column labels
         -> [[(Double, Double)]]  -- ^ size and color
         -> Diagram B
spotPlot rowlab collab = vsep 1 . (header:) . zipWith g rowlab . map draw
  where
    draw dat = hsep 1 $ zipWith mkCircle ranks' expr'
      where
        ranks' = linearMap (2, 15) $ map log ranks
        expr' = linearMap (0, 1) $ map ihs' expr
        (ranks, expr) = unzip dat
    header = alignR $ hsep 1 $ map (\x ->
        (alignB $ center $ scale 6 $ texterific x # rotate (90 @@ deg)) <> box) collab
    mkCircle x y = withEnvelope box $ circle x # lw 0 # fc (blend y red white)
    g lab x = alignR $ center (scale 6 $ texterific lab) ||| strutX 2 ||| x
    box = circle 15 # lw 0 :: Diagram B

readTSV :: B.ByteString -> M.HashMap (B.ByteString, B.ByteString) Double
readTSV input = M.fromList $ concatMap (f . B.split '\t') content
  where
    f (x:xs) = zipWith (\s v -> ((x,s), readDouble v)) samples xs
    (header:content) = B.lines input
    samples = tail $ B.split '\t' header

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
