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
import qualified Data.Matrix as M

data Table a = Table
    { rowNames :: [String]
    , colNames :: [String]
    , matrix :: M.Matrix a
    }

type ReodrderFn a = [(String, V.Vector a)] -> [(String, V.Vector a)]

filterRows :: ((String, V.Vector a) -> Bool) -> Table a -> Table a
filterRows fn table = table
    { rowNames = names
    , matrix = M.fromRows rows }
  where
    (names, rows) = unzip $ filter fn $ zip (rowNames table) $ M.toRows $ matrix table

reorderRows :: ReodrderFn a -> Table a -> Table a
reorderRows fn table = table
    { rowNames = names
    , matrix = M.fromRows rows }
  where
    (names, rows) = unzip $ fn $ zip (rowNames table) $ M.toRows $ matrix table

reorderColumns :: ReodrderFn a -> Table a -> Table a
reorderColumns fn table = table
    { colNames = names
    , matrix = M.fromColumns cols}
  where
    (names, cols) = unzip $ fn $ zip (colNames table) $ M.toColumns $ matrix table

spotPlot :: Table (Double, Double) -> Diagram B
spotPlot (Table rowlab collab xs) = vsep 1 $ (header:) $ zipWith g rowlab $
    map draw $ M.toLists xs
  where
    draw dat = hsep 1 $ zipWith mkCircle ranks' expr'
      where
        ranks' = linearMap (0, 1) ranks
        expr' = linearMap (2, 15) expr
        (ranks, expr) = unzip dat
    header = alignR $ hsep 1 $ map (\x ->
        (alignB $ center $ scale 6 $ texterific x # rotate (90 @@ deg)) <> box) collab
    mkCircle x y = withEnvelope box $ circle y # lw 0 # fc (blend x red white)
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
