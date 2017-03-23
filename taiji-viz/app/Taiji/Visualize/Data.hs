{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs            #-}

module Taiji.Visualize.Data where

import           Bio.Utils.Functions            (ihs')
import           Bio.Utils.Misc                 (readDouble)
import qualified Data.ByteString.Char8          as B
import           Data.ByteString.Lex.Fractional
import qualified Data.CaseInsensitive           as CI
import qualified Data.HashMap.Strict            as M
import qualified Data.Matrix                    as M
import           Data.Maybe
import qualified Data.Vector                    as V
import Statistics.Function (sort)
import Statistics.Sample (meanVarianceUnb)
import Statistics.Distribution.Normal (normalDistr)
import Statistics.Distribution (complCumulative)

data Table a = Table
    { rowNames :: [String]
    , colNames :: [String]
    , matrix   :: M.Matrix a
    }

type ReodrderFn a = [(String, V.Vector a)] -> [(String, V.Vector a)]
type FilterFn a = (String, V.Vector a) -> Bool

filterRows :: FilterFn a -> Table a -> Table a
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

readTSV :: B.ByteString -> M.HashMap (CI.CI B.ByteString, CI.CI B.ByteString) Double
readTSV input = M.fromList $ concatMap (f . B.split '\t') content
  where
    f (x:xs) = zipWith (\s v -> ((CI.mk x, CI.mk s), readDouble v)) samples xs
    (header:content) = B.lines input
    samples = tail $ B.split '\t' header

-- | Convert a list of values to p-values, assuming a Gaussian distribution.
-- The average and variance is calculated from the lower 90% of the data.
pValueGaussian :: V.Vector Double -> V.Vector Double
pValueGaussian xs = V.map (complCumulative distribution) xs
  where
    distribution = normalDistr m $ sqrt v
    (m, v) = meanVarianceUnb $ V.take n $ sort xs
    n = truncate $ fromIntegral (V.length xs) * 0.9
