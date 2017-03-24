{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE GADTs            #-}

module Taiji.Visualize.Data where

import           Bio.Utils.Misc                 (readDouble)
import qualified Data.ByteString.Char8          as B
import qualified Data.CaseInsensitive           as CI
import qualified Data.HashMap.Strict            as M
import qualified Data.Matrix                    as M
import           Data.Maybe
import           Data.List (sortBy, isPrefixOf)
import qualified Data.Vector                    as V
import qualified Data.Vector.Unboxed            as U
import           Statistics.Distribution        (complCumulative)
import           Statistics.Distribution.Normal (normalDistr)
import           Statistics.Function            (sort)
import           Statistics.Sample              (meanVarianceUnb)

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
pValueGaussian :: U.Vector Double -> U.Vector Double
pValueGaussian xs = U.map (complCumulative distribution) xs
  where
    distribution = normalDistr m $ sqrt v
    (m, v) = meanVarianceUnb $ U.take n $ sort xs
    n = truncate $ fromIntegral (U.length xs) * 0.9

orderByName :: [String] -> ReodrderFn a
orderByName prefix = sortBy $ \(a,_) (b,_) ->
    let idx1 = findIdx a
        idx2 = findIdx b
    in case () of
        _ | isJust idx1 && isJust idx2 -> case compare idx1 idx2 of
                LT -> LT
                GT -> GT
                EQ -> compare a b
          | otherwise -> compare a b
  where
    findIdx x = go prefix 0
      where
        go (y:ys) !i | isPrefixOf y x = Just i
                     | otherwise = go ys (i+1)
        go _ _ = Nothing

filterByName :: [String] -> FilterFn a
filterByName xs = \(x, _) -> x `elem` xs
