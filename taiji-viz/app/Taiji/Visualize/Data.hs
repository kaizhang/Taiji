{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs            #-}

module Taiji.Visualize.Data where

import           Bio.Utils.Functions
import           Bio.Utils.Misc                 (readDouble)
import qualified Data.ByteString.Char8          as B
import qualified Data.CaseInsensitive           as CI
import           Data.Function                  (on)
import qualified Data.HashMap.Strict            as HM
import           Data.List                      (groupBy, isPrefixOf, sort,
                                                 sortBy)
import qualified Data.Matrix                    as M
import           Data.Maybe
import qualified Data.Vector                    as V
import qualified Data.Vector.Unboxed            as U
import           Statistics.Distribution        (complCumulative)
import           Statistics.Distribution.Normal (normalDistr, normalFromSample)
import qualified Statistics.Function            as S
import           Statistics.Sample              (meanVarianceUnb, mean)

data Table a = Table
    { rowNames :: [String]
    , colNames :: [String]
    , matrix   :: M.Matrix a
    } deriving (Show)

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

data DataFiltOpts = DataFiltOpts
    { coefficientOfVariance :: Double }

defaultDataFiltOpts :: DataFiltOpts
defaultDataFiltOpts = DataFiltOpts
    { coefficientOfVariance = 1 }

-- | Read data, normalize and calculate p-values.
readData :: FilePath   -- ^ PageRank
         -> FilePath   -- ^ Gene expression
         -> DataFiltOpts
         -> IO (Table (Double, Double))  -- ^ ranks, expression and p-values
readData input1 input2 opts = do
    rank <- (fmap ihs' . readTSV) <$> B.readFile input1
    expr <- (fmap ihs' . readTSV) <$> B.readFile input2

    let (labels, xs) = unzip $ map unzip $ groupBy ((==) `on` (fst.fst)) $ sort $
            HM.toList $ HM.intersectionWith (,) rank expr
        rowlab = map (B.unpack . CI.original) $ fst $ unzip $ map head labels
        collab = map (B.unpack . CI.original) $ snd $ unzip $ head $ labels
    return $ filterRows f $ Table rowlab collab $ M.fromLists xs
  where
    f (_, xs) = V.any (>1e-4) xs' && case () of
        _ | n >= 5 -> let (m, v) = meanVarianceUnb $ fst $ V.unzip xs
                      in sqrt v / m > coefficientOfVariance opts
          | otherwise -> V.maximum xs' / V.minimum xs' >= 2.5
      where
        n = V.length xs
        xs' = fst $ V.unzip xs

readTSV :: B.ByteString -> HM.HashMap (CI.CI B.ByteString, CI.CI B.ByteString) Double
readTSV input = HM.fromList $ concatMap (f . B.split '\t') content
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
    (m, v) = meanVarianceUnb $ U.take n $ S.sort xs
    n = truncate $ fromIntegral (U.length xs) * 0.8

pooledPValue :: [Int] -> U.Vector Double -> U.Vector Double
pooledPValue groups xs
    | sum groups /= U.length xs = error "Group assignment error"
    | otherwise = U.fromList $ concat $ zipWith replicate groups $ U.toList $
        U.map (complCumulative (normalFromSample dat)) dat
  where
    dat = U.fromList $ f groups xs
    f (i:rest) x = mean (U.take i x) : f rest (U.drop i x)
    f _ _ = []


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
