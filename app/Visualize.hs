{-# LANGUAGE OverloadedStrings #-}

module Visualize where

import qualified Data.ByteString.Char8        as B
import qualified Data.Text                    as T
import qualified Data.Vector.Unboxed          as U
import           Data.Double.Conversion.ByteString
import           Statistics.Sample            (meanVarianceUnb)


data Metrics = Metrics
    { averageRank   :: Double
    , variability   :: Double
    , logFoldChange :: Double
    }

getMetrics :: ([T.Text], [B.ByteString], [[Double]]) -> [(B.ByteString, Metrics)]
getMetrics (cts, geneNames, dat) = zip geneNames $ map (f . U.fromList) dat
  where
    f xs = let (m, v) = meanVarianceUnb xs
               fold = log $ max 1 $ U.maximum xs / U.minimum xs
           in Metrics m (sqrt v / m) fold

outputData :: FilePath
           -> ([T.Text], [B.ByteString], [[Double]])
           -> [(B.ByteString, Metrics)] -> IO ()
outputData out (cts, geneNames, dat) metrics = do
    let dat' = map toBS $ filter f $ zip metrics dat
        header = B.pack $ T.unpack $ T.intercalate "\t" $ "Gene" : cts
    B.writeFile out $ B.unlines $ header : dat'
  where
    f ((_, Metrics m v fold), _) = m > 0.0001 && fold > 1.5
    toBS ((nm, _), xs) = B.intercalate "\t" $ nm : map toShortest xs
