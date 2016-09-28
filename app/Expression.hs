-- Gene expression
{-# LANGUAGE GADTs             #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}

module Expression where

import           Bio.Utils.Functions   (scale)
import           Bio.Utils.Misc        (readDouble)
import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict   as M
import           Data.List
import qualified Data.Vector.Unboxed   as U

-- | Read RNA expression data
readExpression :: FilePath
               -> IO [( B.ByteString    -- ^ cell type
                     , M.HashMap B.ByteString (Double, Double)  -- ^ absolute value and z-score
                     )]
readExpression fl = do
    c <- B.readFile fl
    let ((_:header):dat) = map (B.split '\t') $ B.lines c
        rowNames = map head dat
        dataTable = filtering $ map (map readDouble . tail) dat
    return $ zipWith (\a b -> (a, M.fromList $ zip rowNames b)) header $
        transpose $ zipWith zip dataTable $ computeZscore dataTable
  where
    filtering = filter (\x -> not $ all (<1) x || all (==head x) x)
    computeZscore = map (U.toList . scale . U.fromList)
