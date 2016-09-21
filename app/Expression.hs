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
               -> IO [(B.ByteString, M.HashMap B.ByteString Double)]  -- ^ cell type and data
readExpression fl = do
    c <- B.readFile fl
    let ((_:header):dat) = map (B.split '\t') $ B.lines c
        rowNames = map head dat
    return $ zipWith (\a b -> (a, M.fromList $ zip rowNames b)) header $
        transpose $ computeZscore $ filtering $ map (map readDouble . tail) dat
  where
    filtering = filter (\x -> not $ all (<1) x || all (==head x) x)
    computeZscore = map (U.toList . scale . U.fromList)
