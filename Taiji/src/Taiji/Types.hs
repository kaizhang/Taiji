{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE DeriveGeneric #-}

module Taiji.Types where

import GHC.Generics (Generic)
import           Bio.Data.Bed          (BED)
import           Data.Binary           (Binary (..))
import qualified Data.ByteString.Char8 as B
import           Data.CaseInsensitive  (CI, mk, original)
import qualified Data.Map.Strict   as M
import qualified Data.Text             as T

type GeneName = CI B.ByteString

instance Binary (CI B.ByteString) where
    put = put . original
    get = fmap mk get

-- | Gene and its regulators
type Linkage = (GeneName, [(GeneName, [BED])])

data RankTable = RankTable
    { rowNames    :: [B.ByteString]
    , colNames    :: [B.ByteString]
    , ranks       :: [[Double]]
    , expressions :: [[Double]]
    } deriving (Generic)

instance Binary RankTable

data TaijiResults = TaijiResults
    { ranktable :: RankTable
    , nets :: M.Map T.Text (M.Map B.ByteString B.ByteString)
    } deriving (Generic)

instance Binary TaijiResults
