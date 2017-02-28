{-# LANGUAGE FlexibleInstances #-}

module Types where

import           Bio.Data.Bed          (BED)
import           Data.Binary           (Binary (..))
import qualified Data.ByteString.Char8 as B
import           Data.CaseInsensitive  (CI, mk, original)

type GeneName = CI B.ByteString

instance Binary (CI B.ByteString) where
    put = put . original
    get = fmap mk get

-- | Gene and its regulators
type Linkage = (GeneName, [(GeneName, [BED])])
