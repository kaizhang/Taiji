{-# LANGUAGE DeriveGeneric      #-}
{-# LANGUAGE FlexibleInstances  #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TemplateHaskell    #-}

module Bio.Pipeline.Instances where

import           Bio.Data.Bed
import           Bio.Motif             (Motif (..), PWM (..))
import           Control.Monad
import           Data.Aeson
import           Data.Aeson.TH
import           Data.Binary           (Binary)
import qualified Data.ByteString.Char8 as B
import           Data.CaseInsensitive  (CI, mk, original)
import           Data.Hashable
import qualified Data.Matrix.Unboxed   as MU
import           Data.Maybe
import           Data.Serialize        (Serialize(..))
import qualified Data.Text             as T
import           Data.Vector.Binary    ()
import           Data.Vector.Serialize ()
import qualified Data.Vector.Unboxed   as U
import           GHC.Generics          (Generic)

instance FromJSON B.ByteString where
    parseJSON (String x) = return $ B.pack $ T.unpack x
    parseJSON _          = mzero

instance ToJSON B.ByteString where
    toJSON x = String $ T.pack $ B.unpack x

instance FromJSON (CI B.ByteString) where
    parseJSON = fmap mk . parseJSON

instance ToJSON (CI B.ByteString) where
    toJSON = toJSON . original

instance FromJSON (MU.Matrix Double) where
    parseJSON = genericParseJSON defaultOptions

instance ToJSON (MU.Matrix Double) where
    toEncoding = genericToEncoding defaultOptions

$(deriveJSON defaultOptions ''BED)
$(deriveJSON defaultOptions ''BED3)
$(deriveJSON defaultOptions ''NarrowPeak)
$(deriveJSON defaultOptions ''PWM)
$(deriveJSON defaultOptions ''Motif)

deriving instance Generic BED
instance Serialize BED
instance Binary BED

instance Serialize (MU.Matrix Double)
instance Binary (MU.Matrix Double)

deriving instance Generic PWM
instance Serialize PWM
instance Binary PWM

deriving instance Generic Motif
instance Serialize Motif
instance Binary Motif

instance Serialize (CI B.ByteString) where
    put = put . original
    get = fmap mk get
