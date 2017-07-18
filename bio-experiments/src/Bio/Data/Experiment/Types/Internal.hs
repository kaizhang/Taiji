{-# LANGUAGE DeriveGeneric          #-}
{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses  #-}
{-# LANGUAGE OverloadedStrings      #-}
{-# LANGUAGE TemplateHaskell        #-}

module Bio.Data.Experiment.Types.Internal where

import           Control.Lens        hiding ((.=))
import           Data.Aeson.TH       (defaultOptions, deriveJSON)
import qualified Data.Map.Strict     as M
import           Data.Serialize      (Serialize (..))
import           Data.Serialize.Text ()
import qualified Data.Text           as T
import           GHC.Generics        (Generic)

data FileType = BamFile
              | BaiFile
              | BedFile
              | BedGZip
              | FastqFile
              | FastqGZip
              | BedgraphFile
              | BigWigFile
              | NarrowPeakFile
              | BroadPeakFile
              | SRA
              | Other
    deriving (Show, Read, Eq, Ord, Generic)

deriveJSON defaultOptions ''FileType
instance Serialize FileType

data File = File
    { fileLocation :: !FilePath
    , fileFormat   :: !FileType
    , fileInfo     :: !(M.Map T.Text T.Text)
    , fileTags     :: ![T.Text]
    } deriving (Show, Read, Eq, Ord, Generic)

makeFields ''File
deriveJSON defaultOptions ''File
instance Serialize File

emptyFile :: File
emptyFile = File
    { fileLocation = ""
    , fileFormat = Other
    , fileInfo = M.empty
    , fileTags = []
    }


data FileSet = Single File
             | Pair File File
             deriving (Show, Read, Eq, Ord, Generic)

makePrisms ''FileSet
deriveJSON defaultOptions ''FileSet
instance Serialize FileSet


data Replicate = Replicate
    { replicateFiles  :: [FileSet]
    , replicateInfo   :: !(M.Map T.Text T.Text)
    , replicateNumber :: !Int
    } deriving (Show, Read, Eq, Ord, Generic)

emptyReplicate :: Replicate
emptyReplicate = Replicate
    { replicateFiles  = []
    , replicateInfo   = M.empty
    , replicateNumber = 0
    }

makeFields ''Replicate
deriveJSON defaultOptions ''Replicate
instance Serialize Replicate


-- Data types representing different kinds of assays.

-- | A set of fields that exist in all kinds of Assays
data CommonFields = CommonFields
    { _commonEid        :: !T.Text
    , _commonGroupName  :: !(Maybe T.Text)
    , _commonCellType   :: !T.Text
    , _commonReplicates :: [Replicate]
    } deriving (Show, Read, Eq, Ord, Generic)

makeLenses ''CommonFields
deriveJSON defaultOptions ''CommonFields
instance Serialize CommonFields

defaultCommonFields :: CommonFields
defaultCommonFields = CommonFields
    { _commonEid = ""
    , _commonGroupName = Nothing
    , _commonCellType = ""
    , _commonReplicates = []
    }


class Experiment e where
    commonFields :: Lens' e CommonFields

    eid :: Lens' e T.Text
    eid = commonFields . commonEid

    groupName :: Lens' e (Maybe T.Text)
    groupName = commonFields . commonGroupName

    cellType :: Lens' e T.Text
    cellType = commonFields . commonCellType

    replicates :: Lens' e [Replicate]
    replicates = commonFields . commonReplicates

    {-# MINIMAL commonFields #-}

-- | Next generation sequencing experiments.
class NGS e where
    pairedEnd :: e -> Bool

data ChIPSeq = ChIPSeq
    { chipseqCommon    :: CommonFields
    , chipseqTarget    :: !T.Text
    , chipseqPairedEnd :: !Bool
    , chipseqControl   :: Maybe ChIPSeq
    } deriving (Show, Read, Eq, Ord, Generic)

target :: Lens' ChIPSeq T.Text
target f e = (\x -> e{chipseqTarget = x}) <$> f (chipseqTarget e)

control :: Lens' ChIPSeq (Maybe ChIPSeq)
control f e = (\x -> e{chipseqControl = x}) <$> f (chipseqControl e)

deriveJSON defaultOptions ''ChIPSeq

instance Experiment ChIPSeq where
    commonFields f e = (\x -> e{chipseqCommon = x}) <$> f (chipseqCommon e)

instance NGS ChIPSeq where
    pairedEnd = chipseqPairedEnd

instance Serialize ChIPSeq

defaultChIPSeq :: ChIPSeq
defaultChIPSeq = ChIPSeq
    { chipseqCommon = defaultCommonFields
    , chipseqTarget = ""
    , chipseqPairedEnd = False
    , chipseqControl = Nothing
    }


data ATACSeq = ATACSeq
    { atacseqCommon    :: CommonFields
    , atacseqPairedEnd :: !Bool
    } deriving (Show, Read, Eq, Ord, Generic)

deriveJSON defaultOptions ''ATACSeq

instance Experiment ATACSeq where
    commonFields f e = (\x -> e{atacseqCommon = x}) <$> f (atacseqCommon e)

instance NGS ATACSeq where
    pairedEnd = atacseqPairedEnd

instance Serialize ATACSeq


data RNASeq = RNASeq
    { rnaseqCommon    :: CommonFields
    , rnaseqPairedEnd :: !Bool
    } deriving (Show, Read, Eq, Ord, Generic)

deriveJSON defaultOptions ''RNASeq

instance Experiment RNASeq where
    commonFields f e = (\x -> e{rnaseqCommon = x}) <$> f (rnaseqCommon e)

instance NGS RNASeq where
    pairedEnd = rnaseqPairedEnd

instance Serialize RNASeq


data HiC = HiC
    { hicCommon :: CommonFields
    } deriving (Show, Read, Eq, Ord, Generic)

deriveJSON defaultOptions ''HiC

instance Experiment HiC where
    commonFields f e = (\x -> e{hicCommon = x}) <$> f (hicCommon e)

instance NGS HiC where
    pairedEnd _ = True

instance Serialize HiC


class IsDNASeq a
instance IsDNASeq ChIPSeq
instance IsDNASeq ATACSeq
instance IsDNASeq HiC
