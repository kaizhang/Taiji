{-# LANGUAGE OverloadedStrings #-}
module Constants where

import           Scientific.Workflow (ProcState, getConfig')

atacSeqDir :: ProcState FilePath
atacSeqDir = (++ "/ATAC_Seq/") <$> getConfig' "outputDir"

networkDir :: ProcState FilePath
networkDir = (++ "/Network/") <$> getConfig' "outputDir"

tfbsDir :: ProcState FilePath
tfbsDir = (++ "/TFBS/") <$> getConfig' "outputDir"
