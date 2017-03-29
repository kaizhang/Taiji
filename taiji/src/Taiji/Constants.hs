{-# LANGUAGE OverloadedStrings #-}
module Taiji.Constants where

import           Scientific.Workflow (ProcState, getConfig', getConfigMaybe')

atacOutput :: ProcState FilePath
atacOutput = (++ "/ATAC_Seq/") <$> getConfig' "outputDir"

netOutput :: ProcState FilePath
netOutput = (++ "/Network/") <$> getConfig' "outputDir"

tfbsOutput :: ProcState FilePath
tfbsOutput = (++ "/TFBS/") <$> getConfig' "outputDir"

rnaOutput :: ProcState FilePath
rnaOutput = (++ "/RNA_Seq/") <$> getConfig' "outputDir"

rankOutput :: ProcState FilePath
rankOutput = (++ "/Rank/") <$> getConfig' "outputDir"

bwaIndex :: ProcState FilePath
bwaIndex = (++ "/genome.fa") <$> getConfig' "bwaIndex"

rsemIndex :: ProcState (Maybe FilePath)
rsemIndex = fmap (++ "/genome") <$> getConfigMaybe' "rsemIndex"
