--------------------------------------------------------------------------------
-- Create output directories, check dependencies
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Component.Initialization (builder) where

import           Bio.Data.Experiment.Types
import           Bio.Pipeline.NGS
import           Bio.Seq.IO                (mkIndex)
import           Control.Lens
import           Control.Monad.IO.Class    (liftIO)
import           Data.Aeson.Types          (Result (..), Value (..), fromJSON)
import qualified Data.HashMap.Strict       as M
import           Data.Maybe                (isNothing)
import qualified Data.Text                 as T
import           Data.Yaml                 (Object, decodeFile)
import           Scientific.Workflow       hiding (Success)
import           System.IO                 (hPutStrLn, stderr)
import           Turtle                    (fromText, mktree, testfile)

import           Constants

builder :: Builder ()
builder = do
    node "init00" [| \() -> mkAllDirs >> mkIndices >> readData |] $ do
        label .= "Initialization"
        stateful .= True


-- | Create directories
mkAllDirs :: ProcState ()
mkAllDirs = mkdir atacOutput >> mkdir netOutput >> mkdir tfbsOutput >>
    mkdir rnaOutput >> mkdir rankOutput
  where
    mkdir x = x >>= liftIO . mktree . fromText . T.pack

-- | Generate a variety of genome indices when they are absent.
mkIndices :: ProcState ()
mkIndices = do
    fastq <- getConfig' "genome"

    -- generate sequence index
    seqIndex <- getConfig "seqIndex"
    fileExist <- liftIO $ testfile (fromText seqIndex)
    liftIO $ if fileExist
        then hPutStrLn stderr "Sequence index exists. Skipped."
        else do
            hPutStrLn stderr "Generating sequence index"
            mkIndex [fastq] $ T.unpack seqIndex

    -- generate BWA index
    bwaIndex <- getConfig' "bwaIndex"
    liftIO $ bwaMkIndex fastq bwaIndex

    -- generate STAR index
    starIndex <- getConfigMaybe' "starIndex"
    case starIndex of
        Nothing -> return ()
        Just dir -> do
            anno <- getConfig' "annotation"
            liftIO $ starMkIndex "STAR" dir [fastq] anno 100
            return ()

    -- generate RSEM index
    rsemIndex <- getConfigMaybe' "rsemIndex"
    case rsemIndex of
        Nothing -> return ()
        Just prefix -> do
            anno <- getConfig' "annotation"
            liftIO $ rsemMkIndex prefix anno [fastq]
            return ()

-- | Read input data information.
readData ::  ProcState ( [Experiment ATAC_Seq]
                       , [Experiment ChIP_Seq]
                       , [Experiment RNA_Seq]
                       , Maybe FilePath )
readData = do
    inputFl <- getConfig' "input"
    Just dat <- liftIO (decodeFile inputFl  :: IO (Maybe Object))
    return ( parse $ M.lookup "atac-seq" dat
           , parse $ M.lookup "chip-seq" dat
           , parse $ M.lookup "rna-seq" dat
           , case M.lookup "expression_profile" dat of
               Nothing -> Nothing
               Just (String x) -> Just $ T.unpack x
           )
  where
    parse x = case x of
        Nothing -> []
        Just x' -> case fromJSON x' of
            Error msg -> error msg
            Success r -> if any (isNothing . (^.groupName)) r
                then error "Missing \"group\" field in the input file"
                else r
