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
    node "Initialization" [| \() -> mkAllDirs >> mkIndices >> readData |] $ do
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

    -- Generate sequence index
    seqIndex <- getConfig "seqIndex"
    fileExist <- liftIO $ testfile (fromText seqIndex)
    liftIO $ if fileExist
        then hPutStrLn stderr "Sequence index exists. Skipped."
        else do
            hPutStrLn stderr "Generating sequence index"
            mkIndex [fastq] $ T.unpack seqIndex

    -- Generate BWA index
    bwaIndex >>= (liftIO . bwaMkIndex fastq)

    -- Generate STAR index
    starIndex <- getConfigMaybe' "starIndex"
    case starIndex of
        Nothing -> return ()
        Just dir -> do
            anno <- getConfig' "annotation"
            liftIO $ starMkIndex "STAR" dir [fastq] anno 100
            return ()

    -- generate RSEM index
    rsemIdx <- rsemIndex
    case rsemIdx of
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
    dat <- liftIO $ readYml inputFl
    return ( parse $ M.lookup "ATAC-SEQ" dat
           , parse $ M.lookup "CHIP-SEQ" dat
           , parse $ M.lookup "RNA-SEQ" dat
           , case M.lookup "EXPRESSION_PROFILE" dat of
               Nothing -> Nothing
               Just (String x) -> Just $ T.unpack x
           )
  where
    readYml :: FilePath -> IO Object
    readYml fl = do
        result <- decodeFile fl
        case result of
            Nothing -> error "Unable to read input file. Formatting error!"
            Just dat -> return $ M.fromList $
                map (\(k, v) -> (T.toUpper k, v)) $ M.toList dat
    parse x = case x of
        Nothing -> []
        Just x' -> case fromJSON x' of
            Error msg -> error msg
            Success r -> if any (isNothing . (^.groupName)) r
                then error "Missing \"group\" field in the input file"
                else r
