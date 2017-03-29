--------------------------------------------------------------------------------
-- Create output directories, check dependencies
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Taiji.Component.Initialization (builder) where

import           Bio.Data.Experiment.Parser
import           Bio.Data.Experiment.Types
import           Bio.Pipeline.NGS
import           Bio.Seq.IO                 (mkIndex)
import           Control.Arrow              (first)
import           Control.Lens
import           Control.Monad.IO.Class     (liftIO)
import           Data.Aeson.Types           (Value (..), parseMaybe)
import           Data.CaseInsensitive       (CI, mk)
import qualified Data.HashMap.Strict        as M
import           Data.List                  (foldl')
import           Data.Maybe                 (fromMaybe, fromJust)
import qualified Data.Text                  as T
import           Data.Yaml                  (decodeFile)
import           Scientific.Workflow        hiding (Success)
import           Shelly                     (fromText, mkdir_p, shelly, test_f)
import           System.IO                  (hPutStrLn, stderr)
import           Text.Printf                (printf)

import           Taiji.Constants

builder :: Builder ()
builder = do
    node "Initialization" [| \() -> do
        dat <- mkAllDirs >> mkIndices >> readData
        case validateData dat of
            Left x  -> error x
            Right _ -> return dat
        |] $ do label .= "Initialization"
                stateful .= True


-- | Create directories
mkAllDirs :: ProcState ()
mkAllDirs = mkdir atacOutput >> mkdir netOutput >> mkdir tfbsOutput >>
    mkdir rnaOutput >> mkdir rankOutput
  where
    mkdir x = x >>= liftIO . shelly . mkdir_p . fromText . T.pack

-- | Generate a variety of genome indices when they are absent.
mkIndices :: ProcState ()
mkIndices = do
    fastq <- getConfig' "genome"

    -- Generate sequence index
    seqIndex <- getConfig "seqIndex"
    fileExist <- liftIO $ shelly $ test_f (fromText seqIndex)
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

    -- Generate RSEM index
    rsemIdx <- rsemIndex
    case rsemIdx of
        Nothing -> return ()
        Just prefix -> do
            anno <- getConfig' "annotation"
            liftIO $ rsemMkIndex prefix anno [fastq]
            return ()

type InputData = ([ATACSeq], [ChIPSeq], [RNASeq], [HiC])

-- | Read input data information.
readData ::  ProcState InputData
readData = do
    inputFl <- getConfig' "input"
    dat <- liftIO $ readYml inputFl
    let parse p key obj = fromMaybe [] $
            M.lookup key obj >>= parseMaybe (parseList p)
    return ( parse parseATACSeq "ATAC-SEQ" dat
           , parse parseChIPSeq "CHIP-SEQ" dat
           , parse parseRNASeq "RNA-SEQ" dat
           , parse parseHiC "HIC" dat
           )
  where
    readYml :: FilePath -> IO (M.HashMap (CI T.Text) Value)
    readYml fl = do
        result <- decodeFile fl
        case result of
            Nothing  -> error "Unable to read input file. Formatting error!"
            Just dat -> return $ M.fromList $ map (first mk) $ M.toList dat

validateData :: InputData -> Either String ()
validateData (atac, chip, rna, hic) = foldl' f (Right ()) counts
  where
    f (Left x) _ = Left x
    f (Right ()) (x,c)
        | c < nExps = Left $ printf
            "Missing Group %s in some experiments. Please check your inputs." (T.unpack x)
        | otherwise = Right ()
    counts = M.toList $ M.fromListWith (+) $ zip (map fromJust $ concat groupNames) $ repeat 1
    nExps = length groupNames
    groupNames = filter (not . null) [ map (^.groupName) atac
        , map (^.groupName) chip, map (^.groupName) rna, map (^.groupName) hic ]
