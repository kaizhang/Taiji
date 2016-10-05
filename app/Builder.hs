{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Builder (graph) where

import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Utils
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Instances    ()
import           Bio.Pipeline.NGS
import           Bio.Pipeline.ScanMotifs
import           Bio.Seq.IO                (mkIndex)
import           Control.Arrow             (first, second, (&&&), (***))
import           Control.Lens
import           Control.Monad.IO.Class    (liftIO)
import           Data.Aeson.Types          (Result (..), fromJSON)
import qualified Data.ByteString.Char8     as B
import qualified Data.HashMap.Strict       as M
import           Data.Maybe
import           Data.Ord
import qualified Data.Text                 as T
import           Data.Yaml                 (Object, decodeFile)
import           Scientific.Workflow       hiding (Success)
import           System.IO                 (hPutStrLn, stderr)
import           Turtle                    (fromText, mktree, testfile)

import           Assign
import           Expression
import           Network
import           Visualize

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


readData ::  ProcState ( [Experiment ATAC_Seq]
                       , [Experiment ChIP_Seq]
                       , [Experiment RNA_Seq] )
readData = do
    inputFl <- getConfig' "input"
    Just dat <- liftIO (decodeFile inputFl  :: IO (Maybe Object))
    return ( parse $ M.lookup "atac-seq" dat
           , parse $ M.lookup "chip-seq" dat
           , parse $ M.lookup "rna-seq" dat
           )
  where
    parse x = case x of
        Nothing -> []
        Just x' -> case fromJSON x' of
            Error msg -> error msg
            Success r -> r

--------------------------------------------------------------------------------
-- Various output dirs
--------------------------------------------------------------------------------
mkAllDirs :: ProcState ()
mkAllDirs = mkdir atacSeqDir >> mkdir networkDir >> mkdir tfbsDir
  where
    mkdir x = x >>= liftIO . mktree . fromText . T.pack

atacSeqDir :: ProcState FilePath
atacSeqDir = (++ "/ATAC_Seq/") <$> getConfig' "outputDir"

networkDir :: ProcState FilePath
networkDir = (++ "/Network/") <$> getConfig' "outputDir"

tfbsDir :: ProcState FilePath
tfbsDir = (++ "/TFBS/") <$> getConfig' "outputDir"


graph :: Builder ()
graph = do
    node "init00" [| \() -> mkAllDirs >> mkIndices >> readData |] $ do
        label .= "Initialization"
        stateful .= True

    node "atac00" [| return . (^._1) |] $ do
        submitToRemote .= Just False
        label .= "Get ATAC-seq data"
    path ["init00", "atac00"]

    node "align00" [| \x -> bwaAlign <$> atacSeqDir <*>
        (getConfig' "bwaIndex") <*> return (bwaCores .= 4) <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True >> remoteParam .= "-pe smp 4"
    node "align01" [| \x -> filterBam <$> atacSeqDir <*> return x
        >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    node "align02" [| \x -> removeDuplicates <$>
        getConfig' "picard" <*> atacSeqDir  <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    node "align03" [| \x -> bamToBed <$> atacSeqDir <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    path ["atac00", "align00", "align01", "align02", "align03"]

    node "peak00" [| mapM $ \x -> return (x, Nothing) |] $
        batch .= 1
    node "peak01" [| \x -> callPeaks <$> atacSeqDir <*>
        return (return ()) <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    node "peak02" [| \xs -> do
        let f x = map (^.location) $ filter ((==NarrowPeakFile) . (^.format)) $ x^.files
        scanMotifs <$> (getConfig' "seqIndex") <*> (getConfig' "motifFile") <*>
            return 1e-5 <*> ((++ "/TFBS_open_chromatin_union.bed") <$> tfbsDir ) <*>
            return (concatMap f xs) >>= liftIO
        |] $ stateful .= True
    path ["align03", "peak00", "peak01", "peak02"]

    node "ass00" [| \x -> getDomains <$> networkDir <*>
        getConfig' "annotation" <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True
    ["peak01"] ~> "ass00"
    node "ass01" [| \(x,y) -> return $ zip x $ repeat y |] $ do
        label .= "prepare input"
        submitToRemote .= Just False
    ["ass00", "peak02"] ~> "ass01"
    node "ass02" [| \xs -> do
        dir <- networkDir
        liftIO $ mapM (linkGeneToTFs dir) xs
        |] $ batch .= 1 >> stateful .= True
    node "ass03" [| \xs -> do
        dir <- networkDir
        liftIO $ mapM (printEdgeList dir) xs
        |] $ batch .= 1 >> stateful .= True
    path ["ass01", "ass02", "ass03"]

    node "rna00" [| return . (^._3) |] $ do
        submitToRemote .= Just False
        label .= "Get RNA-seq data"
    node "rna01" [| \x -> do
        dir <- getConfig' "outputDir"
        starAlign <$> return (dir++"/RNA_Seq/") <*> getConfig' "starIndex" <*>
            return (starCores .= 4) <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True >> remoteParam .= "-l vmem=10G -pe smp 4"
    node "rna02" [| \x -> do
        dir <- getConfig' "outputDir"
        rsemQuant <$> return (dir++"/RNA_Seq/") <*> getConfig' "rsemIndex" <*>
            return (rsemCores .= 4) <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True >> remoteParam .= "-l vmem=10G -pe smp 4"
    node "rna03" [| \x -> do
        dir <- getConfig' "outputDir"
        combineExpression <$> return (dir++"/RNA_Seq/gene_expression.tsv") <*>
            getConfig' "annotation" <*> return x >>= liftIO
        |] $ stateful .= True
    path ["init00", "rna00", "rna01", "rna02", "rna03"]

    node "net00" [| \x -> do
        expression <- getConfigMaybe' "expression_profile"
        liftIO $ case expression of
            Nothing -> pageRank x
            Just e -> personalizedPageRank (e, x)
        |] $ stateful .= True
    path ["ass02", "net00"]

    node "vis00" [| \x -> outputData "ranks" x (getMetrics x) |] $ submitToRemote .= Just False
    ["net00"] ~> "vis00"
