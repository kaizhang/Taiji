{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses  #-}
{-# LANGUAGE OverloadedLists        #-}
{-# LANGUAGE OverloadedStrings      #-}
{-# LANGUAGE TemplateHaskell        #-}

module Bio.Pipeline.CallPeaks
    ( CallPeakOpts(..)
    , CallPeakOptSetter
    , Cutoff(..)
    , tmpDir
    , cutoff
    , gSize
    , pair
    , defaultCallPeakOpts
    , callPeaks
    , idr
    , idrMultiple
    ) where

import           Bio.Data.Bed
import           Bio.Data.Experiment.Types
import           Conduit
import           Control.Lens
import           Control.Monad.State.Lazy
import qualified Data.ByteString.Char8     as B
import           Data.Conduit.Zlib         (ungzip)
import           Data.Ord
import qualified Data.Text                 as T
import           Shelly                    (fromText, mv, run_, shelly)
import           System.IO
import           System.IO.Temp            (withTempDirectory)

import           Bio.Pipeline.Utils
import           Data.List

type CallPeakOptSetter = State CallPeakOpts ()

data CallPeakOpts = CallPeakOpts
    { callPeakOptsTmpDir :: FilePath
    , callPeakOptsCutoff :: Cutoff
    , callPeakOptsGSize  :: String
    , callPeakOptsPair   :: Bool
    --, callPeakOptsBroad :: Bool
    --, callPeakOptsBroadCutoff :: Double
    }

data Cutoff = PValue Double
            | QValue Double

makeFields ''CallPeakOpts

defaultCallPeakOpts :: CallPeakOpts
defaultCallPeakOpts = CallPeakOpts
    { callPeakOptsTmpDir = "./"
    , callPeakOptsCutoff = QValue 0.01
    , callPeakOptsGSize = "mm"
    , callPeakOptsPair  = False
    --, callPeakOptsBroad = False
    --, callPeakOptsBroadCutoff = 0.05
    }

-- | Call peaks using MACS2.
callPeaks :: FilePath           -- ^ Ouptut file
          -> FileSet            -- ^ Sample
          -> Maybe FileSet      -- ^ Input/control sample
          -> CallPeakOptSetter  -- ^ Options
          -> IO FileSet
callPeaks output target input setter = do
    macs2 output (target^._Single.location) (fmap (^._Single.location) input)
        fileFormat opt
    f <- frip (target^._Single.location) output
    return $ Single $ format .~ NarrowPeakFile $ location .~ output $
        tags .~ ["macs2"] $ info .~ [("FRiP", T.pack $ show f)] $ emptyFile
  where
    opt = execState setter defaultCallPeakOpts
    fileFormat
        | opt^.pair = case target of
            Single fl -> case fl^.format of
                BedFile -> "BEDPE"
                BedGZip -> "BEDPE"
                _ -> error "Only BED input is supported in pairedend mode."
            _ -> error "Paired files detected."
        | otherwise = "AUTO"
{-# INLINE callPeaks #-}

macs2 :: FilePath        -- ^ Output
      -> FilePath        -- ^ Target
      -> Maybe FilePath  -- ^ Input
      -> String          -- ^ File format
      -> CallPeakOpts
      -> IO ()
macs2 output target input fileformat opt = withTempDirectory (opt^.tmpDir)
    "tmp_macs2_dir." $ \tmp -> shelly $ do
        run_ "macs2" $
            [ "callpeak", "-f", T.pack fileformat, "-g", T.pack $ opt^.gSize
            , "--outdir", T.pack tmp, "--tempdir", T.pack tmp, "--keep-dup"
            , "all", "-t", T.pack target
            ] ++ control ++ cut
        mv (fromText $ T.pack $ tmp ++ "/NA_peaks.narrowPeak") $ fromText $
            T.pack output
  where
    control = case input of
        Nothing -> []
        Just x -> ["-c", T.pack x]
    cut = case opt^.cutoff of
        QValue x -> ["--qvalue", T.pack $ show x]
        PValue x -> ["--pvalue", T.pack $ show x]
{-# INLINE macs2 #-}

-- | Fraction of reads in peaks
frip :: FilePath   -- ^ reads, in BedGzip format
     -> FilePath   -- ^ peaks, in bed format
     -> IO Double
frip tags peaks = do
    p <- readBed' peaks :: IO [BED3]
    withFile tags ReadMode $ \h -> do
        (n, m) <- flip execStateT (0::Int, 0::Int) $ 
            sourceHandle h =$= ungzip =$= linesUnboundedAsciiC =$=
            mapC (fromLine :: B.ByteString -> BED3) =$= total =$= intersectBed p $$ count
        return $ fromIntegral m / fromIntegral n
  where
    total = awaitForever $ \i -> do
        (c, x) <- get
        put (c+1,x)
        yield i
    count = awaitForever $ \i -> do
        (x, c) <- get
        put (x, c+1)

idrMultiple :: [FileSet]   -- ^ Peaks
            -> FileSet     -- ^ Merged peaks
            -> Double
            -> FilePath
            -> IO FileSet
idrMultiple [x] _ _ _ = return x
idrMultiple peakFiles merged th output =
    withTempDirectory "./" "tmp_idr_dir." $ \tmp -> do
        peaks <- forM (zip [1..] peakPair) $ \(i, (p1, p2)) -> do
            result <- idr p1 p2 merged th $ tmp ++ "/" ++ show i
            n <- numLine $ result^._Single.location
            return (result^._Single.location, n)
        let final = fst $ maximumBy (comparing snd) peaks
        shelly $ mv (fromText $ T.pack final) $ fromText $ T.pack output
        return $ Single $ format .~ NarrowPeakFile $ location .~ output $
            tags .~ ["IDR"] $ emptyFile
  where
    numLine x = do
        c <- B.readFile x
        return $ length $ B.lines c
    peakPair = comb peakFiles
    comb (x:xs) = zip (repeat x) xs ++ comb xs
    comb _ = []

-- | Perform Irreproducible Discovery Rate (IDR) analysis
idr :: FileSet   -- ^ Peak 1
    -> FileSet   -- ^ Peak 2
    -> FileSet   -- ^ Peaks called from merged replicates (relax threshold)
    -> Double    -- ^ IDR threshold
    -> FilePath  -- ^ Output
    -> IO FileSet
idr peak1 peak2 peakMerged th output = do
    shelly $ run_ "idr" [ "--samples", p1, p2, "--peak-list", pm
        , "--input-file-type", "narrowPeak", "--rank", "signal.value"
        , "--idr-threshold", T.pack $ show th, "-o", T.pack output ]
    return $ Single $ format .~ NarrowPeakFile $ location .~ output $
        tags .~ ["IDR"] $ emptyFile
  where
    p1 = T.pack $ peak1^._Single.location
    p2 = T.pack $ peak2^._Single.location
    pm = T.pack $ peakMerged^._Single.location
