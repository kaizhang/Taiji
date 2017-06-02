{-# LANGUAGE ExtendedDefaultRules #-}
{-# LANGUAGE FlexibleContexts     #-}
{-# LANGUAGE GADTs                #-}
{-# LANGUAGE OverloadedLists      #-}
{-# LANGUAGE OverloadedStrings    #-}
{-# LANGUAGE TemplateHaskell      #-}

module Bio.Pipeline.NGS
    ( sra2fastq
    , BWAOpts
    , BWAOptSetter
    , bwaCores
    , bwaSeedLen
    , bwaMaxMis
    , bwaReadTrim
    , defaultBWAOpts
    , bwaMkIndex
    , bwaAlign
    , filterBam
    , removeDuplicates
    , bam2Bed
    , sortedBam2BedPE
    , mergeReplicatesBed

    , STAROpts
    , STAROptSetter
    , starCmd
    , starCores
    , starSort
    , starTmpDir
    , starMkIndex
    , starAlign

    , RSEMOpts
    , RSEMOptSetter
    , rsemPath
    , rsemCores
    , rsemSeed
    , rsemMkIndex
    , rsemQuant
    ) where

import           Bio.Data.Bam              (bamToBed, readBam, runBam,
                                            sortedBamToBedPE)
import           Bio.Data.Bed              (BED3 (..), BEDLike (..), toLine)
import           Bio.Data.Experiment.Types
import           Conduit
import           Control.Lens
import           Control.Monad             (forM)
import           Control.Monad.State.Lazy
import           Data.Conduit.Zlib         (gzip, ungzip)
import           Data.Conduit.Zlib         (gzip)
import           Data.Int                  (Int32)
import           Data.List.Split           (splitOn)
import           Data.Maybe                (fromJust)
import qualified Data.Text                 as T
import           Shelly                    (cmd, cp, escaping, fromText,
                                            mkdir_p, mv, run_, shelly, silently,
                                            test_d, test_f)
import           System.FilePath
import           System.IO                 (hPutStrLn, stderr)
import           System.IO.Temp            (withTempDirectory)
import           Text.Printf               (printf)

import           Bio.Pipeline.Utils        (mapOfFiles)
default (T.Text)

-- | Convert SRA to gzipped Fastq file
sra2fastq :: (NGS e, Experiment e)
          => FilePath
          -> e
          -> IO e
sra2fastq outDir = mapOfFiles fn
  where
    fn e r (Single fl) = if fl^.format == SRA
        then if pairedEnd e then dumpPair fl else dumpSingle fl
        else return [Single fl]
    fn _ _ x = return [x]
    dumpSingle fl = do
        let f1_name = outDir ++ "/" ++ fl^.location ++ ".fastq.gz"
            f1 = format .~ FastqGZip $ location .~ f1_name $ fl
        case splitOn "+" (fl^.location) of
            [f] -> shelly $ cmd "fastq-dump" "--origfmt" "-I" "--gzip"
                "-O" (T.pack outDir) (T.pack f)
            fs -> do
                forM_ fs $ \f -> shelly $ cmd "fastq-dump" "--origfmt" "-I"
                    "-O" (T.pack outDir) (T.pack f)
                let f1s = map (\x -> T.pack $ outDir ++ "/" ++ x ++ ".fastq") fs
                shelly $ escaping False $ do
                    run_ "cat" $ f1s ++ ["|", "gzip", "-c", ">", T.pack f1_name]
                    run_ "rm" f1s
        return [Single f1]
    dumpPair fl = do
        let f1_name = outDir ++ "/" ++ fl^.location ++ "_1.fastq.gz"
            f2_name = outDir ++ "/" ++ fl^.location ++ "_2.fastq.gz"
            f1 = format .~ FastqGZip $ location .~ f1_name $ fl
            f2 = format .~ FastqGZip $ location .~ f2_name $ fl

        case splitOn "+" (fl^.location) of
            [f] -> shelly $ cmd "fastq-dump" "--origfmt" "-I" "--split-files" "--gzip"
                "-O" (T.pack outDir) (T.pack f)
            fs -> do
                forM_ fs $ \f -> shelly $ cmd "fastq-dump" "--origfmt" "-I"
                    "--split-files" "-O" (T.pack outDir) (T.pack f)
                let f1s = map (\x -> T.pack $ outDir ++ "/" ++ x ++ "_1.fastq") fs
                    f2s = map (\x -> T.pack $ outDir ++ "/" ++ x ++ "_2.fastq") fs
                shelly $ escaping False $ do
                    run_ "cat" $ f1s ++ ["|", "gzip", "-c", ">", T.pack f1_name]
                    run_ "rm" f1s
                    run_ "cat" $ f2s ++ ["|", "gzip", "-c", ">", T.pack f2_name]
                    run_ "rm" f2s
        return [Pair f1 f2]

--------------------------------------------------------------------------------
-- DNA-seq
--------------------------------------------------------------------------------

data BWAOpts = BWAOpts
    { _bwaCores    :: Int          -- ^ number of cpu cores
    , _bwaSeedLen  :: Int     -- ^ seed length, equivalent to -l
    , _bwaMaxMis   :: Int    -- ^ max mismatches in seed, equivalent to -k
    , _bwaReadTrim :: Int       -- ^ dynamic read trimming, equivalent to -q
    , _bwaTmpDir   :: FilePath  -- ^ temp dir
    } deriving (Show)

makeLenses ''BWAOpts

defaultBWAOpts :: BWAOpts
defaultBWAOpts = BWAOpts
    { _bwaCores = 1
    , _bwaSeedLen = 32
    , _bwaMaxMis = 2
    , _bwaReadTrim = 5
    , _bwaTmpDir = "./"
    }

type BWAOptSetter = State BWAOpts ()

{-
-- Determine whether bwa index has been generated
isBWAIndexExist :: FilePath -> IO Bool
isBWAIndexExist dir = do
    fileExist <- testfile (fromText $ T.pack dir ++ "/")
    ls
  where
    filename = "the_suffix_of_this_file_is_the_prefix_of_bwa_indices."
    -}

-- | Generate BWA genome index
bwaMkIndex :: FilePath
           -> FilePath   -- ^ Index prefix, e.g., /path/genome.fa
           -> IO FilePath
bwaMkIndex input prefix = do
    fileExist <- shelly $ test_f (fromText $ T.pack prefix)
    if fileExist
        then hPutStrLn stderr "BWA index exists. Skipped."
        else shelly $ do
            mkdir_p $ fromText $ T.pack $ takeDirectory prefix
            cp (fromText $ T.pack input) $ fromText $ T.pack prefix
            liftIO $ hPutStrLn stderr "Generating BWA index"
            cmd "bwa" "index" "-p" (T.pack prefix) "-a" "bwtsw" $ T.pack input
    return prefix

-- | Tag alignment with BWA aligner.
bwaAlign :: (NGS e, IsDNASeq e, Experiment e)
         => FilePath  -- ^ Directory to save the results
         -> FilePath  -- ^ Genome index
         -> BWAOptSetter
         -> e
         -> IO e
bwaAlign dir index setter = mapOfFiles fn
  where
    fn e r (Single fl) = if fl^.format == FastqFile || fl^.format == FastqGZip
        then do
            shelly $ mkdir_p $ fromText $ T.pack dir
            let output = printf "%s/%s_rep%d.bam" dir (T.unpack $ e^.eid) $ r^.number
                input = T.pack $ fl^.location

            stats <- withTempDirectory (opt^.bwaTmpDir) "bwa_align_tmp_dir." $
                \tmpdir -> shelly $ escaping False $ do
                    let tmp_sai = T.pack $ tmpdir ++ "/tmp.sai"
                    -- Align reads and save the results to tmp_sai.
                    cmd "bwa" "aln" "-q" (T.pack $ show $ opt^.bwaReadTrim) "-l"
                        (T.pack $ show $ opt^.bwaSeedLen) "-k"
                        (T.pack $ show $ opt^.bwaMaxMis) "-t"
                        (T.pack $ show $ opt^.bwaCores) (T.pack index) input ">" tmp_sai
                    -- Convert sai to bam.
                    cmd "bwa" "samse"  (T.pack index) tmp_sai input "|"
                        "samtools" "view" "-Su" "-" ">" (T.pack output)

            return [ Single $ format .~ BamFile $
                              location .~ output $ fl
                   ]
        else return []
    fn e r (Pair f1 f2) =
        if (f1^.format == FastqFile || f1^.format == FastqGZip) &&
           (f2^.format == FastqFile || f2^.format == FastqGZip)
           then do
                shelly $ mkdir_p $ fromText $ T.pack dir
                let output = printf "%s/%s_rep%d.bam" dir (T.unpack $ e^.eid) $
                        r^.number
                    input1 = T.pack $ f1^.location
                    input2 = T.pack $ f2^.location

                stats <- shelly $ escaping False $
                    cmd "bwa" "mem" "-M" "-k" (T.pack $ show $ opt^.bwaSeedLen)
                        "-t" (T.pack $ show $ opt^.bwaCores) (T.pack index)
                        input1 input2 "|" "samtools" "view" "-Su" "-" ">"
                        (T.pack $ output)
                return [ Single $ format .~ BamFile $
                                  location .~ output $ f1
                       ]
            else return []
    opt = execState setter defaultBWAOpts

-- | Remove low quality and redundant tags, fill in mate information.
filterBam :: (NGS e, IsDNASeq e, Experiment e)
          => FilePath  -- ^ directory to save the results
          -> e -> IO e
filterBam dir = mapOfFiles fn
  where
    fn e r (Single fl)
        | fl^.format == BamFile = do
            shelly $ mkdir_p $ fromText $ T.pack dir
            let output = T.pack $ printf "%s/%s_rep%d_filt.bam" dir
                    (T.unpack $ e^.eid) (r^.number)
                input = T.pack $ fl^.location
            bamFilter (pairedEnd e) input output
            return [ Single $ info .~ [] $ format .~ BamFile $
                location .~ T.unpack output $ fl ]
        | otherwise = return []
    fn _ _ _ = return []
    bamFilter isPair input output = withTempDirectory dir "tmp_filt_dir." $ \tmp ->
        shelly $ escaping False $ silently $ do
            let tmp_filt = T.pack $ tmp ++ "/tmp_filt.bam"
                tmp_fixmate = T.pack $ tmp ++ "/tmp_fixmate.bam"
                tmp_sort = T.pack $ tmp ++ "/tmp_sort"
            run_ "samtools" $ ["view"] ++
                (if isPair then ["-f", "2"] else []) ++
                ["-F", "0x70c", "-q", "30", "-u", input] ++
                ( if isPair
                    then [ "|", "samtools", "sort", "-", "-n", "-T", tmp_sort
                        , "-l", "0", "-o", tmp_filt ]
                    else [ "|", "samtools", "sort", "-", "-T", tmp_sort
                        , "-l", "9", "-o", output ] )
            when isPair $ do
                run_ "samtools" ["fixmate", "-r", tmp_filt, tmp_fixmate]
                run_ "samtools" [ "view", "-F", "1804", "-f", "2", "-u"
                    , tmp_fixmate, "|", "samtools", "sort", "-", "-T"
                    , tmp_sort, "-l", "9", "-o", output ]

-- | Remove duplicates
removeDuplicates :: (NGS e, IsDNASeq e, Experiment e)
                 => FilePath -> FilePath -> e -> IO e
removeDuplicates picardPath dir = mapOfFiles fn
  where
    fn e r (Single fl) = if fl^.format == BamFile
        then do
            shelly $ mkdir_p $ fromText $ T.pack dir
            let output = printf "%s/%s_rep%d_filt_mono.bam" dir
                    (T.unpack $ e^.eid) (r^.number)
                input = fl^.location
                qcFile = printf ("%s/%s_rep%d_picard.qc") dir
                    (T.unpack $ e^.eid) (r^.number)

            withTempDirectory "./" "tmp_picard_dir." $ \tmp -> shelly $ do
                let markdupTmp = tmp++"/dup_marked.bam"
                    filtTmp = tmp++"/dup_filt.bam"
                -- Mark duplicates
                run_ "java" ["-Xmx4G", "-jar", T.pack picardPath
                    , "MarkDuplicates", T.pack $ "INPUT=" ++ input
                    , T.pack $ "OUTPUT=" ++ markdupTmp
                    , T.pack $ "TMP_DIR=" ++ tmp
                    , T.pack $ "METRICS_FILE=" ++ qcFile
                    , "VALIDATION_STRINGENCY=LENIENT"
                    , "ASSUME_SORT_ORDER=coordinate", "REMOVE_DUPLICATES=false"]

                -- Remove duplicates.
                escaping False $ run_ "samtools" [ "view", "-F", "0x70c", "-b"
                    , T.pack markdupTmp, ">", T.pack filtTmp ]

                -- Re-sort by names for pairedend sequencing
                if pairedEnd e
                    then run_ "samtools" [ "sort", T.pack filtTmp, "-n", "-T"
                        , T.pack $ tmp ++ "/tmp_sort", "-o", T.pack output ]
                    else mv (fromText $ T.pack filtTmp) $ fromText $ T.pack output

                let finalBam = format .~ BamFile $
                               tags .~ ["processed bam file"] $
                               location .~ output $ fl
                    dupQC = format .~ Other $
                            info .~ [] $
                            tags .~ ["picard qc file"] $
                            location .~ qcFile $ fl

                return [Single finalBam, Single dupQC]
        else return []
    fn _ _ _ = return []

bam2Bed :: String -> FileSet -> IO (Maybe FileSet)
bam2Bed prefix (Single fl)
    | fl^.format == BamFile = do
        shelly $ mkdir_p $ fromText $ T.pack $ takeDirectory prefix
        let output = prefix ++ takeBaseName (fl^.location) ++ ".bed.gz"
            bedFile = format .~ BedGZip $
                      location .~ output $ fl
        runBam $ readBam (fl^.location) =$= bamToBed =$= mapC toLine =$=
            unlinesAsciiC =$= gzip $$ sinkFileBS output
        return $ Just $ Single bedFile
    | otherwise = return Nothing
bam2Bed _ _ = return Nothing
{-# INLINE bam2Bed #-}

-- | Convert name sorted BAM to BEDPE suitable for MACS2.
sortedBam2BedPE :: String -> FileSet -> IO (Maybe FileSet)
sortedBam2BedPE prefix (Single fl)
    | fl^.format == BamFile = do
        shelly $ mkdir_p $ fromText $ T.pack $ takeDirectory prefix
        let output = prefix ++ takeBaseName (fl^.location) ++ ".bed.gz"
            bedFile = format .~ BedGZip $
                      location .~ output $ fl
        runBam $ readBam (fl^.location) =$= sortedBamToBedPE =$=
            concatMapC f =$= mapC toLine =$= unlinesAsciiC =$= gzip $$
            sinkFileBS output
        return $ Just $ Single bedFile
    | otherwise = return Nothing
  where
    f (b1, b2)
        | chrom b1 /= chrom b2 || bedStrand b1 == bedStrand b2 = Nothing
        | otherwise =
            let left = if fromJust (bedStrand b1)
                    then chromStart b1 else chromStart b2
                right = if not (fromJust $ bedStrand b2)
                    then chromEnd b2 else chromEnd b1
            in if left <= right
                  then Just $ BED3 (chrom b1) left right
                  else error "Left coordinate is smaller than right coordinate."
sortedBam2BedPE _ _ = return Nothing
{-# INLINE sortedBam2BedPE #-}

-- | Merge multiple BED files.
mergeReplicatesBed :: Experiment e => FilePath -> e -> IO e
mergeReplicatesBed dir e = do
    shelly $ mkdir_p $ fromText $ T.pack dir
    let fls = e^..replicates.folded.files.folded._Single
        output = printf "%s/%s_rep0.bed.gz" dir (T.unpack $ e^.eid)
        bedFile = Single $ format .~ BedGZip $ location .~ output $ emptyFile
        source = forM_ fls $ \fl -> case fl^.format of
            BedGZip -> sourceFileBS (fl^.location) =$= ungzip
            BedFile -> sourceFileBS (fl^.location)
            _       -> return ()
    runResourceT $ source =$= gzip $$ sinkFile output
    return $ replicates .~ [files .~ [bedFile] $ emptyReplicate] $ e
{-# INLINE mergeReplicatesBed #-}


--------------------------------------------------------------------------------
-- RNA-seq
--------------------------------------------------------------------------------

data STAROpts = STAROpts
    { _starCmd    :: FilePath
    , _starCores  :: Int
    , _starTmpDir :: FilePath
    , _starSort   :: Bool
    }

makeLenses ''STAROpts

type STAROptSetter = State STAROpts ()

defaultSTAROpts :: STAROpts
defaultSTAROpts = STAROpts
    { _starCmd = "STAR"
    , _starCores = 1
    , _starTmpDir = "./"
    , _starSort = False
    }


-- | Create index files for STAR
starMkIndex :: FilePath   -- ^ STAR command path
            -> FilePath   -- ^ Directory used to store genome indices
            -> [FilePath] -- ^ Fastq files
            -> FilePath   -- ^ Annotation file
            -> Int        -- ^ The length of the genomic sequence
                          -- around the annotated junction to be used in
                          -- constructing the splice junctions database. Set it
                          -- to "ReadLength-1" or 100 for general purpose.
            -> IO FilePath
starMkIndex star dir fstqs anno r = do
    dirExist <- shelly $ test_d $ fromText $ T.pack dir
    if dirExist
        then hPutStrLn stderr "STAR index directory exists. Skipped."
        else shelly $ do
            mkdir_p $ fromText $ T.pack dir
            liftIO $ hPutStrLn stderr "Generating STAR indices"
            run_ (fromText $ T.pack star) $
                [ "--runThreadN", "1", "--runMode", "genomeGenerate", "--genomeDir"
                , T.pack dir, "--genomeFastaFiles" ] ++ map T.pack fstqs ++
                ["--sjdbGTFfile", T.pack anno, "--sjdbOverhang", T.pack $ show r]
    return dir

-- | Align RNA-seq raw reads with STAR
starAlign :: FilePath                    -- ^ Output directory
          -> FilePath                    -- ^ STAR genome index
          -> STAROptSetter               -- ^ Options
          -> RNASeq
          -> IO RNASeq
starAlign dir index setter = mapOfFiles fn
  where
    isFastq :: File -> Bool
    isFastq x = x^.format == FastqFile || x^.format == FastqGZip
    opt = execState setter defaultSTAROpts
    fn e r flset = if null fls then return [] else shelly $ do
        mkdir_p $ fromText $ T.pack dir
        withTempDirectory (opt^.starTmpDir) "STAR_align_tmp_dir." $ \tmp_dir -> do
            let outputGenome = printf "%s/%s_rep%d_genome.bam" dir
                    (T.unpack $ e^.eid) $ r^.number
                outputAnno = printf "%s/%s_rep%d_anno.bam" dir
                    (T.unpack $ e^.eid) $ r^.number
                inputs = map (T.pack . (^.location)) fls

            run_ (fromText $ T.pack $ opt^.starCmd) $
                ["--genomeDir", T.pack index, "--readFilesIn"] ++ inputs ++
                [ "--outFileNamePrefix", T.pack $ tmp_dir ++ "/"
                , "--runThreadN",  T.pack $ show $ opt^.starCores ] ++
                ( if (head fls)^.format == FastqGZip
                    then ["--readFilesCommand", "zcat"]
                    else [] ) ++
                [ "--genomeLoad", "NoSharedMemory"
                , "--outFilterType", "BySJout"     -- reduces the number of ”spurious” junctions
                , "--outFilterMultimapNmax", "20"  -- max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
                , "--alignSJoverhangMin", "8"      -- minimum overhang for unannotated junctions
                , "--alignSJDBoverhangMin", "1"    -- minimum overhang for annotated junctions
                , "--outFilterMismatchNmax", "999" -- maximum number of mismatches per pair, large number switches off this filter
                , "--outFilterMismatchNoverReadLmax", "0.04" -- max number of mismatches per pair relative to read length: for 2x100b, max number of mismatches is 0.06*200=8 for the paired read
                , "--alignIntronMin", "20"        -- minimum intron length
                , "--alignIntronMax", "1000000"    -- maximum intron length
                , "--alignMatesGapMax", "1000000"  -- maximum genomic distance between mates
                , "--outSAMunmapped", "Within", "--outSAMattributes"
                , "NH", "HI", "AS", "NM", "MD"
                , "--outSAMheaderCommentFile", "COfile.txt"
                , "--outSAMheaderHD", "@HD", "VN:1.4", "SO:coordinate" ] ++
                ( if opt^.starSort
                    then [ "--outSAMtype", "BAM", "SortedByCoordinate"
                         , "--limitBAMsortRAM", "60000000000" ]
                    else ["--outSAMtype", "BAM", "Unsorted"] ) ++
                ( if pairedEnd e
                    then []
                    else ["--outSAMstrandField", "intronMotif"] ) ++
                ["--quantMode", "TranscriptomeSAM", "--sjdbScore", "1"]

            let starOutput | opt^.starSort = "/Aligned.sortedByCoord.out.bam"
                           | otherwise = "/Aligned.out.bam"
            mv (fromText $ T.pack $ tmp_dir ++ starOutput) $ fromText $
                T.pack $ outputGenome

            -- Sorting annotation bam
            if opt^.starSort
                then cmd "samtools" "sort" "-@" (T.pack $ show $ opt^.starCores)
                        "-T" (T.pack $ tmp_dir ++ "/sort_bam_tmp") "-o"
                        (T.pack outputAnno)
                        (T.pack $ tmp_dir ++ "/Aligned.toTranscriptome.out.bam")
                else mv (fromText $ T.pack $ tmp_dir ++
                        "/Aligned.toTranscriptome.out.bam") $ fromText $
                        T.pack outputAnno

            let genomeAlignFile = Single $
                    format .~ BamFile $
                    location .~ outputGenome $
                    tags .~ ["RNA genome align bam"] $ emptyFile
                annoFile = Single $
                    format .~ BamFile $
                    location .~ outputAnno $
                    tags .~ ["RNA anno align bam"] $ emptyFile
            return [genomeAlignFile, annoFile]
      where
        fls = case flset of
            Single f -> if not (pairedEnd e) && isFastq f then [f] else []
            Pair f1 f2 -> if pairedEnd e && isFastq f1 && isFastq f2
                then [f1,f2] else []


rsemMkIndex :: FilePath   -- ^ Prefix
            -> FilePath   -- ^ annotation file in GFF3 format
            -> [FilePath] -- ^ fastq files
            -> IO FilePath
rsemMkIndex prefix anno fstqs = do
    dirExist <- shelly $ test_d $ fromText $ T.pack dir
    if dirExist
        then hPutStrLn stderr "RSEM index directory exists. Skipped."
        else shelly $ do
            mkdir_p $ fromText $ T.pack dir
            liftIO $ hPutStrLn stderr "Generating RSEM indices"
            cmd "rsem-prepare-reference" "--gtf" (T.pack anno)
                (T.intercalate "," $ map T.pack fstqs) $ T.pack prefix
    return prefix
  where
    dir = takeDirectory prefix


data RSEMOpts = RSEMOpts
    { _rsemPath  :: FilePath
    , _rsemCores :: Int
    , _rsemSeed  :: Int32
    }

makeLenses ''RSEMOpts

type RSEMOptSetter = State RSEMOpts ()

defaultRSEMOpts :: RSEMOpts
defaultRSEMOpts = RSEMOpts
    { _rsemPath = ""
    , _rsemCores = 1
    , _rsemSeed = 12345
    }

-- | Gene and transcript quantification using rsem
rsemQuant :: FilePath
          -> FilePath
          -> RSEMOptSetter
          -> RNASeq
          -> IO RNASeq
rsemQuant dir indexPrefix setter = mapOfFiles fn
  where
    opt = execState setter defaultRSEMOpts
    isAnnoBam f = f^.format == BamFile && "RNA anno align bam" `elem` f^.tags
    fn e r (Single fl) = if not (isAnnoBam fl) then return [] else shelly $ do
        mkdir_p $ fromText $ T.pack dir
        let input = fl^.location
            output = printf "%s/%s_rep%d_rsem" dir (T.unpack $ e^.eid) $ r^.number

        run_ (fromText $ T.pack $ opt^.rsemPath ++ "rsem-calculate-expression") $
            [ "--bam", "--estimate-rspd", "--calc-ci", "--seed"
            , T.pack $ show $ opt^.rsemSeed, "-p", T.pack $ show $ opt^.rsemCores
            , "--no-bam-output", "--ci-memory", "30000" ] ++
            ( if pairedEnd e
                then ["--paired-end", "--forward-prob", "0"]
                else [] ) ++
            [T.pack input, T.pack indexPrefix, T.pack output]

        let geneQuant = Single $ format .~ Other $
                location .~ output ++ ".genes.results" $
                tags .~ ["gene quantification"] $ fl
            transcirptQuant = Single $ format .~ Other $
                location .~ output ++ ".isoforms.results" $
                tags .~ ["transcript quantification"] $ fl
        return [geneQuant, transcirptQuant]
    fn _ _ _ = return []
