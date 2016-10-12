--------------------------------------------------------------------------------
-- Network construction
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Component.Network (builder) where

import           Bio.Data.Bed
import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Utils
import           Bio.GO.GREAT
import           Bio.Pipeline.Instances            ()
import           Bio.Utils.Functions               (scale)
import           Bio.Utils.Misc                    (readDouble, readInt)
import           Conduit
import           Control.Arrow                     (first, second, (&&&), (***))
import           Control.Lens
import           Control.Monad
import           Data.Binary                       (decodeFile, encodeFile)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (mk, original)
import           Data.Double.Conversion.ByteString (toShortest)
import           Data.Function                     (on)
import qualified Data.HashMap.Strict               as M
import           Data.List
import           Data.List.Ordered                 (nubSort)
import           Data.Maybe
import           Data.Ord
import qualified Data.Text                         as T
import qualified Data.Vector.Unboxed               as U
import           IGraph
import           IGraph.Structure                  (pagerank,
                                                    personalizedPagerank)
import           Scientific.Workflow
import           Statistics.Sample                 (meanVarianceUnb)
import           System.IO                         (hPutStrLn, stderr)

import           Constants
import           Types

builder :: Builder ()
builder = do
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

    node "net00" [| \(x, expr) -> do
        expression <- getConfigMaybe' "expression_profile"
        liftIO $ do
            gene_expr <- case () of
                _ | isJust expression -> do
                        hPutStrLn stderr "Use user supplied gene expression profile."
                        return expression
                  | isJust expr -> return expr
                  | otherwise -> return Nothing
            case gene_expr of
                Nothing -> do
                    hPutStrLn stderr "Running PageRank..."
                    pageRank x
                Just e -> do
                    hPutStrLn stderr "Running personalized PageRank..."
                    personalizedPageRank (e, x)
        |] $ stateful .= True
    ["ass02", "rna03"] ~> "net00"

    node "vis00" [| \x -> outputData "ranks" x (getMetrics x) |] $ submitToRemote .= Just False
    ["net00"] ~> "vis00"


getDomains :: FilePath   -- output
           -> FilePath   -- annotation
           -> [Experiment ATAC_Seq] -> IO [Experiment ATAC_Seq]
getDomains dir anno = mapM $ \e -> do
    let [fl] = filter ((==NarrowPeakFile) . (^.format)) $ e^.files
        output = dir ++ "/" ++ T.unpack (e^.eid) ++ "_gene_reg_domains.bed"
        newFile = format .~ BedFile $
                  info .~ [] $
                  keywords .~ ["gene regulatory domain"] $
                  location .~ output $ fl
    peaks <- readBed' $ fl^.location
    activePromoters <- gencodeActiveGenes anno peaks
    writeBed' output $ getGeneDomains GREAT activePromoters
    return $ files %~ (newFile:) $ e

linkGeneToTFs :: FilePath -> (Experiment ATAC_Seq, FilePath) -> IO (Experiment ATAC_Seq)
linkGeneToTFs dir (e, tfbs) = do
    let [peakFl] = filter ((==NarrowPeakFile) . (^.format)) $ e^.files
        [domainFl] = filter ((==["gene regulatory domain"]) . (^.keywords)) $ e^.files
    genes <- readBed' $ domainFl^.location
    peaks <- readBed' $ peakFl^.location

    regulators <- findRegulators genes peaks tfbs
    let result :: [Linkage]
        result = map (second (map ((head *** id) . unzip) . groupBy ((==) `on` fst) .
            sortBy (comparing fst) . map (getTFName &&& id))) $ M.toList $
            M.fromListWith (++) $ map (first (mk . fromJust . bedName)) $
            regulators
        output = dir ++ "/" ++ T.unpack (e^.eid) ++ ".assign"
        newFile = format .~ Other $
                  info .~ [] $
                  keywords .~ ["gene-TF assignment"] $
                  location .~ output $ peakFl
    encodeFile output result
    return $ files .~ [newFile] $ e
  where
    getTFName = mk . head . B.split '+' . fromJust . bedName

findRegulators :: [BED]          -- ^ gene domains
               -> [NarrowPeak]   -- ^ open chromatin
               -> FilePath          -- ^ TFBS
               -> IO [(BED, [BED])]
findRegulators genes peaks tfbs = do
    tfbs' <- readBed tfbs =$= intersectBed peaks' $$ sinkList
    yieldMany genes =$= intersectBedWith id tfbs' =$=
        filterC (not . null .snd) $$ sinkList
  where
    peaks' = flip map peaks $ \p -> let c = chromStart p + (fromJust . _npPeak) p
                                    in BED3 (chrom p) (c - 50) (c+ 50)

data Method = GREAT  -- ^ using GREAT strategy

getGeneDomains :: Method -> [((B.ByteString, Int, Bool), B.ByteString)]
               -> [BED]   -- ^ genes' regulatory domains. A single gene can
                          -- potentially have multiple domains.
getGeneDomains GREAT = map ( \(b, x) ->
    BED (chrom b) (chromStart b) (chromEnd b) (Just x) Nothing Nothing ) .
    getRegulatoryDomains (BasalPlusExtension 5000 1000 1000000)

-- | Identify active genes by overlapping their promoters with activity indicators.
gencodeActiveGenes :: FilePath   -- ^ gencode file in GTF format
                   -> [BED3]     -- ^ feature that is used to determine the activity
                                 -- of promoters, e.g., H3K4me3 peaks or ATAC-seq peaks
                   -> IO [((B.ByteString, Int, Bool), B.ByteString)]  -- ^ chr, tss, strand and gene name.
gencodeActiveGenes input peaks = do
    c <- B.readFile input
    let promoters = map ( \((chr, i, isForward), geneName) ->
            if isForward
                then BED chr (i-5000) (i+1000) (Just geneName)
                     (Just $ fromIntegral i) (Just True)
                else BED chr (i-1000) (i+5000) (Just geneName)
                     (Just $ fromIntegral i) (Just False)
            ) $ map f $ filter g $ map (B.split '\t') $ B.lines c
    return $ map (\x -> ( (chrom x, truncate $ fromJust $ bedScore x,
        fromJust $ bedStrand x), fromJust $ bedName x )) $ fst $ unzip $
        runIdentity $ yieldMany promoters =$= intersectBedWith id peaks =$=
        filterC (not . null . snd) $$ sinkList
  where
    f xs = let name = B.filter (/='\"') $ last $ B.split ' ' $ head $
                      filter (B.isInfixOf "gene_name") $ B.split ';' $ last xs
               start = readInt $ xs !! 3
               end = readInt $ xs !! 4
               strand = if xs !! 6 == "+" then True else False
          in ((xs !! 0, if strand then start else end, strand), name)
    g xs = not $ B.isPrefixOf "#" (head xs) || (xs !! 2 /= "transcript")

-- | Read RNA expression data
readExpression :: FilePath
               -> IO [( B.ByteString    -- ^ cell type
                     , M.HashMap GeneName (Double, Double)  -- ^ absolute value and z-score
                     )]
readExpression fl = do
    c <- B.readFile fl
    let ((_:header):dat) = map (B.split '\t') $ B.lines c
        rowNames = map (mk . head) dat
        dataTable = map (map readDouble . tail) dat
    return $ zipWith (\a b -> (a, M.fromList $ zip rowNames b)) header $
        transpose $ zipWith zip dataTable $ map computeZscore dataTable
  where
    computeZscore xs
        | all (<1) xs || all (==head xs) xs = replicate (length xs) (-10)
        | otherwise = U.toList $ scale $ U.fromList xs

pageRank :: [Experiment a] -> IO ([T.Text], [B.ByteString], [[Double]])
pageRank es = do
    results <- forM es $ \e -> do
        gr <- buildNet e
        let labs = map (nodeLab gr) $ nodes gr
        return $ zip labs $ pagerank gr Nothing 0.85
    let genes = nubSort $ concatMap (fst . unzip) results
        expNames = map (fromJust . (^.groupName)) es
        ranks = flip map results $ \xs ->
            let geneRanks = M.fromList xs
            in flip map genes $ \g -> M.lookupDefault 0 g geneRanks
    return (expNames, map original genes, transpose ranks)

personalizedPageRank :: (FilePath, [Experiment a])
                     -> IO ([T.Text], [B.ByteString], [[Double]])
personalizedPageRank (rnaseq, es) = do
    rnaseqData <- readExpression rnaseq

    results <- forM es $ \e -> do
        gr <- buildNet e
        let lookupExpr x = M.lookupDefault (0.01,-10) x $ fromJust $ lookup
                (B.pack $ T.unpack $ fromJust $ e^.groupName) rnaseqData
            nodeWeights = map (exp . snd . lookupExpr) labs
            edgeWeights = map (sqrt . fst . lookupExpr . nodeLab gr . snd) $ edges gr
            labs = map (nodeLab gr) $ nodes gr
        return $ zip labs $
            personalizedPagerank gr nodeWeights (Just edgeWeights) 0.85
    let genes = nubSort $ concatMap (fst . unzip) results
        expNames = map (fromJust . (^.groupName)) es
        ranks = flip map results $ \xs ->
            let geneRanks = M.fromList xs
            in flip map genes $ \g -> M.lookupDefault 0 g geneRanks
    return (expNames, map original genes, transpose ranks)

buildNet :: Experiment a -> IO (LGraph D GeneName ())
buildNet e = do
    let [fl] = filter (\x -> x^.keywords == ["gene-TF assignment"]) $ e^.files
    result <- decodeFile $ fl^.location :: IO [Linkage]
    return $ fromLabeledEdges $ flip concatMap result $ \(a, b) ->
        zip (zip (repeat a) $ fst $ unzip b) $ repeat ()

printEdgeList :: FilePath -> Experiment a -> IO ()
printEdgeList dir e = do
    let [fl] = filter (\x -> x^.keywords == ["gene-TF assignment"]) $ e^.files
    result <- decodeFile $ fl^.location :: IO [Linkage]
    let output = dir ++ "/" ++ T.unpack (e^.eid) ++ "_network.tsv"
    B.writeFile output $ B.unlines $ flip map result $ \(a,b) ->
        B.intercalate "\t" [original a, B.intercalate "," $ map original $ fst $ unzip b]

writeTSV :: FilePath -> ([T.Text], [B.ByteString], [[Double]]) -> IO ()
writeTSV output (colNames, rowNames, dat) = B.writeFile output $ B.unlines $
    B.intercalate "\t" ("Gene" : map (B.pack . T.unpack) colNames) :
    map (\(a,b) -> B.intercalate "\t" $ a : map toShortest b)
        (zip rowNames dat)

data Metrics = Metrics
    { averageRank   :: Double
    , variability   :: Double
    , logFoldChange :: Double
    }

getMetrics :: ([T.Text], [B.ByteString], [[Double]]) -> [(B.ByteString, Metrics)]
getMetrics (cts, geneNames, dat) = zip geneNames $ map (f . U.fromList) dat
  where
    f xs = let (m, v) = meanVarianceUnb xs
               fold = log $ max 1 $ U.maximum xs / U.minimum xs
           in Metrics m (sqrt v / m) fold

outputData :: FilePath
           -> ([T.Text], [B.ByteString], [[Double]])
           -> [(B.ByteString, Metrics)] -> IO ()
outputData out (cts, geneNames, dat) metrics = do
    let filteredData = map toBS $ filter f $ zip metrics dat
        allData = map toBS $ zip metrics dat
        header = B.pack $ T.unpack $ T.intercalate "\t" $ "Gene" : cts
    B.writeFile (out++"_filtered.tsv") $ B.unlines $ header : filteredData
    B.writeFile (out++"_all.tsv") $ B.unlines $ header : allData
  where
    f ((_, Metrics m v fold), _) = m > 0.0001 && fold > 1.5
    toBS ((nm, _), xs) = B.intercalate "\t" $ nm : map toShortest xs
