--------------------------------------------------------------------------------
-- RNA-seq data processing
--------------------------------------------------------------------------------
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Taiji.Component.RNASeq (builder) where

import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Utils
import           Bio.Pipeline.NGS
import           Bio.RealWorld.GENCODE
import           Bio.Utils.Misc                    (readDouble)
import           Control.Arrow                     (second)
import           Control.Lens
import           Control.Monad                     (forM)
import           Control.Monad.IO.Class            (liftIO)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (CI, mk, original)
import           Data.Double.Conversion.ByteString (toShortest)
import qualified Data.HashMap.Strict               as M
import           Data.List
import           Data.List.Ordered                 (nubSort)
import           Data.Maybe                        (fromJust)
import qualified Data.Text                         as T
import           Scientific.Workflow

import           Taiji.Constants

builder :: Builder ()
builder = do
    node "Get_RNA_data" [| return . (^._3) |] $ do
        submitToRemote .= Just False
        label .= "Get RNA-seq data"
    node "RNA_alignment_prepare" [| \input ->
        return $ concatMap splitExpByFile $ filterExpByFile
            (\x -> formatIs FastqFile x || formatIs FastqGZip x) input
        |] $ submitToRemote .= Just False
    node "RNA_alignment" [| \x -> starAlign <$> rnaOutput <*>
        getConfig' "starIndex" <*> return (starCores .= 4) <*> return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True >> remoteParam .= "-l vmem=10G -pe smp 4"
    node "RNA_quantification" [| \x -> rsemQuant <$> rnaOutput <*>
            fmap fromJust rsemIndex <*> return (rsemCores .= 4) <*>
            return x >>= liftIO
        |] $ batch .= 1 >> stateful .= True >> remoteParam .= "-l vmem=10G -pe smp 4"
    node "RNA_convert_ID_to_name" [| \x -> geneId2Name <$> rnaOutput <*>
            getConfig' "annotation" <*> return x >>= liftIO
        |] $ batch .= 5 >> stateful .= True
    path [ "Initialization", "Get_RNA_data", "RNA_alignment_prepare", "RNA_alignment"
         , "RNA_quantification", "RNA_convert_ID_to_name"]

    -- Gene expression profile can optionally be provided in original input data.
    node "RNA_average_prepare" [| \(oriInput, quant) -> do
        let filt (Single fl) = "gene quantification" `elem` fl^.tags
            filt _ = False
        return $ mergeExps $ filterExpByFile filt oriInput ++ quant
        |] $ submitToRemote .= Just False
    node "RNA_average" [| \x -> averageExpr <$> rnaOutput <*>
        return x >>= liftIO
        |] $ batch .= 5 >> stateful .= True
    node "Output_expression" [| \x -> combineExpression <$>
            fmap (++"/gene_expression.tsv") rnaOutput <*> return x >>= liftIO
        |] $ stateful .= True
    ["Get_RNA_data", "RNA_convert_ID_to_name"] ~> "RNA_average_prepare"
    path ["RNA_average_prepare", "RNA_average", "Output_expression"]


-- | Retrieve gene names
-- Requirements: one replicate should contain only one file with keywords
-- "gene quantification"
geneId2Name :: FilePath   -- ^ Output directory
            -> FilePath   -- ^ Annotation in GTF format
            -> RNASeq
            -> IO RNASeq
geneId2Name outdir anno e = do
    id2Name <- fmap (M.fromList . map (\x -> (geneId x, original $ geneName x))) $
        readGenes' anno

    rs <- forM (e^.replicates) $ \r -> do
        let [fl] = r^..files.folded._Single.
                filtered (elem "gene quantification" . (^.tags))
            output = outdir ++ "/" ++ T.unpack (e^.eid) ++ "_rep" ++
                show (r^.number) ++ "_TPM_by_names.tsv"
            newFile = Single $ location .~ output $
                tags .~ ["gene quantification"] $ emptyFile

        c <- B.readFile $ fl^.location
        B.writeFile output $ B.unlines $ map ( (\xs -> B.intercalate "\t"
            [M.lookupDefault undefined (head xs) id2Name, xs!!4]) . B.split '\t' ) $
            tail $ B.lines c
        return $ files .~ [newFile] $ r

    return $ replicates .~ rs $ e
{-# INLINE geneId2Name #-}


-- | Retrieve gene names and compute the average expression of replicates.
averageExpr :: FilePath   -- ^ Output directory
            -> RNASeq
            -> IO RNASeq
averageExpr outdir e = do
    let fls = e^..replicates.folded.files.folded._Single.
            filtered (elem "gene quantification" . (^.tags))
        output = outdir ++ "/" ++ T.unpack (e^.eid) ++ "_average_gene_quant.tsv"
        newFile = Single $ location .~ output $
            tags .~ ["average gene quantification"] $ emptyFile
        readExpr fl = do
            c <- B.readFile fl
            return $ map (\xs -> let [a,b] = B.split '\t' xs in (mk a, readDouble b)) $
                B.lines c
    expr <- mapM (readExpr . (^.location)) fls
    B.writeFile output $ B.unlines $
        map (\(a,b) -> B.intercalate "\t" [original a, toShortest $ average b]) $
        combine expr
    return $ replicates .~ [files .~ [newFile] $ emptyReplicate] $ e
{-# INLINE averageExpr #-}

combine :: [[(CI B.ByteString, Double)]] -> [(CI B.ByteString, [Double])]
combine xs = flip map names $ \nm -> (nm, map (M.lookupDefault 0.01 nm) xs')
  where
    names = nubSort $ concatMap (fst . unzip) xs
    xs' = map (fmap average . M.fromListWith (++) . map (second return)) xs
{-# INLINE combine #-}

average :: [Double] -> Double
average [a] = a
average [a,b] = (a + b) / 2
average [a,b,c] = (a + b + c) / 3
average xs = foldl1' (+) xs / fromIntegral (length xs)
{-# INLINE average #-}

-- | Combine RNA expression data into a table.
combineExpression :: FilePath
                  -> [RNASeq]
                  -> IO (Maybe FilePath)
combineExpression output es
    | null es = return Nothing
    | otherwise = do
        dat <- forM es $ \e -> do
            let [fl] = e^..replicates.folded.filtered ((==0) . (^.number)).
                    files.folded._Single.
                    filtered (elem "average gene quantification" . (^.tags))
            expr <- readExpr $ fl^.location
            return (fromJust $ e^.groupName, expr)
        let (expNames, values) = unzip dat
        B.writeFile output $ B.unlines $
            B.pack (T.unpack $ T.intercalate "\t" $ "Name" : expNames) :
            map (\(x,xs) -> B.intercalate "\t" $ original x : map toShortest xs)
                (combine values)
        return $ Just output
  where
    readExpr fl = do
        c <- B.readFile fl
        return $ map (\xs -> let [a,b] = B.split '\t' xs in (mk a, readDouble b)) $
            B.lines c
{-# INLINE combineExpression #-}
