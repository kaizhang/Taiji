-- Gene expression
{-# LANGUAGE GADTs             #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}

module Expression where

import           Bio.Data.Experiment.Types
import           Bio.RealWorld.GENCODE
import           Bio.Utils.Functions               (scale)
import           Bio.Utils.Misc                    (readDouble)
import           Control.Arrow                     (second)
import           Control.Lens
import           Control.Monad                     (forM)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (mk, original)
import           Data.Double.Conversion.ByteString (toShortest)
import qualified Data.HashMap.Strict               as M
import           Data.List
import           Data.List.Ordered                 (nubSort)
import           Data.Maybe                        (fromJust)
import qualified Data.Text                         as T
import qualified Data.Vector.Unboxed               as U

import           Types

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

-- | Combine RNA expression data into a table
combineExpression :: FilePath
                  -> FilePath   -- ^ annotation in GTF format
                  -> [Experiment RNA_Seq]
                  -> IO (Maybe FilePath)
combineExpression output anno es
    | null es = return Nothing
    | otherwise = do
        id2Name <- fmap (M.fromList . map (\x -> (geneId x, original $ geneName x))) $
            readGenes' anno
        dat <- forM es $ \e -> do
            let fls = filter (elem "gene quantification" . (^.keywords)) $ e^.files
            expr <- mapM (readExpr . (^.location)) fls
            return (fromJust $ e^.groupName, map (second average) $ combine expr)
        let (expNames, values) = unzip dat
        B.writeFile output $ B.unlines $
            (B.pack $ T.unpack $ T.intercalate "\t" $ "Name" : expNames) :
            (map (\(x,xs) -> B.intercalate "\t" $
            M.lookupDefault (error $ show x) x id2Name :
            map toShortest xs) $ combine values)
        return $ Just output
  where
    combine xs = flip map names $ \nm -> (nm, map (M.lookupDefault 0.01 nm) xs')
      where
        names = nubSort $ concatMap (fst . unzip) xs
        xs' = map M.fromList xs
    readExpr fl = do
        c <- B.readFile fl
        return $ map ((\xs -> (head xs, readDouble $ xs!!4)) . B.split '\t') $
            tail $ B.lines c
    average xs = foldl' (+) 0 xs / fromIntegral (length xs)
