{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module:      Bio.Data.Experiment.Parser
-- Copyright:   (c) 2016 Kai Zhang
-- License:     BSD3
-- Maintainer:  Kai Zhang <kai@kzhang.org>
-- Stability:   experimental
-- Portability: portable
--
-- The default serialization of most types are too verbose and not suitable for
-- direct use by human beings. This module provids a more succint and
-- user-friendly way to write YAML files. Featuring:
-- 1. Most keywords have been renamed.
-- 2. Keywords are case-insensitive.
-- 3. Providing default values for some fields.

module Bio.Data.Experiment.Parser
    ( parseChIPSeq
    , parseATACSeq
    , parseRNASeq
    , parseHiC
    , parseList
    , guessFormat
    ) where

import           Control.Arrow                      (first)
import           Data.Aeson
import           Data.Aeson.Internal                (JSONPathElement (..),
                                                     (<?>))
import           Data.Aeson.Types
import qualified Data.HashMap.Strict                as HM
import qualified Data.Map.Strict                    as M
import qualified Data.Text                          as T
import qualified Data.Vector                        as V

import           Bio.Data.Experiment.Types.Internal

parseFile :: Value -> Parser File
parseFile = withObject "File" $ \obj' -> do
    let obj = toLowerKey obj'
    path <- obj .: "path"
    File <$> return path <*>
             obj .:? "format" .!= guessFormat path <*>
             obj .:? "info" .!= M.empty <*>
             obj .:? "tags" .!= []

guessFormat :: FilePath -> FileType
guessFormat fl = case () of
    _ | ".bam" `T.isSuffixOf` fl' -> BamFile
      | ".bai" `T.isSuffixOf` fl' -> BaiFile
      | ".bed" `T.isSuffixOf` fl' -> BedFile
      | ".bed.gz" `T.isSuffixOf` fl' -> BedGZip
      | ".fastq" `T.isSuffixOf` fl' -> FastqFile
      | ".fastq.gz" `T.isSuffixOf` fl' -> FastqGZip
      | ".bw" `T.isSuffixOf` fl' -> BigWigFile
      | otherwise -> Other
  where
    fl' = T.pack fl

parseFileSet :: Value -> Parser FileSet
parseFileSet = withObject "FileSet" $ \obj' -> do
    let obj = toLowerKey obj'
    fls <- obj .:? "pair"
    case fls of
        Nothing -> Single <$> parseFile (Object obj')
        Just array -> flip (withArray "FileSet") array $ \xs ->
            if V.length xs == 2
                then Pair <$> parseFile (xs `V.unsafeIndex` 0)
                          <*> parseFile (xs `V.unsafeIndex` 1)
                else error "The number of files must be 2."

parseReplicate :: Value -> Parser Replicate
parseReplicate = withObject "Replicate" $ \obj' -> do
    let obj = toLowerKey obj'
    Replicate <$> withParser (parseList parseFileSet) obj "files" <*>
                  obj .:? "info" .!= M.empty <*>
                  obj .:? "rep" .!= 0

parseCommonFields :: Value -> Parser CommonFields
parseCommonFields = withObject "CommonFields" $ \obj' -> do
    let obj = toLowerKey obj'
    CommonFields <$> obj .: "id" <*>
                     obj .:? "group" <*>
                     obj .:? "celltype" .!= "" <*>
                     withParser (parseList parseReplicate) obj "replicates"

parseChIPSeq :: Value -> Parser ChIPSeq
parseChIPSeq = withObject "ChIPSeq" $ \obj' -> do
    let obj = toLowerKey obj'
    ChIPSeq <$> parseCommonFields (Object obj') <*>
                obj .: "target" <*>
                obj .:? "pairedend" .!= False <*>
                obj .:? "control"

parseATACSeq :: Value -> Parser ATACSeq
parseATACSeq = withObject "ATACSeq" $ \obj' -> do
    let obj = toLowerKey obj'
    ATACSeq <$> parseCommonFields (Object obj') <*>
                obj .:? "pairedend" .!= False

parseRNASeq :: Value -> Parser RNASeq
parseRNASeq = withObject "RNASeq" $ \obj' -> do
    let obj = toLowerKey obj'
    RNASeq <$> parseCommonFields (Object obj') <*>
                obj .:? "pairedend" .!= False

parseHiC :: Value -> Parser HiC
parseHiC = withObject "HiC" $ \obj' -> do
    let obj = toLowerKey obj'
    HiC <$> parseCommonFields (Object obj')


--------------------------------------------------------------------------------

toLowerKey :: HM.HashMap T.Text a -> HM.HashMap T.Text a
toLowerKey = HM.fromList . map (first T.toLower) . HM.toList
{-# INLINE toLowerKey #-}

withParser :: (Value -> Parser a) -> Object -> T.Text -> Parser a
withParser p obj key = case HM.lookup key obj of
    Nothing -> fail $ "key " ++ show key ++ " not present"
    Just v  -> p v <?> Key key
{-# INLINE withParser #-}

parseList :: (Value -> Parser a) -> Value -> Parser [a]
parseList p (Array a) = mapM p $ V.toList a
parseList _ _ = error "Not a list"
{-# INLINE parseList #-}
