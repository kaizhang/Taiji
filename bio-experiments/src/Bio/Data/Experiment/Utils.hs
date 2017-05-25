{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
module Bio.Data.Experiment.Utils where

import           Control.Lens              ((.~), (^.), (%~))
import           Data.Function             (on)
import           Data.List                 (foldl1', groupBy, sortBy)
import           Data.List.Ordered         (nubSort)
import qualified Data.Map.Strict           as M
import           Data.Maybe                (mapMaybe)
import           Data.Ord                  (comparing)
import qualified Data.Text                 as T

import           Bio.Data.Experiment.Types

formatIs :: FileType -> FileSet -> Bool
formatIs t (Single f) = f^.format == t
formatIs t (Pair f1 f2) = f1^.format == t && f2^.format == t
{-# INLINE formatIs #-}

-- | Split a single experiment into multiple experiments, each containing a
-- fileset.
splitExpByFile :: Experiment e => e -> [e]
splitExpByFile e = zipWith (\x y -> replicates .~ [x] $ y)
    (concatMap f $ e^.replicates) $ repeat e
  where
    f r = zipWith (\x y -> files .~ [x] $ y) (r^.files) $ repeat r
{-# INLINE splitExpByFile #-}

-- | Merge experiments with same id.
mergeExps :: Experiment e => [e] -> [e]
mergeExps es = map combineExp expGroup
  where
    expGroup = groupBy ((==) `on` (^.eid)) $ sortBy (comparing (^.eid)) es
    combineExp e = if allEqual (map (\x -> (x^.groupName, x^.cellType)) e)
        then replicates .~ mergeReps (concatMap (^.replicates) e) $ head e
        else error "Abort: Found experiments with same id but with different contents"
    allEqual (x:xs) = all (==x) xs
    allEqual _ = True
{-# INLINE mergeExps #-}

-- | Merge replicates with same number.
mergeReps :: [Replicate] -> [Replicate]
mergeReps rs = map (foldl1' combineRep) repGroup
  where
    repGroup = groupBy ((==) `on` (^.number)) $ sortBy (comparing (^.number)) rs
    combineRep r1 r2 = files .~ nubSort (r1^.files ++ r2^.files) $
        info .~ M.unionWith f (r1^.info) (r2^.info) $ r1
    f a b | a == b = a
          | otherwise = a `T.append` " | " `T.append` b
{-# INLINE mergeReps #-}

-- | Keep Experiments, Replicates and FileSets that satisfy the predicate.
filterExp :: Experiment e
          => (Replicate -> Bool)   -- ^ first fiter by replicates
          -> (FileSet -> Bool)     -- ^ then filter by files
          -> [e] -> [e]
filterExp f g = mapMaybe $ \e ->
    let rs = flip mapMaybe (filter f $ e^.replicates) $ \r ->
            let fls = filter g $ r^.files
            in if null fls then Nothing else Just $ files .~ fls $ r
    in if null rs then Nothing else Just $ replicates .~ rs $ e
{-# INLINE filterExp #-}

-- | Keep those experiments having specific files.
filterExpByFile :: Experiment e => (FileSet -> Bool) -> [e] -> [e]
filterExpByFile f = mapMaybe $ \e ->
    let e' = replicates %~ filterRepByFile f $ e
    in if null (e'^.replicates) then Nothing else Just e'
{-# INLINE filterExpByFile #-}

-- | Keep those experiments having specific files.
filterRepByFile :: (FileSet -> Bool) -> [Replicate] -> [Replicate]
filterRepByFile f = mapMaybe $ \r ->
    let r' = files %~ filter f $ r
    in if null (r'^.files) then Nothing else Just r'
{-# INLINE filterRepByFile #-}
