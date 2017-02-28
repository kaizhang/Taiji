{-# LANGUAGE OverloadedStrings #-}

module Bio.Pipeline.Utils where

import           Bio.Data.Experiment.Types
import           Control.Lens
import qualified Data.Text                 as T

-- | Get the prefix of a filename, e.g., "path/xxx.txt" -> "path/xxx"
getPrefix :: T.Text -> T.Text
getPrefix x = if suffix == "gz"
    then T.init $ fst $ T.breakOnEnd "." $ T.init prefix
    else T.init prefix
  where
    (prefix, suffix) = T.breakOnEnd "." x

mapOfFiles :: (Experiment e, Monad m)
           => (e -> Replicate -> FileSet -> m [FileSet])
           -> e -> m e
mapOfFiles fn e = flip (id (replicates.traverse)) e $ \r ->
    id files (fmap concat . mapM (fn e r)) r
{-# INLINE mapOfFiles #-}
