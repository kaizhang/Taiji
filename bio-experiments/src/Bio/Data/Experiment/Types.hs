{-# LANGUAGE DeriveGeneric          #-}
{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses  #-}
{-# LANGUAGE OverloadedStrings      #-}
{-# LANGUAGE TemplateHaskell        #-}

module Bio.Data.Experiment.Types
    ( FileType(..)
    , File
    , location
    , format
    , tags
    , emptyFile

    , FileSet(..)
    , _Single
    , _Pair

    , Replicate
    , emptyReplicate
    , files
    , info
    , number
    , Experiment(..)
    , NGS(..)
    , IsDNASeq

    , ChIPSeq
    , target
    , control
    , defaultChIPSeq

    , ATACSeq
    , RNASeq
    , HiC
    ) where

import Bio.Data.Experiment.Types.Internal
