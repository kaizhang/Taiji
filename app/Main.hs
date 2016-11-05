{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Main where

import           Bio.Pipeline.Instances   ()
import           Data.Version             (showVersion)
import           Paths_Taiji              (version)
import           Scientific.Workflow.Main (MainOpts (..), defaultMainOpts,
                                           mainWith)
import           Text.Printf              (printf)

import qualified Component.ATACSeq        as ATACSeq
import qualified Component.Initialization as Initialization
import qualified Component.Network        as Network
import qualified Component.Rank           as Rank
import qualified Component.RNASeq         as RNASeq
import qualified Component.TFBS           as TFBS

mainWith defaultMainOpts
    { programHeader = printf "Taiji-v%s" (showVersion version)
    } $ Initialization.builder >> ATACSeq.builder >> Network.builder >>
        RNASeq.builder >> TFBS.builder >> Rank.builder
