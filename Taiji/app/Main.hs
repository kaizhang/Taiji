{-# LANGUAGE CPP               #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE ViewPatterns #-}

module Main where

import           Bio.Pipeline.Instances           ()
import           Data.Version                     (showVersion)
import           Paths_Taiji                      (version)
import           Scientific.Workflow.Main         (MainOpts (..),
                                                   defaultMainOpts, mainWith)
import           Text.Printf                      (printf)

import qualified Taiji.Component.ATACSeq                as ATACSeq
import qualified Taiji.Component.Initialization         as Initialization
import qualified Taiji.Component.Network                as Network
import qualified Taiji.Component.Rank                   as Rank
import qualified Taiji.Component.RNASeq                 as RNASeq
import qualified Taiji.Component.TFBS                   as TFBS
import qualified Taiji.Component.Exporter as Exporter

#ifdef IDR_PEAK_CALLER
import qualified Taiji.Component.ATACSeq.CallPeak.IDR   as CallPeak
#else
import qualified Taiji.Component.ATACSeq.CallPeak.MACS2 as CallPeak
#endif


mainWith defaultMainOpts
    { programHeader = printf "Taiji-v%s" (showVersion version)
    } $ Initialization.builder >> ATACSeq.builder >> RNASeq.builder >>
        TFBS.builder >> Network.builder >> Rank.builder >> CallPeak.builder >>
        Exporter.builder
