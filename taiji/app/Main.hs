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

import qualified Component.ATACSeq                as ATACSeq
import qualified Component.Initialization         as Initialization

#ifdef IDR_PEAK_CALLER
import qualified Component.ATACSeq.CallPeak.IDR   as CallPeak
#else
import qualified Component.ATACSeq.CallPeak.MACS2 as CallPeak
#endif

import qualified Component.Network                as Network
import qualified Component.Rank                   as Rank
import qualified Component.RNASeq                 as RNASeq
import qualified Component.TFBS                   as TFBS

mainWith defaultMainOpts
    { programHeader = printf "Taiji-v%s" (showVersion version)
    } $ Initialization.builder >> ATACSeq.builder >> RNASeq.builder >>
        TFBS.builder >> Network.builder >> Rank.builder >> CallPeak.builder
