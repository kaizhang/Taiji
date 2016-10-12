{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Main where

import           Bio.Pipeline.Instances   ()
import           Scientific.Workflow.Main

import qualified Component.Initialization as Initialization
import qualified Component.ATACSeq        as ATACSeq
import qualified Component.Network        as Network
import qualified Component.RNASeq         as RNASeq
import qualified Component.TFBS           as TFBS

defaultMain $ Initialization.builder >> ATACSeq.builder >>
    Network.builder >> RNASeq.builder >> TFBS.builder
