{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
module Main where

import Scientific.Workflow.Main
import Bio.Pipeline.Instances ()

import Pipeline
import qualified Network

defaultMain $ pipeline >> Network.graph
