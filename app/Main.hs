{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
module Main where

import Scientific.Workflow.Main
import Bio.Pipeline.Instances ()

import Builder (graph)

defaultMain graph
