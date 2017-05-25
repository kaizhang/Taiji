{-# LANGUAGE OverloadedStrings #-}

import           Bio.Data.Experiment.Types
import qualified Data.Map                  as M
import           Data.Maybe
import           Data.Yaml
import           Test.Tasty
import           Test.Tasty.HUnit
import Control.Lens
import System.IO.Unsafe

chipData :: [ChIPSeq]
chipData = fromJust $ unsafePerformIO $ decodeFile "tests/test.yml"

{-
chipData :: ChIPSeq
chipData = eid .~ "chip_id" $
           celltype .~ "chip_cell" $
           groupName .~ Just "chip_group" $
           target .~ "chip_target" $
           control .~ Just "chip_control" $ defaultChIPSeq

commonParseTest :: Assertion
commonParseTest = chipData^.commonFields @=? fromEither
    (decodeEither $ encode $ chipData^.commonFields)
    -}

chipParseTest :: Assertion
chipParseTest = chipData @=? fromEither (decodeEither $ encode chipData)

fromEither (Left x) = error x
fromEither (Right x) = x

main :: IO ()
main = defaultMain $ testGroup "Main"
    [ testCase "ChIP-Seq data parsing" chipParseTest
    ]
