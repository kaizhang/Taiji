module Config (config, (!)) where

import Data.Yaml (decodeFile)
import qualified Data.Map as M
import System.IO.Unsafe (unsafePerformIO)

defaultConfig :: M.Map String String
defaultConfig = M.fromList
    [ ("input", "data/input.yaml")
    , ("outputDir", "output/")
    , ("motifFile", "data/cisbp.meme")
    , ("genome", "/home/kai/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa")
    , ("genomeIndex", "/home/kai/data/genome/mm10.index")
    , ("gencode", "/home/kai/project/TRM/data/gencode.vM9.annotation.gff3")
    , ("STARIndex", "/home/kai/Mus_musculus/UCSC/mm10/Sequence/STARIndex/")
    ]

config :: M.Map String String
config = unsafePerformIO $ do
    Just c <- decodeFile "taiji.yml"
    return $ M.union c defaultConfig

(!) :: M.Map String String -> String -> String
(!) m x = M.findWithDefault undefined x m
