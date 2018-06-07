{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
module Taiji.Tools where

import Control.Monad
import Data.Maybe
import           Bio.Pipeline.Utils        (mapOfFiles)
import Bio.Data.Experiment.Types
import Bio.Data.Experiment.Parser (guessFormat)
import Control.Lens
import Conduit
import Data.List (isPrefixOf)
import Network.HTTP.Conduit
import qualified Data.ByteString.Char8 as B
import Control.Monad.IO.Class (liftIO)
import qualified Data.Text as T
import Shelly hiding (FilePath)
import           Data.List.Split           (splitOn)

downloadData :: (NGS e, Experiment e)
             => FilePath
             -> e
             -> IO e
downloadData outDir = mapOfFiles fn
  where
    fn e _ (Single fl) = case () of
        _ | fl^.format == SRA -> if pairedEnd e then dumpPair fl else dumpSingle fl
          | otherwise -> do
              fl' <- downloadfile fl
              return [Single fl']
    fn _ _ (Pair a b) = do
        f1 <- downloadfile a
        f2 <- downloadfile b
        return [Pair f1 f2]
    downloadfile fl
      | "ENCODE:" `isPrefixOf` (fl^.location) = do
          let acc = T.unpack $ T.strip $ snd $ T.breakOnEnd "ENCODE:" $ T.pack $
                fl^.location
          f1_name <- downloadENCODE acc outDir
          let file_format = if fl^.format == Other then guessFormat f1_name else fl^.format
              f1 = location .~ f1_name $ format .~ file_format $ fl
          return f1
      | otherwise = return fl
    dumpSingle fl = do
        let f1_name = outDir ++ "/" ++ fl^.location ++ ".fastq.gz"
            f1 = format .~ FastqGZip $ location .~ f1_name $ fl
        case splitOn "+" (fl^.location) of
            [f] -> shelly $ run_ "fastq-dump" ["--origfmt", "-I", "--gzip",
                "-O" ,T.pack outDir, T.pack f]
            fs -> do
                forM_ fs $ \f -> shelly $ run_ "fastq-dump" ["--origfmt", "-I",
                    "-O", T.pack outDir, T.pack f]
                let f1s = map (\x -> T.pack $ outDir ++ "/" ++ x ++ ".fastq") fs
                shelly $ escaping False $ do
                    run_ "cat" $ f1s ++ ["|", "gzip", "-c", ">", T.pack f1_name]
                    run_ "rm" f1s
        return [Single f1]
    dumpPair fl = do
        let f1_name = outDir ++ "/" ++ fl^.location ++ "_1.fastq.gz"
            f2_name = outDir ++ "/" ++ fl^.location ++ "_2.fastq.gz"
            f1 = format .~ FastqGZip $ location .~ f1_name $ fl
            f2 = format .~ FastqGZip $ location .~ f2_name $ fl

        case splitOn "+" (fl^.location) of
            [f] -> shelly $ run_ "fastq-dump" ["--origfmt", "-I",
                "--split-files", "--gzip", "-O", T.pack outDir, T.pack f]
            fs -> do
                forM_ fs $ \f -> shelly $ run_ "fastq-dump" ["--origfmt", "-I",
                    "--split-files", "-O", T.pack outDir, T.pack f]
                let f1s = map (\x -> T.pack $ outDir ++ "/" ++ x ++ "_1.fastq") fs
                    f2s = map (\x -> T.pack $ outDir ++ "/" ++ x ++ "_2.fastq") fs
                shelly $ escaping False $ do
                    run_ "cat" $ f1s ++ ["|", "gzip", "-c", ">", T.pack f1_name]
                    run_ "rm" f1s
                    run_ "cat" $ f2s ++ ["|", "gzip", "-c", ">", T.pack f2_name]
                    run_ "rm" f2s
        return [Pair f1 f2]


-- | Download data from ENCODE portal to a given directory.
downloadENCODE :: String    -- ^ Accession number
               -> FilePath  -- ^ Output dir
               -> IO FilePath
downloadENCODE acc dir = do
     request <- parseRequest url
     manager <- newManager tlsManagerSettings
     runResourceT $ do
         response <- http request manager
         let filename = T.unpack $ snd $ T.breakOnEnd "filename=" $ T.pack $
                B.unpack $ fromJust $ lookup "Content-Disposition" $
                responseHeaders response
         runConduit $ responseBody response .| sinkFileBS (dir ++ "/" ++ filename)
         return $ dir ++ "/" ++ filename
  where
    url = "https://www.encodeproject.org/files/" ++ acc ++ "/@@download"
