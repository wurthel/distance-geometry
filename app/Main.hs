module Main where

import Data.Either
import DistanceGeometry
import Generators
import IO
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data
import System.Directory

idir = "test/molecules"
odir = "test/results"

main :: IO ()
main = do
  putStrLn ""
  putStrLn "generateOneMolecule: Begin"
  --generateOneMolecule idir odir "example"
  putStrLn "generateOneMolecule: End"
  putStrLn ""
  putStrLn "generateManyMolecules: Begin"
  generateManyMolecules idir odir "C2Br4" 0.01 100
  putStrLn "generateManyMolecules: End"
