module Main where

import           Data.Either
import           DistanceGeometry
import           Generators
import           Numeric.LinearAlgebra
import           Numeric.LinearAlgebra.Data
import           ReadWrite
import           System.Directory

idir = "test/molecules"
odir = "test/results"

main :: IO ()
main = do
  putStrLn ""
  putStrLn "generateOneMolecule: Begin"
  generateOneMolecule idir odir "example" 16
  putStrLn "generateOneMolecule: End"
  putStrLn ""
  putStrLn "generateManyMolecules: Begin"
  generateManyMolecules idir odir "methan" 16 100
  putStrLn "generateManyMolecules: End"
