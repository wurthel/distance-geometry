module Main where

import Data.Either
import DistanceGeometry
import IO
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data

main :: IO ()
main = do
  let (atoms, bonds) = readMoleculeMolV2000 "example/example.mol"
  let 
    s0@(u, l) = (triangleInequalitySmoothingFloyd . generateDistanceBoundsMatrix atoms) bonds
  s1 <- randomDistanceMatrix s0
  let s2 = (generateCoordinFromEigValAndVec . largestEigValAndVec . distanceMatrixToMetricMatrix) s1
  let dm = coordMatrixToDistanceMatrix s2
      errf1 = distanceErrorFunction1 dm u l
      errf2 = distanceErrorFunction2 dm u l
      errf3 = distanceErrorFunction3 dm u l
  print atoms
  writeMoleculeXYZ "test.xyz" atoms
  print s2
  print errf1
  print errf2
  print errf3