module Main where

import Data.Either
import DistanceGeometry
import IO
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data

main :: IO ()
main = do
  let (atoms, bonds) = readMolV2000 "example/methan.mol"
  let s0@(u, l) =
        (triangleInequalitySmoothingFloyd . generateDistanceBoundsMatrix atoms) bonds
  s1 <- randomDistanceMatrix s0
  let s2 =
        (generateCoordinFromEigValAndVec .
         largestEigValAndVec . distanceMatrixToMetricMatrix)
          s1
  let dm = coordMatrixToDistanceMatrix s2
      errf1 = distanceErrorFunction1 dm u l
      errf2 = distanceErrorFunction2 dm u l
      errf3 = distanceErrorFunction3 dm u l
  let newmole = updateCoordinates s2 atoms
  writeXYZ "test.xyz" "Comment" newmole
  --writeMoleculeXYZ "test.xyz" "Comment Line" atoms
  print newmole
  print errf1
  print errf2
  print errf3