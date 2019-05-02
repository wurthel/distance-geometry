module Main where

import Data.Maybe
import DistanceGeometry
import Parser

main :: IO ()
main = do
  let (atoms, bonds) = parserMolV2000 "example/retinal.mol"
  q <-
    (randomDistanceMatrix .
     fromJust .
     triangleInequalitySmoothingFloyd . generateDistanceBoundsMatrix atoms)
      bonds
  let p =
        (generateCoordinFromEigValAndVec .
         largestEigValAndVec . distanceToMetricMatrix)
          q
  print p