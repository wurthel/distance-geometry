module Main where

import Parser
import DistanceGeometry
import Data.Maybe

main :: IO ()
main = do
    let (atoms, bonds) = parserMolV2000 "example/retinal.mol"
    q <- (randomDistanceMatrix . fromJust .  triangleInequalitySmoothingFloyd . generateDistanceBoundsMatrix atoms) bonds
    let p = (generateCoordinFromEigValAndVec . largestEigValAndVec . distanceToMetricMatrix) q
    print p