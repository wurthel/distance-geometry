module Main where

import Parser
import DistanceGeometry
import Data.Maybe

main :: IO ()
main = do
    let (atoms, bonds) = parserMolV2000 "example/example.mol"
        n = (triangleInequalitySmoothingFloyd . generateDistanceBoundsMatrix atoms) bonds
        z = fromJust n
    q <- randomDistanceMatrix z
    print $ distanceToMetricMatrix q
