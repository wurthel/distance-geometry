module Main where

import Parser
import DistanceGeometry

main :: IO ()
main = do
    let (atoms, bonds) = parserMolV2000 "example/retinal.mol"
        n = (triangleInequalitySmoothingFloyd . generateOfDistanceBoundsMatrix atoms) bonds
    print n
