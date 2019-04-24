module DistanceGeometry where

import Types
import Numeric.LinearAlgebra
import Control.Category
import Data.Label
import Prelude hiding ((.), id)

generateOfDistanceBoundsMatrix :: [Atom] -> Matrix Double
generateOfDistanceBoundsMatrix atoms = 
    let n = length atoms
        distance i j =
            let (i', j') = (fromEnum i, fromEnum j)
                (xi, yi, zi) = get coordin (atoms !! i')
                (xj, yj, zj) = get coordin (atoms !! j')
            in sqrt $ (xj - xi)^2 + (yj - yi)^2 + (zj - zi)^2
    in build (n, n) (\i j -> if i >= j then distance i j else distance j i)     

triangleInequalitySmoothingFloyd :: [Atom] -> [Bond] -> [[Double]]
triangleInequalitySmoothingFloyd atoms bonds = undefined