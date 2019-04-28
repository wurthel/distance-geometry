{-# LANGUAGE BangPatterns     #-}

module DistanceGeometry 
    ( generateOfDistanceBoundsMatrix
    , triangleInequalitySmoothingFloyd)
    where

import Types
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data
import Numeric.LinearAlgebra.Devel
import Control.Category
import Data.Label
import Data.List (foldl')
import Prelude hiding ((.), id)

van_der_Waals_radii = 
            [0.001, 31, 28, 128, 96, 84, 73, 71, 66, 57, 58, 166, 141, 121, 111, 107, 105, 102, 106, 203, 176, 170, 160,
            153, 139, 139, 132, 126, 124, 132, 122, 122, 120, 119, 120, 120, 116, 220, 195, 190, 175, 164, 154, 147,
            146, 142, 139, 145, 144, 142, 139, 139, 138, 139, 140, 244, 215, 207, 204, 203, 201, 199, 198, 196, 194,
            192, 192, 189, 190, 187, 187, 175, 170, 162, 151, 144, 141, 136, 136, 132, 145, 146, 148, 140, 150, 150,
            260, 221, 215, 206, 200, 196, 190, 187, 180, 169]

-- | If we know the bound distance between atoms we set it
-- in upper and lower bounds matrix. If we assime that the atoms have 
-- van der Waals radii, then we can set all the other lower
-- bounds as sum of van der Waals radii (the default lower
-- bound between any two atoms is the sum of their
-- van der Waals radii). If the upper bounds are not
-- known then we shall enter a default value of 100.
generateOfDistanceBoundsMatrix :: [Atom] -> [Bond] -> Matrix Double
generateOfDistanceBoundsMatrix atoms bonds = 
    build (n, n) (\i' j' -> let (i, j) = (ID . fromEnum $ i' + 1, ID . fromEnum $ j' + 1)
                            in  if i == j 
                                then 0
                                else if isBonded i j bonds 
                                    then fixedDist i j 
                                    else if i < j 
                                        then upperDist i j
                                        else lowerDist i j)
    where n = length atoms
          upperDist (ID i) (ID j) = 100
          lowerDist (ID i) (ID j) = (vdmr atoms !! (i - 1)) + (vdmr atoms !! (j - 1))
          fixedDist (ID i) (ID j) = sqrt $ (xj - xi)^2 + (yj - yi)^2 + (zj - zi)^2
                where (xi, yi, zi) = get coordin (atoms !! (i - 1))
                      (xj, yj, zj) = get coordin (atoms !! (j - 1))
    
triangleInequalitySmoothingFloyd :: Matrix Double -> Maybe (Matrix Double)
triangleInequalitySmoothingFloyd matrix =
    do
    let n = rows matrix  
        upperBounds = build (n, n) (\i' j' -> let (i, j) = (fromEnum i', fromEnum j')
                                              in if i <= j 
                                                 then matrix `atIndex` (i, j) 
                                                 else matrix `atIndex` (j, i))
        lowerBounds = build (n, n) (\i' j' -> let (i, j) = (fromEnum i', fromEnum j')
                                              in if i <= j 
                                                 then matrix `atIndex` (j, i) 
                                                 else matrix `atIndex` (i, j))
        kij = [(k, i, j) | k <- [0 .. n-1], i <- [0 .. n-2], j <- [i+1 .. n-1]]
        smoothing Nothing _ = Nothing
        smoothing (Just (u0, l0)) (k, i, j) = 
            do
            let (u1, l1) = if (e1 > e2 + e3) then (changeMatrix u0 (i, j) (e2 + e3), l0) else (u0, l0)
                    where e1 = u0 `atIndex` (i,j)
                          e2 = u0 `atIndex` (i,k)
                          e3 = u0 `atIndex` (k,j)
                (u2, l2) = if (e1 < e2 - e3) then (u1, changeMatrix l1 (i, j) (e2 - e3)) else (u1, l1)
                    where e1 = l1 `atIndex` (i,j)
                          e2 = l1 `atIndex` (i,k)
                          e3 = u1 `atIndex` (k,j)
                (u3, l3) = if (e1 < e2 - e3) then (u2, changeMatrix l2 (i, j) (e2 - e3)) else (u2, l2)
                    where e1 = l2 `atIndex` (i,j)
                          e2 = l2 `atIndex` (j,k)
                          e3 = u2 `atIndex` (k,i)
            if (l3 ! i ! j > u3 !i !j) then fail "Erroneous Bounds" else Just (u3, l3)
    (upperBounds', lowerBounds') <- foldl' smoothing (Just (upperBounds, lowerBounds)) kij
    return $ build (n, n) (\i' j' -> let (i, j) = (fromEnum i', fromEnum j')
                                     in if i <= j 
                                        then upperBounds' `atIndex` (i, j)
                                        else lowerBounds' `atIndex` (i, j))

-- | UTILITS. НАЧАЛО
-- | Изменяет значение матрицы @m@ по индексам (@i@,@j@) на указанное @v@
changeMatrix :: Matrix Double -> (Int, Int) -> Double -> Matrix Double
changeMatrix m (i, j) v = runSTMatrix $ do
    m' <- thawMatrix m
    writeMatrix m' i j v
    return m'

-- | Массив Ван-дер-Ваальсовых радиусов
vdmr :: [Atom] -> [Double]
vdmr = map (getVDWR' . get name )
    where getVDWR' a = case a of "C" -> 0.500
                                 "H" -> 0.500
                                 "N" -> 0.500
                                 "O" -> 0.500
                                 "S" -> 0.500

-- | Определяет, связаны ли атомы
isBonded :: ID -> ID -> [Bond]-> Bool
isBonded n m s = or [(n, m) `elem` (bondsID s), (m, n) `elem` (bondsID s)]

-- | Возвращает массив пар связанных атомов, 
bondsID :: [Bond] -> [(ID, ID)]
bondsID = map (\x -> (get fid x, get sid x)) 
-- | UTILITS. КОНЕЦ