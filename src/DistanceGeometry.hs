module DistanceGeometry
  ( 
  -- * Main functions
    generateDistanceBoundsMatrix
  , triangleInequalitySmoothingFloyd
  , randomDistanceMatrix
  , distanceMatrixToMetricMatrix
  , largestEigValAndVec
  , generateCoordinFromEigValAndVec
  , coordMatrixToDistanceMatrix

  -- * Error Functions
  , distanceErrorFunction1
  , distanceErrorFunction2
  , distanceErrorFunction3
  ) where

import Control.Category
import Data.List (foldl', zipWith3)
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data
import Numeric.LinearAlgebra.Devel
import System.Random
import Types

-- * Main functions
-- | Generate of a distance bounds matrix
-- If we know the bound distance between atoms we set it
-- in upper and lower bounds matrix. If we assime that the atoms have 
-- van der Waals radii, then we can set all the other lower
-- bounds as sum of van der Waals radii (the default lower
-- bound between any two atoms is the sum of their
-- van der Waals radii). If the upper bounds are not
-- known then we shall enter a default value of 100.
generateDistanceBoundsMatrix :: [Atom] -> [Bond] -> (Matrix Double, Matrix Double)
generateDistanceBoundsMatrix atoms bonds =
  let n = length atoms
      upperDist (ID i) (ID j) = 100
      lowerDist (ID i) (ID j) = (vdmr atoms !! (i - 1)) + (vdmr atoms !! (j - 1))
      fixedDist (ID i) (ID j) = sqrt $ (xj - xi) ^ 2 + (yj - yi) ^ 2 + (zj - zi) ^ 2
        where
          (xi, yi, zi) = get acoord (atoms !! (i - 1))
          (xj, yj, zj) = get acoord (atoms !! (j - 1))
      upper =
        build
          (n, n)
          (\i' j' ->
             let (i, j) = (ID . fromEnum $ i' + 1, ID . fromEnum $ j' + 1)
              in if i == j
                   then 0
                   else if isBonded i j bonds
                          then fixedDist i j
                          else upperDist i j)
      lower =
        build
          (n, n)
          (\i' j' ->
             let (i, j) = (ID . fromEnum $ i' + 1, ID . fromEnum $ j' + 1)
              in if i == j
                   then 0
                   else if isBonded i j bonds
                          then fixedDist i j
                          else lowerDist i j)
   in (upper, lower)

-- | Triangle Inequality Bounds Smoothing
-- triangleInequalitySmoothingFloyd :: Matrix Double -> Either String (Matrix Double)
triangleInequalitySmoothingFloyd ::
     (Matrix Double, Matrix Double) -> (Matrix Double, Matrix Double)
triangleInequalitySmoothingFloyd (upper, lower) =
  let n = rows upper
      kij = [(k, i, j) | k <- [0 .. n - 1], i <- [0 .. n - 2], j <- [i + 1 .. n - 1]]
      smoothing (u0, l0) (k, i, j) = do
        let (u1, l1) =
              if e1 > e2 + e3
                then (changeMatrix u0 (i, j) (e2 + e3), l0)
                else (u0, l0)
              where
                e1 = u0 `atIndex` (i, j)
                e2 = u0 `atIndex` (i, k)
                e3 = u0 `atIndex` (k, j)
            (u2, l2) =
              if e1 < e2 - e3
                then (u1, changeMatrix l1 (i, j) (e2 - e3))
                else (u1, l1)
              where
                e1 = l1 `atIndex` (i, j)
                e2 = l1 `atIndex` (i, k)
                e3 = u1 `atIndex` (k, j)
            (u3, l3) =
              if e1 < e2 - e3
                then (u2, changeMatrix l2 (i, j) (e2 - e3))
                else (u2, l2)
              where
                e1 = l2 `atIndex` (i, j)
                e2 = l2 `atIndex` (j, k)
                e3 = u2 `atIndex` (k, i)
        if l3 ! i ! j > u3 ! i ! j
          then error "Erroneous Bounds"
          else (u3, l3)
          --then Left "Erroneous Bounds"
          --else Right (u3, l3)
      (upper', lower') = foldl' smoothing (upper, lower) kij
      upper'' =
        build
          (n, n)
          (\i' j' ->
             let (i, j) = (fromEnum i', fromEnum j')
              in if i < j
                   then upper' `atIndex` (i, j)
                   else upper' `atIndex` (j, i))
      lower'' =
        build
          (n, n)
          (\i' j' ->
             let (i, j) = (fromEnum i', fromEnum j')
              in if i < j
                   then lower' `atIndex` (i, j)
                   else lower' `atIndex` (j, i))
   in (upper'', lower'')

-- | Generation of a distance matrix by random selection 
-- of distances between the bounds.
randomDistanceMatrix :: (Matrix Double, Matrix Double) -> IO (Matrix Double)
randomDistanceMatrix (upper, lower) = do
  let n = rows upper
  distanceVector <-
    sequence
      [ r
      | i <- [0 .. n - 1]
      , j <- [0 .. n - 1]
      , let a = lower `atIndex` (i, j)
            b = upper `atIndex` (i, j)
            r =
              if j <= i
                then return 0
                else randomRIO (a, b)
      ]
  let distanceMatrix = matrix n distanceVector
  return $ distanceMatrix `add` (tr' distanceMatrix)

-- | Improving Random Sampling: Metrization
metrization :: Matrix Double -> Matrix Double
metrization matrix = undefined

-- | Improving Random Sampling: Partial Metrization
partialMetrization :: Matrix Double -> Matrix Double
partialMetrization matrix = undefined

-- | Conversion of the distance matrix to a metric matrix
distanceMatrixToMetricMatrix :: Matrix Double -> Matrix Double
distanceMatrixToMetricMatrix matr =
  let n = rows matr
      m = fromIntegral n
      uTriangleMatr =
        build
          (n, n)
          (\i' j' ->
             let (i, j) = (fromEnum i', fromEnum j')
              in if j <= i
                   then 0
                   else matr `atIndex` (i, j))
      d0 i =
        (1 / m) * (sumElements $ (matr ! i) ^ 2) -
        (1 / m ^ 2) * (sumElements $ uTriangleMatr ^ 2)
   in build
        (n, n)
        (\i' j' ->
           let (i, j) = (fromEnum i', fromEnum j')
               a = d0 i
               b = d0 j
               c = matr `atIndex` (i, j)
            in (a + b - c ^ 2) / 2)

-- | distanceToMetricMatrix is symmetric matrix => eigenvalues is real
-- Function return @k@ pairs of eigenvalues and eigenvectors (as columns)    
-- in descending order
largestEigValAndVec :: Matrix Double -> (Vector Double, Matrix Double)
largestEigValAndVec matr =
  if (not . isSymmetric) matr
    then error "The matrix is not symmetric"
    else (subVector 0 k eigval, takeColumns k eigvec)
  where
    (eigval, eigvec) = (eigSH . trustSym) matr
    k = min 3 (size eigval)

-- | Generation of three-dimensional coordinates from 
-- eigenvalues and eigenvectors.
generateCoordinFromEigValAndVec :: (Vector Double, Matrix Double) -> Matrix Double
generateCoordinFromEigValAndVec (val, vec) = vec Numeric.LinearAlgebra.<> diag val

-- | Convert coordinate matrix to 
-- distance matrix.
coordMatrixToDistanceMatrix :: Matrix Double -> Matrix Double
coordMatrixToDistanceMatrix matr =
  let n = rows matr
      dist i j = sqrt $ (xj - xi) ^ 2 + (yj - yi) ^ 2 + (zj - zi) ^ 2
        where
          [xi, yi, zi] = toList $ matr ! i
          [xj, yj, zj] = toList $ matr ! j
      dmatr =
        build
          (n, n)
          (\i' j' ->
             let (i, j) = (fromEnum i', fromEnum j')
              in if j <= i
                   then 0
                   else dist i j)
   in dmatr `add` (tr' dmatr)

-- * Error Functions
-- | Distance Error Function, type 1
distanceErrorFunction1 :: Matrix Double -> Matrix Double -> Matrix Double -> Double
distanceErrorFunction1 dist upper lower =
  let n = rows dist
      d2 = dist ^ 2
      u2 = upper ^ 2
      l2 = lower ^ 2
      ij = [(i, j) | i <- [0 .. n - 2], j <- [i + 1 .. n - 1], i /= j]
      f (i, j) =
        (+) $
        max 0 (d2 ! i ! j - u2 ! i ! j) ^ 2 + max 0 (l2 ! i ! j ^ 2 - d2 ! i ! j) ^ 2
   in foldr f 0 ij

-- | Distance Error Function, type 3
distanceErrorFunction2 :: Matrix Double -> Matrix Double -> Matrix Double -> Double
distanceErrorFunction2 dist upper lower =
  let n = rows dist
      d2 = dist ^ 2
      u2 = upper ^ 2
      l2 = lower ^ 2
      ij = [(i, j) | i <- [0 .. n - 2], j <- [i + 1 .. n - 1], i /= j]
      f (i, j) =
        (+) $
        max 0 (d2 ! i ! j / u2 ! i ! j - 1) ^ 2 + max 0 (l2 ! i ! j / d2 ! i ! j - 1) ^ 2
   in foldr f 0 ij

-- | Distance Error Function, type 3
distanceErrorFunction3 :: Matrix Double -> Matrix Double -> Matrix Double -> Double
distanceErrorFunction3 dist upper lower =
  let n = rows dist
      d2 = dist ^ 2
      u2 = upper ^ 2
      l2 = lower ^ 2
      ij = [(i, j) | i <- [0 .. n - 2], j <- [i + 1 .. n - 1], i /= j]
      f (i, j) =
        (+) $
        max 0 (d2 ! i ! j / u2 ! i ! j - 1) ^ 2 +
        max 0 (2 * l2 ! i ! j / (l2 ! i ! j + d2 ! i ! j) - 1) ^ 2
   in foldr f 0 ij

-- * Utils.
-- | Изменяет значение матрицы @m@ по индексам (@i@,@j@) на указанное @v@
changeMatrix :: Matrix Double -> (Int, Int) -> Double -> Matrix Double
changeMatrix m (i, j) v =
  runSTMatrix $ do
    m' <- thawMatrix m
    writeMatrix m' i j v
    return m'

isSymmetric :: Matrix Double -> Bool
isSymmetric matr =
  let r = rows matr
      c = cols matr
      s = [matr ! i ! j == matr ! j ! i | i <- [0 .. r - 1], j <- [i .. c - 1]]
   in r == c && and s

-- | Массив Ван-дер-Ваальсовых радиусов
vdmr :: [Atom] -> [Double]
vdmr = map (getVDWR' . get aelem)
  where
    getVDWR' (Element a) =
      case a of
        "H" -> 0.5 -- 1.000
        "O" -> 0.5 -- 1.300
        "N" -> 0.5 -- 1.400
        "C" -> 0.5 -- 1.500
        "S" -> 0.5 -- 1.900
        othrewise -> error $ "vdmr not found for: " ++ show a

-- | Определяет, связаны ли атомы
isBonded :: ID -> ID -> [Bond] -> Bool
isBonded n m s = (n, m) `elem` bondsID s || (m, n) `elem` bondsID s

-- | Возвращает массив пар связанных атомов, 
bondsID :: [Bond] -> [(ID, ID)]
bondsID = map (\x -> (get bfid x, get bsid x))