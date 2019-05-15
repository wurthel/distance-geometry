module DistanceGeometry
    -- * Steps of algorithm
  ( generateDistanceBoundsMatrix
  , triangleInequalitySmoothingFloyd
  , randomDistanceMatrix
  , distanceMatrixToMetricMatrix
  , largestEigValAndVec
  , generateCoordinFromEigValAndVec
  , coordMatrixToDistanceMatrix
  , updateCoordinates
  , isBonded
  -- * Error Functions
  , distanceErrorFunction
  , chiralErrorFunction
  ) where

import Control.Lens
import Control.Monad.State
import Data.List (foldl', zipWith3)
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data
import Numeric.LinearAlgebra.Devel
import System.Random
import Types

-- * Steps of algorithm
-- | Generate of a distance bounds matrix
generateDistanceBoundsMatrix ::
     Molecule -> [Bond] -> (Matrix Double, Matrix Double)
generateDistanceBoundsMatrix molecule bonds =
  let n = views atoms length molecule
      upperDist _ _ = 100
      lowerDist i j =
        getAtom i molecule ^. avdwrad + getAtom j molecule ^. avdwrad
      fixedDist i j = sqrt $ (xj - xi) ^ 2 + (yj - yi) ^ 2 + (zj - zi) ^ 2
        where
          (Point xi yi zi) = getAtom i molecule ^. acoordin
          (Point xj yj zj) = getAtom j molecule ^. acoordin
      upper =
        build
          (n, n)
          (\i' j' ->
             let (i, j) = (fromEnum i', fromEnum j')
              in if i == j
                   then 0
                   else if isBonded i j bonds
                          then fixedDist i j
                          else upperDist i j)
      lower =
        build
          (n, n)
          (\i' j' ->
             let (i, j) = (fromEnum i', fromEnum j')
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
      kij =
        [ (k, i, j)
        | k <- [0 .. n - 1]
        , i <- [0 .. n - 2]
        , j <- [i + 1 .. n - 1]
        ]
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
      (upper', lower') = foldl' smoothing (upper, lower) kij
      newUpper =
        build
          (n, n)
          (\i' j' ->
             let (i, j) = (fromEnum i', fromEnum j')
              in if i < j
                   then upper' `atIndex` (i, j)
                   else upper' `atIndex` (j, i))
      newLower =
        build
          (n, n)
          (\i' j' ->
             let (i, j) = (fromEnum i', fromEnum j')
              in if i < j
                   then lower' `atIndex` (i, j)
                   else lower' `atIndex` (j, i))
   in (newUpper, newLower)

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
        (1 / m) * (sumElements $ (matr ! i) ^ 2) - (1 / m ^ 2) * 
        (sumElements $ uTriangleMatr ^ 2)
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
generateCoordinFromEigValAndVec ::
     (Vector Double, Matrix Double) -> Matrix Double
generateCoordinFromEigValAndVec (val, vec) =
  vec Numeric.LinearAlgebra.<> (sqrt . diag) val

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

-- | Update coordinates
updateCoordinates :: Matrix Double -> Molecule -> Molecule
updateCoordinates coord mol = set atoms newatoms mol
  where
    newatoms = zipWith updatec (toRows coord) (view atoms mol)
    updatec coord' =
      execState
        (do acoordin . x .= coord' ! 0
            acoordin . y .= coord' ! 1
            acoordin . z .= coord' ! 2)

-- * Error Functions
-- | Distance error function.
distanceErrorFunction ::
     Matrix Double -> Matrix Double -> Matrix Double -> [Bond] -> Double
distanceErrorFunction dist upper lower bonds =
  let n = rows dist
      [d2, u2, l2] = map (^2) [dist, upper, lower]
      ij' = [(i, j)| i <- [0 .. n - 2], j <- [i + 1 .. n - 1]]
      ij = filter (\(i,j) -> not $ isBonded i j bonds) ij'
      f1 (i, j) =
        (+) $ max 0 (d2 ! i ! j - u2 ! i ! j) ^ 2 +
        max 0 (l2 ! i ! j ^ 2 - d2 ! i ! j) ^ 2
      f2 (i, j) =
        (+) $ max 0 (d2 ! i ! j / u2 ! i ! j - 1) ^ 2 +
        max 0 (l2 ! i ! j / d2 ! i ! j - 1) ^ 2
      f3 (i, j) =
        (+) $ max 0 (d2 ! i ! j / u2 ! i ! j - 1) ^ 2 +
        max 0 (2 * l2 ! i ! j / (l2 ! i ! j + d2 ! i ! j) - 1) ^ 2
   in foldr f1 0 ij + foldr f2 0 ij + foldr f3 0 ij
-- | Chiral error function.
-- 
chiralErrorFunction :: Matrix Double -> Matrix Double -> Matrix Double -> Double
chiralErrorFunction coord upper lower = undefined

-- * Utils.
-- | Изменяет значение матрицы @m@ по индексам (@i@,@j@) на указанное @v@
changeMatrix :: Matrix Double -> (Int, Int) -> Double -> Matrix Double
changeMatrix m (i, j) v =
  runSTMatrix $ do
    m' <- thawMatrix m
    writeMatrix m' i j v
    return m'

-- | Check the symmetric matrix
isSymmetric :: Matrix Double -> Bool
isSymmetric matr =
  let r = rows matr
      c = cols matr
      s = [matr ! i ! j == matr ! j ! i | i <- [0 .. r - 1], j <- [i .. c - 1]]
   in r == c && and s

-- | Определяет, связаны ли атомы
isBonded :: Serial -> Serial -> [Bond] -> Bool
isBonded n m s = (n, m) `elem` bonds s || (m, n) `elem` bonds s
  where bonds = map (\x -> (view bfid x, view bsid x))

-- | Get atom from molecule
getAtom :: Serial -> Molecule -> Atom
getAtom i = views atoms (!! i)