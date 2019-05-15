module Generators
    -- * Molecular Generators
  ( generateOneMolecule
  , generateManyMolecules
  ) where

import Control.Lens
import System.Directory
import DistanceGeometry
import IO
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data
import Types

generateOneMolecule ::
     FilePath -- ^ Input directory
  -> FilePath -- ^ Output dorectory
  -> FilePath -- ^ Molecule
  -> IO ()
generateOneMolecule idir odir mol = do
  let (atoms, bonds) = readMolV2000 (idir ++ "/" ++ mol ++ ".mol")
  let s0@(u, l) =
        (triangleInequalitySmoothingFloyd . generateDistanceBoundsMatrix atoms) bonds
  s1 <- randomDistanceMatrix s0
  let s2 =
        (generateCoordinFromEigValAndVec .
         largestEigValAndVec . distanceMatrixToMetricMatrix)
          s1
      newmole = updateCoordinates s2 atoms
  -- Write result molecule in output files
  writeXYZ (odir ++ "/" ++ mol ++ ".xyz") "Comment" newmole
  -- Write logs in output files
  let dm = coordMatrixToDistanceMatrix s2
      errf1 = distanceErrorFunction1 dm u l
      errf2 = distanceErrorFunction2 dm u l
      errf3 = distanceErrorFunction3 dm u l
  let logs = (show newmole ++ "\n") ++
          ("Errf1 = " ++ show errf1 ++ "\n") ++
          ("Errf2 = " ++ show errf2 ++ "\n") ++
          ("Errf3 = " ++ show errf3 ++ "\n")
  writeFile (odir ++ "/" ++ mol ++ ".log") logs
  -- Print to console
  putStrLn ("generateOneMolecule: " ++ mol ++ ": " ++ " OK!")

generateManyMolecules ::
     FilePath -- ^ Input directory
  -> FilePath -- ^ Output directory
  -> FilePath -- ^ Molecule
  -> Double -- ^ Error
  -> Int -- ^ The number of molecules generated
  -> IO ()
generateManyMolecules idir odir mol err n = do
  let odir' = odir ++ "/" ++ mol
  doesPathExist odir' >>= (`createDirectoryIfMissing` odir')
  step 1
  where
    step s
      | s > n = return ()
      | otherwise = do
        let (atoms, bonds) = readMolV2000 (idir ++ "/" ++ mol ++ ".mol")
        let s0@(u, l) =
              (triangleInequalitySmoothingFloyd . generateDistanceBoundsMatrix atoms)
                bonds
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
        let logs =
              (show newmole ++ "\n") ++
                ("Errf1 = " ++ show errf1 ++ "\n") ++
                ("Errf2 = " ++ show errf2 ++ "\n") ++
                ("Errf3 = " ++ show errf3 ++ "\n")
        if errf1 + errf2 + errf3 > err
          then step s
          else do
            -- Write result molecule in output files
            writeXYZ (odir ++ "/" ++ mol ++ "-" ++ show s ++ ".xyz") "Comment" newmole
            -- Write logs in output files
            writeFile (odir ++ "/" ++ mol ++ "-" ++ show s ++ ".log") logs
            -- Print to console
            putStrLn ("generateManyMolecules: " ++ mol ++ ": step " ++ show s ++ ": OK!")
            step (s+1)
