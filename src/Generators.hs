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
  -> Double -- ^ Error
  -> IO ()
generateOneMolecule idir odir mol err = generateManyMolecules idir odir mol err 1

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
  step odir' 1
  where
    step odir s
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
            derr = distanceErrorFunction dm u l bonds
            cerr = 0
        let newmole = updateCoordinates s2 atoms
        let logs =
              (show newmole ++ "\n") ++
                ("Derr = " ++ show derr ++ "\n") ++
                ("Cerr = " ++ show cerr ++ "\n")
        if derr + cerr > err
          then step odir s
          else do
            print dm
            print u
            print l
            -- Write result molecule in output files
            writeXYZ (odir ++ "/" ++ mol ++ "-" ++ show s ++ ".xyz") "Comment" newmole
            -- Write logs in output files
            writeFile (odir ++ "/" ++ mol ++ "-" ++ show s ++ ".log") logs
            -- Print to console
            putStrLn ("generateManyMolecules: " ++ mol ++ ": step " ++ show s ++ ": OK!")
            step odir (s+1)
