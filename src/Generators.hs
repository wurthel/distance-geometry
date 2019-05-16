module Generators
  ( 
  -- * Molecular Generators
    generateOneMolecule
  , generateManyMolecules
  ) where

import Control.Lens
import System.Directory
import DistanceGeometry
import ReadWrite
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
              (triangleSmooth . generateDistBoundsMatr atoms)
                bonds
        s1 <- randomDistMatr s0
        let s2 =
              (generateCoorFromEigValVec .
               largestEigValVec . distMatrToMetricMatr)
                s1
        let dm = coordMatrToDistMatr s2
            derr = distanceErrorFunction dm u l bonds
            cerr = 0
        if derr + cerr > err
          then step odir s
          else do
            -- Write result molecule in output files
            let newmole = updateCoord s2 atoms
            writeXYZ (odir ++ "/" ++ mol ++ "-" ++ show s ++ ".xyz") "Comment" newmole
            -- Write logs in output files
            let logs =
                 show newmole ++ "\n" ++
                 "Derr = " ++ show derr ++ "\n" ++
                 "Cerr = " ++ show cerr ++ "\n"
            writeFile (odir ++ "/" ++ mol ++ "-" ++ show s ++ ".log") logs
            -- Print to console
            putStrLn ("generateManyMolecules: " ++ mol ++ ": step " ++ show s ++ ": OK!")
            step odir (s+1)
