import Data.Either
import DistanceGeometry
import IO
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data
import System.Directory

idir = "test/molecules"
odir = "test/results"

main :: IO ()
main = do
  putStrLn ""
  putStrLn "generateOneMolecule: Begin"
  generateOneMolecule idir odir "example"
  generateOneMolecule idir odir "methan"
  putStrLn "generateOneMolecule: End"
  putStrLn ""
  putStrLn "generateManyMolecules: Begin"
  generateManyMolecules idir odir "example" 60
  generateManyMolecules idir odir "methan" 100
  putStrLn "generateManyMolecules: End"

generateOneMolecule idir odir mol
  -- Computing
 = do
  let (atoms, bonds) = readMolV2000 (idir ++ "/" ++ mol ++ ".mol")
  let s0@(u, l) = (triangleInequalitySmoothingFloyd . generateDistanceBoundsMatrix atoms) bonds
  s1 <- randomDistanceMatrix s0
  let s2 = (generateCoordinFromEigValAndVec . largestEigValAndVec . distanceMatrixToMetricMatrix) s1
  let dm = coordMatrixToDistanceMatrix s2
      errf1 = distanceErrorFunction1 dm u l
      errf2 = distanceErrorFunction2 dm u l
      errf3 = distanceErrorFunction3 dm u l
  let newmole = updateCoordinates s2 atoms
  -- Write result molecule in output files
  writeXYZ (odir ++ "/" ++ mol ++ ".xyz") "Comment" newmole
  -- Write logs in output files
  let logs =
        show newmole ++
        "\n" ++
        "Errf1 = " ++
        show errf1 ++ "\n" ++ "Errf2 = " ++ show errf2 ++ "\n" ++ "Errf3 = " ++ show errf3 ++ "\n"
  writeFile (odir ++ "/" ++ mol ++ ".log") logs
  -- Print to console
  putStrLn ("generateOneMolecule: " ++ mol ++ ": " ++ " OK!")

generateManyMolecules idir odir mol n = do
  let odir' = odir ++ "/" ++ mol
  doesPathExist odir' >>= (`createDirectoryIfMissing` odir')
  mapM_ (step idir odir' mol) [1 .. n]
  where
    step idir odir mol n
    -- Computing
     = do
      let (atoms, bonds) = readMolV2000 (idir ++ "/" ++ mol ++ ".mol")
      let s0@(u, l) = (triangleInequalitySmoothingFloyd . generateDistanceBoundsMatrix atoms) bonds
      s1 <- randomDistanceMatrix s0
      let s2 =
            (generateCoordinFromEigValAndVec . largestEigValAndVec . distanceMatrixToMetricMatrix)
              s1
      let dm = coordMatrixToDistanceMatrix s2
          errf1 = distanceErrorFunction1 dm u l
          errf2 = distanceErrorFunction2 dm u l
          errf3 = distanceErrorFunction3 dm u l
      let newmole = updateCoordinates s2 atoms
    -- Write result molecule in output files
      writeXYZ (odir ++ "/" ++ mol ++ "-" ++ show n ++ ".xyz") "Comment" newmole
    -- Write logs in output files
      let logs =
            show newmole ++
            "\n" ++
            "Errf1 = " ++
            show errf1 ++
            "\n" ++ "Errf2 = " ++ show errf2 ++ "\n" ++ "Errf3 = " ++ show errf3 ++ "\n"
      writeFile (odir ++ "/" ++ mol ++ "-" ++ show n ++ ".log") logs
    -- Print to console
      putStrLn ("generateManyMolecules: " ++ mol ++ ": step " ++ show n ++ ": OK!")