import Data.Either
import DistanceGeometry
import IO
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data

main :: IO ()
main = do
  let idir = "test/molecules"
  let odir = "test/results"
  putStrLn ""
  putStrLn "test: Begin"
  test idir odir "example"
  test idir odir "arg"
  test idir odir "methan"
  test idir odir "retinal"
  putStrLn "test: End"

test idir odir mol = do
  -- Computing
  let (atoms, bonds) = readMolV2000 (idir ++ "/" ++ mol ++ ".mol")
  let s0@(u, l) =
        (triangleInequalitySmoothingFloyd . generateDistanceBoundsMatrix atoms) bonds
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

  -- Write result molecule in output files
  writeXYZ (odir ++ "/" ++ mol ++ ".xyz") "Comment" newmole

  -- Write logs in output files
  let logs = show newmole ++ "\n" ++
             "Errf1 = " ++ show errf1 ++ "\n" ++ 
             "Errf2 = " ++ show errf2 ++ "\n" ++ 
             "Errf3 = " ++ show errf3 ++ "\n"
  writeFile (odir ++ "/" ++ mol ++ ".log") logs

  -- Print to console
  putStrLn ("test: " ++ mol ++ ": "++ " OK!")