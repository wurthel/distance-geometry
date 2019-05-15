# distance-geometry

Installation with Stack
-----------------------
Stackage is a stable package archive. Stackage builds are supposed to
be reproducible. Stackage also provides Long Term Support releases.
To build `distance-geometry` with Stackage dependencies, use the `stack` tool:

  * install [`stack`](https://docs.haskellstack.org/)
  * if necessary, install GHC: run `stack setup`
  * run: `stack update`
  * in the project source directory run: `stack build`

### Build Status

[![Build Status](https://travis-ci.org/wurthel/distance-geometry.svg?branch=master)](https://travis-ci.org/wurthel/distance-geometry)

An example of test
----------
To run test: `stack test`

See [test/Spec.hs](test/Spec.hs)
```haskell
idir = "test/molecules"
odir = "test/results"

main :: IO ()
main = do
  putStrLn ""
  putStrLn "generateOneMolecule: Begin"
  generateOneMolecule idir odir "example"
  putStrLn "generateOneMolecule: End"
  putStrLn ""
  putStrLn "generateManyMolecules: Begin"
  generateManyMolecules idir odir "C2Br4" 0.01 100
  putStrLn "generateManyMolecules: End"
```
