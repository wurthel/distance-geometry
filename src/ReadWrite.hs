module ReadWrite   
  ( -- * Read
    readMolV2000
    -- * Write
  , writeXYZ
  ) where

import Prelude hiding (readFile)

import Control.Lens
import Control.Monad.State
import Numeric.LinearAlgebra.Data
import System.Directory
import System.IO
import System.IO.Unsafe
import Text.Printf (hPrintf)

import Types

-- * Read
-- | Read molecule at *.pdb format (V2000)
readMolV2000 :: FilePath -> (Molecule, [Bond])
readMolV2000 inf =
  let txt = (lines . unsafePerformIO . readFile) inf
      counts_line = words (txt !! 3)
      count_atoms = read (counts_line !! 0)
      count_bonds = read (counts_line !! 1)
      atoms_lines = (take count_atoms . drop 4) txt
      bonds_lines = (take count_bonds . drop (4 + count_atoms)) txt
   in (foldr addAtom molecule atoms_lines, foldr addBond [] bonds_lines)
  where
    addAtom a = over atoms ((:) (readline' a))
    addBond b = (:) (readline'' b)
    readline' l =
      execState
        (do let w = words l
            acoordin . x .= read (w !! 0)
            acoordin . y .= read (w !! 1)
            acoordin . z .= read (w !! 2)
            aelement .= w !! 3
            avdwrad .= vdwr (w !! 3))
        atom
    readline'' l =
      execState
        (do let w = words l
            bfid .= read (w !! 0) - 1
            bsid .= read (w !! 1) - 1
            btype .= read (w !! 2)
            bster .= read (w !! 3)
            btop .= read (w !! 4))
        bond

-- | Get VDW radius
vdwr :: Element -> Double
vdwr a =
  case a of
    "H" -> 0.500 -- 0.5 -- 1.000
    "O" -> 0.650 -- 0.5 -- 1.300
    "N" -> 0.700 -- 0.5 -- 1.400
    "C" -> 0.750 -- 0.5 -- 1.500
    "S" -> 0.950 -- 0.5 -- 1.900
    "Br"-> 0.950
    othrewise -> error $ "vdwr not found for: " ++ show a

-- * Write
-- | Write molecule in *.xyz format
writeXYZ :: FilePath -> String -> Molecule -> IO ()
writeXYZ ouf comment molecule = do
  (tmp_name, tmp_handle) <- openTempFile "." "temp"
  hPrint tmp_handle (views atoms length molecule)
  hPutStrLn tmp_handle comment
  mapM_ (writeData tmp_handle) (view atoms molecule)
  hClose tmp_handle
  renameFile tmp_name ouf
  where
    writeData hdl atom = do
      let e = view aelement atom
          (Point x y z) = view acoordin atom
      hPrintf hdl "%s\t%8.6f\t%8.6f\t%8.6f\n" e x y z