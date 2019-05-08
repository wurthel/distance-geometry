module IO
  ( readMoleculeMolV2000
  , writeMoleculeXYZ
  ) where

import Prelude hiding (readFile)

import Numeric.LinearAlgebra.Data
import System.Directory (renameFile)
import System.IO
import System.IO.Unsafe
import Text.Printf (hPrintf)
import Control.Lens
import Control.Monad.State

import Types

-- | Я рапсарсиваю файл неполность. 
-- Смотри реализации addAtoms и addBonds:
-- В .mol файле есть поля, предназачение которых мне неизвестно
readMoleculeMolV2000 :: FilePath -> (Molecule, [Bond])
readMoleculeMolV2000 inf =
  let txt = (lines . unsafePerformIO . readFile) inf
      counts_line = words $ txt !! 3
      count_atoms = read $ counts_line !! 0
      count_bonds = read $ counts_line !! 1
      atoms_line = (take count_atoms . drop 4) txt
      bonds_line = (take count_bonds . drop (4 + count_atoms)) txt
   in (foldr addAtom molecule atoms_line, addBond bonds_line)
    --addIDs 
  where
    addAtom a = over atoms ((:) (reada a))
    reada l = execState (do
      let w = words l
      acoord.x += read $ w !! 0
      scribe (acoord.y) (read $ w !! 1)
      scribe (acoord.z) (read $ w !! 2)
      scribe aelement (read $ w !! 3))
      atom

    addBond ::
         [String] -- ^ Lines containing information about bonds 
      -> [Bond] -- ^ Bonds
    addBond [] = []
    addBond (b:bs) =
      let w = words b
          i0 = read $ w !! 0
          i1 = read $ w !! 1
          i2 = read $ w !! 2
          i3 = read $ w !! 3
          i4 = read $ w !! 4
          --bond = Bond i0 i1 i2 i3 i4
       in bond : addBond bs

updateMoleculeXYZ :: [Atom] -> Matrix Double -> [Atom]
updateMoleculeXYZ atomx matr = undefined

-- writeMoleculeXYZ :: FilePath -> String -> [Atom] -> IO ()
-- writeMoleculeXYZ ouf comment atoms = do
--   (tmp_name, tmp_handle) <- openTempFile "." "temp"
--   (hPrint tmp_handle . length) atoms
--   hPutStrLn tmp_handle comment
--   mapM_ (writeData tmp_handle) atoms
--   hClose tmp_handle
--   renameFile tmp_name ouf
--   where
--     writeData hdl atom = do
--       let (e) = get aelem atom
--           (x, y, z) = get acoord atom
--       hPrintf hdl "%s\t%8.6f\t%8.6f\t%8.6f\n" e x y z