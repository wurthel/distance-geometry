module Parser
  ( parserMolV2000
  ) where

import Prelude hiding (readFile)

import System.IO.Strict
import System.IO.Unsafe
import Types

-- | Я рапсарсиваю файл неполность. 
-- Смотри реализации addAtoms и addBonds:
-- В .mol файле есть поля, предназачение которых мне неизвестно
parserMolV2000 :: FilePath -> ([Atom], [Bond])
parserMolV2000 inf =
  let txt = lines . unsafePerformIO . readFile $ inf
      counts_line = words $ txt !! 3
      count_atoms = read $ counts_line !! 0
      count_bonds = read $ counts_line !! 1
      atoms_line = take count_atoms . drop 4 $ txt
      bonds_line = take count_bonds . drop (4 + count_atoms) $ txt
   in (addAtom atoms_line 1, addBond bonds_line)
  where
    addAtom :: [String] -> Int -> [Atom]
    addAtom [] _ = []
    addAtom (l:ls) n =
      let w = words l
          i = ID n
          x = read $ w !! 0
          y = read $ w !! 1
          z = read $ w !! 2
          c = (x, y, z)
          s = Element $ w !! 3
          atom = Atom i s c
       in atom : addAtom ls (n + 1)
    addBond :: [String] -> [Bond]
    addBond [] = []
    addBond (b:bs) =
      let w = words b
          i0 = ID . read $ w !! 0
          i1 = ID . read $ w !! 1
          i2 = BondType . read $ w !! 2
          i3 = BondSter . read $ w !! 3
          i4 = BondTop . read $ w !! 4
          bond = Bond i0 i1 i2 i3 i4
       in bond : addBond bs
