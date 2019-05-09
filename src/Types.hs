{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Types where

import Control.Lens

type Serial = Int

type Element = String

data Point =
  Point
    { _x :: Double
    , _y :: Double
    , _z :: Double
    }
  deriving (Show)

data Atom =
  Atom
    { _aelement :: Element
    , _acoordin :: Point
    , _avdwrad :: Double
    }
  deriving (Show)

data Bond =
  Bond
    { _bfid :: Serial
    , _bsid :: Serial
    , _btype :: Int
    , _bster :: Int
    , _btop :: Int
    }
  deriving (Show)

data Molecule =
  Molecule
    { _atoms :: [Atom]
    }
  deriving (Show)

molecule = Molecule {_atoms = []}

point = Point {_x = 0, _y = 0, _z = 0}

atom = Atom {_aelement = "", _acoordin = point, _avdwrad = 0}

bond = Bond {_bfid = 0, _bsid = 0, _btype = 0, _bster = 0, _btop = 0}

makeLenses ''Point

makeLenses ''Atom

makeLenses ''Molecule

makeLenses ''Bond
