{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Types where
  
import Control.Lens

data Point = Point
  { _x :: Double
  , _y :: Double
  , _z :: Double
  } deriving (Show)

data Atom = Atom
  { _aserial :: Int
  , _aelement :: String
  , _acoord :: Point
  } deriving (Show)

data Bond = Bond
  { _bfid :: Int
  , _bsid :: Int
  , _btype :: Int
  , _bster :: Int
  , _btop :: Int
  } deriving (Show)

data Molecule = Molecule
  { _atoms :: [Atom]
  } deriving (Show)

molecule = Molecule { _atoms = []}

point = Point {_x = 1, _y = 0, _z = 0}

atom = Atom {_aserial = 0, _aelement = "", _acoord = point}

bond = Bond {_bfid = 0, _bsid = 0, _btype = 0, _bster = 0, _btop = 0}

makeLenses ''Point

makeLenses ''Atom

makeLenses ''Molecule

makeLenses ''Bond
