{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Types where

import Control.Category
import Data.Label
import Prelude hiding ((.), id)

newtype ID =
  ID Int
  deriving (Num, Enum, Eq, Ord, Show)

newtype BondType =
  BondType Int
  deriving (Num, Enum, Eq, Ord, Show)

newtype BondSter =
  BondSter Int
  deriving (Num, Enum, Eq, Ord, Show)

newtype BondTop =
  BondTop Int
  deriving (Num, Enum, Eq, Ord, Show)

type Point = (Double, Double, Double)

type Name = String

data Atom =
  Atom
    { _aid :: ID
    , _coordin :: Point
    , _name :: Name
    }
  deriving (Show)

data Bond =
  Bond
    { _fid :: ID
    , _sid :: ID
    , _btype :: BondType
    , _bster :: BondSter
    , _btop :: BondTop
    }
  deriving (Show)

mkLabels [''Atom, ''Bond]
