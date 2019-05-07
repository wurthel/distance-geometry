{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Types where

import Control.Category
import Data.Label
import Prelude hiding ((.), id)

newtype ID =
  ID Int
  deriving (Num, Enum, Eq, Ord, Show)

newtype Element = 
  Element String
  deriving Show

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

data Atom =
  Atom
    { _aid :: ID
    , _aelem :: Element
    , _acoord :: Point
    }
  deriving (Show)

data Bond =
  Bond
    { _bfid :: ID
    , _bsid :: ID
    , _btype :: BondType
    , _bster :: BondSter
    , _btop :: BondTop
    }
  deriving (Show)

mkLabels [''Atom, ''Bond]
