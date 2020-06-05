{-# language BangPatterns        #-}
{-# language NoImplicitPrelude   #-}
{-# language ScopedTypeVariables #-}

-- | This module exposes functions for performing
--   Fast Fourier Transform (FFT) and Inverse Fast Fourier Transform (IFFT)
--   over 'Contiguous' data structures.
module Data.Primitive.Contiguous.FFT
  ( fft
  , ifft
  , mfft
  ) where

import qualified Prelude

import Control.Applicative (pure)
import Control.Monad (when,unless)
import Control.Monad.Primitive (PrimMonad(..))
import Control.Monad.ST (runST)
import Data.Bits (shiftR,shiftL,(.&.))
import Data.Bool (Bool,(&&))
import Data.Complex (Complex(..),conjugate)
import Data.Eq (Eq(..))
import Data.Function (($))
import Data.Ord (Ord(..))
import Data.Primitive.Contiguous (Contiguous,Element,Mutable)
import Data.Semiring (negate,(+),(*),(-))
import GHC.Exts
import GHC.Real ((/))
import qualified Data.Primitive.Contiguous as Contiguous

{-# RULES
"fft/ifft" forall x. fft (ifft x) = x
"ifft/fft" forall x. ifft (fft x) = x
  #-}

-- | Radix-2 decimation-in-time fast Fourier Transform.
--   The given array must have a length that is a power of two.
fft :: forall arr. (Contiguous arr, Element arr (Complex Double))
  => arr (Complex Double)
  -> arr (Complex Double)
{-# inlinable [1] fft #-}
fft arr = if arrOK arr
  then runST $ do
    marr <- copyWhole arr
    mfft marr
    Contiguous.unsafeFreeze marr
  else Prelude.error "Data.Primitive.Contiguous.FFT.fft: bad array length"

-- | Inverse fast Fourier transform.
ifft :: forall arr. (Contiguous arr, Element arr (Complex Double))
  => arr (Complex Double)
  -> arr (Complex Double)
{-# inlinable [1] ifft #-}
ifft arr = if arrOK arr
  then runST $ do
    marr <- copyWhole arr
    mifft marr
    Contiguous.unsafeFreeze marr
  else Prelude.error "Data.Primitive.Contiguous.FFT.ifft: bad vector length"

copyWhole :: forall arr m a. (PrimMonad m, Contiguous arr, Element arr a)
  => arr a
  -> m (Mutable arr (PrimState m) a)
{-# inline copyWhole #-}
copyWhole arr = do
  let len = Contiguous.size arr
  marr <- Contiguous.new len
  Contiguous.copy marr 0 arr 0 len
  pure marr

arrOK :: forall arr a. (Contiguous arr, Element arr a)
  => arr a
  -> Bool
{-# inline arrOK #-}
arrOK arr = (n > 0) && (n .&. n - 1 == 0)
  where n = Contiguous.size arr

-- | Radix-2 decimation-in-time fast Fourier Transform.
--   The given array must have a length that is a power of two,
--   though this property is not checked.
mfft :: forall arr m. (PrimMonad m, Contiguous arr, Element arr (Complex Double))
  => Mutable arr (PrimState m) (Complex Double)
  -> m ()
mfft mut = do
  len <- Contiguous.sizeMutable mut
  let
    bitReverse !i !j = if i == len - 1
      then stage 0 1
      else do
        when (i < j) $ swap mut i j
        let inner k l = if k <= l
              then inner (k `shiftR` 1) (l - k)
              else bitReverse (i + 1) (l + k)
        inner (len `shiftR` 1) j
    stage l l1 = unless (1 `shiftL` l == len) $ do
      let
        !l2 = l1 `shiftL` 1
        !e = (negate twoPi) / (intToDouble l2)
        flight j !a = if j == l1
          then stage (l + 1) l2
          else do
            let
              butterfly i = if i >= len
                then flight (j + 1) (a + e)
                else do
                  let i1 = i + l1
                  xi1 :+ yi1 <- Contiguous.read mut i1
                  let !c = Prelude.cos a
                      !s' = Prelude.sin a
                      d = (c*xi1 - s'*yi1) :+ (s'*xi1 + c*yi1)
                  ci <- Contiguous.read mut i
                  Contiguous.write mut i1 (ci - d)
                  Contiguous.write mut i (ci + d)
                  butterfly (i + l2)
            butterfly j
      flight 0 0
  bitReverse 0 0

-- | Radix-2 decimation-in-time inverse fast Fourier Transform.
--   The given array must have a length that is a power of two,
--   though this property is not checked.
mifft :: forall arr m. (PrimMonad m, Contiguous arr, Element arr (Complex Double))
  => Mutable arr (PrimState m) (Complex Double)
  -> m ()
mifft mut = do
  lenComplex <- Prelude.fmap intToComplexDouble $ Contiguous.sizeMutable mut
  Contiguous.mapMutable conjugate mut
  mfft mut
  Contiguous.mapMutable (\x -> conjugate x / lenComplex) mut

swap :: forall arr m x. (PrimMonad m, Contiguous arr, Element arr x)
  => Mutable arr (PrimState m) x
  -> Int
  -> Int
  -> m ()
{-# inline swap #-}
swap mut i j = do
  atI <- Contiguous.read mut i
  atJ <- Contiguous.read mut j
  Contiguous.write mut i atJ
  Contiguous.write mut j atI

twoPi :: Double
{-# inline twoPi #-}
twoPi = 6.283185307179586

intToDouble :: Int -> Double
{-# inline intToDouble #-}
intToDouble = Prelude.fromIntegral

intToComplexDouble :: Int -> Complex Double
{-# inline intToComplexDouble #-}
intToComplexDouble = Prelude.fromIntegral
