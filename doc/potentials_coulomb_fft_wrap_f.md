# src/potentials/coulomb/fft_wrap.f

## Overview

This Fortran 77 (fixed-form) file, `fft_wrap.f`, serves as a wrapper layer for performing 3D Fast Fourier Transforms (FFTs). Its purpose is to provide a consistent API for FFT operations (`fft_setup`, `fft_forward`, `fft_back`) while allowing the underlying FFT implementation to be chosen at compile time via preprocessor directives. The code is stated to be adopted from the ORAC package.

It supports:
*   SGI's proprietary FFT library (`#ifdef SGIFFT`).
*   Cray's FFT library (`#ifdef CRAY`).
*   A "public domain" implementation (`pubz3d` and `pubz3di`) provided within this file, which constructs 3D FFTs by applying 1D FFT routines (from `fft3-public.f`) sequentially along each dimension.

The file also includes `get_fftdims` to help determine the necessary array dimensions and work space sizes based on the selected FFT backend. These routines are critical for methods like Particle Mesh Ewald (PME) that rely on efficient computation of discrete Fourier transforms.

## Key Components

### Subroutines

*   `get_fftdims(nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sizfftab, sizffwrk)`
    *   **Description**: Calculates array dimensions and work array sizes required for the FFT routines.
    *   `nfft1, nfft2, nfft3 :: INTEGER, intent(in)`: Dimensions of the 3D grid for FFT.
    *   `nfftdim1, nfftdim2, nfftdim3 :: INTEGER, intent(out)`: Actual first dimensions of arrays to be used (may include padding, e.g., `nfft1+1` for even `nfft1` in some RFFT libraries, though this implementation seems to use it for complex-to-complex).
    *   `nfftable, nffwork :: INTEGER, intent(out)`: Recommended sizes for `fftable` and `ffwork` arrays based on the FFT backend.
    *   `sizfftab, sizffwrk :: INTEGER, intent(out)`: Sizes used for allocation (can be same as `nfftable`/`nffwork` or different, e.g. `fftable` for `pubz3d` is `(ntable,3)`).
    *   The calculated sizes vary significantly depending on whether `SGIFFT`, `CRAY`, or the default public domain routines are used.

*   `fft_setup(array, fftable, ffwork, nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork)`
    *   **Description**: Initializes the FFT plan or pre-computes tables (twiddle factors, etc.) required by the chosen FFT library. This routine must be called before `fft_forward` or `fft_back`.
    *   `array :: REAL*8 (*), intent(inout)`: The data array (or a pointer to it). For Cray's `CCFFT3D`, this is used in initialization mode.
    *   `fftable :: REAL*8 (*), intent(out)`: Array to store FFT tables/plans.
    *   `ffwork :: REAL*8 (*), intent(out)`: Work array, may be used by some initialization routines.
    *   Other arguments define grid dimensions and table/work array sizes.
    *   **Backend calls**:
        *   `SGIFFT`: `ZFFT3DI(nfft1, nfft2, nfft3, fftable)`
        *   `CRAY`: `CCFFT3D(isign=0, ...)`
        *   Default: `pubz3di(nfft1, nfft2, nfft3, fftable, nfftable)`

*   `fft_forward(array, fftable, ffwork, ...)`
    *   **Description**: Performs a forward 3D FFT on the input `array`.
    *   `array :: COMPLEX*16 (*), intent(inout)`: The 3D data array (passed as 1D). On input, contains the data in spatial domain; on output, data in reciprocal (frequency) domain.
    *   `fftable :: REAL*8 (*), intent(in)`: Pre-initialized FFT tables.
    *   `ffwork :: COMPLEX*16 (*), intent(inout)`: Work array.
    *   **Backend calls**:
        *   `SGIFFT`: `ZFFT3D(isign=1, ...)`
        *   `CRAY`: `CCFFT3D(isign=1, ...)`
        *   Default: `pubz3d(isign=-1, ...)` (Note: `isign=-1` for forward in `pubz3d`).

*   `fft_back(array, fftable, ffwork, ...)`
    *   **Description**: Performs a backward (inverse) 3D FFT on the input `array`.
    *   Arguments are similar to `fft_forward`.
    *   **Backend calls**:
        *   `SGIFFT`: `ZFFT3D(isign=-1, ...)`
        *   `CRAY`: `CCFFT3D(isign=-1, ...)`
        *   Default: `pubz3d(isign=1, ...)` (Note: `isign=1` for backward in `pubz3d`).

*   `pubz3di(n1, n2, n3, table, ntable)`
    *   **Description**: Initializes the `table` array for the public domain 3D FFT. It does this by calling `cffti` (from `fft3-public.f`) for each of the three dimensions (n1, n2, n3). The `table` array is dimensioned `(ntable, 3)` to store these three sets of 1D FFT initialization data.
    *   `ntable` should be `4*max(n1,n2,n3) + 15`.

*   `pubz3d(isign, n1, n2, n3, w, ld1, ld2, table, ntable, work, nwork)`
    *   **Description**: Implements a 3D complex-to-complex FFT by performing a sequence of 1D FFTs.
        1.  For each plane in Z, and for each row in Y, it performs a 1D FFT along X on `w(i,j,k)`.
        2.  For each plane in Z, and for each column in X, it performs a 1D FFT along Y.
        3.  For each column in X, and for each row in Y, it performs a 1D FFT along Z.
        It uses `cfftf` if `isign == -1` (forward) and `cfftb` if `isign == 1` (backward), which are routines from `fft3-public.f`.
    *   `w(ld1, ld2, n3) :: COMPLEX*16, intent(inout)`: The 3D data array.
    *   `work(nwork) :: COMPLEX*16, intent(out)`: A 1D work array, size `max(n1,n2,n3)`.
    *   `table(ntable,3) :: REAL*8, intent(in)`: Initialization tables from `pubz3di`.

## Important Variables/Constants

*   **Data Types**: `REAL*8` for most table/work arrays, `COMPLEX*16` for the data being transformed by `fft_forward`/`fft_back` and the `work` array in `pubz3d`.
*   **Preprocessor Flags**: `SGIFFT`, `CRAY` control which library is used. If neither is set, the `pubz3d/i` routines (public domain fallback) are used.
*   **`isign` convention for `pubz3d`**: `isign = -1` for forward FFT, `isign = 1` for backward FFT. This is consistent with `cfftf` being forward and `cfftb` being backward if `cfftf` performs `exp(-j*...*)` and `cfftb` performs `exp(+j*...*)`.

## Usage Examples

These routines provide a common interface for FFTs. A typical workflow in a PME implementation would be:
1. Call `get_fftdims` to determine required array sizes.
2. Allocate `fftable` and `ffwork` arrays.
3. Call `fft_setup` to prepare `fftable`.
4. When needed, call `fft_forward` to transform data to reciprocal space.
5. Perform calculations in reciprocal space.
6. Call `fft_back` to transform results back to real space.

```fortran
      INTEGER N1, N2, N3, NFFDIM1, NFFDIM2, NFFDIM3
      INTEGER NFFTABLE_SIZE, NFFWORK_SIZE, SIZFFTAB, SIZFFWRK
      REAL*8, ALLOCATABLE :: FFTABLE(:)
      COMPLEX*16, ALLOCATABLE :: MY_3D_ARRAY(:,:,:), FFWORK(:)

      ! Define grid dimensions N1, N2, N3
      N1 = 32; N2 = 32; N3 = 32

      CALL GET_FFTDIMS(N1, N2, N3, NFFDIM1, NFFDIM2, NFFDIM3, &
     &                 NFFTABLE_SIZE, NFFWORK_SIZE, SIZFFTAB, SIZFFWRK)

      ALLOCATE(FFTABLE(SIZFFTAB))
      IF (SIZFFWRK .GT. 0) ALLOCATE(FFWORK(SIZFFWRK)) ! FFWORK might be zero size for some backends
      ALLOCATE(MY_3D_ARRAY(NFFDIM1, NFFDIM2, NFFDIM3))

      CALL FFT_SETUP(MY_3D_ARRAY, FFTABLE, FFWORK, N1, N2, N3, &
     &               NFFDIM1, NFFDIM2, NFFDIM3, SIZFFTAB, SIZFFWRK)

      ! ... Fill MY_3D_ARRAY with data ...

      CALL FFT_FORWARD(MY_3D_ARRAY, FFTABLE, FFWORK, N1, N2, N3, &
     &                 NFFDIM1, NFFDIM2, NFFDIM3, SIZFFTAB, SIZFFWRK)

      ! ... Process data in MY_3D_ARRAY (now in frequency space) ...

      CALL FFT_BACK(MY_3D_ARRAY, FFTABLE, FFWORK, N1, N2, N3, &
     &              NFFDIM1, NFFDIM2, NFFDIM3, SIZFFTAB, SIZFFWRK)

      ! ... MY_3D_ARRAY is back in real space (may need scaling) ...
```

## Dependencies and Interactions

*   **External FFT Libraries**: Conditionally depends on SGI's `ZFFT3DI`/`ZFFT3D` or Cray's `CCFFT3D` if the corresponding preprocessor flags are set.
*   **`fft3-public.f`**: If no external library is specified, the `pubz3di` and `pubz3d` routines are used, which in turn call `cffti`, `cfftf`, and `cfftb` from `fft3-public.f`.
*   This module acts as an abstraction layer, allowing the main Atomistica code (e.g., PME implementation) to use a consistent set of calls for 3D FFTs regardless of the available system-specific FFT libraries.
```
