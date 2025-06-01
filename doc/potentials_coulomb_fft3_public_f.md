# src/potentials/coulomb/fft3-public.f

## Overview

This file, `fft3-public.f`, provides a collection of Fortran 77 (fixed-form) subroutines for performing Fast Fourier Transforms (FFTs). The code is stated to be an adaptation of routines from ORAC, which in turn is a `REAL*8` (double precision) version of the well-known FFTPACK library from Netlib, developed by Paul N. Swartztrauber at NCAR (National Center for Atmospheric Research).

While the filename might suggest 3D FFT capabilities, the core routines provided (`CFFTF`, `CFFTB`, `CFFTI`) are designed for 1D complex-to-complex transforms. A comment at the end of the file indicates that 3D FFTs would be achieved by applying these 1D transforms sequentially along each dimension of a 3D data grid. This is a common "brute force" approach for constructing multi-dimensional FFTs from 1D routines.

These FFT routines are fundamental for advanced electrostatic calculations, particularly for methods like Particle Mesh Ewald (PME), where solving Poisson's equation in reciprocal space is a key step.

## Key Components

The file consists of several groups of subroutines:

### Main User-Callable FFT Routines

These are the primary routines a user would typically interact with:

*   `CFFTI(N, WSAVE)`
    *   **Description**: Initializes the work array `WSAVE` for a given transform length `N`. This routine must be called before `CFFTF` or `CFFTB` can be used for that length. It computes the prime factorization of `N` and pre-calculates necessary trigonometric factors (twiddle factors).
    *   **Arguments**:
        *   `N :: INTEGER, intent(in)`: The length of the complex data sequence to be transformed.
        *   `WSAVE :: REAL*8 array, intent(out)`: A work and save array. Its required size depends on `N` (at least `4*N + 15` based on typical FFTPACK requirements, though the exact size partitioning is internal).

*   `CFFTF(N, C, WSAVE)`
    *   **Description**: Computes the forward complex discrete Fourier Transform (DFT) of a sequence `C`.
    *   **Arguments**:
        *   `N :: INTEGER, intent(in)`: Length of the transform.
        *   `C :: REAL*8 array, intent(inout)`: On input, the complex sequence to be transformed. The real and imaginary parts are stored contiguously (e.g., `C(1)=Re0, C(2)=Im0, C(3)=Re1, C(4)=Im1, ...`). On output, `C` contains the transformed sequence.
        *   `WSAVE :: REAL*8 array, intent(in)`: The work array initialized by `CFFTI`.

*   `CFFTB(N, C, WSAVE)`
    *   **Description**: Computes the backward complex discrete Fourier Transform (DFT) of a sequence `C`. This is typically the inverse transform.
    *   **Arguments**: Same as `CFFTF`.

### Core 1D FFT Worker Routines

These routines are called by `CFFTF`, `CFFTB`, and `CFFTI` and contain the main logic:

*   `CFFTI1(N, WA, IFAC)`: Called by `CFFTI`. `WA` and `IFAC` are parts of `WSAVE`. `IFAC` stores the prime factors of `N`, and `WA` stores the twiddle factors.
*   `CFFTF1(N, C, CH, WA, IFAC)`: Called by `CFFTF`. Performs the actual forward transform logic using various `PASSFx` routines. `CH` is a work array.
*   `CFFTB1(N, C, CH, WA, IFAC)`: Called by `CFFTB`. Performs the actual backward transform logic using various `PASSBx` routines.

### `PASSFx` and `PASSBx` Subroutines (Computational Kernels)

These are the core computational stages of the FFT algorithm, performing passes for different prime factors of `N`.
*   `PASSF(NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)`: Generic forward pass routine.
*   `PASSF2, PASSF3, PASSF4, PASSF5`: Specialized and optimized forward pass routines for factors 2, 3, 4, and 5.
*   `PASSB(NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)`: Generic backward pass routine.
*   `PASSB2, PASSB3, PASSB4, PASSB5`: Specialized and optimized backward pass routines for factors 2, 3, 4, and 5.

The arguments to these `PASSx` routines involve various indexing and length parameters (`IDO`, `IP`, `L1`, `IDL1`), data arrays (`CC`, `C1`, `CH`, `CH2`), and twiddle factor arrays (`WA`).

## Important Variables/Constants

*   **Implicit Typing**: The code uses `implicit REAL*8 (a-h,o-z)`, meaning variables starting with these letters are double precision real numbers by default.
*   `NTRYH(4) = (/3,4,2,5/)`: An array in `CFFTI1` defining the preferred order of prime factors (3, 4, 2, 5) for decomposing `N`. This is typical for FFTPACK to optimize for these common factors.
*   `TPI = 6.28318530717959d0`: Represents 2*pi.
*   Trigonometric constants in `PASSx3` and `PASSx5` routines (e.g., `TAUR, TAUI` in `PASSF3/B3` are related to `cos(2*pi/3)` and `sin(2*pi/3)`).

## Usage Examples

These subroutines form a library. A typical usage pattern would be:

```fortran
      INTEGER N, IER
      REAL*8 C(2*N), WSAVE(4*N+15) ! Ensure WSAVE is sufficiently large

      ! 1. Initialize WSAVE for a given transform length N
      CALL CFFTI(N, WSAVE)

      ! 2. Populate C with complex data (Re, Im, Re, Im, ...)
      ! ...

      ! 3. Perform forward FFT
      CALL CFFTF(N, C, WSAVE)

      ! C now contains the transformed data

      ! 4. Perform backward FFT (e.g., after some operation in frequency space)
      CALL CFFTB(N, C, WSAVE)

      ! C now contains the inverse transformed data (may need scaling by 1/N)
```
To perform a 3D FFT on a 3D array `data(NX,NY,NZ)`, one would typically apply these 1D routines:
1.  Transform along NX for each (y,z) pair.
2.  Transform the result along NY for each (x,z) pair.
3.  Transform the result along NZ for each (x,y) pair.

## Dependencies and Interactions

*   The routines are largely self-contained but rely on the correct initialization of the `WSAVE` array by `CFFTI` before `CFFTF` or `CFFTB` are called.
*   They are fundamental building blocks for more complex algorithms that require spectral methods, such as the Particle Mesh Ewald (PME) summation technique for calculating long-range Coulomb interactions in periodic systems.
*   The comment "3D (slow) Fourier Transform ... brute force approach" suggests that while these 1D routines can be used for 3D transforms, more optimized 3D-specific FFT libraries might offer better performance for that specific task.
```
