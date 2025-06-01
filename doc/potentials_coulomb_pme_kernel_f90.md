# src/potentials/coulomb/pme_kernel.f90

## Overview

The `pme_kernel` module provides the core computational routines and data structures for the smooth Particle-Mesh Ewald (PME) method. It is responsible for managing the PME grid, performing charge assignment using B-splines, carrying out Fast Fourier Transforms (FFTs), calculating the reciprocal space energy and field contributions, and interpolating these fields back to particle positions. This module is a crucial backend for the higher-level `pme.f90` module. The code is stated to be adopted from the ORAC package.

It can be compiled to use the FFTW3 library (`#ifdef HAVE_FFTW3`) for improved performance, or it falls back to using FFT routines provided by `fft_wrap.f` (which itself can use `fft3-public.f`).

## Key Components

### Modules

*   `pme_kernel`
    *   **Uses**: `supplib` (for `DP`, `PI`, logging, error handling, etc.). Conditionally uses `omp_lib` (for OpenMP support with FFTW3) and FFTW3 specific modules (`iso_c_binding`, `fftw3.f03`).
    *   **Provides**: The `pme_grid_t` derived data type and key subroutines for PME calculations.

### Data Types

*   `pme_grid_t` (Public)
    *   **Description**: Stores all data associated with the PME grid and calculation setup.
    *   **Fields**:
        *   `numatoms :: integer`: Number of atoms.
        *   `order :: integer`: Order of B-spline interpolation.
        *   `nfft1, nfft2, nfft3 :: integer`: Dimensions of the FFT grid.
        *   `nfftdim1, nfftdim2, nfftdim3 :: integer`: Actual dimensions of the complex grid array `Q` (may differ from `nfft` for some FFT libraries or real-to-complex transforms, though PME typically uses complex-to-complex for `Q`).
        *   `ntheta :: integer`: Size of `theta` arrays (`numatoms * order`).
        *   `bsp_mod1(:), bsp_mod2(:), bsp_mod3(:) :: real(DP), allocatable`: Stores the squared modulus of the Fourier transform of the B-splines along each grid dimension. Used in the reciprocal space sum.
        *   `theta1(:), theta2(:), theta3(:) :: real(DP), allocatable`: Stores B-spline coefficients for each particle for charge spreading. Dimension `(order, numatoms)`.
        *   `dtheta1(:), dtheta2(:), dtheta3(:) :: real(DP), allocatable`: Stores derivatives of B-spline coefficients for force calculation. Dimension `(order, numatoms)`.
        *   `fr1(:), fr2(:), fr3(:) :: real(DP), allocatable`: Scaled fractional coordinates of particles on the grid. Dimension `(numatoms)`.
        *   **FFT Data**:
            *   If `HAVE_FFTW3`:
                *   `fftw_is_initialized :: logical`: Flag for FFTW plan status.
                *   `plan_forward, plan_backward :: type(C_PTR)`: FFTW plans for forward and backward transforms.
                *   `Q_ptr :: type(C_PTR)`: C pointer to the allocated memory for the charge grid.
                *   `Q(:,:,:) :: complex(C_DOUBLE_COMPLEX), pointer`: Fortran pointer to `Q_ptr`, representing the 3D charge/potential grid.
            *   Else (fallback FFT):
                *   `nfftable, nffwork, sizfftable, sizffwork :: integer`: Sizes for FFT work/table arrays.
                *   `fftable(:), ffwork(:) :: real(DP), allocatable`: Work/table arrays for fallback FFT (`fft_wrap.f`).
                *   `Q(:,:,:) :: complex(DP), allocatable`: 3D charge/potential grid.

### Public Subroutines & Interfaces

*   `init(this, grid, numatoms, order, error)` (maps to `pme_grid_init`)
    *   **Description**: Initializes the `pme_grid_t` object. Allocates all necessary arrays based on grid dimensions, number of atoms, and spline order. It pre-computes B-spline moduli by calling `load_bsp_moduli`. It also prepares FFT plans if using FFTW3 or initializes tables for the fallback FFT routines via `fft_setup` (from `fft_wrap.f`).
*   `del(this)` (maps to `pme_grid_del`)
    *   **Description**: Destructor. Deallocates all allocatable arrays within `pme_grid_t` and destroys FFTW plans if they were created.
*   `potential_and_field(this, x,y,z, charge, recip, volume, ewald_coeff, eer, virial, phi, dx,dy,dz)` (maps to `pme_grid_potential_and_field`)
    *   **Description**: This is the main computational routine for the reciprocal space part of PME.
        1.  `get_scaled_fractionals`: Converts Cartesian coordinates `x,y,z` to scaled fractional coordinates `fr1,fr2,fr3` relative to the PME grid dimensions.
        2.  `get_bspline_coeffs`: Calculates B-spline values (`theta1`, etc.) and their derivatives (`dtheta1`, etc.) for each particle at its scaled fractional coordinates.
        3.  `fill_charge_grid`: Spreads particle charges onto the complex grid `this%Q` using the B-spline coefficients.
        4.  **Backward FFT**: Performs a backward FFT on `this%Q` (transforming charge density to electrostatic potential in k-space). Uses FFTW or fallback `fft_back`.
        5.  `energy_and_virial_sum`: Calculates the reciprocal space energy (`eer`) and virial contribution (`virial`) by summing over the grid points in k-space. It applies the PME structure factor (related to `ewald_coeff` and `bsp_moduli`) to the transformed charges. The grid `this%Q` is modified in place to store `StructureFactor * Q_k`.
        6.  **Forward FFT**: Performs a forward FFT on the modified `this%Q` (transforming potential/field in k-space back to real space).
        7.  **Interpolation**:
            *   If `dx,dy,dz` are present, calls `potential_and_field_sum` to interpolate the potential `phi` and electric field components `dx,dy,dz` at particle positions from the real-space grid `this%Q` using B-spline values and their derivatives.
            *   Else, calls `potential_sum` to interpolate only the potential `phi`.
    *   **Arguments**: `this (pme_grid_t)`, particle coordinates `x,y,z`, `charge`, reciprocal cell vectors `recip`, cell `volume`, `ewald_coeff` (Ewald splitting parameter \(\alpha\)), outputs `eer` (reciprocal energy), `virial`, and updates `phi` and optionally `dx,dy,dz` (electric field components).

### Private Helper Subroutines

*   `get_scaled_fractionals`: Converts Cartesian to scaled fractional coordinates for grid mapping.
*   `load_bsp_moduli`: Calls `fill_bspline` then `dftmod` to compute and store Fourier transforms of B-spline basis functions.
*   `dftmod`: Computes the modulus of the Discrete Fourier Transform (DFT) of the B-spline array.
*   `fill_charge_grid`: Implements the charge spreading (assignment) from particle positions to the grid `Q` using B-spline interpolation.
*   `get_bspline_coeffs`: Calls `fill_bspline` for each particle to get B-spline values and derivatives.
*   `fill_bspline`: Core routine for B-spline calculation, uses `bsp_init`, `bsp_one_pass`, `bsp_diff`.
*   `bsp_init, bsp_one_pass, bsp_diff`: Recursive B-spline basis function evaluation and differentiation.
*   `energy_and_virial_sum`: Performs the k-space sum for energy and virial. Iterates over grid points `k1,k2,k3`, calculates \(m^2\) (squared magnitude of reciprocal lattice vector), the structure factor `exp(-pi^2*m^2/alpha^2) / (pi*V*m^2 * |BS(m)|^2)`, and sums up contributions. Modifies `Q(k)` by multiplying with the structure factor term relevant for potential/field.
*   `potential_and_field_sum`: Interpolates potential and electric field from the grid `Q` back to particle positions using B-spline coefficients and their derivatives.
*   `potential_sum`: Similar to above, but only for the potential.

## Important Variables/Constants
*   `order`: Order of B-spline interpolation (e.g., 4 for cubic, 6 for quintic).
*   `nfft1, nfft2, nfft3`: Dimensions of the FFT grid.
*   `ewald_coeff`: The Ewald splitting parameter \(\alpha\).
*   `bsp_mod1, bsp_mod2, bsp_mod3`: Moduli of the DFT of the B-splines, used in the denominator of the reciprocal space sum.
*   `Q`: The 3D complex array representing charges/potentials/fields on the grid.

## Usage Examples
This module's routines are primarily called by the `pme.f90` module. `pme_grid_init` is called once, `pme_grid_del` at the end, and `pme_grid_potential_and_field` at each step where energy/forces/potential are needed.

## Dependencies and Interactions
*   **`supplib`**: For basic utilities and constants.
*   **FFT Implementation**:
    *   If `HAVE_FFTW3` is defined: Uses the FFTW3 library for FFT operations. OpenMP support can be enabled for FFTW3.
    *   Else: Uses routines from `fft_wrap.f` (i.e., `fft_setup`, `fft_back`, `fft_forward`), which in turn may use `fft3-public.f` for the actual 1D FFT calculations if no other system FFT library like SGIFFT or CRAYFFT is selected within `fft_wrap.f`.
*   The module performs the computationally intensive parts of the PME algorithm, particularly the grid assignments, FFTs, and k-space summation.
*   It is designed to be used by `pme.f90`, which handles the real-space part of the Ewald sum and combines all contributions.
```
