# src/potentials/coulomb/pme.f90

## Overview

The `pme` module implements the smooth Particle-Mesh Ewald (PME) method, a widely used and efficient algorithm for calculating long-range electrostatic interactions in periodic systems. The PME method splits the Coulomb sum into two parts:
1.  A short-range direct sum performed in real space, typically using a complementary error function (`erfc`) to screen interactions.
2.  A long-range sum performed in reciprocal (Fourier) space, calculated efficiently using Fast Fourier Transforms (FFTs).
A self-energy correction term is also included.

This implementation references the foundational PME papers by Darden et al. (J. Chem. Phys. 98, 10089 (1993)) and Essmann et al. (J. Chem. Phys. 103, 8577 (1995)). The code is stated to be adopted from the ORAC package. It relies on a `pme_kernel` module for core PME grid operations and reciprocal space calculations.

## Key Components

### Modules

*   `pme`
    *   **Uses**: `supplib`, `particles`, `neighbors`, `pme_kernel`. Conditionally uses `mpi`.
    *   **Provides**: The `pme_t` type and the standard set of potential interface routines for PME calculations.

### Data Types

*   `pme_t` (Public)
    *   **Description**: Holds parameters and state for the PME calculation.
    *   **Fields**:
        *   `alpha :: real(DP)`: The Ewald splitting parameter, determining the range separation between direct and reciprocal sums. Default: `0.4_DP`. (Note: This is often automatically determined from `cutoff` in `bind_to`).
        *   `cutoff :: real(DP)`: Cutoff radius for the real-space (direct sum) part of the calculation. Default: `10.0_DP`.
        *   `grid(3) :: integer`: Dimensions of the FFT grid used for the reciprocal space calculation. Default: `(/ 64, 64, 64 /)`.
        *   `order :: integer`: Order of the B-spline interpolation used for spreading charges onto the grid and interpolating forces/potentials from the grid. Default: `8`.
        *   `pme_grid :: type(pme_grid_t)`: An object (from `pme_kernel`) that stores data related to the PME grid, FFT plans, and precomputed factors for the reciprocal sum.
        *   `sqrt_alpha :: real(DP)`: Square root of `alpha`.
        *   `sqrt_alpha_pi :: real(DP)`: `sqrt(alpha / PI)`.
        *   `cutoff_sq :: real(DP)`: Square of `cutoff`.

### Public Subroutines & Interfaces

*   `init(this, cutoff, order, grid, error)`: Constructor. Sets user-specified `cutoff`, `order`, and `grid` dimensions.
*   `del(this)`: Destructor. Calls `del(this%pme_grid)` to deallocate resources held by the PME grid object in the `pme_kernel`.
*   `bind_to(this, p, nl, ierror)`: Binds the PME potential to a particle system `p` and neighbor list handler `nl`.
    1.  Checks that the system `p` is 3D periodic, as PME is designed for such systems.
    2.  Calculates an optimal `this%alpha` based on `this%cutoff` using the formula: `alpha = log(10.0d0)*12 / cutoff_sq`. This aims to balance the computational load between real and reciprocal space.
    3.  Computes derived constants `sqrt_alpha`, `sqrt_alpha_pi`.
    4.  Calls `request_interaction_range(nl, this%cutoff)` to set the real-space cutoff for neighbor finding.
    5.  Initializes the `this%pme_grid` object by calling `init(this%pme_grid, this%grid, p%nat, this%order)` from the `pme_kernel` module.
*   `potential(this, p, nl, q, phi, ierror)`: Calculates the electrostatic potential `phi` at each atom.
    1.  **Direct Space**: Computes the short-range contribution to `phi` using `erfc(sqrt_alpha * rij) / rij` for pairs within `cutoff`. Includes a self-term `-2*q(i)*sqrt_alpha_pi`. OpenMP parallelized.
    2.  **Reciprocal Space**: Calls `potential_and_field` from `pme_kernel` (passing `this%pme_grid`, particle positions, charges, cell information, and `sqrt_alpha`) to get the long-range contribution to `phi`.
*   `energy_and_forces(this, p, nl, q, epot, f, wpot, error)`: Calculates the total PME energy, forces, and virial.
    1.  **Direct Space**: Computes short-range energy (`epot_dir`), forces, and virial (`wpot_dir`) using `erfc` damping. OpenMP parallelized.
    2.  **Reciprocal Space**: Calls `potential_and_field` from `pme_kernel` to get reciprocal space energy (`epot_rec`), forces (via electric fields `Ex, Ey, Ez`), and virial (`wpot_rec`).
    3.  **Self-Energy Correction**: Calculates `epot_self = sum(q_i^2) * sqrt_alpha_pi`.
    4.  **Total**: `epot = epot_dir + epot_rec - epot_self`. Forces and virials are summed.
*   `register(this, cfg, m)`: Registers `alpha`, `cutoff`, `grid`, and `order` parameters with `ptrdict`.

### Private Parameters
*   `EPS = 1d-10`: A small constant used for floating-point comparisons (e.g., to check if a charge `q(i)` is non-zero).

## Important Variables/Constants
*   `alpha`: Ewald splitting parameter.
*   `cutoff`: Real-space cutoff distance.
*   `grid`: FFT grid dimensions (e.g., `64x64x64`).
*   `order`: Interpolation order for charge spreading and force/potential evaluation on the grid.

## Usage Examples
PME is a standard method for full electrostatic calculations in periodic molecular simulations.

```fortran
! Conceptual usage:
! USE pme_module
! TYPE(pme_t) :: pme_calculator
!
! CALL pme_calculator%init(cutoff=12.0_DP, order=6, grid=(/48,48,48/))
! ! ... setup particles (p_data), neighbor_list_handler (nl_data), charges (q_array) ...
! CALL pme_calculator%bind_to(p_data, nl_data)
!
! ! Calculate energy, forces, virial
! CALL pme_calculator%energy_and_forces(p_data, nl_data, q_array, E_total, F_array, V_tensor)
!
! CALL pme_calculator%del()
```

## Dependencies and Interactions
*   **`supplib`**: For utilities and constants.
*   **`particles`, `neighbors`**: For core simulation data structures.
*   **`mpi`**: Conditionally used, likely within the `pme_kernel` for distributed FFTs if Atomistica is compiled with MPI support.
*   **`pme_kernel`**: This is a crucial dependency. The `pme` module delegates all reciprocal space calculations (including FFTs, grid operations, and solving Poisson's equation in reciprocal space) to routines within `pme_kernel`, particularly `potential_and_field`. The `pme_grid_t` type, which encapsulates the grid and FFT plans, is also defined in `pme_kernel`.
*   The `pme_kernel` likely uses the FFT routines from `fft_wrap.f` and `fft3-public.f`.
*   The PME method correctly handles periodic boundary conditions for Coulomb interactions, offering O(N log N) scaling, which is much more efficient than direct summation for large periodic systems.
```
