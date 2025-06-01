# src/potentials/coulomb/coulomb_short_gamma.f90

## Overview

The `coulomb_short_gamma` module provides functions for calculating specific terms related to the short-range gamma function, which is a component of the third-order Density Functional Tight Binding (DFTB3) method. These terms are essential for computing the DFTB3 Hamiltonian and associated forces, particularly the contributions arising from fluctuations in Mulliken charges. The formulation is based on the work by Gaus et al., J. Chem. Theory Comput. 7, 931 (2001).

This module works in conjunction with `damp_short_gamma`, which provides damping functions (`hij`) and their derivatives, optionally controlled by a `zeta` parameter. The functions herein often distinguish between cases where Hubbard-like parameters `U_i` and `U_j` for two interacting atoms are similar or different.

## Key Components

### Modules

*   `coulomb_short_gamma`
    *   **Uses**: `supplib`, `damp_short_gamma`
    *   **Provides**: Functions for calculating "capital short gamma" terms and their derivatives, which are essentially weighted derivatives of the underlying short-range gamma interaction.

### Public Functions

*   `capital_short_gamma(abs_rij, dU_i, U_i, U_j, zeta)`
    *   **Description**: Calculates a term contributing to the energy or forces, proportional to the derivative of a Hubbard-like parameter `U_i`. It's computed as `part_deriv_sgamma_wrt_Ui * dU_i`.
    *   **Arguments**:
        *   `abs_rij :: real(DP), intent(in)`: The absolute distance between two atoms i and j.
        *   `dU_i :: real(DP), intent(in)`: The derivative of the Hubbard-like parameter `U_i` (e.g., with respect to an atomic coordinate or external field).
        *   `U_i, U_j :: real(DP), intent(in)`: Hubbard-like parameters for atoms i and j.
        *   `zeta :: real(DP), intent(in), optional`: An optional parameter, likely for controlling damping via the `hij` function from `damp_short_gamma`.
    *   **Returns**: `real(DP)` - The calculated "capital short gamma" value.

*   `derivative_capital_short_gamma(abs_rij, dU_i, U_i, U_j, zeta)`
    *   **Description**: Calculates the derivative of the `capital_short_gamma` term with respect to interatomic distance `abs_rij`. This is computed as `second_part_deriv_csgamma_wrt_Ui_and_r * dU_i`.
    *   **Arguments**: Same as `capital_short_gamma`.
    *   **Returns**: `real(DP)` - The calculated derivative.

*   `Sfij(abs_rij, U_i, U_j)`
    *   **Description**: An intermediate function used when `U_i` and `U_j` are different. It calculates `exp(-U_i*abs_rij)*fij(...) + exp(-U_j*abs_rij)*fij(...)`. While not explicitly private, its primary use is within `short_gamma`.
    *   **Arguments**: `abs_rij`, `U_i`, `U_j`.
    *   **Returns**: `real(DP)`.

*   `Sgij(abs_rij, U_i)`
    *   **Description**: An intermediate function used when `U_i` and `U_j` are very similar. It calculates `exp(-U_i*abs_rij)*gij(...)`. While not explicitly private, its primary use is within `short_gamma`.
    *   **Arguments**: `abs_rij`, `U_i`.
    *   **Returns**: `real(DP)`.

### Core Internal Functions (Used by Public Interface)

*   `short_gamma(abs_rij, U_i, U_j, zeta)`: The fundamental short-range gamma interaction.
    *   If `abs(U_i - U_j)` is small, it uses `-Sgij(...)`.
    *   Otherwise, it uses `-Sfij(...)`.
    *   If `zeta` is present, the result is multiplied by `hij(abs_rij, U_i, U_j, zeta)` (from `damp_short_gamma`).
*   `part_deriv_sgamma_wrt_Ui(abs_rij, U_i, U_j, zeta)`: Calculates the partial derivative of `short_gamma` with respect to `U_i`. Includes a factor of `3.20_DP`.
*   `part_deriv_sgamma_wrt_r(abs_rij, U_i, U_j, zeta)`: Calculates the partial derivative of `short_gamma` with respect to `abs_rij`.
*   `second_part_deriv_csgamma_wrt_Ui_and_r(abs_rij, U_i, U_j, zeta)`: Calculates the mixed second partial derivative of `short_gamma` with respect to `U_i` and `abs_rij`.
*   `fij(abs_rij, U_i, U_j)` and `gij(abs_rij, U_i)`: These are the base mathematical functions whose forms depend on `U_i`, `U_j`, and `abs_rij`. `gij` is a polynomial in `U_i*abs_rij` divided by `abs_rij`. `fij` has a more complex rational form.
*   A suite of partial derivative functions for `Sfij`, `Sgij`, `fij`, and `gij` with respect to `U_i`, `U_j`, and `abs_rij` (e.g., `part_deriv_Sfij_wrt_Ui`, `second_part_deriv_gij_wrt_Ui_and_r`).

## Important Variables/Constants

*   `U_i, U_j`: Hubbard-like parameters, central to DFTB methods.
*   `dU_i`: Derivative of `U_i`, indicating these terms are likely part of a force calculation or response property.
*   `abs_rij`: Interatomic distance.
*   `zeta`: Optional parameter for damping.
*   `3.20_DP`: A numerical constant appearing in some derivative calculations.

## Usage Examples

These functions are specialized mathematical components for DFTB3 calculations. They would be called by higher-level routines that assemble the DFTB3 energy expression or its derivatives for forces.

```fortran
! Conceptual usage within a DFTB3 code:
! USE coulomb_short_gamma
! REAL(DP) :: distance, U_atom_i, U_atom_j, deriv_U_atom_i, zeta_val, gamma_term, force_contrib
!
! ! ... obtain distance, U values, derivative of U, zeta ...
!
! gamma_term = capital_short_gamma(distance, deriv_U_atom_i, U_atom_i, U_atom_j, zeta_val)
! force_contrib = derivative_capital_short_gamma(distance, deriv_U_atom_i, U_atom_i, U_atom_j, zeta_val)
!
! ! ... use these terms in total energy / force accumulation ...
```

## Dependencies and Interactions

*   **`supplib`**: Likely used for general utility functions or constants (e.g., `DP`).
*   **`damp_short_gamma`**: This is a critical dependency. The `coulomb_short_gamma` module uses:
    *   `hij(abs_rij, U_i, U_j, zeta)`: A damping function.
    *   `part_deriv_hij_wrt_Ui(abs_rij, U_i, U_j, zeta)`
    *   `part_deriv_hij_wrt_r(abs_rij, U_i, U_j, zeta)`
    *   `second_part_deriv_hij_wrt_Ui_and_r(abs_rij, U_i, U_j, zeta)`
    These functions from `damp_short_gamma` are combined with the terms calculated in the present module, especially when `zeta` is provided.
*   The internal logic frequently branches based on the presence of the `zeta` parameter and whether `U_i` and `U_j` are numerically close, selecting different functional forms (`fij` vs. `gij`) accordingly.
```
