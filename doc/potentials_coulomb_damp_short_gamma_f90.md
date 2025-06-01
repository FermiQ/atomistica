# src/potentials/coulomb/damp_short_gamma.f90

## Overview

The `damp_short_gamma` module provides a damping function, denoted `hij`, and its various partial derivatives. This damping function is specifically designed for use in the third-order Density Functional Tight Binding (DFTB3) method, where it modulates the short-range gamma function (representing electron-electron interaction terms). The formulation and application are based on the work by Gaus et al., J. Chem. Theory Comput. 7, 931 (2011).

The functions within this module take input parameters `tau_i` and `tau_j` (related to atomic Hubbard U-like parameters), an interatomic distance `rij`, and a damping exponent `zeta`. Internally, `tau` values are scaled to effective `U` values, and distances are converted to Bohr radii for calculations.

## Key Components

### Modules

*   `damp_short_gamma`
    *   **Uses**: `supplib` (primarily for the `DP` kind parameter and the `Bohr` constant).
    *   **Provides**: The damping function `hij` and its partial derivatives with respect to `U_i` (derived from `tau_i`) and `rij`.

### Public Functions

All functions take `rij`, `tau_i`, `tau_j`, and `zeta` as input. Internally, they convert:
*   `U_k = (5.0_DP/16.0_DP) * tau_k`
*   `abs_rij_internal = rij / Bohr`
*   `U_k_internal = U_k * Bohr` (Effective U parameters in atomic units of energy, assuming tau is also in energy units that need conversion or U is dimensionless before this scaling)

*   `hij(rij, tau_i, tau_j, zeta) result(res)`
    *   **Description**: Calculates the value of the damping function `h(rij, Ui, Uj)`.
    *   **Arguments**:
        *   `rij :: real(DP), intent(in)`: The interatomic distance.
        *   `tau_i, tau_j :: real(DP), intent(in)`: Input parameters related to atomic hardness or Hubbard U values for atoms i and j.
        *   `zeta :: real(DP), intent(in)`: The damping exponent.
    *   **Returns**: `res :: real(DP)` - The value of the damping function.
    *   **Form**: `res = exp(-((0.5*(U_i_internal + U_j_internal))^zeta) * abs_rij_internal^2)`.

*   `part_deriv_hij_wrt_Ui(rij, tau_i, tau_j, zeta) result(res)`
    *   **Description**: Calculates the partial derivative of the `hij` function with respect to the internal Hubbard-like parameter `U_i` (derived from `tau_i`). The final result is scaled by `Bohr` to reflect the derivative with respect to `tau_i` if `tau_i` were in energy units and `U_i_internal` is effectively dimensionless or if other unit convention is followed.
    *   **Arguments**: Same as `hij`.
    *   **Returns**: `res :: real(DP)` - The partial derivative `dh/dU_i`.

*   `part_deriv_hij_wrt_r(rij, tau_i, tau_j, zeta) result(res)`
    *   **Description**: Calculates the partial derivative of the `hij` function with respect to the input interatomic distance `rij`. The internal calculation uses `abs_rij_internal`, and the result is scaled by `1.0_DP/Bohr` to reflect the derivative with respect to `rij` in original distance units.
    *   **Arguments**: Same as `hij`.
    *   **Returns**: `res :: real(DP)` - The partial derivative `dh/drij`.

*   `second_part_deriv_hij_wrt_Ui_and_r(rij, tau_i, tau_j, zeta) result(res)`
    *   **Description**: Calculates the mixed second partial derivative of the `hij` function, first with respect to `U_i` (derived from `tau_i`) and then with respect to `rij`. The scaling of the result is consistent with the first derivatives.
    *   **Arguments**: Same as `hij`.
    *   **Returns**: `res :: real(DP)` - The mixed second partial derivative `d^2h/(dU_i dr_ij)`.

## Important Variables/Constants

*   `rij`: Interatomic distance.
*   `tau_i`, `tau_j`: Input parameters related to atomic properties (e.g., Hubbard U).
*   `zeta`: Damping exponent controlling the behavior of the damping function.
*   `Bohr`: The Bohr radius, used for internal unit conversion of distances.
*   The scaling factor `5.0_DP/16.0_DP` is used to convert input `tau` parameters to internal `U` parameters before they are used in the exponential damping term.

## Usage Examples

These functions are not typically called directly by end-users but are essential components for DFTB3 calculations. They are used by the `coulomb_short_gamma` module, which combines this damping with the core short-range gamma interaction terms.

```fortran
! Conceptual usage within coulomb_short_gamma module:
! USE damp_short_gamma
! REAL(DP) :: distance, tau1, tau2, zeta_exp, damping_factor, d_damping_dU1
!
! ! ... obtain distance, tau values, zeta ...
!
! damping_factor = hij(distance, tau1, tau2, zeta_exp)
! d_damping_dU1  = part_deriv_hij_wrt_Ui(distance, tau1, tau2, zeta_exp)
!
! ! core_gamma_term = ...
! ! damped_gamma_term = core_gamma_term * damping_factor
```

## Dependencies and Interactions

*   **`supplib`**: Used for the `DP` double precision kind parameter and the `Bohr` constant.
*   **`coulomb_short_gamma`**: This module is the primary consumer of the `damp_short_gamma` functions. The `hij` function and its derivatives are multiplied with the corresponding terms from `coulomb_short_gamma` to produce the final damped short-range gamma contributions for DFTB3.
*   The internal unit conversions (scaling of `tau` to `U`, and `rij` to Bohr units) are important for the correct application of the damping formula as specified in the DFTB3 literature.
```
