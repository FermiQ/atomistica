# src/potentials/bop/tersoff/tersoff_func.f90

## Overview

This file provides the specific mathematical functions that define the interactions for the Tersoff Bond-Order Potential (BOP). These functions are designed to be called by the generic `BOP_KERNEL` when the Tersoff potential is active. The file includes implementations for:
*   Attractive (`VA`) and repulsive (`VR`) pair potential terms.
*   The angular contribution (`g`) to the bond order, which depends on the angle formed by three atoms.
*   The bond order function (`bo`) itself, which modulates the strength of the attractive interaction.
*   A length-dependent contribution (`h`) to the bond order, affecting interactions based on relative bond lengths.
Standard cutoff function implementations (like `fCin`) are made available by including `../default_cutoff.f90`.

The functional forms used are characteristic of the Tersoff potential model.

## Key Components

### Included Files
*   `../default_cutoff.f90`: This file is included to provide standard implementations for cutoff functions (e.g., `fCin` for inner cutoff, and conditionally `fCar`, `fCbo` if screening were enabled for a Tersoff variant, though Tersoff typically doesn't use these extra screening cutoffs). These ensure smooth truncation of interactions.

### Pair Potentials

*   `VA(this, ijpot, dr, val, dval)`
    *   **Type**: Subroutine
    *   **Description**: Calculates the attractive part (`val`) of the pair potential and its derivative (`dval`) with respect to distance `dr`.
    *   **Implementation**: `val = -B(ijpot) * exp(-mu(ijpot) * dr)`. The parameters `B` (pre-exponential factor) and `mu` (decay factor) are specific to the atom pair type `ijpot` and are retrieved from `this%db`.
*   `VR(this, ijpot, dr, val, dval)`
    *   **Type**: Subroutine
    *   **Description**: Calculates the repulsive part (`val`) of the pair potential and its derivative (`dval`).
    *   **Implementation**: `val = A(ijpot) * exp(-lambda(ijpot) * dr)`. The parameters `A` (pre-exponential factor) and `lambda` (decay factor) are specific to the atom pair type `ijpot` and are retrieved from `this%db`.

### Angular Contribution to Bond Order

*   `g(this, ktypj, ktypi, ktypk, ijpot, ikpot, costh, val, dval_dcosth)`
    *   **Type**: Subroutine
    *   **Description**: Computes the angular contribution (`val`) to the bond order term for an angle formed by atoms j-i-k, where `i` is the central atom. It also returns the derivative of this contribution with respect to `costh` (`dval_dcosth`).
    *   **Implementation**: `val = omega(ikpot) * (1.0 + c(ktypi)^2/d(ktypi)^2 - c(ktypi)^2 / (d(ktypi)^2 + (h(ktypi) - costh)^2))`.
        *   `omega` is a parameter for the i-k bond type (`ikpot`).
        *   `c`, `d`, and `h` are parameters specific to the central atom type `ktypi`.
        All parameters are retrieved from `this%db`.

### Bond Order Function

*   `bo(this, ktypi, ijpot, zij, fcij, faij, bij, dfbij)`
    *   **Type**: Subroutine
    *   **Description**: Calculates the bond order parameter (`bij`) for a bond of type `ijpot`. `zij` is the sum of angular (`g`) and length-dependent (`h`) contributions from its environment. `fcij` and `faij` are cutoff function values.
    *   **Implementation**: If `zij > 0`, `bij = xi(ijpot) * (1.0 + beta(ktypi)^n(ktypi) * zij^n(ktypi))^(-0.5 / n(ktypi))`. Otherwise, `bij = 1.0`.
        *   `xi` is a parameter for the i-j pair type (`ijpot`).
        *   `beta` and `n` are parameters for the central atom type `ktypi`.
        `dfbij` is the derivative of `bij`. All parameters from `this%db`.

### Length-Dependent Contribution to Bond Order

*   `h(this, ktypj, ktypi, ktypk, ijpot, ikpot, dr, val, dval)`
    *   **Type**: Subroutine
    *   **Description**: Computes the length-dependent contribution (`val`) to the bond order term for the i-k bond (where `i` is the central atom of the j-i-k angle), based on `dr` (typically `rij - rik`). It also returns its derivative (`dval`).
    *   **Implementation**: `val = exp((2*mu*dr)^m)` or variants depending on the integer parameter `m(ikpot)`. If `mu(ikpot)` (named `mubo` in parameters) is zero, `val = 1.0`. Parameters `mubo` and `m` are for the i-k pair type (`ikpot`) from `this%db`. This functional form is similar to that used in some Brenner potentials.

### Pair Indexing

*   `Z2pair(this, ktypi, ktypj)`
    *   **Type**: Function
    *   **Description**: Generates a unique integer index for a given pair of (internal) element types `ktypi` and `ktypj`.
    *   **Implementation**: Uses the `PAIR_INDEX(ktypi, ktypj, this%db%nel)` macro. `this%db%nel` provides the number of elements defined in the parameter set.

## Important Variables/Constants

These functions are parameterized by values stored in the `this` argument (an instance of `tersoff_t`, via the `BOP_TYPE` macro), specifically within its `db` component:
*   Pair potential parameters: `A, B, lambda, mu` (per `ijpot` pair type).
*   Angular function `g` parameters: `omega` (per `ikpot` pair type), and `c, d, h` (per central atom type `ktypi`).
*   Bond order function `bo` parameters: `xi` (per `ijpot` pair type), and `beta, n` (per central atom type `ktypi`).
*   Length contribution `h` parameters: `mubo, m` (per `ikpot` pair type).

## Usage Examples

These functions are low-level routines called by the `BOP_KERNEL` (specifically `tersoff_kernel`, which is an alias for the generic BOP kernel) during the computation of energies and forces. They are not intended for direct use by end-users.

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on the `BOP_TYPE` (`tersoff_t`) structure for parameters.
    *   Uses the `PAIR_INDEX` macro.
    *   Includes `../default_cutoff.f90` for standard cutoff function implementations (`fCin`, etc.), which are then used by the `BOP_KERNEL`.
*   **External Libraries:** None.
*   **Interactions:**
    *   These functions collectively define the energy landscape according to the Tersoff potential model.
    *   The `BOP_KERNEL` orchestrates calls to these functions based on atomic geometries and neighbor lists.
```
