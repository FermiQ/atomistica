# src/potentials/bop/juslin/juslin_func.f90

## Overview

This file provides the specific mathematical function implementations that define the Juslin W-C-H (Tungsten-Carbon-Hydrogen) Bond-Order Potential. These functions are designed to be called by the generic `BOP_KERNEL` when the Juslin potential is selected for a simulation. The file includes routines for calculating cutoff function values, attractive and repulsive pair potential terms, the angular contribution to bond order, the bond order function itself, and the length-dependent contribution to the bond order.

A key difference from some other BOP implementations in Atomistica is that the cutoff functions (`fCin`, `fCar`, `fCbo`) are implemented directly here using a cosine-based form, rather than delegating to a method of a generic `CUTOFF_T` object. The length-dependent term `h` also has a functional form specific to the Juslin potential.

## Key Components

### Cutoff Functions

These functions implement a cosine-based smoothing for interactions.

*   `fCin(this, ijpot, dr, val, dval)`
    *   **Type**: Elemental Subroutine
    *   **Description**: Calculates the value (`val`) and derivative (`dval`) of the **inner cutoff function**.
    *   **Implementation**: If `dr` (distance) is greater than `this%cut_in_h(ijpot)` (high cutoff limit for pair type `ijpot`), `val` is 0. If `dr` is less than `this%cut_in_l(ijpot)` (low cutoff limit), `val` is 1. Between these limits, `val = 0.5 * (1.0 + cos(arg))`, where `arg = this%cut_in_fca(ijpot) * (dr - this%cut_in_l(ijpot))`. `dval` is derived from this. `cut_in_fca` and `cut_in_fc` are pre-calculated factors.

*   `fCar(this, ijpot, dr, val, dval)` (Compiled only if `#ifdef SCREENING` is defined)
    *   **Type**: Subroutine
    *   **Description**: Similar cosine-based cutoff, typically for the attractive/repulsive pair potential range when screening is active. Uses parameters `this%cut_out_h`, `this%cut_out_l`, `this%cut_out_fca`, and `this%cut_out_fc`.

*   `fCbo(this, ijpot, dr, val, dval)` (Compiled only if `#ifdef SCREENING` is defined)
    *   **Type**: Subroutine
    *   **Description**: Similar cosine-based cutoff, applied specifically to bond-order terms when screening is active. Uses parameters `this%cut_bo_h`, `this%cut_bo_l`, `this%cut_bo_fca`, and `this%cut_bo_fc`.

### Pair Potentials

*   `VA(this, ijpot, dr, val, dval)`
    *   **Type**: Elemental Subroutine
    *   **Description**: Calculates the attractive part (`val`) of the pair potential and its derivative (`dval`). The functional form is Morse-like: `val = -VA_f * exp(-expA * (dr - r0))`. Parameters `VA_f`, `expA`, `r0` are from `this` and `this%db`.
*   `VR(this, ijpot, dr, val, dval)`
    *   **Type**: Elemental Subroutine
    *   **Description**: Calculates the repulsive part (`val`) of the pair potential and its derivative (`dval`). The functional form is Morse-like: `val = VR_f * exp(-expR * (dr - r0))`. Parameters `VR_f`, `expR`, `r0` are from `this` and `this%db`.

### Angular Contribution to Bond Order

*   `g(this, ktypj, ktypi, ktypk, ijpot, ikpot, costh, val, dval_dcosth)`
    *   **Type**: Elemental Subroutine
    *   **Description**: Computes the angular contribution (`val`) for an angle j-i-k based on `costh` (cosine of the angle), and its derivative (`dval_dcosth`). The functional form `val = gamma * (1 + c_d - c_sq / (d_sq + (h_param + costh)^2))` is similar to that used in other Brenner-type potentials. Parameters are from `this` and `this%db`.

### Bond Order Function

*   `bo(this, ktypi, ijpot, zij, fcij, faij, bij, dfbij)`
    *   **Type**: Subroutine
    *   **Description**: Calculates the bond order parameter (`bij`) and its derivative (`dfbij`). `zij` is the sum of angular and length contributions. `fcij` and `faij` are cutoff values. The form is `bij = (1 + zij^n)^bo_exp` or similar, depending on `n`. Parameters are from `this` and `this%db`.

### Length-Dependent Contribution to Bond Order

*   `h(this, ktypj, ktypi, ktypk, ijpot, ikpot, dr, val, dval)`
    *   **Type**: Elemental Subroutine
    *   **Description**: Computes the length-dependent contribution (`val`) to the bond order term for the i-j bond, considering neighbor k. This function is specific to the Juslin potential.
    *   **Implementation**: It uses `TRIPLET_INDEX_NS(ktypi, ktypj, ktypk, this%db%nel)` to get parameters `alpha`, `omega`, and `m` which are specific to the triplet of element types i, j, k (where i is the central atom of the j-i-k angle, and the h-function is for the i-k bond usually). The functional form is `val = omega * exp((alpha*dr)^m)` or variants based on `m`. `dr` here usually refers to the length of the i-k bond.

### Pair Indexing

*   `Z2pair(this, i, j)`
    *   **Type**: Elemental Function
    *   **Description**: Generates a unique integer index for a given pair of (internal) element types `i` and `j`.
    *   **Implementation**: Uses the `PAIR_INDEX_NS` macro (likely a non-symmetric version, meaning i-j can be different from j-i if needed, though parameters are usually symmetrized). `this%db%nel` provides the number of elements.

## Important Variables/Constants

These functions are parameterized by values stored in the `this` argument (an instance of `juslin_t`, via the `BOP_TYPE` macro), particularly:
*   `this%db`: Contains the raw Juslin potential parameters (e.g., `r0`, `D0`, `S`, `beta`, `gamma`, `c`, `d`, `h_param_angular`, `alpha`, `omega`, `m_exp_h`).
*   Pre-calculated factors stored in `this` for efficiency: `expA`, `expR`, `VA_f`, `VR_f`, `c_sq`, `d_sq`, `c_d`, `bo_exp`, `bo_fac`, `bo_exp1`.
*   Cutoff limit parameters also stored in `this`: `cut_in_h`, `cut_in_l`, `cut_in_fca`, `cut_in_fc` (and similar for `cut_out_` and `cut_bo_` if screening is enabled).
*   The `h` function is notable for using triplet-specific parameters (`alpha`, `omega`, `m`) obtained via `TRIPLET_INDEX_NS`.

## Usage Examples

These are low-level routines called by the `BOP_KERNEL` (specifically `juslin_kernel`) during the computation of energies and forces. They are not intended for direct use by end-users.

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on the `BOP_TYPE` (`juslin_t`) structure for parameters.
    *   Uses macros `TRIPLET_INDEX_NS` and `PAIR_INDEX_NS` for indexing into parameter arrays.
*   **External Libraries:** None.
*   **Interactions:**
    *   These functions collectively define the energy landscape according to the Juslin W-C-H potential.
    *   The `BOP_KERNEL` orchestrates calls to these functions based on atomic geometries and neighbor lists.
    *   The direct implementation of cosine-based cutoff functions within this file is a specific choice for this potential, differing from a more object-oriented delegation to `CUTOFF_T` methods seen in `default_cutoff.f90`.
```
