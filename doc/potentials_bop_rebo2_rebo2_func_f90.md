# src/potentials/bop/rebo2/rebo2_func.f90

## Overview

This file provides the specific mathematical functions that define the second-generation REBO (REBO2) potential by Brenner et al. These functions are called by the specialized `bop_kernel_rebo2.f90` to calculate various components of the potential energy and forces. The REBO2 potential is known for its complexity, including terms for conjugation, sophisticated cutoff schemes, and detailed angular dependencies.

The functions cover:
*   A conjugation switching function (`fconj`).
*   Cutoff functions (`fCin`, and conditionally `fCar`, `fCbo`, `fCnc` for screening) that rely on spline objects for interpolation between hard limits.
*   Attractive (`VA`) and repulsive (`VR`) pair potentials, which can use either analytical forms or spline interpolation based on the `SPLINE_POTENTIAL` macro.
*   A highly complex angular contribution function (`g`) that behaves differently for Carbon (using splines) versus other elements (like Hydrogen, using direct polynomial evaluation).
*   The bond order function (`bo`).
*   The length-dependent contribution to the bond order (`h`).
*   A REBO2-specific pair indexing function (`Z2pair`).

## Key Components

### Conjugation Function

*   `fconj(this, x, fx, dfx)`
    *   **Type**: Elemental Subroutine
    *   **Description**: A switching function that typically depends on a measure of conjugation `x`.
    *   **Implementation**: Returns `fx=1` (derivative `dfx=0`) if `x <= 2.0`. Returns `fx=0` (`dfx=0`) if `x >= 3.0`. Between `x=2.0` and `x=3.0`, it provides a smooth transition using `fx = 0.5 * (1.0 + cos(pi * (x - 2.0)))`.

### Cutoff Functions

These subroutines provide smoothly truncated interactions. They check for distances outside hard lower (`_l`) and upper (`_h`) limits (returning 1/0 and 0/0 for value/derivative respectively). If within the active range, they delegate to spline interpolation.

*   `fCin(this, ijpot, dr, val, dval)`: For inner cutoff. Uses `this%spl_fCin(ijpot)`.
*   `fCar(this, ijpot, dr, val, dval)` (`#ifdef SCREENING`): For attractive/repulsive range cutoff. Uses `this%spl_fCar(ijpot)`.
*   `fCbo(this, ijpot, dr, val, dval)` (`#ifdef SCREENING`): For bond-order specific cutoff. Uses `this%spl_fCbo(ijpot)`.
*   `fCnc(this, ijpot, dr, val, dval)` (`#ifdef SCREENING`): For neighbor counting/conjugation cutoff. Uses `this%spl_fCnc(ijpot)`.
    *   **Arguments**: `this (BOP_TYPE)`, `ijpot (pair type)`, `dr (distance)`, `val (out)`, `dval (out)`.
    *   **Implementation**: `call f_and_df(spline_object, dr, val, dval)`.

### Pair Potentials

*   `VA(this, ijpot, dr, val, dval)` and `VR(this, ijpot, dr, val, dval)`
    *   **Type**: Subroutines
    *   **Description**: Calculate attractive (`VA`) and repulsive (`VR`) pair potential terms and their derivatives.
    *   **Implementation**:
        *   If `#ifdef SPLINE_POTENTIAL` is defined and `dr` is within the spline's range, `call f_and_df(this%spl_VA(ijpot), ...)` (or `spl_VR`) is used.
        *   Otherwise, analytical forms are used:
            *   `VA` for C-C: Sum of three exponentials: `- (B1*exp(-beta1*dr) + B2*exp(-beta2*dr) + B3*exp(-beta3*dr))`. Parameters from `this%cc_B1` etc.
            *   `VA` for C-H or H-H: Single exponential: `- B1*exp(-beta1*dr)`.
            *   `VR` for all types: `(1 + Q/dr) * A*exp(-alpha*dr)`. Parameters from `this%cc_A`, etc.

### Angular Contribution (`g` function)

*   `g(this, ktyp, costh, n, val, dval_dcosth, dval_dN)`
    *   **Type**: Elemental Subroutine
    *   **Description**: Calculates the complex angular contribution to the bond order. Its behavior depends on the central atom type `ktyp` and coordination `n`.
    *   **Implementation**:
        *   If `ktyp == rebo2_C_` (Carbon):
            *   If coordination `n < 3.2`, uses spline `this%cc_g2_coeff` via `cc_g_from_spline`.
            *   If `n > 3.7`, uses spline `this%cc_g1_coeff` via `cc_g_from_spline`.
            *   If `3.2 <= n <= 3.7`, interpolates between values from `cc_g1_coeff` and `cc_g2_coeff` using a cosine switching function dependent on `n`. `dval_dN` captures the derivative w.r.t this switching.
        *   If `ktyp` is not Carbon (e.g., Hydrogen): Uses a direct polynomial evaluation: `val = spgh(1,ig) + spgh(2,ig)*costh + ...` where `ig` is an index derived from `costh` and `this%spgh` contains coefficients. `dval_dN` is not explicitly set here (likely zero).
*   `cc_g_from_spline(this, g_coeff, costh, val, dval)`
    *   **Type**: Elemental Subroutine
    *   **Description**: Helper to evaluate the C-C angular spline (`g_coeff` being `this%cc_g1_coeff` or `this%cc_g2_coeff`) at `costh`. It selects the correct polynomial segment based on `costh` relative to `this%cc_g_theta` node points.

### Bond Order Function

*   `bo(this, ktypi, ijpot, zij, fcij, faij, bij, dfbij)`
    *   **Type**: Elemental Subroutine
    *   **Description**: Calculates the bond order parameter `bij`.
    *   **Implementation**: `arg = 1.0 + zij; bij = arg ^ this%conpe(ktypi)`. `dfbij` is its derivative. `this%conpe(ktypi)` and related `this%conan(ktypi)`, `this%conpf(ktypi)` are pre-calculated exponents/factors.

### Length-Dependent Contribution (`h` function)

*   `h(this, ktypj, ktypi, ktypk, ijpot, ikpot, dr, val, dval)`
    *   **Type**: Elemental Subroutine
    *   **Description**: Calculates the length-dependent factor in the bond order.
    *   **Implementation**:
        *   If `ijpot + ikpot <= 4` (logic based on C-C, C-H, H-H pair type indices): `val = 1.0`, `dval = 0.0`.
        *   Otherwise: `val = this%conear(ijpot, ikpot) * exp(this%conalp * dr)`. `this%conear` and `this%conalp` are parameters.

### Pair Indexing

*   `Z2pair(this, ktypi, ktypj)`
    *   **Type**: Elemental Function
    *   **Description**: REBO2-specific mapping of two element types (`ktypi`, `ktypj`) to a single pair index.
    *   **Implementation**: If `ktypi` is Carbon (`rebo2_C_`), index is `ktypj`. If `ktypj` is Carbon, index is `ktypi`. Otherwise (H-H), index is `ktypi + ktypj`. (This assumes specific integer values for `rebo2_C_` and `rebo2_H_`).

## Important Variables/Constants
*   `rebo2_C_`: An integer parameter representing the type index for Carbon.
*   The functions extensively use parameters stored in `this (type BOP_TYPE)`, including simple scalar/array parameters (e.g., `this%cc_B1`, `this%conpe`, `this%conear`, `this%conalp`), spline objects (e.g., `this%spl_fCin`, `this%spl_VA`), and spline coefficients (e.g., `this%cc_g1_coeff`, `this%spgh`).

## Usage Examples
These are low-level routines called by `bop_kernel_rebo2.f90` as part of the overall energy and force computation for the REBO2 potential.

## Dependencies and Interactions
*   Relies on the `BOP_TYPE` (`rebo2_t`) structure for parameters and spline objects.
*   The behavior is significantly modified by preprocessor flags: `#ifdef SPLINE_POTENTIAL` (switches between analytical and spline pair potentials) and `#ifdef SCREENING` (enables additional cutoff functions).
*   The `g` function's complexity, with its dual approach for C-C versus C-H/H-H angles and spline interpolation, is a hallmark of REBO2.
```
