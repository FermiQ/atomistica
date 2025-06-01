# src/potentials/bop/kumagai/kumagai_params.f90

## Overview

This file defines the parameter sets for the Kumagai Bond-Order Potential (BOP), as used within the Atomistica simulation package. It specifies the `BOP_DB_TYPE` (typically aliased as `kumagai_db_t`) derived data type for storing these parameters. Currently, it provides one specific parameter set for Silicon (Si), sourced from Kumagai, Izumi, Hara, & Sakai, Comp. Mater. Sci. 39, 457 (2007).

The parameterization in this file is configured for a single element system (`KUMAGAI_MAX_EL = 1`).

## Key Components

### Constants

*   `KUMAGAI_MAX_REF :: integer, parameter`: Defines the maximum length of the reference string (e.g., citation) for a parameter set. Value: 1000.
*   `KUMAGAI_MAX_EL :: integer, parameter`: Maximum number of distinct chemical elements this specific parameter definition file handles. Value: 1 (implying parameters for a single element like Si).
*   `KUMAGAI_MAX_PAIRS :: integer, parameter`: Maximum number of unique chemical element pairs. Given `KUMAGAI_MAX_EL = 1`, this evaluates to 1.

### Data Types

*   `BOP_DB_TYPE` (Public, typically aliased to `kumagai_db_t`)
    *   **Description**: A Fortran derived data type designed to store a complete set of parameters for a Kumagai BOP.
    *   **Fields**:
        *   `nel :: integer`: The number of chemical elements defined (1 for the provided Si set).
        *   `nA, nB, nlambda1, ... :: integer`: Integer counts for each parameter array (e.g., `nA` would be 1 for a single element system).
        *   `ref :: character(KUMAGAI_MAX_REF)`: A character string containing the bibliographic reference.
        *   `el(2, KUMAGAI_MAX_EL) :: character`: Array storing the chemical symbol of the element.
        *   **Pair Interaction Parameters** (arrays of size `KUMAGAI_MAX_PAIRS`=1):
            *   `A, B :: real(DP)`: Pre-exponential factors for repulsive and attractive pair terms.
            *   `lambda1, lambda2 :: real(DP)`: Exponential decay factors for repulsive and attractive pair terms.
        *   **Bond Order Parameters**:
            *   `eta(KUMAGAI_MAX_EL) :: real(DP)`: Parameter for the `zij^eta` term in the bond order function.
            *   `delta(KUMAGAI_MAX_EL) :: real(DP)`: Exponent parameter in the bond order function `(1 + zij^eta)^(-delta)`.
            *   `alpha(KUMAGAI_MAX_PAIRS) :: real(DP)`: Parameter for the length-dependent `h` function.
            *   `beta(KUMAGAI_MAX_PAIRS) :: integer`: Integer exponent parameter for the length-dependent `h` function (`dr^beta`).
        *   **Angular `g` Function Parameters** (arrays of size `KUMAGAI_MAX_EL`=1):
            *   `c1, c2, c3, c4, c5 :: real(DP)`: Parameters for the complex angular function.
            *   `h :: real(DP)`: Parameter `h` used in the `(h - cos(theta))` terms of the angular function.
        *   **Cutoff Parameters** (arrays of size `KUMAGAI_MAX_PAIRS`=1):
            *   `r1, r2 :: real(DP)`: Cutoff start and end distances.
        *   **Screening Parameters** (conditional on `#ifdef SCREENING`, arrays of size `KUMAGAI_MAX_PAIRS`=1):
            *   `or1, or2 :: real(DP)`: Outer cutoff start/end.
            *   `bor1, bor2 :: real(DP)`: Bond-order specific cutoff start/end.
            *   `Cmin, Cmax :: real(DP)`: Parameters for the screening function S(C).
    *   *Note on parameter naming*: The comments in the `Kumagai_CompMaterSci_39_457_Si` data block sometimes use different local names for parameters (e.g., `n, c, d, h`). The field names listed above are taken from the `BOP_DB_TYPE` definition itself.

### Predefined Parameter Sets

*   `Kumagai_CompMaterSci_39_457_Si :: type(BOP_DB_TYPE), parameter`:
    *   **Description**: Parameters for pure Silicon (Si) from Kumagai, T., Izumi, S., Hara, S., & Sakai, S. (2007). Computational Materials Science, 39(2), 457-464.
    *   This instance is initialized with specific numerical values for all the fields in `BOP_DB_TYPE` relevant to a single-element Si system.

### Global Parameter Database

*   `kumagai_db(1) :: type(BOP_DB_TYPE), parameter, private`: A private array containing the single defined parameter set (`Kumagai_CompMaterSci_39_457_Si`). This allows the `INIT_FUNC` in `kumagai_module.f90` to access this default parameter set.

## Important Variables/Constants

*   The file's primary content is the definition of `BOP_DB_TYPE` for the Kumagai potential and the `Kumagai_CompMaterSci_39_457_Si` parameter set.
*   `KUMAGAI_MAX_EL = 1` significantly simplifies the array dimensions compared to multi-element potentials.

## Usage Examples

This file provides data that is consumed by other parts of the Kumagai potential module, primarily `INIT_FUNC` in `kumagai_module.f90`.

```fortran
! Conceptual lookup in INIT_FUNC (kumagai_module.f90):
! ! If no 'db' argument is passed, 'this%db' often defaults to the first
! ! (and in this case, only) entry in the internal 'kumagai_db' array.
! this%db = kumagai_db(1) ! Which is Kumagai_CompMaterSci_39_457_Si
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Uses the `PAIR_INDEX` macro (likely from `macros.inc`) for calculating array sizes, though with `MAX_EL=1` its effect is trivial.
*   **External Libraries:** None.
*   **Interactions:**
    *   Provides the fundamental numerical parameters that define the Kumagai potential for Silicon.
    *   `kumagai_module.f90::INIT_FUNC` uses the `kumagai_db` array (or a passed-in `BOP_DB_TYPE` instance) to configure a `kumagai_t` potential object.
    *   The structure of `BOP_DB_TYPE` and its parameters are used by `kumagai_func.f90` to calculate energies and forces.
```
