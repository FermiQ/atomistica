# src/potentials/bop/juslin/juslin_params.f90

## Overview

This file defines the parameter sets for Juslin-type Bond-Order Potentials (BOPs), specifically tailored for systems like W-C-H (Tungsten-Carbon-Hydrogen) as described in the referenced literature. It introduces a `BOP_DB_TYPE` (typically aliased to `juslin_db_t` in the Juslin module context) derived data type to store these parameters. The file then declares specific, constant instances of this type, each corresponding to a published Juslin potential parameterization. These instances are collected into a private array `BOP_DB` for runtime selection.

A key distinction for the Juslin potential parameters compared to some other BOPs (like the basic Brenner) is the inclusion of parameters (`alpha`, `omega`, `m`) for the length-dependent `h` function that depend on a **triplet** of interacting atoms (j-i-k). Pair parameters are also stored in arrays sized for potentially non-symmetric interactions before explicit symmetrization in `juslin_module.f90`.

## Key Components

### Constants

*   `JUSLIN_MAX_REF :: integer, parameter`: Defines the maximum length of the reference string (e.g., citation) for a parameter set. Value: 1000.
*   `JUSLIN_MAX_EL :: integer, parameter`: Maximum number of distinct chemical elements a single Juslin parameter set can handle (typically 3 for W-C-H). Value: 3.
*   `JUSLIN_MAX_PAIRS :: integer, parameter`: Maximum number of unique ordered element pairs. Calculated using `PAIR_INDEX_NS(JUSLIN_MAX_EL, JUSLIN_MAX_EL, JUSLIN_MAX_EL)`. This is `JUSLIN_MAX_EL**2` if `PAIR_INDEX_NS` maps (i,j) to a flat array.

### Data Types

*   `BOP_DB_TYPE` (Public, typically aliased to `juslin_db_t`)
    *   **Description**: A Fortran derived data type designed to store a complete set of parameters for a Juslin interatomic potential.
    *   **Fields**:
        *   `nel :: integer`: The number of chemical elements defined in this parameter set (e.g., 3 for W-C-H).
        *   `nD0, nr0, ..., nalpha, nomega, nm, ... :: integer`: Counts for each parameter array, storing the actual number of defined values for validation.
        *   `ref :: character(JUSLIN_MAX_REF)`: Bibliographic reference string for the parameter set.
        *   `el(2, JUSLIN_MAX_EL) :: character`: Array storing chemical symbols of the elements.
        *   **Pair Interaction Parameters** (arrays typically dimensioned to `JUSLIN_MAX_EL**2`, allowing for non-symmetric i-j vs j-i parameters initially):
            *   `D0, r0, S, beta :: real(DP)`: Standard parameters for dimer binding energy, equilibrium distance, Pauling plot slope, and dimer stiffness.
            *   `gamma, c, d, h :: real(DP)`: Parameters for the angular `g` function. Note: `h` here is the angular parameter, distinct from the length-dependent `h` function.
            *   `n :: real(DP)`: Exponent in the bond order function `(1 + zij^n)`.
            *   `r1, r2 :: real(DP)`: Cutoff start/end distances.
        *   **Triplet Parameters for Length-Dependent `h` Function** (arrays typically dimensioned to `JUSLIN_MAX_EL**3` for j-i-k triplets):
            *   `alpha :: real(DP)`: Exponential factor in the `h` function.
            *   `omega :: real(DP)`: Prefactor for the `h` function.
            *   `m :: integer`: Exponent in the `h` function: `omega * exp((alpha*dr)^m)`.
        *   **Screening Parameters** (conditional on `#ifdef SCREENING`, arrays dimensioned to `JUSLIN_MAX_EL**2`):
            *   `or1, or2 :: real(DP)`: Outer cutoff start/end.
            *   `bor1, bor2 :: real(DP)`: Bond-order specific cutoff start/end.
            *   `Cmin, Cmax :: real(DP)`: Parameters for the screening function S(C).

### Predefined Parameter Sets

The file defines the following `BOP_DB_TYPE` instances as `parameter` constants:

*   `Juslin_J_Appl_Phys_98_123520_WCH`: Parameters for W-C-H systems from Juslin, N. et al. (2005). J. Appl. Phys., 98(12), 123520.
*   `Kuopanportti_P_Comp_Mat_Sci_111_525_FeCH`: Parameters for Fe-C-H systems from Kuopanportti, P. et al. (2016). Comp. Mat. Sci., 111, 525. (Note: This parameter set is for FeCH, also using the Juslin functional form).

These parameter sets are initialized with specific values for all fields, including the triplet parameters `alpha`, `omega`, and `m`. The pair parameter arrays (e.g., `D0`, `r0`) are sized for `JUSLIN_MAX_EL**2` entries, and the triplet parameters for `JUSLIN_MAX_EL**3` entries. A value of `r0 = -1.0_DP` is used as a flag in the initialization data to indicate that parameters for a given pair (i,j) should be copied from its symmetric counterpart (j,i) during the `BIND_TO_FUNC` stage in `juslin_module.f90`.

### Global Parameter Database

*   `BOP_DB(2) :: type(BOP_DB_TYPE), parameter, private`: A private array that aggregates the defined Juslin (and Juslin-like) parameter sets. This array is used by `INIT_FUNC` in `juslin_module.f90` to load parameters by reference string.

## Important Variables/Constants

*   The file's main purpose is to define the `BOP_DB_TYPE` for Juslin potentials and provide concrete parameter sets.
*   `PAIR_INDEX_NS` and `TRIPLET_INDEX_NS` (implicitly used for indexing `alpha`, `omega`, `m`) macros are important for mapping element indices to flat array indices.
*   The sizing of parameter arrays (`JUSLIN_MAX_EL**2` for pairs, `JUSLIN_MAX_EL**3` for triplets) reflects the potential's complexity.

## Usage Examples

This file provides data consumed by other parts of the Juslin potential module, primarily `INIT_FUNC` in `juslin_module.f90`.

```fortran
! Conceptual lookup in INIT_FUNC (juslin_module.f90):
! IF (TRIM(this%ref) == TRIM(Juslin_J_Appl_Phys_98_123520_WCH%ref)) THEN
!   this%db = Juslin_J_Appl_Phys_98_123520_WCH
! END IF
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Uses macros like `PAIR_INDEX_NS` (and implicitly `TRIPLET_INDEX_NS` for indexing `alpha`, `omega`, `m`).
*   **External Libraries:** None.
*   **Interactions:**
    *   Provides the fundamental numerical parameters that define the Juslin W-C-H potential and similar potentials.
    *   `juslin_module.f90::INIT_FUNC` uses the `BOP_DB` array to initialize `juslin_t` objects.
    *   The structure of `BOP_DB_TYPE`, especially the inclusion of triplet parameters and the non-symmetric pair parameter storage, is specific to the needs of the Juslin potential's functional form and setup.
```
