# src/potentials/bop/tersoff/tersoff_params.f90

## Overview

This file defines the parameter sets for the Tersoff Bond-Order Potential (BOP) as used within the Atomistica simulation package. It specifies the `BOP_DB_TYPE` (typically aliased to `tersoff_db_t` in the Tersoff module context) derived data type for storing these parameters. A notable feature of this type definition is that most parameters have default values assigned directly.

The file provides one specific parameter set: `Tersoff_PRB_39_5566_SiC` for Silicon Carbide (SiC), based on the work by J. Tersoff in Phys. Rev. B 39, 5566 (1989). This parameter set is stored in a private global array `tersoff_db`, making it available for selection at runtime, typically by a reference string.

## Key Components

### Constants

*   `TERSOFF_MAX_REF :: integer, parameter`: Defines the maximum length of the reference string (e.g., citation) for a parameter set. Value: 1000.
*   `TERSOFF_MAX_EL :: integer, parameter`: Maximum number of distinct chemical elements a single Tersoff parameter set can handle. Value: 3.
*   `TERSOFF_MAX_PAIRS :: integer, parameter`: Maximum number of unique chemical element pairs. Calculated using `PAIR_INDEX(TERSOFF_MAX_EL, TERSOFF_MAX_EL, TERSOFF_MAX_EL)`.

### Data Types

*   `BOP_DB_TYPE` (Public, typically aliased to `tersoff_db_t`)
    *   **Description**: A Fortran derived data type designed to store a complete set of parameters defining a specific Tersoff interatomic potential. Many fields have default initializations.
    *   **Fields**:
        *   `nel :: integer`: The number of chemical elements defined in this parameter set (e.g., 2 for SiC). Initialized to -1.
        *   `nA, nB, nxi, ... :: integer`: Integer counts for each parameter array (e.g., `nA` would store the actual number of `A` parameters, up to `TERSOFF_MAX_PAIRS`).
        *   `ref :: character(TERSOFF_MAX_REF)`: A character string containing the bibliographic reference for the parameter set.
        *   `el(2, TERSOFF_MAX_EL) :: character`: Array storing the chemical symbols of the elements.
        *   **Pair Interaction Parameters** (arrays dimensioned to `TERSOFF_MAX_PAIRS`):
            *   `A, B :: real(DP)`: Pre-exponential factors for repulsive and attractive pair terms. Default: 1.0.
            *   `xi :: real(DP)`: Prefactor for the bond order term `b_ij`. Default: 1.0.
            *   `lambda, mu :: real(DP)`: Exponential decay factors for repulsive and attractive pair terms. Default: 1.0.
            *   `omega :: real(DP)`: Prefactor in the angular function `g`. Default: 1.0.
            *   `mubo :: real(DP)`: Decay factor for the length-dependent `h` function. Default: 0.0.
            *   `m :: integer`: Integer exponent in the length-dependent `h` function. Default: 1.
            *   `r1, r2 :: real(DP)`: Cutoff start and end distances. Defaults: 1.0 and 2.0.
        *   **Single-Element Parameters** (arrays dimensioned to `TERSOFF_MAX_EL`):
            *   `beta :: real(DP)`: Parameter in the bond order function. Default: 1.0.
            *   `n :: real(DP)`: Exponent `n` in the bond order function. Default: 1.0.
            *   `c, d, h :: real(DP)`: Parameters for the angular function `g`. Default: 1.0.
        *   **Screening Parameters** (conditional on `#ifdef SCREENING`, arrays dimensioned to `TERSOFF_MAX_PAIRS`):
            *   `or1, or2 :: real(DP)`: Outer cutoff start/end.
            *   `bor1, bor2 :: real(DP)`: Bond-order specific cutoff start/end.
            *   `Cmin, Cmax :: real(DP)`: Parameters for the screening function S(C).

### Predefined Parameter Sets

*   `Tersoff_PRB_39_5566_SiC :: type(BOP_DB_TYPE), parameter`:
    *   **Description**: Parameters for Silicon Carbide (SiC) from Tersoff, J. (1989). Phys. Rev. B, 39(8), 5566–5568.
    *   This instance is initialized with specific numerical values for Si-Si, C-C, and Si-C interactions. For mixed Si-C terms, parameters are often constructed from the pure element parameters (e.g., `A_SiC = sqrt(A_Si * A_C)` or `lambda_SiC = (lambda_Si + lambda_C)/2`).
    *   The `mubo` and `m` parameters for the length-dependent `h` function are explicitly set for the screened version, otherwise `mubo` is 0.0, effectively disabling this term for the non-screened version in the 1989 paper.

### Global Parameter Database

*   `tersoff_db(1) :: type(BOP_DB_TYPE), parameter, private`: A private array that holds the `Tersoff_PRB_39_5566_SiC` parameter set. This allows `INIT_FUNC` in `tersoff_module.f90` to access this default set if no other is specified.

## Important Variables/Constants

*   The file's main purpose is to define `BOP_DB_TYPE` for Tersoff potentials and provide the SiC parameter set.
*   The default initialization of many parameters directly within the `BOP_DB_TYPE` definition is a key characteristic.
*   `PAIR_INDEX` macro is used for calculating array sizes.
*   The `#ifdef SCREENING` directive conditionally includes parameters for screening and also affects the default values for `mubo` and `m` in the `Tersoff_PRB_39_5566_SiC` set.

## Usage Examples

This file provides data that is consumed by other parts of the Tersoff potential module, primarily `INIT_FUNC` in `tersoff_module.f90`.

```fortran
! Conceptual lookup in INIT_FUNC (tersoff_module.f90):
! IF (TRIM(this%ref) == TRIM(Tersoff_PRB_39_5566_SiC%ref)) THEN
!   this%db = Tersoff_PRB_39_5566_SiC
! ELSE ! Or if no ref provided, might default to the first entry
!   this%db = tersoff_db(1) ! which is Tersoff_PRB_39_5566_SiC
! END IF
```

## Dependencies and Interactions

*   **Internal Dependencies:** Uses the `PAIR_INDEX` macro.
*   **External Libraries:** None.
*   **Interactions:**
    *   Provides the fundamental numerical parameters for the Tersoff potential, particularly for SiC.
    *   `tersoff_module.f90::INIT_FUNC` uses the `tersoff_db` array or a passed `BOP_DB_TYPE` instance to configure a `tersoff_t` potential object.
    *   The structure of `BOP_DB_TYPE` and its parameters are used by `tersoff_func.f90` to calculate energies and forces.
```
