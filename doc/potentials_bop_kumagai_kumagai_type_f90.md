# src/potentials/bop/kumagai/kumagai_type.f90

## Overview

This file defines the `BOP_TYPE` derived data type, specifically tailored for Kumagai-type Bond-Order Potentials (BOPs). When included by `kumagai.f90` or `kumagai_scr.f90`, `BOP_TYPE` is effectively aliased to `kumagai_t` or `kumagai_scr_t`, respectively. This data type serves as the primary container for all data associated with a Kumagai potential instance. This includes:
*   The base parameter set (defined in `kumagai_params.f90`).
*   Mapping information for element types.
*   Cutoff function objects (of type `CUTOFF_T`, which is externally defined by the including module) and associated cutoff distances.
*   Allocatable arrays for storing internal neighbor lists and related bond data generated during calculations.

The structure of this `BOP_TYPE` is conditional on the `SCREENING` preprocessor macro to include fields necessary for screened interaction calculations. Notably, unlike some other BOP type definitions (e.g., for Brenner or Juslin), this version does not explicitly list fields for pre-storing many derived constants (like `expA`, `VA_f`); these are likely computed on-the-fly within the Kumagai-specific functions in `kumagai_func.f90`.

The file also declares the public interfaces for standard potential operations (`init`, `del`, `bind_to`, `energy_and_forces`, `register`), mapping them to macro-defined implementation functions.

## Key Components

### Data Types

*   `BOP_TYPE` (Public)
    *   **Description**: The central derived data type for Kumagai-like potentials.
    *   **Core Parameter Storage and Mapping**:
        *   `db :: type(BOP_DB_TYPE)`: An instance of `BOP_DB_TYPE` (which should be `kumagai_db_t`, defined in `kumagai_params.f90`). Stores the fundamental parameter set for the potential. Defaults to `Kumagai_CompMaterSci_39_457_Si`.
        *   `Z2db(MAX_Z) :: integer`: An array mapping atomic numbers (Z) to the internal indices used within the `db` parameter set.
    *   **State Variables**:
        *   `neighbor_list_allocated :: logical`: A flag, initialized to `.false.`, indicating if internal neighbor list arrays are allocated.
        *   `it :: integer`: An iteration counter, initialized to 0.
    *   **Cutoff Information** (Arrays typically dimensioned to `KUMAGAI_MAX_PAIRS`):
        *   `cut_in(:) :: type(CUTOFF_T)`: Array of cutoff function objects for the primary (inner) interaction range. The actual type of `CUTOFF_T` (e.g., `trig_off_t` or `exp_cutoff_t`) is determined by the including module (`kumagai.f90` or `kumagai_scr.f90`).
        *   `cut_in_l, cut_in_h, cut_in_h2 :: real(DP)`: Low cutoff distance, high cutoff distance, and its square, for the inner cutoff.
    *   **Screening-Specific Cutoff Fields** (conditional on `#ifdef SCREENING`):
        *   `cut_out(:), cut_bo(:) :: type(CUTOFF_T)`: Cutoff function objects for outer and bond-order specific ranges.
        *   `cut_out_l, cut_out_h, cut_bo_l, cut_bo_h :: real(DP)`: Corresponding low and high cutoff distances.
        *   `max_cut_sq :: real(DP)`: The square of the largest cutoff distance considering all active cutoffs.
        *   `Cmin, Cmax, dC :: real(DP)`: Parameters defining the S(C) screening function.
        *   `C_dr_cut :: real(DP)`: An effective cutoff distance used to identify potential screening neighbors.
        *   `screening_threshold, dot_threshold :: real(DP)`: Numerical thresholds for screening logic.
    *   **Internal Allocatable Arrays for Neighbor Lists and Bond Data**:
        *   Core lists: `neb(:)` (neighbor atom indices), `nbb(:)` (bond indices).
        *   Periodic boundary conditions: `dcell(:)` (displacement cell indices, not used if compiled for LAMMPS via `#ifndef LAMMPS`).
        *   Bond properties: `bndtyp(:)` (bond type index), `bndlen(:)` (bond length), `bndnm(:, :)` (normalized bond vector).
        *   Cutoff values: `cutfcnar(:)` (cutoff function value), `cutdrvar(:)` (cutoff function derivative).
    *   **Screening-Specific Allocatable Arrays** (conditional on `#ifdef SCREENING`):
        *   `cutfcnbo(:), cutdrvbo(:)`: Cutoff values/derivatives for bond-order terms.
        *   Screening neighbor lists: `sneb_seed(:)`, `sneb_last(:)`, `sneb(:)`, `sbnd(:)`. The type of `sbnd` is conditional on `LAMMPS`.
        *   Screening force components: `sfacbo(:)`.
        *   Screening derivative terms: `cutdrarik(:)`, `cutdrarjk(:)`, `cutdrboik(:)`, `cutdrbojk(:)`.

### Public Interfaces

Standard potential operation interfaces are declared, mapping generic names to macro-defined specific implementations from `kumagai_module.f90` (e.g., `init` maps to `INIT_FUNC` which becomes `kumagai_init` or `kumagai_scr_init`):
*   `interface init`
*   `interface del`
*   `interface bind_to`
*   `interface energy_and_forces`
*   `interface register` (also exports `REGISTER_FUNC` directly for some reason, possibly legacy)

## Important Variables/Constants

*   `KUMAGAI_MAX_PAIRS`, `MAX_Z`: Constants (defined in `kumagai_params.f90` or via macros) that dictate static array sizes.
*   `CUTOFF_T`: This is a macro whose definition (e.g., `trig_off_t` or `exp_cutoff_t`) is provided by the including file (`kumagai.f90` or `kumagai_scr.f90`). This determines the actual type of the `cut_in`, `cut_out`, and `cut_bo` objects.
*   The default value for the `db` component is `Kumagai_CompMaterSci_39_457_Si`.

## Usage Examples

This file defines the `BOP_TYPE` for Kumagai potentials. Other files, particularly `kumagai_module.f90`, declare variables of this type and use the associated procedures (`init`, `del`, etc.) to manage and execute the potential calculations.

```fortran
! In kumagai_module.f90 or user code:
! TYPE(BOP_TYPE) :: my_kumagai_potential ! Will be kumagai_t or kumagai_scr_t
!
! CALL my_kumagai_potential%init() ! Using the interface, db can be optional
! ! ... further operations ...
! CALL my_kumagai_potential%del()
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on `BOP_DB_TYPE` (expected to be `kumagai_db_t` from `kumagai_params.f90`).
    *   Depends on constants like `KUMAGAI_MAX_PAIRS`, `MAX_Z`.
    *   The `CUTOFF_T` macro is crucial and is defined externally by the module that includes this type definition.
*   **External Libraries:** None directly used in this file.
*   **Interactions:**
    *   `BOP_TYPE` is the central data structure for the Kumagai potential.
    *   It is instantiated and populated by `INIT_FUNC` (e.g., `kumagai_init`).
    *   The cutoff function objects (`cut_in`, etc.) are initialized by `BIND_TO_FUNC` (likely the default implementation) using the type `CUTOFF_T`. These objects' methods are then called by the actual cutoff functions (e.g., from `default_cutoff.f90` as included by `kumagai_func.f90`).
    *   The `#ifdef SCREENING` directive alters the type's structure to include fields for screened interactions.
    *   Allocatable arrays are managed by other parts of the module (typically allocated in the kernel or `BIND_TO_FUNC`, deallocated in `DEL_FUNC`).
```
