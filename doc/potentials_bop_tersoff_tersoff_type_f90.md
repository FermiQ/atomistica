# src/potentials/bop/tersoff/tersoff_type.f90

## Overview

This file defines the `BOP_TYPE` derived data type, specifically tailored for Tersoff-type Bond-Order Potentials (BOPs). When included by `tersoff.f90` or `tersoff_scr.f90`, `BOP_TYPE` is effectively aliased to `tersoff_t` or `tersoff_scr_t`, respectively. This data type serves as the primary container for all data associated with a Tersoff potential instance. This includes:
*   The base parameter set (defined in `tersoff_params.f90`).
*   Mapping information for element types.
*   Cutoff function objects (of type `CUTOFF_T`). The actual type of `CUTOFF_T` (e.g., `trig_off_t` or `exp_cutoff_t`) is defined by the module that includes this type definition file (`tersoff.f90` or `tersoff_scr.f90`).
*   Associated cutoff distances.
*   Allocatable arrays for storing internal neighbor lists and related bond data generated during calculations.

The structure of this `BOP_TYPE` is conditional on the `SCREENING` preprocessor macro to include fields necessary for screened interaction calculations. Similar to the Kumagai potential's type definition, this version does not explicitly list fields for pre-storing many derived constants from the base parameters; these are typically calculated on-the-fly within the Tersoff-specific functions in `tersoff_func.f90`.

The file also declares the public interfaces for standard potential operations (`init`, `del`, `bind_to`, `energy_and_forces`, `register`), mapping them to macro-defined implementation functions.

## Key Components

### Data Types

*   `BOP_TYPE` (Public)
    *   **Description**: The central derived data type for Tersoff-like potentials.
    *   **Core Parameter Storage and Mapping**:
        *   `db :: type(BOP_DB_TYPE)`: An instance of `BOP_DB_TYPE` (which should be `tersoff_db_t`, defined in `tersoff_params.f90`). Stores the fundamental parameter set. Defaults to `Tersoff_PRB_39_5566_SiC`.
        *   `Z2db(MAX_Z) :: integer`: An array mapping atomic numbers (Z) to the internal indices used within the `db` parameter set.
    *   **State Variables**:
        *   `neighbor_list_allocated :: logical`: A flag, initialized to `.false.`, indicating if internal neighbor list arrays are allocated.
        *   `it :: integer`: An iteration counter, initialized to 0.
    *   **Cutoff Information** (Arrays typically dimensioned to `TERSOFF_MAX_PAIRS`):
        *   `cut_in(:) :: type(CUTOFF_T)`: Array of cutoff function objects for the primary (inner) interaction range. The specific type of `CUTOFF_T` is determined by the including module.
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
        *   Cutoff values: `cutfcnar(:)` (cutoff function value for pair interaction), `cutdrvar(:)` (its derivative).
    *   **Screening-Specific Allocatable Arrays** (conditional on `#ifdef SCREENING`):
        *   `cutfcnbo(:), cutdrvbo(:)`: Cutoff values/derivatives for bond-order terms.
        *   Screening neighbor lists: `sneb_seed(:)`, `sneb_last(:)`, `sneb(:)`, `sbnd(:)`.
        *   Screening force components: `sfacbo(:)`.
        *   Screening derivative terms: `cutdrarik(:)`, `cutdrarjk(:)`, `cutdrboik(:)`, `cutdrbojk(:)`.

### Public Interfaces

Standard potential operation interfaces are declared, mapping generic names to macro-defined specific implementations from `tersoff_module.f90` (e.g., `init` maps to `INIT_FUNC` which becomes `tersoff_init` or `tersoff_scr_init`):
*   `interface init`
*   `interface del`
*   `interface bind_to`
*   `interface energy_and_forces`
*   `interface register` (also exports `REGISTER_FUNC` directly)

## Important Variables/Constants

*   `TERSOFF_MAX_PAIRS`, `MAX_Z`: Constants (defined in `tersoff_params.f90` or via macros) that dictate static array sizes.
*   `CUTOFF_T`: This is a macro whose definition (e.g., `trig_off_t` or `exp_cutoff_t`) is provided by the including file (`tersoff.f90` or `tersoff_scr.f90`). This determines the actual type of the `cut_in`, `cut_out`, and `cut_bo` cutoff function objects.
*   The default value for the `db` component is `Tersoff_PRB_39_5566_SiC`.

## Usage Examples

This file defines the `BOP_TYPE` for Tersoff potentials. Other files, particularly `tersoff_module.f90`, declare variables of this type and use the associated procedures (`init`, `del`, etc.) to manage and execute the potential calculations.

```fortran
! In tersoff_module.f90 or user code:
! TYPE(BOP_TYPE) :: my_tersoff_potential ! Will be tersoff_t or tersoff_scr_t
!
! CALL my_tersoff_potential%init(db=some_tersoff_db_parameters) ! Using the interface
! ! ... further operations ...
! CALL my_tersoff_potential%del()
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on `BOP_DB_TYPE` (expected to be `tersoff_db_t` from `tersoff_params.f90`).
    *   Depends on constants like `TERSOFF_MAX_PAIRS`, `MAX_Z`.
    *   The `CUTOFF_T` macro is crucial and is defined externally by the module that includes this type definition.
*   **External Libraries:** None directly used in this file.
*   **Interactions:**
    *   `BOP_TYPE` is the central data structure for the Tersoff potential.
    *   It is instantiated and populated by `INIT_FUNC` (e.g., `tersoff_init`).
    *   The cutoff function objects (`cut_in`, etc.) are initialized by `BIND_TO_FUNC` (likely the default implementation from `default_bind_to_func.f90`) using the type `CUTOFF_T`. These objects' methods are then called by the actual cutoff functions (e.g., from `default_cutoff.f90` as included by `tersoff_func.f90`).
    *   The `#ifdef SCREENING` directive alters the type's structure to include fields for screened interactions.
    *   Allocatable arrays for neighbor lists are managed by other parts of the module.
```
