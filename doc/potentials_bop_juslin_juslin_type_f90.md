# src/potentials/bop/juslin/juslin_type.f90

## Overview

This file defines the `BOP_TYPE` derived data type, which is specifically tailored for Juslin-type Bond-Order Potentials (BOPs). When included by `juslin.f90` or `juslin_scr.f90`, `BOP_TYPE` is effectively aliased to `juslin_t` or `juslin_scr_t`, respectively. This data type is the central container for all information related to a Juslin potential instance. This includes:
*   The base parameters read from a database (defined in `juslin_params.f90`).
*   Pre-calculated constants derived from these base parameters for computational efficiency.
*   Specific parameters for the cosine-based cutoff functions that are implemented directly in `juslin_func.f90`.
*   Allocatable arrays for storing internal neighbor lists and associated bond data during calculations.

The structure of `BOP_TYPE` is conditional on the `SCREENING` preprocessor macro to include fields necessary for screened interaction calculations. Unlike some other BOP type definitions (e.g., for Brenner if it were to use generic `CUTOFF_T` objects), this version stores pre-calculated factors (`_fca`, `_fc`) for its specific cosine cutoff functions directly.

The file also declares the public interfaces for standard potential operations (`init`, `del`, `bind_to`, `energy_and_forces`, `register`), which map to macro-defined implementation functions.

## Key Components

### Data Types

*   `BOP_TYPE` (Public)
    *   **Description**: The main derived data type for Juslin-like potentials, holding the complete state and parameterization.
    *   **Core Parameter Storage**:
        *   `db :: type(BOP_DB_TYPE)`: An instance of `BOP_DB_TYPE` (which should be `juslin_db_t`, defined in `juslin_params.f90`). Stores the fundamental parameter set. Defaults to `Juslin_J_Appl_Phys_98_123520_WCH`.
        *   `ref(JUSLIN_MAX_REF) :: character`: A character string for selecting a parameter set by reference. Defaults to "*".
        *   `Z2db(MAX_Z) :: integer`: Maps atomic numbers (Z) to internal element indices used by the `db`.
    *   **Pre-calculated Constants** (Arrays typically dimensioned to `JUSLIN_MAX_PAIRS`):
        *   `c_sq, d_sq, c_d :: real(DP)`: Derived parameters for angular terms.
        *   `VR_f, expR, VA_f, expA :: real(DP)`: Pre-calculated components for repulsive/attractive pair potentials.
        *   `bo_exp, bo_fac, bo_exp1 :: real(DP)`: Pre-calculated terms for the bond order calculation.
    *   **Cosine Cutoff Function Parameters** (Arrays typically dimensioned to `JUSLIN_MAX_PAIRS`):
        *   `cut_in_h, cut_in_h2, cut_in_l :: real(DP)`: High cutoff distance, its square, and low cutoff distance for the inner cutoff range.
        *   `cut_in_fca, cut_in_fc :: real(DP)`: Pre-calculated factors (`PI / (r2 - r1)` and `-0.5 * PI / (r2 - r1)`) for the cosine-based inner cutoff function implemented in `juslin_func.f90`.
    *   **Screening-Specific Fields** (conditional on `#ifdef SCREENING`):
        *   `cut_out_h, cut_out_l, cut_out_fca, cut_out_fc :: real(DP)`: Parameters for the cosine-based outer cutoff.
        *   `cut_bo_h, cut_bo_l, cut_bo_fca, cut_bo_fc :: real(DP)`: Parameters for the cosine-based bond-order cutoff.
        *   `max_cut_sq :: real(DP)`: Square of the largest cutoff distance.
        *   `Cmin, Cmax, dC :: real(DP)`: Parameters for the S(C) screening function.
        *   `C_dr_cut :: real(DP)`: Effective cutoff for screening neighbors.
        *   `screening_threshold, dot_threshold :: real(DP)`: Numerical thresholds for screening logic.
    *   **State Variables**:
        *   `neighbor_list_allocated :: logical`: Flag, initialized to `.false.`, indicating if internal neighbor list arrays are allocated.
        *   `it :: integer`: Iteration counter, initialized to 0.
    *   **Internal Allocatable Arrays for Neighbor Lists and Bond Data**:
        *   These are identical to those in `brenner_type.f90`, including `neb`, `nbb`, `dcell` (conditional on `LAMMPS`), `bndtyp`, `bndlen`, `bndnm`, `cutfcnar`, `cutdrvar`.
        *   If screening is enabled: `cutfcnbo`, `cutdrvbo`, `sneb_seed`, `sneb_last`, `sneb`, `sbnd`, `sfacbo`, `cutdrarik`, `cutdrarjk`, `cutdrboik`, `cutdrbojk`.

### Public Interfaces

Standard potential operation interfaces are declared, mapping generic names to macro-defined specific implementations (e.g., `init` maps to `INIT_FUNC` which becomes `juslin_init` or `juslin_scr_init`):
*   `interface init`
*   `interface del`
*   `interface bind_to`
*   `interface energy_and_forces`
*   `interface register`

## Important Variables/Constants

*   `JUSLIN_MAX_PAIRS`, `MAX_Z`, `JUSLIN_MAX_REF`: Constants (defined in `juslin_params.f90` or via macros) dictating static array sizes.
*   The direct storage of pre-calculated factors for cosine cutoff functions (`cut_in_fca`, `cut_in_fc`, etc.) is a distinguishing feature of this type, corresponding to the direct implementation of these cutoffs in `juslin_func.f90`.
*   The default value for the `db` component is `Juslin_J_Appl_Phys_98_123520_WCH`.

## Usage Examples

This file defines the `BOP_TYPE` for Juslin potentials. Other files, particularly `juslin_module.f90`, declare variables of this type and use the associated procedures (`init`, `del`, etc.) to manage and execute the potential calculations.

```fortran
! In another file (e.g., juslin_module.f90 or user code):
! TYPE(BOP_TYPE) :: my_juslin_potential ! Will be juslin_t or juslin_scr_t
!
! CALL my_juslin_potential%init(ref="Juslin_WCH_2005") ! Using the interface
! ! ... further operations ...
! CALL my_juslin_potential%del()
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on `BOP_DB_TYPE` (expected to be `juslin_db_t` from `juslin_params.f90`).
    *   Depends on constants like `JUSLIN_MAX_PAIRS`, `MAX_Z`, `JUSLIN_MAX_REF`.
*   **External Libraries:** None directly used in this file.
*   **Interactions:**
    *   `BOP_TYPE` is the central data structure for the Juslin potential.
    *   It is instantiated and populated by `INIT_FUNC` (e.g., `juslin_init`).
    *   The cutoff parameter fields (`_fca`, `_fc`) are specifically calculated and stored by `BIND_TO_FUNC` in `juslin_module.f90` for use by the cosine cutoff functions in `juslin_func.f90`.
    *   The `#ifdef SCREENING` directive alters the type's structure to include fields for screened interactions.
    *   Allocatable arrays for neighbor lists are managed by other parts of the module (typically allocated in the kernel or `BIND_TO_FUNC`, deallocated in `DEL_FUNC`).
```
