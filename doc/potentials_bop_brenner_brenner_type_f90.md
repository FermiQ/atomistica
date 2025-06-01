# src/potentials/bop/brenner/brenner_type.f90

## Overview

This file defines the core Fortran derived data type, `BOP_TYPE`, used to represent a Brenner-type Bond-Order Potential (BOP) within the Atomistica simulation package. When `brenner.f90` or `brenner_scr.f90` includes this file, `BOP_TYPE` is effectively aliased to `brenner_t` or `brenner_scr_t`, respectively. This type serves as a container for all parameters, pre-calculated values, cutoff function objects, and allocatable arrays for internal neighbor lists required during energy and force calculations. The structure of `BOP_TYPE` can be modified by preprocessor directives, particularly `#ifdef SCREENING`, to include fields necessary for screened interaction calculations.

The file also declares the public interfaces for the standard operations associated with a potential: `init`, `del`, `bind_to`, `energy_and_forces`, and `register`. These generic interface names are mapped to specific implementation subroutines (like `brenner_init`, `brenner_energy_and_forces`, etc.) via preprocessor macros defined in the main module files (`brenner.f90`, `brenner_scr.f90`).

## Key Components

### Data Types

*   `BOP_TYPE` (Public)
    *   **Description**: This is the central derived data type for Brenner-like potentials. It holds the complete state and parameterization of a potential instance.
    *   **Core Parameter Storage**:
        *   `db :: type(BOP_DB_TYPE)`: An instance of `BOP_DB_TYPE` (defined in `brenner_params.f90`), which stores the fundamental parameter set for the potential (e.g., interaction strengths, equilibrium distances, angular terms). Defaults to `Erhart_PRB_71_035211_SiC`.
        *   `ref(BRENNER_MAX_REF) :: character`: A character array to store a reference string, used to select a specific parameter set from the global `BOP_DB` array at initialization. Defaults to "*".
        *   `Z2db(MAX_Z) :: integer`: An array mapping atomic numbers (Z) to the internal indices used within the `db` parameter set.
    *   **Pre-calculated Constants** (Arrays typically dimensioned to `BRENNER_MAX_PAIRS` for different pair types):
        *   `c_sq, d_sq, c_d :: real(DP)`: Derived parameters for angular terms, pre-calculated for efficiency.
        *   `VR_f, expR, VA_f, expA :: real(DP)`: Pre-calculated components of the repulsive and attractive Morse-style pair potentials.
        *   `bo_exp, bo_fac, bo_exp1 :: real(DP)`: Pre-calculated exponents and factors for the bond order calculation.
    *   **Cutoff Parameters**:
        *   `cut_in(:) :: type(CUTOFF_T)`: Array of cutoff function objects (type determined by `CUTOFF_T` macro, e.g., `trig_off_t` or `exp_cutoff_t`) for the primary (inner) interaction range.
        *   `cut_in_h, cut_in_h2, cut_in_l :: real(DP)`: High cutoff distance, its square, and low cutoff distance for the inner cutoff.
    *   **Screening-Specific Fields** (conditional on `#ifdef SCREENING`):
        *   `cut_out(:), cut_bo(:) :: type(CUTOFF_T)`: Cutoff function objects for outer and bond-order specific ranges.
        *   `cut_out_h, cut_out_l, cut_bo_h, cut_bo_l :: real(DP)`: Corresponding high and low cutoff distances.
        *   `max_cut_sq :: real(DP)`: The square of the largest cutoff distance considering all active cutoffs.
        *   `Cmin, Cmax, dC :: real(DP)`: Parameters defining the screening function S(C).
        *   `C_dr_cut :: real(DP)`: An effective cutoff distance used to identify potential screening neighbors.
        *   `screening_threshold, dot_threshold :: real(DP)`: Numerical thresholds for screening logic.
    *   **State Variables**:
        *   `neighbor_list_allocated :: logical`: A flag, initialized to `.false.`, indicating whether the internal neighbor list arrays have been allocated.
        *   `it :: integer`: An iteration counter, possibly for debugging or tracking usage, initialized to 0.
    *   **Internal Allocatable Arrays for Neighbor Lists and Bond Data**:
        *   `neb(:), nbb(:)`: Store indices of neighboring atoms and corresponding bond indices.
        *   `dcell(:)`: Stores displacement cell indices for periodic boundary conditions (not used if compiled for LAMMPS).
        *   `bndtyp(:)`: Stores the type of each bond (an integer index).
        *   `bndlen(:)`: Stores the length of each bond.
        *   `bndnm(:, :)`: Stores the normalized direction vector of each bond.
        *   `cutfcnar(:), cutdrvar(:)`: Store the value of the primary cutoff function and its derivative for each bond.
    *   **Screening-Specific Allocatable Arrays** (conditional on `#ifdef SCREENING`):
        *   `cutfcnbo(:), cutdrvbo(:)`: Values and derivatives of the bond-order cutoff function.
        *   `sneb_seed(:), sneb_last(:)`: Pointers defining segments in `sneb` and `sbnd` for each primary bond's screening environment.
        *   `sneb(:)`: Stores indices of atoms that screen a particular bond.
        *   `sbnd(:)`: Stores bond indices related to the screening atoms.
        *   `sfacbo(:)`: Stores a screening-related factor for force calculations.
        *   `cutdrarik(:), cutdrarjk(:), cutdrboik(:), cutdrbojk(:)`: Store derivatives of cutoff functions with respect to the positions of screening atoms.

### Public Interfaces

The file declares generic interfaces for standard potential operations. These allow other parts of the code to call, for example, `call potential%init(...)` regardless of the specific potential type. The actual implementations (`brenner_init`, `brenner_scr_init`, etc.) are provided by `INIT_FUNC` (and similarly for `DEL_FUNC`, `BIND_TO_FUNC`, `COMPUTE_FUNC`, `REGISTER_FUNC`) which are macros defined in the main module files (`brenner.f90`, `brenner_scr.f90`).
*   `interface init`: Maps to `INIT_FUNC`.
*   `interface del`: Maps to `DEL_FUNC`.
*   `interface bind_to`: Maps to `BIND_TO_FUNC`.
*   `interface energy_and_forces`: Maps to `COMPUTE_FUNC`.
*   `interface register`: Maps to `REGISTER_FUNC`.

## Important Variables/Constants

*   `CUTOFF_T`: A preprocessor macro that determines the data type for cutoff function objects (e.g., `trig_off_t`, `exp_cutoff_t`). This affects the type of `cut_in`, `cut_out`, and `cut_bo`.
*   `BRENNER_MAX_PAIRS`, `MAX_Z`, `BRENNER_MAX_REF`: Constants (typically defined in `brenner_params.f90` or via macros) that dictate the size of various static arrays within `BOP_TYPE`.
*   The default value for the `db` component is `Erhart_PRB_71_035211_SiC`, which is one of the predefined parameter sets in `brenner_params.f90`.

## Usage Examples

This file defines the `BOP_TYPE`. Code in other files, particularly `brenner_module.f90`, will declare variables of this type and use the associated procedures (`init`, `del`, etc.) to manage and utilize the potential.

```fortran
! In another file (e.g., brenner_module.f90 or a user's simulation setup):
! TYPE(BOP_TYPE) :: my_brenner_potential
!
! CALL my_brenner_potential%init(ref="SomeRefString") ! Using the interface
! ! ... further operations ...
! CALL my_brenner_potential%del()
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on `BOP_DB_TYPE` (expected to be available from `brenner_params.f90`).
    *   Depends on constants like `BRENNER_MAX_PAIRS`, `MAX_Z`, `BRENNER_MAX_REF`.
    *   The `CUTOFF_T` macro points to specific cutoff type definitions (e.g., `trig_off_t` from `cutoff_types.f90`).
*   **External Libraries:** None directly used in this file.
*   **Interactions:**
    *   `BOP_TYPE` is the cornerstone of the Brenner potential implementation, holding all its data.
    *   It is instantiated and managed by the subroutines exposed through the public interfaces (`init`, `del`, etc.).
    *   The `#ifdef SCREENING` directive significantly alters the fields included in the type, tailoring it to either standard or screened Brenner potential calculations.
    *   The allocatable arrays are managed internally by the module procedures, typically allocated in `BIND_TO_FUNC` or the `BOP_KERNEL` and deallocated in `DEL_FUNC`.
```
