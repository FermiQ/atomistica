# src/potentials/bop/kumagai/kumagai_scr.f90

## Overview

This file defines the `kumagai_scr` Fortran module, which implements a **screened** version of the Kumagai-Izumi-Hara-Sakai Bond-Order Potential (BOP). The base potential formulation is from Kumagai et al., Comp. Mater. Sci. 39, 457 (2007), with screening modifications potentially drawing from Pastewka et al., arXiv:1301.2142.

The `kumagai_scr` module, like its non-screened counterpart, is constructed by including several shared source files from the `kumagai` subdirectory and the generic `bop_kernel.f90`. The key differences in this version are the definition of the `SCREENING` preprocessor macro, which activates environmental screening logic within the included components, and the change of the default cutoff function type (`CUTOFF_T`) to `exp_cutoff_t`.

## Key Components

### Modules

*   `kumagai_scr`: This module encapsulates the screened Kumagai BOP.
    *   **Functionality**: It handles the initialization, binding to particle data, and computation of energy/forces, with screening effects enabled and potentially different cutoff behavior due to `exp_cutoff_t`.
    *   **Construction**: It includes the following files, processed with `SCREENING` defined and `CUTOFF_T` as `exp_cutoff_t`:
        *   `kumagai_params.f90`: Defines parameters for the Kumagai potential.
        *   `kumagai_type.f90`: Defines the `kumagai_scr_t` (via `BOP_TYPE`) derived data type. The `CUTOFF_T` members within this type will now be `exp_cutoff_t`.
        *   `kumagai_module.f90`: Contains high-level module procedures (e.g., `kumagai_scr_init`).
        *   `../bop_kernel.f90`: The generic BOP kernel, compiled with `SCREENING` active.
        *   `kumagai_func.f90`: Provides Kumagai-specific functions. If this file includes `default_cutoff.f90`, those default cutoff wrappers will now operate on `exp_cutoff_t` objects.
        *   `kumagai_registry.f90`: Contains routines to register the screened Kumagai potential.

### Classes/Types

*   `kumagai_scr_t`: (Defined via `BOP_TYPE` in `kumagai_type.f90`) The primary derived data type, now using `exp_cutoff_t` for its cutoff members.
*   `kumagai_scr_db_t`: (Defined via `BOP_DB_TYPE` in `kumagai_type.f90`) Type for Kumagai parameter database.

### Functions/Subroutines

The public interface is defined in `kumagai_module.f90` and aliased to `_scr_` versions (e.g., `kumagai_scr_init`).

### Preprocessor Definitions

Key definitions that modify the behavior of included files:
*   `SCREENING`: Enables environmental screening effects.
*   `CUTOFF_T`: Explicitly defined as `exp_cutoff_t` (exponential cutoff type). This will change the type of cutoff objects (e.g., `this%cut_in`) stored within the `kumagai_scr_t` type and thus alter the cutoff function's mathematical form.
*   Aliases for constants and function names (e.g., `KUMAGAI_SCR_MAX_EL`, `BOP_NAME` set to `kumagai_scr`).

## Important Variables/Constants

Parameters are sourced from `kumagai_params.f90`. The `SCREENING` macro and the specific definition of `CUTOFF_T` as `exp_cutoff_t` are the most critical distinguishing features of this module, influencing both the core interaction calculations and how interactions are truncated.

## Usage Examples

Usage is analogous to the non-screened Kumagai potential, but this version would be selected if screened interactions with exponential cutoffs are desired.

```fortran
! Conceptual usage:
! USE kumagai_scr_module
! TYPE(kumagai_scr_t) :: my_screened_kumagai_potential
!
! CALL kumagai_scr_init(my_screened_kumagai_potential, ...)
! CALL kumagai_scr_bind_to(my_screened_kumagai_potential, particles, neighbor_list)
! CALL kumagai_scr_energy_and_forces(my_screened_kumagai_potential, particles, neighbor_list, E, F, V)
! CALL kumagai_scr_del(my_screened_kumagai_potential)
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Same base dependencies as `kumagai.f90` (Atomistica core modules, and the set of included `kumagai_*.f90` and `bop_kernel.f90` files).
    *   The behavior of included files, especially `kumagai_type.f90` (regarding the type of its `cut_in`, etc. members) and any files that use `CUTOFF_T` or call cutoff functions (like `default_cutoff.f90` included in `kumagai_func.f90`), will be modified by `CUTOFF_T = exp_cutoff_t`.
*   **External Libraries:** None explicitly listed.
*   **Interactions:**
    *   The `kumagai_scr` module provides a screened Kumagai BOP, featuring exponential cutoffs.
    *   The `SCREENING` macro enables specific code paths in shared components for screened interactions.
    *   The redefinition of `CUTOFF_T` to `exp_cutoff_t` means that the cutoff function objects stored and used within this potential will be of the exponential type, leading to a different falloff behavior of interactions compared to the `trig_off_t` used in the non-screened version.
```
