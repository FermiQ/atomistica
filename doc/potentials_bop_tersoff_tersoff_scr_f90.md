# src/potentials/bop/tersoff/tersoff_scr.f90

## Overview

This file defines the `tersoff_scr` Fortran module, which implements a **screened** version of the Tersoff Bond-Order Potential (BOP). The base potential formulation is derived from the works of J. Tersoff, with screening modifications potentially drawing from methodologies described by Pastewka et al. (arXiv:1301.2142).

The `tersoff_scr` module leverages many of the same underlying component files as the non-screened Tersoff potential but compiles them with the `SCREENING` preprocessor macro defined. A significant difference in this screened version is the explicit definition of `CUTOFF_T` as `exp_cutoff_t`, indicating the use of exponential cutoff functions instead of the trigonometric ones used in the non-screened `tersoff.f90`.

## Key Components

### Modules

*   `tersoff_scr`: This module encapsulates the screened Tersoff BOP.
    *   **Functionality**: It handles the initialization, binding to particle systems, computation of energy and forces (including screening effects), and parameter registration, using exponential cutoff functions.
    *   **Construction**: It is built by including the following files, which are processed with `SCREENING` defined and `CUTOFF_T` set to `exp_cutoff_t`:
        *   `tersoff_params.f90`: Defines parameters for various Tersoff parameterizations.
        *   `tersoff_type.f90`: Defines `tersoff_scr_t` (via `BOP_TYPE`). The `CUTOFF_T` members within this type will now be `exp_cutoff_t`.
        *   `tersoff_module.f90`: Contains high-level module procedures (e.g., `tersoff_scr_init`).
        *   `../bop_kernel.f90`: The generic Bond-Order Potential kernel, compiled with `SCREENING` active.
        *   `tersoff_func.f90`: Provides Tersoff-specific functions. If this includes `default_cutoff.f90`, those wrappers will now operate on `exp_cutoff_t` objects.
        *   `tersoff_registry.f90`: Routines for registering the screened Tersoff potential.

### Classes/Types

*   `tersoff_scr_t`: (Defined via `BOP_TYPE` in `tersoff_type.f90`) The primary derived data type for the screened Tersoff potential, configured to use `exp_cutoff_t` for its cutoff function objects.
*   `tersoff_scr_db_t`: (Defined via `BOP_DB_TYPE` in `tersoff_type.f90`) Type for Tersoff parameter database entries.

### Functions/Subroutines

The public interface is largely defined in `tersoff_module.f90` and aliased to `_scr_` versions (e.g., `tersoff_scr_init`, `tersoff_scr_get_cutoff`).

### Preprocessor Definitions

Key definitions that modify the behavior of included files:
*   `SCREENING`: Enables environmental screening effects in the potential calculations.
*   `CUTOFF_T`: Explicitly defined as `exp_cutoff_t`. This changes the type of cutoff objects (e.g., `this%cut_in`) stored within the `tersoff_scr_t` type and thus alters the mathematical form of the cutoff functions used.
*   Aliases for constants (`TERSOFF_SCR_MAX_REF`, etc.) and function names (`BOP_NAME` set to `tersoff_scr`, etc.).

## Important Variables/Constants

Parameters are sourced from `tersoff_params.f90`. The `SCREENING` macro and the specific definition of `CUTOFF_T` as `exp_cutoff_t` are the most critical distinguishing features of this module, influencing both the core interaction calculations and how interactions are truncated.

## Usage Examples

Usage is analogous to the non-screened Tersoff potential, but this version would be selected if screened interactions with exponential cutoffs are desired.

```fortran
! Conceptual usage:
! USE tersoff_scr_module
! TYPE(tersoff_scr_t) :: my_screened_tersoff_potential
!
! CALL tersoff_scr_init(my_screened_tersoff_potential, ...)
! CALL tersoff_scr_bind_to(my_screened_tersoff_potential, particles, neighbor_list)
! CALL tersoff_scr_energy_and_forces(my_screened_tersoff_potential, particles, neighbor_list, E, F, V)
! CALL tersoff_scr_del(my_screened_tersoff_potential)
```

## Dependencies and Interactions

*   **Internal Dependencies (within Atomistica project):**
    *   `supplib`, `particles`, `neighbors`.
    *   The included `tersoff_*.f90` files and `../bop_kernel.f90`. The behavior of these included files, especially `tersoff_type.f90` (regarding the type of its `cut_in`, etc. members) and any files using `CUTOFF_T` or calling cutoff functions (like `default_cutoff.f90` included by `tersoff_func.f90`), will be modified by `CUTOFF_T = exp_cutoff_t`.
*   **External Libraries:** None explicitly listed.
*   **Interactions:**
    *   The `tersoff_scr` module provides a screened Tersoff BOP with exponential cutoffs.
    *   The `SCREENING` macro enables specific code paths in shared components for screened interactions.
    *   The redefinition of `CUTOFF_T` to `exp_cutoff_t` means that the cutoff function objects stored and used within this potential will be of the exponential type, leading to a different falloff behavior of interactions compared to the `trig_off_t` used in the non-screened version.
```
