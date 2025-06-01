# src/potentials/bop/juslin/juslin_scr.f90

## Overview

This file defines the `juslin_scr` Fortran module, which implements a **screened** version of the Juslin W-C-H (Tungsten-Carbon-Hydrogen) Bond-Order Potential. The base functional forms are taken from the work by Juslin et al. (J. Appl. Phys. 98, 123520 (2005)). This version explicitly includes screening functions, which modify interatomic interactions based on the local atomic environment.

The `juslin_scr` module is constructed by including several other source files. Crucially, the `SCREENING` and `EXP_BOP` preprocessor macros are defined before these includes. The `SCREENING` macro activates logic for environmental screening effects within the generic BOP kernel and other components. The `EXP_BOP` macro may further influence the functional forms used, potentially related to cutoff functions or other parts of the bond order expression, possibly favoring exponential forms.

## Key Components

### Modules

*   `juslin_scr`: This module encapsulates the screened Juslin W-C-H potential. It provides the necessary types, parameters, and subroutines to calculate energies and forces for atomic systems where screening effects are important.
    *   **Functionality**: The module handles initialization, binding to particle data, and computation of energy/forces, with screening and potentially other `EXP_BOP`-related modifications enabled.
    *   **Construction**: It includes the following files, processed with `SCREENING` and `EXP_BOP` defined:
        *   `juslin_params.f90`: Defines parameters for the Juslin potential.
        *   `juslin_type.f90`: Defines the `juslin_scr_t` (via `BOP_TYPE`) derived data type. This type may include additional fields or altered interpretations due to the defined macros.
        *   `juslin_module.f90`: Contains high-level module procedures (e.g., `juslin_scr_init`, `juslin_scr_energy_and_forces`).
        *   `../bop_kernel.f90`: The generic Bond-Order Potential kernel, now compiled with screening logic.
        *   `juslin_func.f90`: Provides Juslin-specific functions, whose behavior or usage by the kernel might be affected by `SCREENING` or `EXP_BOP`.
        *   `juslin_registry.f90`: Contains routines to register the screened Juslin potential.

### Classes/Types

*   `juslin_scr_t`: (Defined via `BOP_TYPE` in `juslin_type.f90`) The primary derived data type for the screened Juslin potential.
*   `juslin_db_scr_t`: (Defined via `BOP_DB_TYPE` in `juslin_type.f90` using parameters from `juslin_params.f90`) Type for Juslin parameter database.

### Functions/Subroutines

The public interface, defined in `juslin_module.f90` and aliased, includes:
*   `juslin_scr_init(...)`: Initializes the `juslin_scr_t` object.
*   `juslin_scr_del(...)`: Destroys the `juslin_scr_t` object.
*   `juslin_scr_bind_to(...)`: Binds the potential to particle data and neighbor lists, accounting for screening.
*   `juslin_scr_energy_and_forces(...)`: Computes energy and forces with screening.
*   `juslin_scr_force(...)`: Potentially a variant force computation.
*   `juslin_scr_register(...)`: Registers the screened potential.

### Preprocessor Definitions

Key definitions that modify the behavior of included files:
*   `SCREENING`: Enables environmental screening effects in the potential calculations.
*   `EXP_BOP`: May select exponential forms for cutoffs or other parts of the potential. The exact impact depends on how this macro is used within the included files (e.g., `juslin_type.f90` or `bop_kernel.f90`).
*   Aliases for constants and function names (e.g., `JUSLIN_SCR_MAX_EL`, `BOP_NAME` set to `juslin_scr`, `BOP_TYPE` set to `juslin_scr_t`).

## Important Variables/Constants

Parameters are sourced from `juslin_params.f90`. The `SCREENING` and `EXP_BOP` macros are the most significant aspect of this file, tailoring the included generic components to this specific potential variant.

## Usage Examples

Usage is analogous to the non-screened Juslin potential, but this version would be selected if screened interactions are desired.

```fortran
! Conceptual usage:
! USE juslin_scr_module
! TYPE(juslin_scr_t) :: my_screened_juslin_potential
!
! CALL juslin_scr_init(my_screened_juslin_potential, ...)
! CALL juslin_scr_bind_to(my_screened_juslin_potential, particles, neighbor_list)
! CALL juslin_scr_energy_and_forces(my_screened_juslin_potential, particles, neighbor_list, E, F, V)
! CALL juslin_scr_del(my_screened_juslin_potential)
```

## Dependencies and Interactions

*   **Internal Dependencies:** Same base dependencies as `juslin.f90` (core Atomistica modules, and the set of included `juslin_*.f90` and `bop_kernel.f90` files).
*   **External Libraries:** None explicitly listed.
*   **Interactions:**
    *   The `juslin_scr` module provides a screened Juslin W-C-H potential.
    *   The `SCREENING` macro activates specific code paths in shared components like `bop_kernel.f90` to calculate and apply screening modifications to bond strengths.
    *   The `EXP_BOP` macro might influence details like the choice of cutoff function types within `juslin_type.f90` if `CUTOFF_T` is defined conditionally based on `EXP_BOP`, or it could affect other functional forms. If `juslin_func.f90` already implements specific cutoffs (as it does), `EXP_BOP` might be used by the kernel in how it interprets or applies those cutoffs or other terms.
```
