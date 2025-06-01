# src/potentials/bop/rebo2/rebo2_scr.f90

## Overview

This file defines the `rebo2_scr` Fortran module, which implements the REBO2+S potential. This is a **screened** version of the second-generation Reactive Empirical Bond-Order (REBO2) potential. The base REBO2 formalism is from Brenner et al. (J. Phys.: Condens. Matter 14, 783 (2002)), while the screening modifications are based on work by Pastewka, Pou, Perez, Gumbsch, & Moseler (Phys. Rev. B 78, 161402(R) (2008)).

The `rebo2_scr` module integrates various components by including shared REBO2 source files. The behavior of these components is tailored for the REBO2+S variant through several key preprocessor definitions active during their compilation:
*   `SCREENING`: Enables environmental screening functions.
*   `ALT_DIHEDRAL`: Selects an alternative formulation for dihedral (torsional) interaction terms.
*   `NUM_NEIGHBORS`: Enables terms dependent on coordination numbers.

The cutoff functions are based on trigonometric forms (`CUTOFF_T` is `trig_off_t`).

## Key Components

### Modules

*   `rebo2_scr`: Encapsulates the REBO2+S interatomic potential.
    *   **Functionality**: Provides initialization, parameter management, binding to atomic systems, and computation of energy, forces, and virial according to the REBO2+S formalism.
    *   **Dependencies**: Uses standard Atomistica modules (`iso_c_binding`, `supplib`, `particles`, `filter`, `neighbors`) and REBO2-specific modules (`table2d`, `table3d`, `rebo2_default_tables`).
    *   **Construction**: Assembled via `#include` directives for its specialized components, compiled with `SCREENING`, `ALT_DIHEDRAL`, and `NUM_NEIGHBORS` active.

### Preprocessor Definitions
This module defines several critical preprocessor macros:
*   `SCREENING`: Enables the calculation of screening effects from the local atomic environment.
*   `ALT_DIHEDRAL`: Activates an alternative functional form or handling of dihedral interactions compared to the standard `DIHEDRAL` macro possibly used in the non-screened `rebo2.f90`.
*   `NUM_NEIGHBORS`: Includes terms in the bond order and energy calculations that depend on the number and type of neighboring atoms (e.g., Pij, Fij correction factors).
*   `CUTOFF_T`: Set to `trig_off_t`, indicating trigonometric cutoff functions.
*   `BOP_NAME`: `rebo2_scr`
*   `BOP_NAME_STR`: `"rebo2_scr"`
*   `BOP_STR`: `"Rebo2Scr"`
*   `BOP_KERNEL`: `rebo2_scr_kernel` (points to the specialized kernel in `bop_kernel_rebo2.f90`, now compiled with these defines)
*   `BOP_TYPE`: `rebo2_scr_t`
*   Function aliases (`REGISTER_FUNC`, `INIT_FUNC`, etc.) are mapped to `rebo2_scr_` prefixed routines. An `INTERNAL_INIT_FUNC` is also defined.

### Included Files
The module structure is built by including:
*   `rebo2_type.f90`: Defines `rebo2_scr_t` (via `BOP_TYPE`).
*   `rebo2_db.f90`: Handles parameter database structures and initialization logic.
*   `rebo2_module.f90`: Contains high-level module procedures (constructor, destructor, bind, compute).
*   `bop_kernel_rebo2.f90`: The specialized REBO2 computational kernel, now compiled with `SCREENING`, `ALT_DIHEDRAL`, and `NUM_NEIGHBORS`.
*   `rebo2_func.f90`: REBO2-specific mathematical functions.
*   `rebo2_registry.f90`: Routines for registering the potential with `ptrdict`.

### Classes/Types
*   `rebo2_scr_t`: (Defined via `BOP_TYPE` in `rebo2_type.f90`) The primary derived data type for the REBO2+S potential.
*   Associated database types (e.g., from `rebo2_db.f90`).

### Functions/Subroutines
The public interface is largely defined in `rebo2_module.f90` (and aliased here). Key procedures include:
*   `rebo2_scr_init(...)` and potentially `rebo2_scr_internal_init(...)`: Initialize the `rebo2_scr_t` potential.
*   `rebo2_scr_del(...)`: Destructor.
*   `rebo2_scr_bind_to(...)`: Binds the potential to a particle system.
*   `rebo2_scr_energy_and_forces(...)`: Computes energy and forces.
*   `rebo2_scr_register(...)`: Registers the potential.

## Important Variables/Constants
Parameters are defined in `rebo2_db.f90` and sourced from `rebo2_default_tables.f90`. The compile-time definitions of `SCREENING`, `ALT_DIHEDRAL`, and `NUM_NEIGHBORS` are crucial for the behavior of this potential.

## Usage Examples
The REBO2+S potential would be used similarly to other BOPs in Atomistica, by creating an instance of `rebo2_scr_t` and using its interface functions.

```fortran
! Conceptual usage:
! USE rebo2_scr_module
! TYPE(rebo2_scr_t) :: my_rebo2s_potential
!
! CALL rebo2_scr_init(my_rebo2s_potential, ...) ! Or an init_default variant
! CALL rebo2_scr_bind_to(my_rebo2s_potential, particles, neighbor_list)
! CALL rebo2_scr_energy_and_forces(my_rebo2s_potential, particles, neighbor_list, E, F, V)
! CALL rebo2_scr_del(my_rebo2s_potential)
```

## Dependencies and Interactions
*   **Internal Dependencies:** Relies on its constituent included files, which are compiled with the specific preprocessor flags of this module.
*   **External Libraries:** Standard Atomistica core modules.
*   **Interactions:**
    *   The `rebo2_scr` module provides a sophisticated version of the REBO2 potential.
    *   The `SCREENING` macro enables environmental screening effects.
    *   The `ALT_DIHEDRAL` macro selects a specific formulation for torsional interactions, which might differ from a version compiled with the standard `DIHEDRAL` flag.
    *   `NUM_NEIGHBORS` ensures that coordination-dependent terms (like Pij, Fij corrections) are active.
    *   The choice of `CUTOFF_T` as `trig_off_t` means trigonometric cutoff functions will be used for truncating interactions.
```
