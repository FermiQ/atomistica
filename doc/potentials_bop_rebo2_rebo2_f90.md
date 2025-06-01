# src/potentials/bop/rebo2/rebo2.f90

## Overview

This file defines the `rebo2` Fortran module, which implements the second-generation Reactive Empirical Bond-Order (REBO2) potential. This potential, developed by Brenner et al. (J. Phys.: Condens. Matter 14, 783 (2002)), is widely used for simulating hydrocarbons and other carbon-based materials. It is known for its complexity, incorporating terms that depend on local coordination, bond angles, and dihedral (torsional) angles to more accurately model chemical bonding and reactivity.

The `rebo2` module is constructed by including several specialized source files:
*   Type definitions (`rebo2_type.f90`)
*   Parameter database structures and default values (`rebo2_db.f90`, `rebo2_default_tables.f90`)
*   High-level module procedures (`rebo2_module.f90`)
*   A specialized computational kernel (`bop_kernel_rebo2.f90`)
*   REBO2-specific mathematical functions (`rebo2_func.f90`)
*   Registration routines (`rebo2_registry.f90`)

Key features enabled in this implementation include `DIHEDRAL` interactions and `NUM_NEIGHBORS` dependent terms. Options for spline-interpolated cutoffs and potentials (`SPLINE_CUTOFF`, `SPLINE_POTENTIAL`) are present but commented out.

## Key Components

### Modules

*   `rebo2`: Encapsulates the REBO2 interatomic potential.
    *   **Functionality**: Provides initialization, parameter management, binding to atomic systems, and computation of energy, forces, and virial according to the REBO2 formalism.
    *   **Dependencies**: Uses several other Atomistica modules:
        *   `iso_c_binding` (Fortran standard)
        *   `supplib` (Supplementary utilities)
        *   `particles` (Particle data structures)
        *   `filter` (Atom filtering capabilities)
        *   `neighbors` (Neighbor list management)
        *   `table2d`, `table3d`: For handling 2D and 3D tabulated functions, which are part of the REBO2 potential (e.g., correction terms).
        *   `rebo2_default_tables`: Provides the actual data for these tables.
    *   **Construction**: Assembled via `#include` directives for its various specialized components.

### Preprocessor Definitions
This module defines several preprocessor macros to control compilation and behavior:
*   `DIHEDRAL`: Enables the calculation of dihedral torsional interaction terms.
*   `NUM_NEIGHBORS`: Enables terms in the bond order and energy that depend on the number and type of neighbors (e.g., Pij, Fij correction factors).
*   `CUTOFF_T`: Set to `trig_off_t`, indicating the use of trigonometric cutoff functions.
*   `BOP_NAME`: `rebo2`
*   `BOP_NAME_STR`: `"rebo2"`
*   `BOP_STR`: `"Rebo2"`
*   `BOP_KERNEL`: `rebo2_kernel` (points to the specialized kernel in `bop_kernel_rebo2.f90`)
*   `BOP_TYPE`: `rebo2_t`
*   Standard function aliases (`REGISTER_FUNC`, `INIT_FUNC`, etc.) are mapped to `rebo2_` prefixed routines. An additional `INIT_DEFAULT_FUNC` is also defined, likely for initializing with default embedded parameters.

### Classes/Types
*   `rebo2_t`: (Defined via `BOP_TYPE` in `rebo2_type.f90`) The primary derived data type for the REBO2 potential.
*   Associated database types (e.g., `rebo2_db_t`) are expected to be defined in `rebo2_db.f90` or `rebo2_type.f90`.

### Functions/Subroutines
The public interface is largely defined in `rebo2_module.f90`. Key procedures include:
*   `rebo2_init(...)` and `rebo2_init_default(...)`: Initialize the `rebo2_t` potential.
*   `rebo2_del(...)`: Destructor.
*   `rebo2_bind_to(...)`: Binds the potential to a particle system.
*   `rebo2_energy_and_forces(...)`: Computes energy and forces.
*   `rebo2_register(...)`: Registers the potential.

## Important Variables/Constants
Specific parameters for the REBO2 potential are defined in `rebo2_db.f90` and `rebo2_default_tables.f90`. The compile-time defines `DIHEDRAL` and `NUM_NEIGHBORS` significantly influence the calculations.

## Usage Examples
The REBO2 potential is used by creating an instance of `rebo2_t`, initializing it (e.g., with default parameters or from a file), and then using its methods to compute system properties.

```fortran
! Conceptual usage:
! USE rebo2_module ! (Assuming rebo2.f90 defines this module)
! TYPE(rebo2_t) :: my_rebo2_potential
!
! CALL rebo2_init_default(my_rebo2_potential) ! Initialize with default parameters
! ! or CALL rebo2_init(my_rebo2_potential, file="params.rebo2")
!
! CALL rebo2_bind_to(my_rebo2_potential, particles, neighbor_list)
! CALL rebo2_energy_and_forces(my_rebo2_potential, particles, neighbor_list, E, F, V)
! CALL rebo2_del(my_rebo2_potential)
```

## Dependencies and Interactions
*   **Internal Dependencies:** Relies on its constituent included files (`rebo2_type.f90`, `rebo2_db.f90`, `rebo2_module.f90`, `bop_kernel_rebo2.f90`, `rebo2_func.f90`, `rebo2_registry.f90`). The `rebo2_default_tables` module is crucial for providing tabulated data.
*   **External Libraries:** None explicitly listed beyond standard Fortran and Atomistica core modules.
*   **Interactions:**
    *   The `rebo2` module provides a comprehensive implementation of the REBO2 potential.
    *   It uses a specialized kernel (`bop_kernel_rebo2.f90`) to handle the complex terms.
    *   The use of `table2d` and `table3d` indicates that parts of the potential (likely correction terms like Pij, Fij, or torsional components) are evaluated using interpolation from pre-calculated tables.
    *   The `DIHEDRAL` and `NUM_NEIGHBORS` flags enable physics critical for accurate hydrocarbon simulations.
```
