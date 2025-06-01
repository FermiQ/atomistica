# src/potentials/bop/tersoff/tersoff.f90

## Overview

This file defines the `tersoff` Fortran module, which implements the Tersoff interatomic potential. The Tersoff potential is a widely used Bond-Order Potential (BOP) for covalent systems, particularly semiconductors like Silicon, Germanium, and Carbon. This implementation is based on the formulation described in a series of papers by J. Tersoff:
*   Phys. Rev. Lett. 56, 632 (1986)
*   Phys. Rev. Lett. 61, 2879 (1988)
*   Phys. Rev. B 37, 6991 (1988)
*   Phys. Rev. B 38, 9902 (1988)
*   Phys. Rev. B 39, 5566 (1989)

The `tersoff` module is constructed by including several specialized source files that provide parameter definitions, the derived data type for the Tersoff potential, specific functional forms for interactions, high-level module procedures, the generic BOP kernel, and registration routines.

## Key Components

### Modules

*   `tersoff`: This module encapsulates the Tersoff BOP. It provides the necessary types, parameters, and subroutines to calculate energies and forces for atomic systems interacting via this potential.
    *   **Functionality**: Handles initialization, binding to particle systems, computation of energy and forces, and parameter registration.
    *   **Construction**: Built by including:
        *   `tersoff_params.f90`: Defines parameters for various Tersoff parameterizations.
        *   `tersoff_type.f90`: Defines `tersoff_t` (via `BOP_TYPE`) and `tersoff_db_t` (via `BOP_DB_TYPE`).
        *   `tersoff_module.f90`: Contains high-level module procedures (e.g., `tersoff_init`, `tersoff_energy_and_forces`).
        *   `../bop_kernel.f90`: The generic Bond-Order Potential kernel, used for the main computation loop.
        *   `tersoff_func.f90`: Provides Tersoff-specific implementations of functions like pair interactions, bond order, and angular terms.
        *   `tersoff_registry.f90`: Routines for registering the Tersoff potential with `ptrdict`.

### Classes/Types

*   `tersoff_t`: (Defined via `BOP_TYPE` in `tersoff_type.f90`) The primary derived data type for the Tersoff potential.
*   `tersoff_db_t`: (Defined via `BOP_DB_TYPE` in `tersoff_type.f90` using parameters from `tersoff_params.f90`) Type for Tersoff parameter database entries.

### Functions/Subroutines

The public interface is largely defined in `tersoff_module.f90` and aliased using preprocessor defines. Key procedures include:
*   `tersoff_init(...)`: Initializes the `tersoff_t` potential.
*   `tersoff_del(...)`: Destroys a `tersoff_t` object.
*   `tersoff_get_cutoff(...)`: A function to retrieve cutoff information, specific to this potential's interface.
*   `tersoff_bind_to(...)`: Binds the potential to a particle system.
*   `tersoff_energy_and_forces(...)`: Computes energies and forces.
*   `tersoff_register(...)`: Registers the potential.

### Preprocessor Definitions
*   `CUTOFF_T`: Set to `trig_off_t` (trigonometric cutoff functions).
*   `BOP_NAME`: `tersoff`
*   `BOP_NAME_STR`: `"tersoff"`
*   `BOP_STR`: `"Tersoff"`
*   `BOP_KERNEL`: `tersoff_kernel` (This implies the generic `bop_kernel.f90` is used, likely aliased to `tersoff_kernel`).
*   `BOP_TYPE`: `tersoff_t`
*   `BOP_DB_TYPE`: `tersoff_db_t`
*   `REGISTER_FUNC, INIT_FUNC, DEL_FUNC, GET_CUTOFF_FUNC, BIND_TO_FUNC, COMPUTE_FUNC`: Standardized names for interface functions.

## Important Variables/Constants
Parameters are defined in `tersoff_params.f90`. The `CUTOFF_T` definition specifies trigonometric cutoffs.

## Usage Examples
The Tersoff potential is used by creating an instance of `tersoff_t`, initializing it (often with parameters for a specific element like Si or Ge from `tersoff_params.f90`), and then using its interface functions.

```fortran
! Conceptual usage:
! USE tersoff_module ! (Assuming tersoff.f90 defines this module)
! TYPE(tersoff_t) :: my_tersoff_potential
! REAL(DP) :: cutoff_val
!
! CALL tersoff_init(my_tersoff_potential, ref_string="Si_Tersoff_1989")
! cutoff_val = tersoff_get_cutoff(my_tersoff_potential, ...) ! Usage depends on actual args
! CALL tersoff_bind_to(my_tersoff_potential, particles, neighbor_list)
! CALL tersoff_energy_and_forces(my_tersoff_potential, particles, neighbor_list, E, F, V)
! CALL tersoff_del(my_tersoff_potential)
```

## Dependencies and Interactions
*   **Internal Dependencies (within Atomistica project):**
    *   `supplib`: Supplementary utilities.
    *   `particles`: Provides `particles_t`.
    *   `neighbors`: Provides `neighbors_t`.
    *   The included `tersoff_*.f90` files and `../bop_kernel.f90`.
*   **External Libraries:** None explicitly listed.
*   **Interactions:**
    *   The `tersoff` module provides a complete implementation of the Tersoff potential.
    *   It uses the generic `bop_kernel.f90` for computations, specialized with Tersoff-specific functions from `tersoff_func.f90`.
    *   The `GET_CUTOFF_FUNC` suggests a more direct way to query cutoff information than some other potentials might offer.
```
