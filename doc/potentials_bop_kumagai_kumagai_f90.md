# src/potentials/bop/kumagai/kumagai.f90

## Overview

This file defines the `kumagai` Fortran module, which implements the Kumagai-Izumi-Hara-Sakai Bond-Order Potential (BOP). The specific functional forms and parameters are based on the work by Kumagai, T., Izumi, S., Hara, S., & Sakai, S. (2007), published in Computational Materials Science, 39(2), 457-464.

The `kumagai` module is structured by including several other source files using preprocessor directives. These includes bring together parameter definitions, the derived data type for the Kumagai potential, specific functional forms for interactions, high-level module procedures, the generic BOP kernel, and registration routines.

## Key Components

### Modules

*   `kumagai`: This module encapsulates the Kumagai BOP. It provides the necessary types, parameters, and subroutines to calculate energies and forces for atomic systems interacting via this specific BOP formulation.
    *   **Functionality**: The module handles the initialization of the potential (typically from parameter data), binding it to a particle system, and computing energy and forces.
    *   **Construction**: It is built by including the following files:
        *   `kumagai_params.f90`: Defines the parameters specific to the Kumagai potential.
        *   `kumagai_type.f90`: Defines the `kumagai_t` derived data type (via the `BOP_TYPE` macro) used to store the potential's state and parameters. It will also likely define `kumagai_db_t` (via `BOP_DB_TYPE`).
        *   `kumagai_module.f90`: Contains the high-level module procedures such as `kumagai_init`, `kumagai_del`, `kumagai_energy_and_forces`, etc., forming the public interface.
        *   `../bop_kernel.f90`: The generic Bond-Order Potential kernel, which performs the core energy and force calculations, specialized by Kumagai-specific functions.
        *   `kumagai_func.f90`: Provides the Kumagai-specific implementations of functions required by the BOP kernel (e.g., pair potential terms, bond order function, angular functions).
        *   `kumagai_registry.f90`: Contains subroutines to register the Kumagai potential with a configuration system (e.g., `ptrdict`).

### Classes/Types

*   `kumagai_t`: (Defined via `BOP_TYPE` in `kumagai_type.f90`) The primary derived data type for the Kumagai potential.
*   `kumagai_db_t`: (Expected to be defined via `BOP_DB_TYPE` in `kumagai_type.f90` or `kumagai_params.f90`) A type to hold database parameters for Kumagai potential sets.

### Functions/Subroutines

The public interface is largely defined in `kumagai_module.f90` and aliased using preprocessor defines. Key procedures typically include:
*   `kumagai_init(...)`: Initializes the `kumagai_t` potential type.
*   `kumagai_del(...)`: Destroys/deallocates a `kumagai_t` potential object.
*   `kumagai_bind_to(...)`: Binds the potential to a particle system and neighbor list handler.
*   `kumagai_energy_and_forces(...)`: Computes energies and forces.
*   `kumagai_register(...)`: Registers the potential.

### Preprocessor Definitions

This file uses several preprocessor definitions:
*   `CUTOFF_T`: Set to `trig_off_t`, indicating a trigonometric cutoff function scheme.
*   `BOP_NAME`: `kumagai`
*   `BOP_NAME_STR`: `"kumagai"` (String representation)
*   `BOP_STR`: `"Kumagai"` (User-facing string)
*   `BOP_KERNEL`: `kumagai_kernel`
*   `BOP_TYPE`: `kumagai_t`
*   `BOP_DB_TYPE`: `kumagai_db_t`
*   `REGISTER_FUNC, INIT_FUNC, DEL_FUNC, BIND_TO_FUNC, COMPUTE_FUNC`: Standardized names for the main interface functions, mapped to their Kumagai-specific implementations.

## Important Variables/Constants

Specific constants and parameters for the Kumagai potential are defined within the included `kumagai_params.f90` file.

## Usage Examples

This file itself does not contain direct usage examples. The Kumagai potential would be used by creating an instance of `kumagai_t`, initializing it, and then passing it to routines that calculate system properties.

```fortran
! Conceptual usage:
! USE kumagai_module ! (Assuming kumagai.f90 defines this module)
! TYPE(kumagai_t) :: my_kumagai_potential
!
! CALL kumagai_init(my_kumagai_potential, ...)
! CALL kumagai_bind_to(my_kumagai_potential, particles, neighbor_list)
! CALL kumagai_energy_and_forces(my_kumagai_potential, particles, neighbor_list, E, F, V)
! CALL kumagai_del(my_kumagai_potential)
```

## Dependencies and Interactions

*   **Internal Dependencies (within Atomistica project):**
    *   `supplib`: Supplementary library, likely for utilities.
    *   `particles`: Provides `particles_t` for atomic data.
    *   `neighbors`: Provides `neighbors_t` for neighbor lists.
    *   The included `.f90` files for `kumagai_*` and `../bop_kernel.f90` are essential components.
*   **External Libraries:** None explicitly listed.
*   **Interactions:**
    *   The `kumagai` module provides a complete implementation of the Kumagai BOP.
    *   It interacts with the simulation core by registering itself and providing the standard potential interface functions.
    *   It utilizes the generic `bop_kernel.f90` by specializing it with Kumagai-specific functions and parameters.
```
