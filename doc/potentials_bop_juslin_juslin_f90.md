# src/potentials/bop/juslin/juslin.f90

## Overview

This file defines the `juslin` Fortran module, which implements the W-C-H (Tungsten-Carbon-Hydrogen) Bond-Order Potential (BOP) as described by Juslin, Erhart, Traskelin, Nord, Henriksson, Nordlund, Salonen, Albe, in J. Appl. Phys. 98, 123520 (2005). This potential is designed to model interactions in systems containing tungsten, carbon, and hydrogen.

The `juslin` module is constructed by including several other source files using preprocessor directives. These includes bring in parameter definitions, the specific derived data type for the Juslin potential, functional forms for interactions, module procedures for the public interface, the generic BOP kernel, and registration routines.

## Key Components

### Modules

*   `juslin`: This module encapsulates the Juslin W-C-H interatomic potential. It provides the necessary types, parameters, and subroutines to calculate energies and forces for atomic systems interacting via this potential.
    *   **Functionality**: The module is responsible for initializing the potential (often from parameter files or internal databases), binding it to a set of particles, and computing energy and forces.
    *   **Construction**: It is built by including the following files:
        *   `juslin_params.f90`: Defines the parameters specific to the Juslin W-C-H potential for W, C, H, and their cross-interactions.
        *   `juslin_type.f90`: Defines the `juslin_t` derived data type (via `BOP_TYPE` macro) used to store the potential's state and parameters.
        *   `juslin_module.f90`: Contains the high-level module procedures such as `juslin_init`, `juslin_del`, `juslin_energy_and_forces`, etc. These form the primary public interface.
        *   `../bop_kernel.f90`: The generic Bond-Order Potential kernel, which performs the core energy and force calculations, adapted for the Juslin potential through its specific functions.
        *   `juslin_func.f90`: Provides the Juslin-specific implementations of functions required by the BOP kernel (e.g., pair potential terms, bond order function, angular function).
        *   `juslin_registry.f90`: Contains subroutines to register the Juslin potential with a configuration system (e.g., `ptrdict`), allowing it to be selected and configured via input files.

### Classes/Types

*   `juslin_t`: (Defined via `BOP_TYPE` in `juslin_type.f90`) The primary derived data type for the Juslin potential.
*   `juslin_db_t`: (Defined via `BOP_DB_TYPE` in `juslin_type.f90` using parameters from `juslin_params.f90`) A type to hold database parameters for Juslin potential sets.

### Functions/Subroutines

The public interface is largely defined in `juslin_module.f90` and aliased using preprocessor defines. Key procedures typically include:
*   `juslin_init(...)`: Initializes the `juslin_t` potential type.
*   `juslin_del(...)`: Destroys/deallocates a `juslin_t` potential object.
*   `juslin_bind_to(...)`: Binds the potential to a particle system and neighbor list handler.
*   `juslin_energy_and_forces(...)`: Computes energies and forces.
*   `juslin_force(...)`: Potentially a variant or component of the force calculation.
*   `juslin_register(...)`: Registers the potential for use by the simulation engine.

### Preprocessor Definitions

This file uses several preprocessor definitions to configure the included files and define common names:
*   `BOP_NAME`: `juslin_bop`
*   `BOP_NAME_STR`: `"juslin"` (String representation)
*   `BOP_STR`: `"Juslin"` (User-facing string)
*   `BOP_KERNEL`: `juslin_kernel`
*   `BOP_TYPE`: `juslin_t`
*   `BOP_DB`: `juslin_db`
*   `BOP_DB_TYPE`: `juslin_db_t`
*   `REGISTER_FUNC, INIT_FUNC, DEL_FUNC, BIND_TO_FUNC, COMPUTE_FUNC, FORCE_FUNC`: Standardized names for the main interface functions, mapped to their Juslin-specific implementations.

## Important Variables/Constants

Specific constants and parameters for the Juslin potential are defined within the included `juslin_params.f90` file.

## Usage Examples

This file itself does not contain direct usage examples. The Juslin potential would be used by creating an instance of `juslin_t`, initializing it (often with data from a parameter file or the internal database), and then passing it to routines that calculate system properties like energy and forces.

```fortran
! Conceptual usage:
! USE juslin_module ! (Assuming juslin.f90 defines a module making these available)
! TYPE(juslin_t) :: my_juslin_potential
!
! CALL juslin_init(my_juslin_potential, parameter_file_or_db = ...)
! ! ... set up particles (p) and neighbor lists (nl) ...
! CALL juslin_bind_to(my_juslin_potential, p, nl)
! ! ...
! CALL juslin_energy_and_forces(my_juslin_potential, p, nl, E_pot, forces_array, virial_tensor)
! ! ...
! CALL juslin_del(my_juslin_potential)
```

## Dependencies and Interactions

*   **Internal Dependencies (within Atomistica project):**
    *   `libAtoms_module`: Core library functionalities.
    *   `ptrdict`: For dictionary-based configuration and registration.
    *   `logging`: For logging messages.
    *   `timer`: For performance timing.
    *   `particles`: Provides `particles_t` for atomic data.
    *   `neighbors`: Provides `neighbors_t` for neighbor lists.
    *   The included `.f90` files are essential components:
        *   `juslin_params.f90`
        *   `juslin_type.f90`
        *   `juslin_module.f90`
        *   `../bop_kernel.f90`
        *   `juslin_func.f90`
        *   `juslin_registry.f90`
*   **External Libraries:** None explicitly listed in this file.
*   **Interactions:**
    *   The `juslin` module provides a complete implementation of the Juslin W-C-H potential.
    *   It interacts with the simulation core by registering itself and providing the standard potential interface functions.
    *   It utilizes the generic `bop_kernel.f90` by specializing it with Juslin-specific functions (pair interactions, bond order terms, angular corrections) and parameters.
```
