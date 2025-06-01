# src/potentials/bop/brenner/brenner.f90

## Overview

This file defines the `brenner` Fortran module, which implements the Abell-Tersoff-Brenner type interatomic potential. It specifically uses Morse-style pair terms as originally proposed by Brenner. The documentation notes that this implementation does not include the more complex correction tables for treating pi-orbitals and other phenomena found in later Brenner potentials (e.g., REBO2).

The `brenner` module is constructed by including several other source files using preprocessor directives. These includes bring in parameter definitions, type definitions, specific potential functions, module procedures, and the generic bond-order potential (BOP) kernel.

## Key Components

### Modules

*   `brenner`: This module encapsulates the Brenner potential. It provides the necessary types, parameters, and subroutines to calculate energies and forces for atomic systems interacting via this potential.
    *   **Functionality**: The module is responsible for initializing the potential (often from parameter files), binding it to a set of particles, and computing energy and forces.
    *   **Construction**: It is built by including the following files:
        *   `brenner_params.f90`: Defines the parameters specific to the Brenner potential.
        *   `brenner_type.f90`: Defines the `brenner_t` derived data type used to store the potential's state and parameters.
        *   `brenner_module.f90`: Contains the high-level module procedures such as `brenner_init`, `brenner_del`, `brenner_energy_and_forces`, etc. These are the primary public interface of the module.
        *   `../bop_kernel.f90`: The generic Bond-Order Potential kernel, which performs the core energy and force calculations.
        *   `brenner_func.f90`: Provides the Brenner-specific implementations of functions required by the BOP kernel (e.g., pair potential terms `VA`, `VR`, bond order function `bo`, angular dependent function `g`).
        *   `brenner_registry.f90`: Contains subroutines to register the Brenner potential with a configuration system (e.g., `ptrdict`), allowing it to be selected and configured via input files.

### Classes/Types

*   `brenner_t`: (Defined in `brenner_type.f90`) The primary derived data type for the Brenner potential. It holds all parameters, spline coefficients, cutoff information, and pointers to potential-specific functions.
*   `brenner_db_t`: (Likely defined in `brenner_params.f90` or `brenner_type.f90`) A type to hold database parameters for different Brenner potential sets.

### Functions/Subroutines

The public interface is largely defined in `brenner_module.f90` and aliased using preprocessor defines. Key procedures typically include:
*   `brenner_init(...)`: Initializes the `brenner_t` potential type.
*   `brenner_del(...)`: Destroys/deallocates a `brenner_t` potential object.
*   `brenner_bind_to(...)`: Binds the potential to a particle system and neighbor list handler.
*   `brenner_energy_and_forces(...)`: Computes energies and forces.
*   `brenner_register(...)`: Registers the potential for use by the simulation engine.

### Preprocessor Definitions

This file uses several preprocessor definitions to configure the included files and define common names:
*   `CUTOFF_T`: Set to `trig_off_t`, indicating the type of cutoff function scheme used.
*   `BOP_NAME`: `brenner`
*   `BOP_NAME_STR`: `"brenner"` (String representation)
*   `BOP_STR`: `"Brenner"` (User-facing string)
*   `BOP_KERNEL`: `brenner_kernel` (Likely the name given to the specialized version of the generic BOP kernel for Brenner).
*   `BOP_TYPE`: `brenner_t`
*   `BOP_DB`: `brenner_db`
*   `BOP_DB_TYPE`: `brenner_db_t`
*   `REGISTER_FUNC`, `INIT_FUNC`, `DEL_FUNC`, `BIND_TO_FUNC`, `COMPUTE_FUNC`: Standardized names for the main interface functions, mapped to their Brenner-specific implementations.

## Important Variables/Constants

Specific constants and parameters for the Brenner potential are defined within the included `brenner_params.f90` file.

## Usage Examples

This file itself does not contain direct usage examples. The Brenner potential would be used by creating an instance of `brenner_t`, initializing it (often with data from a parameter file), and then passing it to routines that calculate system properties like energy and forces.

```fortran
! Conceptual usage:
! USE brenner
! TYPE(brenner_t) :: my_brenner_potential
!
! CALL brenner_init(my_brenner_potential, parameter_file_or_db = ...)
! ! ... set up particles (p) and neighbor lists (nl) ...
! CALL brenner_bind_to(my_brenner_potential, p, nl)
! ! ...
! CALL brenner_energy_and_forces(my_brenner_potential, p, nl, E_pot, forces_array, virial_tensor)
! ! ...
! CALL brenner_del(my_brenner_potential)
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
        *   `brenner_params.f90`
        *   `brenner_type.f90`
        *   `brenner_module.f90`
        *   `../bop_kernel.f90`
        *   `brenner_func.f90`
        *   `brenner_registry.f90`
*   **External Libraries:** None explicitly listed in this file.
*   **Interactions:**
    *   The `brenner` module provides a complete implementation of a specific Brenner potential.
    *   It interacts with the simulation core by registering itself and providing the standard potential interface functions (initialize, calculate energy/forces, cleanup).
    *   It utilizes the generic `bop_kernel.f90` by specializing it with Brenner-specific functions (pair interactions, bond order terms, angular corrections) and parameters.
```
