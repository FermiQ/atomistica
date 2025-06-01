# src/potentials/bop/brenner/brenner_scr.f90

## Overview

This file defines the `brenner_scr` Fortran module, which implements a **screened** version of the Abell-Tersoff-Brenner type interatomic potential. Like the standard Brenner potential in `brenner.f90`, it uses Morse-style pair terms. The key difference is the activation of screening functions, which modify interactions based on the local environment. The documentation notes that this implementation, similar to the non-screened version, does not include the more complex correction tables for treating pi-orbitals.

The `brenner_scr` module is constructed by including several other source files using preprocessor directives. The crucial `SCREENING` macro is defined before these includes, which alters the behavior of the included components, particularly the `bop_kernel.f90`.

## Key Components

### Modules

*   `brenner_scr`: This module encapsulates the screened Brenner potential. It provides the necessary types, parameters, and subroutines to calculate energies and forces for atomic systems where screening effects are considered.
    *   **Functionality**: The module is responsible for initializing the potential, binding it to a set of particles, and computing energy and forces, with screening effects enabled.
    *   **Construction**: It is built by including the following files, processed with the `SCREENING` macro defined:
        *   `brenner_params.f90`: Defines the parameters specific to the Brenner potential. These parameters are used for both screened and non-screened versions.
        *   `brenner_type.f90`: Defines the `brenner_scr_t` (aliased from `BOP_TYPE`) derived data type. This type might have additional fields or different interpretations of fields when `SCREENING` is active.
        *   `brenner_module.f90`: Contains the high-level module procedures such as `brenner_scr_init`, `brenner_scr_del`, `brenner_scr_energy_and_forces`, etc.
        *   `../bop_kernel.f90`: The generic Bond-Order Potential kernel. With `SCREENING` defined, it includes logic to compute and apply screening functions.
        *   `brenner_func.f90`: Provides the Brenner-specific implementations of functions required by the BOP kernel. Some of these functions or their usage by the kernel might be modified by the `SCREENING` macro.
        *   `brenner_registry.f90`: Contains subroutines to register the screened Brenner potential with a configuration system.

### Classes/Types

*   `brenner_scr_t`: (Defined via `BOP_TYPE` in `brenner_type.f90`) The primary derived data type for the screened Brenner potential.
*   `brenner_db_scr_t`: (Defined via `BOP_DB_TYPE` in `brenner_type.f90` using parameters from `brenner_params.f90`) A type to hold database parameters.

### Functions/Subroutines

The public interface is largely defined in `brenner_module.f90` and aliased using preprocessor defines to `_scr_` versions. Key procedures include:
*   `brenner_scr_init(...)`: Initializes the `brenner_scr_t` potential type.
*   `brenner_scr_del(...)`: Destroys/deallocates a `brenner_scr_t` potential object.
*   `brenner_scr_bind_to(...)`: Binds the potential to a particle system and neighbor list handler, accounting for screening requirements (e.g., potentially larger cutoff distances).
*   `brenner_scr_energy_and_forces(...)`: Computes energies and forces, including screening effects.
*   `brenner_scr_register(...)`: Registers the screened potential for use by the simulation engine.

### Preprocessor Definitions

This file uses several preprocessor definitions:
*   `SCREENING`: This macro is defined to enable screening logic within the included files.
*   `CUTOFF_T`: Set to `exp_cutoff_t`, possibly an exponential cutoff scheme.
*   `BRENNER_MAX_REF, BRENNER_MAX_EL, BRENNER_MAX_PAIRS`: These are aliased (e.g., `BRENNER_SCR_MAX_REF`) before including `brenner_params.f90`. This might be for clarity or to allow `brenner_type.f90` to use different constant names if its structure changes significantly for the screened version, though typically `brenner_params.f90` provides the canonical definitions.
*   `BOP_NAME`: `brenner_scr`
*   `BOP_NAME_STR`: `"brenner_scr"`
*   `BOP_STR`: `"BrennerScr"` (User-facing string for this version)
*   `BOP_KERNEL`: `brenner_scr_kernel`
*   `BOP_TYPE`: `brenner_scr_t`
*   `BOP_DB`: `brenner_db_scr`
*   `BOP_DB_TYPE`: `brenner_db_scr_t`
*   `REGISTER_FUNC, INIT_FUNC, DEL_FUNC, BIND_TO_FUNC, COMPUTE_FUNC`: Standardized names mapped to their screened Brenner-specific implementations.

## Important Variables/Constants

Parameters are sourced from the included `brenner_params.f90`. The `SCREENING` definition is the most critical aspect of this file, fundamentally altering how interactions are computed by the kernel and other components.

## Usage Examples

This file itself does not contain direct usage examples. The screened Brenner potential would be used similarly to the non-screened version but would be selected using its specific name (e.g., "BrennerScr") if registered with the configuration system.

```fortran
! Conceptual usage:
! USE brenner_scr_module ! (Assuming brenner_scr.f90 defines a module making these available)
! TYPE(brenner_scr_t) :: my_screened_brenner_potential
!
! CALL brenner_scr_init(my_screened_brenner_potential, parameter_file_or_db = ...)
! ! ... setup particles (p) and neighbor lists (nl) ...
! CALL brenner_scr_bind_to(my_screened_brenner_potential, p, nl)
! ! ...
! CALL brenner_scr_energy_and_forces(my_screened_brenner_potential, p, nl, E_pot, forces_array, virial_tensor)
! ! ...
! CALL brenner_scr_del(my_screened_brenner_potential)
```

## Dependencies and Interactions

*   **Internal Dependencies (within Atomistica project):**
    *   Same base dependencies as `brenner.f90`: `libAtoms_module`, `ptrdict`, `logging`, `timer`, `particles`, `neighbors`.
    *   The included `.f90` files are compiled with `SCREENING` defined, altering their behavior.
*   **External Libraries:** None explicitly listed in this file.
*   **Interactions:**
    *   The `brenner_scr` module provides a complete implementation of a screened Brenner potential.
    *   The `SCREENING` macro enables specific code paths in `bop_kernel.f90` and potentially other included files, which handle the calculation and application of screening functions to the bond order and interaction energies. This typically involves considering more distant atoms or specific geometric configurations that can modify the strength of a bond.
    *   The `BIND_TO_FUNC` might request larger cutoff distances from the neighbor list handler to ensure all potentially screening atoms are found.
```
