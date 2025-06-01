# src/potentials/bop/rebo2/rebo2_module.f90

## Overview

This file provides the core module procedures that define the public Application Programming Interface (API) for the second-generation REBO (REBO2) potential within the Atomistica framework. These subroutines manage the lifecycle (initialization, destruction), system binding, and computation calls for the REBO2 potential.

The procedures included are:
*   `INIT_FUNC`: A basic constructor. For REBO2, the more substantial parameter setup (including table processing) is often deferred to `BIND_TO_FUNC` or an `INIT_DEFAULT_FUNC`.
*   `DEL_FUNC`: The destructor, responsible for deallocating memory.
*   `BIND_TO_FUNC`: Connects the potential to the particle system and neighbor lists, and triggers the detailed parameter processing and setup of internal tables and splines.
*   `COMPUTE_FUNC`: The main routine that initiates the calculation of potential energy, forces, and virial by calling the specialized REBO2 kernel.

The names `INIT_FUNC`, `DEL_FUNC`, etc., are preprocessor macros defined in `rebo2.f90` (or `rebo2_scr.f90`) to point to REBO2-specific implementations (e.g., `rebo2_init`, `rebo2_bind_to`).

## Key Components

### Functions/Subroutines

*   `INIT_FUNC(this)`
    *   **Description**: This constructor for the `rebo2_t` (BOP_TYPE) potential object has a minimal role in this specific REBO2 implementation. Its primary function here is to check if `this%zero_tables` is true. If so, it initializes the raw input table arrays within the `this` object (e.g., `this%in_Fcc`, `this%in_Pcc`) to zero. The comprehensive setup of parameters, including processing these tables into functional objects (like splines or `table3d_type`), is typically handled later, either by `INIT_DEFAULT_FUNC` (which calls `rebo2_db_init` to load default tables) or within `BIND_TO_FUNC` when `rebo2_db_init_with_parameters` is invoked.
    *   **Key Arguments**: `this :: type(BOP_TYPE), intent(inout)`: The REBO2 potential object.

*   `DEL_FUNC(this)`
    *   **Description**: This is the destructor for the `rebo2_t` object.
        1.  It calls `rebo2_db_del(this)`, which is a specialized routine (from `rebo2_db.f90`) to deallocate all data structures related to the REBO2 parameterization. This includes 2D/3D lookup tables, spline objects, and also the internal neighbor list arrays managed by `BOP_TYPE`.
        2.  It additionally deallocates `this%internal_el` if it was allocated.
    *   **Key Arguments**: `this :: type(BOP_TYPE), intent(inout)`.

*   `BIND_TO_FUNC(this, p, nl, ierror)`
    *   **Description**: This crucial subroutine "binds" an initialized REBO2 potential object (`this`) to particle data (`p`) and a neighbor list handler (`nl`).
        1.  **Filtering**: Sets `this%els` based on `this%elements` string and particle data `p` using `filter_from_string`.
        2.  **Parameter Initialization**: Calls `rebo2_db_init_with_parameters(this, this%in_Fcc, ...)` using the raw table data stored in `this%in_Fcc`, etc. (which should have been populated by a prior call to `rebo2_db_init` or set up externally). This call performs the main processing of REBO2 parameters, including setting up derived constants, initializing cutoff function parameters, internal 2D/3D table objects, and splines, as detailed in `rebo2_db.f90`.
        3.  **Cutoff Setup for Neighbor Lists**: Determines the specific cutoff distances (`c_cc`, `c_ch`, `c_hh`) for C-C, C-H, and H-H interactions based on the base cutoff parameters from `this` and screening factors if `#ifdef SCREENING` is active.
        4.  **Request Interaction Ranges**: Calls `request_interaction_range(nl, ...)` for each relevant element pair (C-C, C-H, H-H) to inform the neighbor list handler.
        5.  **Request Border**: Calls `request_border(p, ...)` to set an appropriate simulation cell border/skin size, likely based on the maximum C-C interaction range.
        6.  **Allocate Internal Element Array**: Allocates `this%internal_el` to store mapped element types for local atoms.
    *   **Key Arguments**: `this`, `p`, `nl`, `ierror`.

*   `COMPUTE_FUNC(this, p, nl, epot, f, wpot, ..., ierror)`
    *   **Description**: This is the main routine for calculating energy, forces, and virial for the REBO2 potential.
        1.  **Timer & Neighbor Update**: Starts a timer and calls `update(nl, p, ...)` for the neighbor lists.
        2.  **Internal Element Mapping**: Reallocates `this%internal_el` if the number of local atoms `p%maxnatloc` has changed. It then populates `this%internal_el` by mapping the element types from `p%el` (filtered by `this%els`) to internal REBO2 integer types (`rebo2_C_`, `rebo2_H_`). Atoms not matching the filter or not being C/H get a 0 or negative type.
        3.  **Neighbor Statistics**: Calculates `nebmax` (maximum neighbors) and `nebavg` (average neighbors).
        4.  **Kernel Call**: Calls the specialized `BOP_KERNEL` (which is `rebo2_kernel` from `bop_kernel_rebo2.f90`), passing the prepared arguments including `this%internal_el`.
        5.  **Timer Stop**: Stops the timer.
    *   **Key Arguments**: `this`, `p`, `nl`, `epot`, `f`, `wpot`, optional per-atom/per-bond terms, `ierror`.

## Important Variables/Constants
*   `C_`, `H_`: Integer parameters (likely from `particles_module` or similar) representing atomic numbers or global types for Carbon and Hydrogen, used for checking `p%el2Z`.
*   `rebo2_C_`, `rebo2_H_`: Integer parameters (likely from `rebo2_type.f90`) representing the internal type indices for Carbon and Hydrogen within the REBO2 formalism.

## Usage Examples
The subroutines in this file constitute the public API of the REBO2 potential module.

```fortran
! Conceptual sequence for using the REBO2 potential:
! USE rebo2_module ! (Assuming rebo2.f90 defines this module)
! TYPE(rebo2_t) :: pot_rebo2
!
! ! 1. Initialize the REBO2 potential (often calls rebo2_db_init internally)
! CALL rebo2_init_default(pot_rebo2) ! Or a more specific init
!
! ! 2. Setup atomic configuration in 'atoms_config'
! ! ...
!
! ! 3. Initialize neighbor list handler 'neighbor_handler_config'
! ! ...
!
! ! 4. Bind the potential to the system
! CALL rebo2_bind_to(pot_rebo2, atoms_config, neighbor_handler_config)
!
! ! 5. Compute energy and forces
! CALL rebo2_energy_and_forces(pot_rebo2, atoms_config, neighbor_handler_config, &
!                               total_potential_energy, forces, virial_tensor)
!
! ! 6. Clean up
! CALL rebo2_del(pot_rebo2)
```

## Dependencies and Interactions
*   **Internal Dependencies:**
    *   Types: `BOP_TYPE` (`rebo2_t`), `particles_t`, `neighbors_t`.
    *   Relies on `rebo2_db_init_with_parameters` and `rebo2_db_del` from `rebo2_db.f90` for detailed parameter setup and teardown.
    *   Uses `filter_from_string` (from `filter` module), `request_interaction_range` and `request_border` (from `particles` or `neighbors` module).
    *   The `COMPUTE_FUNC` calls `BOP_KERNEL` (defined as `rebo2_kernel` in `rebo2.f90`, implemented in `bop_kernel_rebo2.f90`).
*   **External Libraries:** None explicitly used in this file beyond Atomistica core modules.
*   **Interactions:**
    *   `INIT_FUNC` performs minimal setup; the bulk of parameter processing and internal data structure initialization (tables, splines) is triggered by `BIND_TO_FUNC` via `rebo2_db_init_with_parameters`.
    *   `DEL_FUNC` ensures comprehensive cleanup of REBO2 specific data and general BOP allocated arrays.
    *   `COMPUTE_FUNC` prepares REBO2-specific element type mapping before calling the specialized kernel.
```
