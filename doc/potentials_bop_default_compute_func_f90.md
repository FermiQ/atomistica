# src/potentials/bop/default_compute_func.f90

## Overview

This file provides a default implementation for the `COMPUTE_FUNC` subroutine, a standardized interface function used by various Bond-Order Potentials (BOPs) within the Atomistica simulation package. When a specific BOP module (e.g., `brenner.f90`, `tersoff.f90`) is compiled, its `COMPUTE_FUNC` macro (which defines the public energy/force calculation routine like `brenner_energy_and_forces`) typically points to this default implementation.

The primary responsibility of this subroutine is to orchestrate the calculation of potential energy, atomic forces, and the virial stress tensor for a given atomic configuration using the specified BOP. It performs necessary preparatory steps and then calls the appropriate low-level `BOP_KERNEL` to carry out the detailed computations.

## Key Components

### Functions/Subroutines

*   `COMPUTE_FUNC(this, p, nl, epot, f, wpot, mask, epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, ierror)`
    *   **Description**: This subroutine serves as the main entry point for computing system properties using a BOP.
        1.  **Timer Start**: It initiates a timer using `call timer_start(BOP_NAME_STR // "_force")` for performance profiling. `BOP_NAME_STR` is a macro providing the name of the specific BOP (e.g., "brenner").
        2.  **Neighbor List Update**: It ensures the neighbor list (`nl`) is up-to-date for the current particle configuration (`p`) by calling `call update(nl, p, ierror)`.
        3.  **Element Index Mapping**: It translates the global element types stored in `p%el` into internal element indices recognized by the potential. This mapping is stored in `this%Z2db` (populated during `BIND_TO_FUNC`) and the resulting internal indices are placed in a local array `el`.
        4.  **Neighbor Statistics**: It calculates `nebmax` (the maximum number of neighbors any single atom has) and `nebavg` (the average number of neighbors per atom). These statistics are passed to the `BOP_KERNEL`.
        5.  **Kernel Call**: It calls the `BOP_KERNEL` (a macro that resolves to the specific kernel implementation, e.g., `brenner_kernel`). The arguments passed to the kernel are adapted based on preprocessor flags:
            *   If `LAMMPS` is defined, a specific set of arguments is used.
            *   Otherwise, a different set is used, which may include `p%Abox` (simulation cell vectors) and `p%shear_dx` (if `PYTHON` is not defined).
        6.  **Error Check**: It propagates any error (`ierror`) from the kernel call.
        7.  **Timer Stop**: It stops the timer using `call timer_stop(BOP_NAME_STR // "_force")`.
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE), intent(inout)`: The BOP object containing parameters and internal state.
        *   `p :: type(particles_t), intent(inout)`: Structure containing particle data (positions, types, simulation cell information).
        *   `nl :: type(neighbors_t), intent(inout)`: Structure containing neighbor list information.
        *   `epot :: real(DP), intent(inout)`: Variable to accumulate the total potential energy of the system.
        *   `f(3, p%maxnatloc) :: real(DP), intent(inout)`: Array to store the calculated forces on each atom.
        *   `wpot(3, 3) :: real(DP), intent(inout)`: Matrix to store the calculated virial stress tensor.
        *   `mask(p%maxnatloc) :: integer, optional, intent(in)`: An optional array to mask out certain atoms from the calculation.
        *   `epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond`: Optional arrays for storing per-atom or per-bond contributions to energy, force, and virial. The shape of `wpot_per_at` and `wpot_per_bond` is conditional on `#ifdef LAMMPS`.
        *   `ierror :: integer, optional, intent(out)`: Error status flag.

## Important Variables/Constants

*   `BOP_NAME_STR`: A preprocessor macro defined in the specific BOP module (e.g., `brenner.f90`) that provides a string name for the potential (e.g., "brenner", "tersoff"). Used for timer labels.
*   `BOP_KERNEL`: A preprocessor macro that resolves to the name of the actual low-level computational kernel subroutine to be called (e.g., `brenner_kernel`).

## Usage Examples

This `COMPUTE_FUNC` is not typically called directly by an end-user. It is the subroutine that gets executed when a user calls the public energy/force computation routine of a specific BOP module (e.g., `call brenner_energy_and_forces(...)`).

```fortran
! Conceptual usage within a specific potential module (e.g., brenner.f90)
! MODULE brenner_module
!   ! ... other definitions ...
!   PUBLIC :: brenner_energy_and_forces
!   INTERFACE brenner_energy_and_forces
!     MODULE PROCEDURE COMPUTE_FUNC ! COMPUTE_FUNC is this default implementation
!   END INTERFACE
!   ! ...
! END MODULE brenner_module

! User code:
! USE brenner_module
! ! ... setup this_potential, particles, neighbor_list ...
! CALL brenner_energy_and_forces(this_potential, particles, neighbor_list, energy, forces, virial)
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on the data structures `BOP_TYPE`, `particles_t`, and `neighbors_t`.
    *   Calls `update(nl, p, ...)` to ensure neighbor lists are current before computation.
    *   The core of its work is dispatching the call to the `BOP_KERNEL`.
    *   Uses `timer_start` and `timer_stop` from a timing utility.
*   **External Libraries:** None explicitly used.
*   **Interactions:**
    *   This subroutine acts as an intermediary between the high-level request for an energy/force calculation and the low-level BOP kernel.
    *   It handles common preparatory tasks like updating neighbor lists and preparing element indices.
    *   The specific `BOP_KERNEL` it calls depends on how the `BOP_KERNEL` macro is defined in the context of the including BOP module (e.g., in `brenner.f90`, `BOP_KERNEL` is defined as `brenner_kernel`).
    *   Conditional compilation (`#ifdef LAMMPS`, `#ifndef PYTHON`) changes the argument list for the `BOP_KERNEL` call, allowing adaptation to different simulation environments or versions.
```
