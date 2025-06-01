# src/potentials/bop/default_bind_to_func.f90

## Overview

This file provides a default implementation of the `BIND_TO_FUNC` subroutine, which is a common interface function used by various Bond-Order Potentials (BOPs) within the Atomistica framework (e.g., Brenner, Tersoff). The `BIND_TO_FUNC` macro in specific potential modules (like `brenner.f90`) typically points to this default implementation if a more specialized one is not needed.

The primary role of this subroutine is to "bind" an initialized BOP object (of `BOP_TYPE`) to a specific set of atomic particles (`particles_t`) and a neighbor list handler (`neighbors_t`). This binding process involves several crucial setup steps that prepare the potential for subsequent energy and force calculations.

## Key Components

### Functions/Subroutines

*   `BIND_TO_FUNC(this, p, nl, ierror)`
    *   **Description**: This subroutine performs a series of setup and initialization tasks on the `BOP_TYPE` object `this` based on the provided particle data `p` and neighbor list configuration `nl`.
        1.  **Cleanup**: Calls `del(this)` to deallocate any existing internal data structures (like previously built neighbor lists) within the `this` object. This ensures a clean state, especially if the potential is being rebound to a new system or if particle/neighbor configurations change.
        2.  **Screening Parameter Initialization** (Conditional on `#ifdef SCREENING`):
            *   If screening is enabled during compilation, it copies screening parameters (`Cmin`, `Cmax`) from the potential's database (`this%db`) into the `this` object and calculates `dC = Cmax - Cmin`.
            *   It computes `this%C_dr_cut`, a factor used to extend the cutoff distance for finding atoms that might screen interactions. This factor depends on `this%Cmax`.
        3.  **Element Mapping**:
            *   Initializes `this%Z2db` (a mapping array from atomic number Z to the potential's internal element index) to -1.
            *   Iterates through the elements defined in the potential's database (`this%db%el`). For each element, it retrieves its atomic number using `atomic_number(a2s(element_symbol))` and stores the mapping in `this%Z2db(Z)`.
            *   Raises an error if any element symbol in the database is not recognized.
        4.  **Cutoff Function Initialization**:
            *   For each pair type defined by the potential, it initializes the primary cutoff function object `this%cut_in(i)` using the `r1` and `r2` parameters from `this%db`.
            *   It stores these cutoff limits directly in `this%cut_in_l(i)` (lower bound), `this%cut_in_h(i)` (upper bound), and `this%cut_in_h2(i)` (upper bound squared).
        5.  **Screening Cutoff Initialization** (Conditional on `#ifdef SCREENING`):
            *   Initializes additional cutoff function objects: `this%cut_out(i)` (for outer cutoff) and `this%cut_bo(i)` (for bond-order specific cutoff) using `or1/or2` and `bor1/bor2` parameters from `this%db`.
            *   Stores their respective lower and upper bounds.
            *   Calculates `this%max_cut_sq(i)`, the square of the maximum cutoff distance relevant for pair type `i`, considering `cut_in_h`, `cut_out_h`, and `cut_bo_h`.
        6.  **Request Interaction Ranges from Neighbor List Handler**:
            *   For every unique pair of element types present in the input particle data `p`:
                *   It determines the necessary cutoff distance for building neighbor lists.
                *   If screening is enabled: `cutoff = sqrt(this%C_dr_cut(pair_idx)) * sqrt(maxval(this%max_cut_sq))`. (Note: `maxval(this%max_cut_sq)` appears to use the global maximum cutoff, which might be a simplification).
                *   If screening is disabled: `cutoff = this%cut_in_h(pair_idx)` for that specific pair.
                *   Calls `request_interaction_range(nl, cutoff, element_type_i, element_type_j)` to inform the neighbor list handler `nl` of this required range.
                *   If compiled for LAMMPS (`#ifdef LAMMPS`), it also calls `set_interaction_range(p, 2*cutoff, element_type_i, element_type_j)`.
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE), intent(inout)`: The Bond-Order Potential object being configured.
        *   `p :: type(particles_t), intent(inout)`: The particle data (providing element types).
        *   `nl :: type(neighbors_t), intent(inout)`: The neighbor list handler object.
        *   `ierror :: integer, optional, intent(out)`: An error status flag.

## Important Variables/Constants

*   The subroutine extensively modifies the `this` object of `BOP_TYPE`, populating its fields related to element mapping, cutoff distances, and pre-calculated screening parameters.
*   It reads parameters from `this%db`, which is the database component of the `BOP_TYPE`.
*   The behavior is significantly altered by the presence of `#ifdef SCREENING` and `#ifdef LAMMPS` preprocessor directives.

## Usage Examples

This `BIND_TO_FUNC` is not called directly by a user setting up a simulation. Instead, it's a component used by specific BOP implementations. For example, when `brenner_bind_to` is called, it might internally call this default implementation if it's suitable for the Brenner potential.

```fortran
! Conceptual usage within a specific potential module (e.g., brenner_module.f90)
! SUBROUTINE brenner_bind_to(this_brenner, particles_data, neighbor_list_handler, err_code)
!   TYPE(brenner_t), INTENT(INOUT) :: this_brenner
!   TYPE(particles_t), INTENT(INOUT) :: particles_data
!   TYPE(neighbors_t), INTENT(INOUT) :: neighbor_list_handler
!   INTEGER, OPTIONAL, INTENT(OUT) :: err_code
!
!   ! ... Brenner-specific pre-processing ...
!
!   ! Call the default BIND_TO_FUNC, which is this implementation
!   CALL BIND_TO_FUNC(this_brenner, particles_data, neighbor_list_handler, err_code)
!
!   ! ... Brenner-specific post-processing ...
! END SUBROUTINE brenner_bind_to
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on the structure of `BOP_TYPE`, `particles_t`, and `neighbors_t`.
    *   Calls other procedures: `del(this)` (destructor of `BOP_TYPE`), `atomic_number`, `a2s` (array to string), `Z2pair` (element types to pair index), `request_interaction_range`, and `set_interaction_range` (if for LAMMPS).
    *   The `init` procedure of the specific `CUTOFF_T` type is called (e.g., `trig_off_init`).
*   **External Libraries:** None explicitly used.
*   **Interactions:**
    *   This routine is a critical step in preparing a BOP for calculations. It ensures that the potential object is correctly configured for the given atomic system and that the neighbor list handler is aware of the required interaction ranges.
    *   Errors during element mapping (e.g., unknown element symbols in the database) will result in `RAISE_ERROR` being called.
```
