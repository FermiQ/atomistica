# src/potentials/pair_potentials/harmonic.f90

## Overview

The `harmonic` module implements a simple harmonic pair potential, often referred to as a spring-like interaction. This potential is commonly used to model bonds between atoms or effective interactions in coarse-grained models. The functional form of the harmonic potential is:

\[ V(r_{ij}) = \frac{1}{2} k (r_0 - r_{ij})^2 \]

where \(r_{ij}\) is the distance between atoms \(i\) and \(j\), \(k\) is the spring constant, and \(r_0\) is the equilibrium bond distance.

This implementation allows the potential to be applied between two specified element types. It includes a cutoff distance and an option to shift the entire potential such that its value is zero at the cutoff distance.

## Key Components

### Modules

*   `harmonic`
    *   **Uses**: `libAtoms_module`, `ptrdict`, `logging`, `timer`, `neighbors`, `particles`, `filter`.

### Data Types

*   `harmonic_t` (Public)
    *   **Description**: Holds the parameters for the harmonic potential.
    *   **Fields**:
        *   `element1 :: character(MAX_EL_STR)`: The chemical symbol of the first element type for this interaction. Default: "*".
        *   `element2 :: character(MAX_EL_STR)`: The chemical symbol of the second element type. Default: "*".
        *   `el1 :: integer`: Internal filter integer corresponding to `element1`.
        *   `el2 :: integer`: Internal filter integer corresponding to `element2`.
        *   `k :: real(DP)`: The spring constant of the harmonic interaction. Default: `1.0`.
        *   `r0 :: real(DP)`: The equilibrium distance for the harmonic spring. Default: `1.0`.
        *   `cutoff :: real(DP)`: The cutoff distance beyond which the interaction is zero. Default: `1.5`.
        *   `shift :: logical(BOOL)`: A boolean flag indicating whether to shift the potential energy to be zero at the `cutoff`. If true, an `offset` is subtracted. Default: `.false.`.
        *   `offset :: real(DP)`: The energy offset value, calculated as \(0.5 k (\text{cutoff} - r_0)^2\) if `shift` is true, otherwise 0.

### Public Subroutines & Interfaces

*   **`init`**: No explicit `harmonic_init` subroutine is defined in this file. Parameters are typically set directly into the `harmonic_t` variable (e.g., via `ptrdict` after registration, or direct assignment) before `bind_to` is called.
*   **`del`**: No explicit `harmonic_del` subroutine. The type does not contain allocatable members that it owns and would need to deallocate.
*   `bind_to(this, p, nl, ierror)` (maps to `harmonic_bind_to`)
    *   **Description**: Initializes internal settings based on particle data and parameters.
        1.  Converts `element1` and `element2` character strings into integer filter IDs (`this%el1`, `this%el2`) using `filter_from_string`.
        2.  Logs the parameters (`k`, `r0`, `cutoff`, `shift`).
        3.  Calls `request_interaction_range(nl, this%cutoff)` for all pairs of element types matching the filters.
        4.  If `this%shift` is true, calculates `this%offset = 0.5_DP * this%k * (this%cutoff - this%r0)**2`. Otherwise, `this%offset` is `0.0_DP`.
*   `energy_and_forces(this, p, nl, epot, f, wpot, epot_per_at, wpot_per_at, ierror)` (maps to `harmonic_energy_and_forces`)
    *   **Description**: Calculates the potential energy, forces, and virial contributions from the harmonic interactions.
        1.  Iterates through each local atom `i` and its neighbors `j`, ensuring each pair is counted once (`i > j`).
        2.  Checks if the atom pair `(i,j)` matches the specified `element1` and `element2` filters.
        3.  If the distance `abs_dr` is less than `this%cutoff` (squared comparison `abs_dr < cut_sq`):
            *   Calculates the energy contribution: `en = 0.5 * k * (r0 - abs_dr)^2 - this%offset`.
            *   Calculates the magnitude of the force: `for_mag = k * (r0 - abs_dr)`.
            *   Adds `en` to the total potential energy `epot`.
            *   Forces `df = for_mag * dr_vec/abs_dr` are applied to atoms `i` and `j`.
            *   The virial tensor `wpot` is updated: `wpot = wpot - outer_product(dr_vec, df)`.
            *   Per-atom energy and virial are also updated if corresponding arrays are present (energy contribution is halved for per-atom accumulation).
*   `register(this, cfg, m)` (maps to `harmonic_register`)
    *   **Description**: Registers the parameters `element1`, `element2`, `k`, `r0`, `cutoff`, and `shift` with the `ptrdict` configuration system.

## Important Variables/Constants

*   `k`: The spring constant.
*   `r0`: The equilibrium bond distance.
*   `cutoff`: The distance at which the interaction is truncated.
*   `shift`: Boolean flag to control whether the potential is shifted to zero at the cutoff.
*   `offset`: The calculated energy shift if `shift` is true.
*   `element1`, `element2`: Define the specific pair of element types to which this potential instance applies.

## Usage Examples

The harmonic potential is a basic model for bonded interactions or effective spring-like forces.

```fortran
! Conceptual usage:
! USE harmonic_module
! TYPE(harmonic_t) :: harm_potential
!
! ! Set parameters (e.g., directly or via ptrdict through registration)
! harm_potential%element1 = "Si"
! harm_potential%element2 = "Si"
! harm_potential%k = 100.0_DP ! Spring constant
! harm_potential%r0 = 2.35_DP  ! Equilibrium distance
! harm_potential%cutoff = 3.0_DP
! harm_potential%shift = .true. ! Shift energy to zero at cutoff
!
! ! ... setup particles (p_data), neighbor_list_handler (nl_data) ...
! CALL harm_potential%bind_to(p_data, nl_data) ! This calculates the offset
!
! ! Calculate energy and forces
! CALL harm_potential%energy_and_forces(p_data, nl_data, E_total, F_array, V_tensor)
```

## Dependencies and Interactions

*   Uses standard Atomistica modules: `libAtoms_module`, `ptrdict`, `logging`, `timer`, `neighbors`, `particles`, `filter`.
*   The interaction is pairwise and applies only between atoms matching the specified `element1` and `element2`.
*   If `shift` is true, the potential energy is guaranteed to be zero at `r = cutoff`.
```
