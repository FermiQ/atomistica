# src/potentials/pair_potentials/r6.f90

## Overview

The `r6` module implements a pair potential with an \(r^{-6}\) dependence, characteristic of attractive London dispersion forces. The specific functional form used is:

\[ V(r_{ij}) = \frac{A}{(r_0 + r_{ij})^6} \]

where:
*   \(r_{ij}\) is the distance between atoms \(i\) and \(j\).
*   \(A\) is the interaction strength prefactor.
*   \(r_0\) is an offset parameter added to the interatomic distance. If \(r_0 = 0\), this reduces to the standard \(A/r^6\) form for dispersion.

This implementation allows the potential to be applied between two specified element types and includes a cutoff distance beyond which the interaction is zero. No energy shifting is explicitly mentioned, so the potential truncates directly to zero at the cutoff.

## Key Components

### Modules

*   `r6`
    *   **Uses**: `libAtoms_module`, `ptrdict`, `logging`, `timer`, `neighbors`, `particles`, `filter`.

### Data Types

*   `r6_t` (Public)
    *   **Description**: Holds the parameters for the \( (r_0+r)^{-6} \) potential.
    *   **Fields**:
        *   `element1 :: character(MAX_EL_STR)`: The chemical symbol of the first element type for this interaction. Default: "*".
        *   `element2 :: character(MAX_EL_STR)`: The chemical symbol of the second element type. Default: "*".
        *   `el1 :: integer`: Internal filter integer corresponding to `element1`.
        *   `el2 :: integer`: Internal filter integer corresponding to `element2`.
        *   `A :: real(DP)`: The prefactor \(A\) determining the strength of the interaction. Default: `1.0`.
        *   `r0 :: real(DP)`: The offset \(r_0\) added to the distance in the denominator. Default: `0.0`.
        *   `cutoff :: real(DP)`: The cutoff distance for the interaction. Default: `1.0`.

### Public Subroutines & Interfaces

*   **`init`**: No explicit `r6_init` subroutine is defined. Parameters are expected to be set directly into the `r6_t` variable (e.g., via `ptrdict` or direct assignment) before `bind_to` is called.
*   **`del`**: No explicit `r6_del` subroutine. The type does not contain allocatable members that it owns.
*   `bind_to(this, p, nl, ierror)` (maps to `r6_bind_to`)
    *   **Description**: Initializes internal settings based on particle data and parameters.
        1.  Converts `element1` and `element2` character strings into integer filter IDs (`this%el1`, `this%el2`).
        2.  Logs the parameters (`A`, `r0`, `cutoff`).
        3.  Calls `request_interaction_range(nl, this%cutoff)` for all pairs of element types matching the filters.
*   `energy_and_forces(this, p, nl, epot, f, wpot, epot_per_at, wpot_per_at, ierror)` (maps to `r6_energy_and_forces`)
    *   **Description**: Calculates the potential energy, forces, and virial contributions.
        1.  Iterates through each local atom `i` and its neighbors `j`, ensuring each pair is counted once (`i > j`).
        2.  Checks if the atom pair `(i,j)` matches the specified `element1` and `element2` filters.
        3.  If the distance `abs_dr` is less than `this%cutoff` (squared comparison `abs_dr < cut_sq`):
            *   Calculates the energy contribution: `en = this%A / (this%r0 + abs_dr)**6`.
            *   Calculates the magnitude of the force: `for_mag = 6 * en / (this%r0 + abs_dr)`.
            *   Adds `en` to the total potential energy `epot`.
            *   Forces `df = for_mag * dr_vec/abs_dr` are applied to atoms `i` and `j`.
            *   The virial tensor `wpot` is updated: `wpot = wpot - outer_product(dr_vec, df)`.
            *   Per-atom energy and virial are also updated if corresponding arrays are present (energy contribution is halved for per-atom accumulation).
*   `register(this, cfg, m)` (maps to `r6_register`)
    *   **Description**: Registers the parameters `element1`, `element2`, `A`, `r0`, and `cutoff` with the `ptrdict` configuration system.

## Important Variables/Constants

*   `A`: The interaction strength prefactor.
*   `r0`: The offset distance in the denominator. Setting \(r_0=0\) recovers the standard \(A/r^6\) form.
*   `cutoff`: The distance at which the interaction is truncated.
*   `element1`, `element2`: Define the specific pair of element types to which this potential instance applies.

## Usage Examples

This potential is typically used to model attractive London dispersion forces, often as part of a more comprehensive force field (e.g., the attractive term in a Lennard-Jones potential).

```fortran
! Conceptual usage:
! USE r6_module
! TYPE(r6_t) :: r6_potential
!
! ! Set parameters
! r6_potential%element1 = "Ar"
! r6_potential%element2 = "Ar"
! r6_potential%A = 100.0_DP ! Example C6 coefficient
! r6_potential%r0 = 0.0_DP    ! Standard r^-6
! r6_potential%cutoff = 10.0_DP
!
! ! ... setup particles (p_data), neighbor_list_handler (nl_data) ...
! CALL r6_potential%bind_to(p_data, nl_data)
!
! ! Calculate energy and forces
! CALL r6_potential%energy_and_forces(p_data, nl_data, E_total, F_array, V_tensor)
```

## Dependencies and Interactions

*   Uses standard Atomistica modules: `libAtoms_module`, `ptrdict`, `logging`, `timer`, `neighbors`, `particles`, `filter`.
*   The interaction is pairwise and applies only between atoms matching the specified `element1` and `element2`.
*   The potential is sharply truncated at the `cutoff` distance (no energy shifting is applied by this module).
```
