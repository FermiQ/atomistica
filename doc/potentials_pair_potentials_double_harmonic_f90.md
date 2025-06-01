# src/potentials/pair_potentials/double_harmonic.f90

## Overview

The `double_harmonic` module implements a pair potential that combines two distinct harmonic (spring-like) interactions. The potential switches from one harmonic form to the other at a midpoint distance \(r_m\). This allows for modeling more complex bond behavior than a single harmonic well, such as having different stiffness constants for compression versus moderate extension, or for different interaction regimes before a final cutoff.

The potential energy \(V(r_{ij})\) for a pair of atoms \(i\) and \(j\) at distance \(r_{ij}\) is defined as:
*   If \(r_{ij} < r_m\): \(V(r_{ij}) = \frac{1}{2} k_1 (r_1 - r_{ij})^2\)
*   If \(r_m \le r_{ij} < \text{cutoff}\): \(V(r_{ij}) = \frac{1}{2} k_2 (r_2 - r_{ij})^2\)
*   If \(r_{ij} \ge \text{cutoff}\): \(V(r_{ij}) = 0\)

where \(k_1, r_1\) are the spring constant and equilibrium distance for the first regime, \(k_2, r_2\) are for the second regime, and \(r_m = (r_1 + r_2) / 2\) is the automatically calculated transition point.

## Key Components

### Modules

*   `double_harmonic`
    *   **Uses**: `libAtoms_module`, `ptrdict`, `logging`, `timer`, `neighbors`, `particles`, `filter`.

### Data Types

*   `double_harmonic_t` (Public)
    *   **Description**: Holds the parameters for the double harmonic potential.
    *   **Fields**:
        *   `element1 :: character(MAX_EL_STR)`: The chemical symbol of the first element type for this interaction. Default: "*".
        *   `element2 :: character(MAX_EL_STR)`: The chemical symbol of the second element type. Default: "*".
        *   `el1 :: integer`: Internal filter integer for `element1`.
        *   `el2 :: integer`: Internal filter integer for `element2`.
        *   `k1 :: real(DP)`: Spring constant for the first harmonic well (for \(r < r_m\)). Default: `1.0`.
        *   `r1 :: real(DP)`: Equilibrium distance for the first harmonic well. Default: `1.0`.
        *   `k2 :: real(DP)`: Spring constant for the second harmonic well (for \(r_m \le r < \text{cutoff}\)). Default: `1.0`.
        *   `r2 :: real(DP)`: Equilibrium distance for the second harmonic well. Default: `1.2`.
        *   `cutoff :: real(DP)`: The cutoff distance beyond which the interaction is zero. Default: `1.5`.
        *   `rm :: real(DP)`: The midpoint distance \((r_1 + r_2) / 2\), automatically calculated in `bind_to`, where the potential switches from the (k1,r1) regime to the (k2,r2) regime.

### Public Subroutines & Interfaces

*   **`init`**: No explicit `double_harmonic_init` subroutine is defined in the provided file. Parameters are expected to be set directly into the `double_harmonic_t` variable, often after it's created and before `bind_to` is called (e.g., via `ptrdict` if registered).
*   **`del`**: No explicit `double_harmonic_del` subroutine. The type does not contain allocatable members that it owns and needs to deallocate.
*   `bind_to(this, p, nl, ierror)` (maps to `double_harmonic_bind_to`)
    *   **Description**: Initializes internal settings based on particle data and parameters.
        1.  Converts `element1` and `element2` strings into integer filter IDs (`this%el1`, `this%el2`).
        2.  Logs the parameters (`k1, r1, k2, r2`).
        3.  Calls `request_interaction_range(nl, this%cutoff)` for all pairs of element types matching the filters.
        4.  Calculates the midpoint distance `this%rm = (this%r1 + this%r2) / 2.0_DP`.
*   `energy_and_forces(this, p, nl, epot, f, wpot, epot_per_at, wpot_per_at, ierror)` (maps to `double_harmonic_energy_and_forces`)
    *   **Description**: Calculates the potential energy, forces, and virial contributions.
        1.  Iterates through each local atom `i` and its neighbors `j`.
        2.  Checks if the atom pair `(i,j)` matches the specified `element1` and `element2` filters.
        3.  If the distance `abs_dr` is less than `this%cutoff`:
            *   If `abs_dr < this%rm`, the potential and force are calculated using `k1` and `r1`. Energy \(en = 0.5 k_1 (r_1 - \text{abs_dr})^2\), force magnitude \(F = k_1 (r_1 - \text{abs_dr})\).
            *   Else (if `abs_dr >= this%rm`), uses `k2` and `r2`. Energy \(en = 0.5 k_2 (r_2 - \text{abs_dr})^2\), force magnitude \(F = k_2 (r_2 - \text{abs_dr})\).
            *   The energy contribution `0.5 * en` is added to `epot` (factor of 0.5 for pair counting).
            *   Forces are applied to atoms `i` and `j`.
            *   The virial tensor `wpot` is updated: `wpot = wpot - outer_product(dr, 0.5 * F * dr_vec/abs_dr)`.
            *   Per-atom energy and virial are also updated if corresponding arrays are present.
*   `register(this, cfg, m)` (maps to `double_harmonic_register`)
    *   **Description**: Registers the parameters `element1`, `element2`, `k1`, `r1`, `k2`, `r2`, and `cutoff` with the `ptrdict` configuration system.

## Important Variables/Constants

*   `k1, r1`: Spring constant and equilibrium distance for the first (typically inner) harmonic well.
*   `k2, r2`: Spring constant and equilibrium distance for the second (typically outer) harmonic well.
*   `rm`: The distance at which the potential switches from the first to the second harmonic regime.
*   `cutoff`: The distance beyond which the potential is zero.
*   `element1`, `element2`: Specify the pair of element types to which this potential applies.

## Usage Examples

This potential can be used to model bonds where the stiffness or equilibrium length changes, or to create a more complex well shape than a single harmonic potential allows.

```fortran
! Conceptual usage:
! USE double_harmonic_module
! TYPE(double_harmonic_t) :: dh_potential
!
! ! Set parameters (e.g., directly or via ptrdict through registration)
! dh_potential%element1 = "A"
! dh_potential%element2 = "B"
! dh_potential%k1 = 10.0_DP ; dh_potential%r1 = 1.5_DP
! dh_potential%k2 = 5.0_DP  ; dh_potential%r2 = 2.0_DP
! dh_potential%cutoff = 2.5_DP
!
! ! ... setup particles (p_data), neighbor_list_handler (nl_data) ...
! CALL dh_potential%bind_to(p_data, nl_data) ! This calculates rm
!
! ! Calculate energy and forces
! CALL dh_potential%energy_and_forces(p_data, nl_data, E_total, F_array, V_tensor)
```

## Dependencies and Interactions

*   Uses standard Atomistica modules: `libAtoms_module`, `ptrdict`, `logging`, `timer`, `neighbors`, `particles`, `filter`.
*   The interaction is pairwise and only applies between atoms matching the specified `element1` and `element2`.
*   The potential is continuous at `rm` only if \(\frac{1}{2} k_1 (r_1 - r_m)^2 = \frac{1}{2} k_2 (r_2 - r_m)^2\). The force may be discontinuous at `rm` if \(k_1 (r_1 - r_m) \neq k_2 (r_2 - r_m)\).
```
