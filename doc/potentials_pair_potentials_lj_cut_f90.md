# src/potentials/pair_potentials/lj_cut.f90

## Overview

The `lj_cut` module implements the classic 12-6 Lennard-Jones (LJ) pair potential, a fundamental model for describing van der Waals forces (attraction and repulsion) between non-bonded atoms or molecules. The functional form of the Lennard-Jones potential is:

\[ V(r_{ij}) = 4\epsilon \left[ \left(\frac{\sigma}{r_{ij}}\right)^{12} - \left(\frac{\sigma}{r_{ij}}\right)^6 \right] \]

where:
*   \(r_{ij}\) is the distance between atoms \(i\) and \(j\).
*   \(\epsilon\) (epsilon) is the depth of the potential well, representing the strength of the attraction.
*   \(\sigma\) (sigma) is the finite distance at which the inter-particle potential is zero.

This implementation allows the potential to be applied between two specified element types. It includes a cutoff distance (`cutoff`) beyond which the interaction is considered zero. Additionally, an option to `shift` the potential is provided, which subtracts the value of the potential at the cutoff distance to ensure \(V(\text{cutoff}) = 0\).

## Key Components

### Modules

*   `lj_cut`
    *   **Uses**: `libAtoms_module`, `ptrdict`, `logging`, `timer`, `neighbors`, `particles`, `filter`.

### Data Types

*   `lj_cut_t` (Public)
    *   **Description**: Holds the parameters for the Lennard-Jones potential with a cutoff.
    *   **Fields**:
        *   `element1 :: character(MAX_EL_STR)`: The chemical symbol of the first element type for this interaction. Default: "*".
        *   `element2 :: character(MAX_EL_STR)`: The chemical symbol of the second element type. Default: "*".
        *   `el1 :: integer`: Internal filter integer corresponding to `element1`.
        *   `el2 :: integer`: Internal filter integer corresponding to `element2`.
        *   `epsilon :: real(DP)`: The depth of the potential well (\(\epsilon\)). Default: `1.0`.
        *   `sigma :: real(DP)`: The distance at which the potential is zero (\(\sigma\)). Default: `1.0`.
        *   `cutoff :: real(DP)`: The cutoff distance for the interaction. Default: `1.0` (Note: This default is very short for typical LJ interactions, which usually use cutoffs like \(2.5\sigma\) or larger).
        *   `shift :: logical(BOOL)`: A boolean flag to indicate whether the potential should be shifted so that its value is zero at the `cutoff` distance. Default: `.false.`.
        *   `offset :: real(DP)`: If `shift` is true, this stores the energy value at the cutoff distance: \(4\epsilon [(\sigma/\text{cutoff})^{12} - (\sigma/\text{cutoff})^6]\). This value is subtracted from the potential if shifting is enabled.

### Public Subroutines & Interfaces

*   `init(this)` (maps to `lj_cut_init`)
    *   **Description**: Constructor. This is an empty subroutine in the current implementation, as parameters are typically set via direct assignment or through the `ptrdict` registration mechanism before `bind_to` is called.
*   `del(this)` (maps to `lj_cut_del`)
    *   **Description**: Destructor. This is an empty subroutine as `lj_cut_t` does not dynamically allocate memory it owns.
*   `bind_to(this, p, nl, ierror)` (maps to `lj_cut_bind_to`)
    *   **Description**: Initializes internal settings based on particle data and parameters.
        1.  Converts `element1` and `element2` character strings into integer filter IDs (`this%el1`, `this%el2`).
        2.  Logs the parameters (`epsilon`, `sigma`, `cutoff`, `shift`).
        3.  Calls `request_interaction_range(nl, this%cutoff)` for all pairs of element types matching the filters.
        4.  If `this%shift` is true, calculates `this%offset = 4*this%epsilon*((this%sigma/this%cutoff)**12 - (this%sigma/this%cutoff)**6)`. Otherwise, `this%offset` is `0.0_DP`.
        5.  Logs the calculated `offset`.
*   `energy_and_forces(this, p, nl, epot, f, wpot, mask, epot_per_at, wpot_per_at, ierror)` (maps to `lj_cut_energy_and_forces`)
    *   **Description**: Calculates the potential energy, forces, and virial contributions from the Lennard-Jones interactions. This routine is OpenMP parallelized using thread-local storage (`tls`) for accumulation.
        1.  Iterates through each local atom `i` and its neighbors `j`, ensuring each pair is counted once (effectively `i < j` due to `i <= j` and `i == j` handling, though the primary loop is over `i` and its full neighbor list).
        2.  Applies element filters (`el1`, `el2`) and an optional `mask` to determine if the pair interacts. A `weight` factor (1 or 2) is determined based on whether atom `j` is also local and not masked, to correctly distribute energy for per-atom quantities.
        3.  If the distance `abs_dr` is less than `this%cutoff` (squared comparison `abs_dr < cut_sq`):
            *   Calculates \(s_6 = (\sigma/\text{abs_dr})^6\) and \(s_{12} = s_6^2\).
            *   Energy contribution: `en = 0.5*weight*(4*\epsilon*(s12 - s6) - this%offset)`. The outer `0.5*weight` is for distributing energy for `epot_per_at` like sums (via `tls_sca1`). The total `epot` is formed by summing `tls_sca1`.
            *   Force magnitude: `for_mag = 0.5*weight*24*\epsilon*(2*s12 - s6)/abs_dr`.
            *   Forces `df = for_mag * dr_vec/abs_dr` are applied to atoms `i` and `j`.
            *   The virial tensor `wpot` is updated.
            *   Per-atom energy and virial are also updated if corresponding arrays are present.
*   `register(this, cfg, m)` (maps to `lj_cut_register`)
    *   **Description**: Registers `element1`, `element2`, `epsilon`, `sigma`, `cutoff`, and `shift` with `ptrdict`.

## Important Variables/Constants

*   `epsilon`: Depth of the potential well (energy scale).
*   `sigma`: Finite distance where the unshifted potential is zero (length scale).
*   `cutoff`: Distance at which the interaction is truncated.
*   `shift`: Boolean flag to control energy shifting at the cutoff.
*   `offset`: The calculated energy shift if `shift` is true.
*   `element1`, `element2`: Define the specific pair of element types for this potential instance.

## Usage Examples

The Lennard-Jones potential is a widely used model for noble gases and as a non-bonded component in many force fields.

```fortran
! Conceptual usage:
! USE lj_cut_module
! TYPE(lj_cut_t) :: lj_potential
!
! ! Set parameters (e.g., for Argon-Argon interaction)
! lj_potential%element1 = "Ar"
! lj_potential%element2 = "Ar"
! lj_potential%epsilon = 0.0103_DP ! eV
! lj_potential%sigma = 3.40_DP    ! Angstrom
! lj_potential%cutoff = 8.5_DP    ! (approx 2.5*sigma)
! lj_potential%shift = .true.
!
! ! ... setup particles (p_data), neighbor_list_handler (nl_data) ...
! CALL lj_potential%bind_to(p_data, nl_data) ! Calculates offset
!
! ! Calculate energy and forces
! CALL lj_potential%energy_and_forces(p_data, nl_data, E_total, F_array, V_tensor)
```

## Dependencies and Interactions

*   Uses standard Atomistica modules: `libAtoms_module`, `ptrdict`, `logging`, `timer`, `neighbors`, `particles`, `filter`.
*   The interaction is pairwise.
*   The `shift` option ensures energy continuity at the cutoff. Force continuity is not guaranteed by this simple shift if the force is non-zero at cutoff.
```
