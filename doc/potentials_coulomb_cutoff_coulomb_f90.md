# src/potentials/coulomb/cutoff_coulomb.f90

## Overview

The `cutoff_coulomb` module provides a straightforward method for calculating Coulomb (electrostatic) interactions by applying a sharp cutoff at a specified distance. Interactions beyond this cutoff radius are abruptly set to zero, and no smoothing or tapering functions are applied at the cutoff boundary.

**Caution**: The module authors explicitly warn that this method should be used with care duep to the discontinuity introduced by the sharp cutoff, which can lead to energy conservation issues in molecular dynamics simulations and other artifacts. It's generally recommended for situations where this approximation is acceptable or as a baseline for comparison. The module also notes that it does not support MPI for parallel execution.

## Key Components

### Modules

*   `cutoff_coulomb`
    *   **Uses**: `supplib` (supplementary library), `particles` (for particle data), `neighbors` (for neighbor lists). It also implicitly uses `ptrdict` via the registration routine and `tls` utilities if OpenMP is active.

### Data Types

*   `cutoff_coulomb_t` (Public)
    *   **Description**: Holds the parameters for the cutoff Coulomb potential.
    *   **Fields**:
        *   `epsilon_r :: real(DP)`: The relative dielectric constant of the medium. Default: `1.0_DP`.
        *   `cutoff :: real(DP)`: The distance beyond which Coulomb interactions are truncated. Default: `10.0_DP`.

### Public Subroutines & Interfaces

The module provides a standard set of interface procedures for an Atomistica potential:

*   `init(this, cutoff, epsilon_r, error)` (maps to `cutoff_coulomb_init`)
    *   **Description**: Constructor. Initializes the `cutoff_coulomb_t` object with the specified cutoff radius and relative dielectric constant. It issues a warning about the sharp cutoff and raises an error if MPI is enabled.
*   `del(this)` (maps to `cutoff_coulomb_del`)
    *   **Description**: Destructor. This is currently an empty subroutine as `cutoff_coulomb_t` does not allocate memory that it needs to deallocate itself.
*   `bind_to(this, p, nl, ierror)` (maps to `cutoff_coulomb_bind_to`)
    *   **Description**: Binds the potential to a particle set (`p`) and a neighbor list handler (`nl`). It calls `request_interaction_range(nl, this%cutoff)` to inform the neighbor list system of the required interaction range.
*   `potential(this, p, nl, q, phi, ierror)` (maps to `cutoff_coulomb_potential`)
    *   **Description**: Calculates the electrostatic potential `phi` at each atomic site `i` due to all other charges `q(j)` within the cutoff radius. This routine is parallelized using OpenMP and employs thread-local storage (`tls`) for accumulating `phi`. It iterates `i < j` to sum contributions.
*   `energy_and_forces(this, p, nl, q, epot, f, wpot, error)` (maps to `cutoff_coulomb_energy_and_forces`)
    *   **Description**: Calculates the total electrostatic potential energy (`epot`), the forces (`f`) on each atom, and the virial tensor (`wpot`) due to Coulomb interactions within the cutoff radius. The charges `q` are provided as input. The loop structure iterates for `i <= j`. The `i == j` case contributes to `epot` and `wpot` with a `0.5_DP` factor, which is unusual for direct pair sums and might be related to a specific energy expression or self-interaction correction not fully detailed. Standard pair forces (`i /= j`) are accumulated symmetrically.
*   `register(this, cfg, m)` (maps to `cutoff_coulomb_register`)
    *   **Description**: Registers the `epsilon_r` and `cutoff` parameters with the `ptrdict` configuration system, allowing them to be set from an input file.

## Important Variables/Constants

*   `this%epsilon_r`: The relative dielectric constant.
*   `this%cutoff`: The cutoff radius for interactions.

## Usage Examples

This potential is a basic implementation for Coulomb interactions and should be used with an understanding of its limitations (sharp cutoff).

```fortran
! Conceptual usage:
! USE cutoff_coulomb_module ! Assuming cutoff_coulomb.f90 defines this
! TYPE(cutoff_coulomb_t) :: coulomb_pot
! REAL(DP), ALLOCATABLE :: charges(:), potential_phi(:), forces(:,:)
! REAL(DP) :: total_energy, virial_tensor(3,3)
!
! ! Initialize potential
! CALL coulomb_pot%init(cutoff=12.0_DP, epsilon_r=1.0_DP)
!
! ! ... setup particles (p_data), neighbor_list_handler (nl_data), and charges ...
!
! ! Bind potential
! CALL coulomb_pot%bind_to(p_data, nl_data)
!
! ! Calculate potential
! CALL coulomb_pot%potential(p_data, nl_data, charges, potential_phi)
!
! ! Calculate energy and forces
! total_energy = 0.0_DP
! forces = 0.0_DP
! virial_tensor = 0.0_DP
! CALL coulomb_pot%energy_and_forces(p_data, nl_data, charges, total_energy, forces, virial_tensor)
!
! ! Clean up
! CALL coulomb_pot%del()
```

## Dependencies and Interactions

*   **`supplib`**: Used for general utilities like `INIT_ERROR`, `PASS_ERROR`, `prscrlog`, `timer_start`, `timer_stop`, `CSTR`, `c_loc`.
*   **`particles`**: Provides `particles_t` for atomic data and related macros like `DISTJ_SQ`, `VEC3`.
*   **`neighbors`**: Provides `neighbors_t` for neighbor lists and `GET_NEIGHBOR`.
*   **`ptrdict`**: Used by `cutoff_coulomb_register` for parameter registration.
*   **OpenMP/TLS**: The `cutoff_coulomb_potential` routine uses OpenMP directives and thread-local storage (`tls_init`, `tls_sca1`, `tls_reduce`) for parallel computation of the electrostatic potential.
*   The potential explicitly checks for and disallows MPI usage in the `init` routine.
```
