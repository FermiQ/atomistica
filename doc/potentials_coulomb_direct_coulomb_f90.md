# src/potentials/coulomb/direct_coulomb.f90

## Overview

The `direct_coulomb` module provides a straightforward implementation for calculating Coulomb (electrostatic) interactions via direct summation over all unique pairs of atoms. This method has a computational complexity that scales as O(N^2), where N is the number of atoms. It is explicitly intended for **non-periodic systems** only.

**Important Considerations**:
*   **Performance**: Due to its N^2 scaling, this method can become computationally expensive for large systems.
*   **Boundary Conditions**: It does not support periodic boundary conditions. For periodic systems, methods like Ewald summation (e.g., PME) are typically required.
*   **MPI Incompatibility**: The module explicitly raises an error if used in an MPI (Message Passing Interface) parallel environment, indicating it's designed for serial or shared-memory (OpenMP) execution.

## Key Components

### Modules

*   `direct_coulomb`
    *   **Uses**: `supplib` (supplementary library), `particles` (for particle data), `neighbors` (for neighbor list type, though not used for pair iteration in computations). It also implicitly uses `ptrdict` via the registration routine and `tls` utilities if OpenMP is active.

### Data Types

*   `direct_coulomb_t` (Public)
    *   **Description**: Holds the parameters for the direct Coulomb potential.
    *   **Fields**:
        *   `epsilon_r :: real(DP)`: The relative dielectric constant of the medium. Default: `1.0_DP`.

### Public Subroutines & Interfaces

The module provides a standard set of interface procedures for an Atomistica potential:

*   `init(this, epsilon_r, error)` (maps to `direct_coulomb_init`)
    *   **Description**: Constructor. Initializes the `direct_coulomb_t` object with the specified relative dielectric constant. Raises an error if MPI is enabled.
*   `del(this)` (maps to `direct_coulomb_del`)
    *   **Description**: Destructor. This is currently an empty subroutine.
*   `bind_to(this, p, nl, ierror)` (maps to `direct_coulomb_bind_to`)
    *   **Description**: Binds the potential to a particle set (`p`) and a neighbor list handler (`nl`). However, this specific implementation is empty and does not utilize the neighbor list `nl` for its direct summation calculations.
*   `potential(this, p, nl, q, phi, ierror)` (maps to `direct_coulomb_potential`)
    *   **Description**: Calculates the electrostatic potential `phi` at each atomic site `i` due to all other charges `q(j)`. It performs a direct summation over all unique pairs (`j > i`). The neighbor list `nl` is passed as an argument but not used in this routine.
*   `energy_and_forces(this, p, nl, q, epot, f, wpot, error)` (maps to `direct_coulomb_energy_and_forces`)
    *   **Description**: Calculates the total electrostatic potential energy (`epot`), the forces (`f`) on each atom, and the virial tensor (`wpot`). It performs a direct summation over all unique pairs (`j > i`). The charges `q` are provided as input. This routine is parallelized using OpenMP, and thread-local storage (`tls`) is used for accumulating forces before reduction. The neighbor list `nl` is passed but not used for pair iteration.
*   `register(this, cfg, m)` (maps to `direct_coulomb_register`)
    *   **Description**: Registers the `epsilon_r` parameter with the `ptrdict` configuration system.

## Important Variables/Constants

*   `this%epsilon_r`: The relative dielectric constant, which scales the Coulomb interaction.

## Usage Examples

This potential is suitable for finite, non-periodic systems where the N^2 cost is manageable.

```fortran
! Conceptual usage:
! USE direct_coulomb_module ! Assuming direct_coulomb.f90 defines this
! TYPE(direct_coulomb_t) :: coulomb_pot
! REAL(DP), ALLOCATABLE :: charges(:), potential_phi(:), forces(:,:)
! REAL(DP) :: total_energy, virial_tensor(3,3)
!
! ! Initialize potential
! CALL coulomb_pot%init(epsilon_r=1.0_DP)
!
! ! ... setup particles (p_data), neighbor_list_handler (nl_data - though not used by direct_coulomb), and charges ...
!
! ! Bind potential (does nothing significant for direct_coulomb)
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

*   **`supplib`**: Used for general utilities like `INIT_ERROR`, `POS3`, `VEC3`, `outer_product`, `CSTR`, `c_loc`.
*   **`particles`**: Provides `particles_t` for atomic data.
*   **`neighbors`**: Provides `neighbors_t` type, though it's not used by this module for pair selection in its core computational routines.
*   **`ptrdict`**: Used by `direct_coulomb_register` for parameter registration.
*   **OpenMP/TLS**: The `direct_coulomb_energy_and_forces` routine uses OpenMP directives and thread-local storage (`tls_init`, `tls_reduce`) for parallel computation.
*   The potential explicitly checks for and disallows MPI usage in the `init` routine.
*   This module calculates interactions between all pairs of atoms directly, without using a cutoff (other than the implicit cutoff of a finite system size).
```
