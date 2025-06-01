# src/potentials/bop/bop_kernel.f90

## Overview

This file contains the `BOP_KERNEL` subroutine, which serves as the computational heart for Bond-Order Potentials (BOPs) of the Tersoff-Brenner type within the Atomistica code. It is designed to work with various specific BOP parameterizations such as Erhart-Albe, Tersoff, and Brenner potentials. The kernel calculates potential energy, atomic forces, and virial stress tensor for a given atomic configuration and potential definition.

The implementation incorporates several advanced features and optimizations, including:
*   Algorithm based on D.W. Brenner's work (Phys. Rev. B 42, 9458 (1990) and Phys. Rev. B 46, 1948 (1990)).
*   Use of linked lists and pointers for efficient neighbor handling.
*   Pre-calculation of pairwise terms.
*   Screening functions for interactions (e.g., M. I. Baskes et al., Modelling Simul. Mater. Sci. Eng. 2, 505 (1994); L. Pastewka et al., Phys. Rev. B 78, 161402(R) (2008)).
*   Fortran 90 compliance.
*   OpenMP directives for parallel execution.
*   Adaptations for use within LAMMPS.

## Key Components

### Functions/Subroutines

*   `BOP_KERNEL(this, cell, maxnat, natloc, nat, r, el, nebmax, nebavg, aptr, a2ptr, bptr, ptrmax, dc, epot, f_inout, wpot_inout, mask, epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, ierror)`
    *   **Description**: This recursive subroutine is the core computational engine for bond-order potentials. It takes the potential definition (`this`), atomic configuration (`r`, `el`), neighbor lists (`aptr`, `bptr`, etc.), and computes the total potential energy (`epot`), forces on each atom (`f_inout`), and the system virial (`wpot_inout`). Optional arguments allow for calculation of per-atom and per-bond contributions. The behavior of the kernel is determined by the specific BOP implementation passed via the `this` argument (of `BOP_TYPE`).
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE), intent(inout)`: Structure containing the specific BOP parameters and functions.
        *   `cell :: real(DP), intent(in)`: Simulation cell vectors (not used in LAMMPS version).
        *   `natloc :: integer, intent(in)`: Number of atoms processed locally.
        *   `nat :: integer, intent(in)`: Total number of atoms in the system.
        *   `r :: real(DP), intent(inout)`: Atomic positions (input), potentially modified if integrating.
        *   `el :: integer, intent(in)`: Array of element types for each atom.
        *   `nebmax :: integer, intent(in)`: Maximum number of neighbors an atom can have.
        *   `aptr, a2ptr, bptr :: integer(NEIGHPTR_T)/integer, intent(in)`: Pointers and arrays defining neighbor lists.
        *   `epot :: real(DP), intent(inout)`: Accumulates the total potential energy.
        *   `f_inout :: real(DP), intent(inout)`: Accumulates forces on atoms.
        *   `wpot_inout :: real(DP), intent(inout)`: Accumulates the virial tensor.
        *   `ierror :: integer, optional, intent(inout)`: Error status flag.
    *   Note: Many arguments related to neighbor lists, per-atom/per-bond quantities, and LAMMPS-specific handling are present. The actual signature can vary based on preprocessor flags like `LAMMPS` or `PYTHON`.

### Modules

*   The file itself does not define a Fortran `MODULE`. It contains the `BOP_KERNEL` subroutine directly.
*   **Uses**: `tls` (Thread Local Storage Utilities), `omp_lib` (if compiled with OpenMP support).

### Classes/Types

*   The file does not define any new Fortran `TYPE`s. It uses `BOP_TYPE`, which is expected to be defined in other modules that provide specific BOP implementations (e.g., Brenner, Tersoff).

## Important Variables/Constants

*   `typemax :: integer, parameter`: Typically set to 3, likely indicating the maximum number of distinct element types the kernel is hardcoded to handle in some older parts of its logic or array dimensions.
*   Preprocessor Macros: The code extensively uses preprocessor macros (e.g., `SCREENING`, `LAMMPS`, `_OPENMP`, `NEB_TOO_SMALL`, `SNEB_TOO_SMALL`, `DCELL_INDEX`, `PARTIAL_SCREENING`) to control compilation and behavior for different environments and features. These are not variables but significantly affect the compiled code.

## Usage Examples

This file provides a low-level computational kernel. It is not typically called directly by a user but by other modules that implement specific bond-order potentials. A higher-level module would initialize a `BOP_TYPE` variable with the parameters and functions for, say, a Tersoff potential, and then pass this variable along with atom and neighbor data to `BOP_KERNEL`.

```fortran
! Conceptual example (actual usage is more complex and part of larger system)
! TYPE(BOP_TYPE) :: my_tersoff_potential
! ! ... initialize my_tersoff_potential with Tersoff parameters and functions ...
!
! ! ... setup atomic coordinates (r), element types (el), neighbor lists etc. ...
!
! CALL BOP_KERNEL(my_tersoff_potential, cell_vectors, &
!                 max_atoms, local_atoms, total_atoms, r, el, &
!                 max_neighbors, avg_neighbors, neighbor_ptr1, neighbor_ptr2, bond_ptr, max_bonds, &
!                 displacement_vectors, &
!                 potential_energy, forces, virial_tensor, &
!                 ! ... other optional arguments ...
!                 error_status)
```

## Dependencies and Interactions

*   **Internal Dependencies (within Atomistica project):**
    *   `BOP_TYPE`: Relies heavily on the structure and content of `BOP_TYPE` which encapsulates the specific potential's parameters (e.g., cutoff radii, interaction strengths) and function pointers (e.g., for pair energies, bond order calculations, angular functions). These are provided by modules like `brenner_type.f90`, `tersoff_type.f90`, etc.
    *   Various potential-specific functions (e.g., `VA`, `VR`, `bo`, `g`, `h`, `fCin`, `fCar`, `fCbo`, `Z2pair`) that are part of the `BOP_TYPE` interface.
    *   `tls` module for thread-local storage management, crucial for OpenMP parallelism.
*   **External Libraries:**
    *   `omp_lib`: If compiled with OpenMP, it uses the OpenMP library for parallel execution.
*   **Interactions:**
    *   `BOP_KERNEL` is called by higher-level potential routines.
    *   It reads atomic positions, element types, and neighbor list information.
    *   It computes energies, forces, and virials and updates the corresponding output variables.
    *   The kernel's internal logic for neighbor list management and memory allocation for these lists is significant.
    *   Its behavior is highly conditional based on preprocessor flags, allowing it to be adapted for different potential features (like screening) and simulation environments (like LAMMPS).
```
