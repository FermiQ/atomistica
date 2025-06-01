# src/potentials/eam/tabulated_alloy_eam.f90

## Overview

The `tabulated_alloy_eam` module implements the Embedded Atom Method (EAM) for multi-component systems (alloys). This formulation is based on the work by Foiles, Baskes, and Daw (Phys. Rev. B 33, 7983 (1986)). The potential relies on tabulated functions for:
1.  The embedding energy of an atom `i` as a function of the host electron density at its site: \(F_i(\rho)\).
2.  The contribution of an atom `j` to the electron density at the site of atom `i`: \(\rho_j(r_{ij})\).
3.  A pair potential interaction between atoms `i` and `j`: \(\phi_{ij}(r_{ij})\).

These functions are read from an EAM potential file, typically in the "setfl" (single element EAM) or "fs" (Finnis-Sinclair style, often for alloys) format. The module uses spline interpolation to evaluate these functions and their derivatives.

## Key Components

### Modules

*   `tabulated_alloy_eam`
    *   **Uses**: `supplib` (for file I/O, logging, constants), `particles`, `neighbors`, `filter`.
    *   **Provides**: The `tabulated_alloy_eam_t` type and standard potential interface routines.

### Constants
*   `MAX_EAM_ELS :: integer, parameter`: Maximum number of distinct element types that can be defined in the EAM potential file. Value: 10.

### Data Types

*   `tabulated_alloy_eam_t` (Public)
    *   **Description**: Holds parameters, spline objects, and state for the tabulated alloy EAM potential.
    *   **Fields**:
        *   `elements :: character(MAX_EL_STR)`: A string to filter which simulation elements this potential applies to (default "*").
        *   `els :: integer`: The compiled element filter.
        *   `fn :: character(100)`: Filename of the EAM potential file (e.g., "alloy.eam.fs"). Default: "default.in".
        *   `dump :: logical(BOOL)`: If true, dumps the parsed spline functions to output files. Default: `.false.`.
        *   `db_nel :: integer`: Number of element types specified in the EAM file header.
        *   `db_elements(MAX_EAM_ELS) :: character(MAX_EL_STR)`: Names of the elements as read from the EAM file header.
        *   `el2db(:) :: integer, allocatable`: Maps the internal element indices of the simulation (from `particles_t`) to the element indices used in the EAM file (0 to `db_nel-1` or 1 to `db_nel`).
        *   `fF(:) :: type(simple_spline_t), allocatable`: Array of spline objects for the embedding functions \(F_i(\rho)\) for each element type `i` defined in the EAM file.
        *   `fphi(:,:) :: type(simple_spline_t), allocatable`: 2D array of spline objects for the pair potential functions \(\phi_{ij}(r)\) for each pair of element types `(i,j)`.
        *   `frho(:) :: type(simple_spline_t), allocatable`: Array of spline objects for the atomic electron density functions \(\rho_j(r)\) for each element type `j`.
        *   `cutoff :: real(DP)`: The cutoff radius for all interactions, read from the EAM file.

### Public Subroutines & Interfaces

*   `init(this, elements, fn, ierror)` (maps to `tabulated_alloy_eam_init`):
    *   **Description**: Constructor. Reads the EAM potential file specified by `fn`.
        1.  Parses the header of the EAM file to determine `db_nel`, `db_elements`, and global parameters like `nF` (number of points for F), `dF` (delta rho for F), `nr` (number of points for \(\rho\) and \(\phi\)), `dr` (delta r), and `cutoff`.
        2.  Allocates `fF`, `frho` (size `db_nel`) and `fphi` (size `db_nel x db_nel`).
        3.  For each element type in the EAM file: reads its atomic number, mass, lattice constant, lattice type, then reads and creates splines for its embedding function \(F_i(\rho)\) and electron density function \(\rho_j(r)\).
        4.  For each unique pair of elements `(i,j)` in the EAM file: reads and creates a spline for the pair potential \(\phi_{ij}(r)\). The pair potential spline `this%fphi(i,j)` is scaled by `0.5` after reading (this usually means the EAM file contains \(2 \phi_{ij}(r)\) or \(r \phi_{ij}(r)\) and the factor of 0.5 is to correct for double counting in the pairwise sum later).
        5.  If `this%dump` is true, writes the parsed splines to `.out` files.
        6.  Handles `#ifdef AVOID_SQRT` by calling `square_x_axis` on density and pair potential splines if defined (modifies splines to be indexed by \(r^2\) instead of \(r\)).
*   `del(this)` (maps to `tabulated_alloy_eam_del`): Destructor. Deallocates all spline arrays (`fF`, `fphi`, `frho`).
*   `bind_to(this, p, nl, ierror)`: Binds the potential to particle data `p` and neighbor list `nl`.
    1.  Sets the element filter `this%els`.
    2.  Allocates `this%el2db` and creates a mapping from simulation element types (via `p%el2Z` and `ElementName`) to the indices of elements found in the EAM file (`this%db_elements`).
    3.  Calls `request_interaction_range(nl, this%cutoff)` for all interacting pairs of types.
*   `energy_and_forces(this, p, nl, epot, f, wpot, mask, epot_per_at, wpot_per_at, ierror)`: Wrapper for `energy_and_forces_kernel`. Updates neighbor lists, determines `maxneb` (maximum number of neighbors for buffer allocation in kernel).
*   `energy_and_forces_kernel(this, p, nl, epot, f, wpot, maxneb, ...)`: The core computational routine.
    1.  **Density Calculation**: For each atom `i`, calculates the total host electron density \(\rho_i\) by summing contributions \(\rho_j(r_{ij})\) from its neighbors `j`. `rho_j(r_{ij})` is obtained from `func(this%frho(dbj), abs_dr)`.
    2.  **Embedding Energy**: Adds \(F_{dbi}(\rho_i)\) to the total energy for atom `i`. `F` is from `f_and_df(this%fF(dbi), rho, ...)`.
    3.  **Pair Potential Energy**: For each neighbor `j` of atom `i`:
        *   Evaluates the pair term \(\Phi = \text{func(this%fphi(dbi, dbj), abs_dr)}\).
        *   Adds \(\Phi / \text{abs_dr}\) to the per-atom energy `tls_sca1(i)`. (Note: This implies the EAM file's pair potential section, after the 0.5 scaling in `init`, effectively stores \(r \phi(r)\) or that the spline itself represents \(r \phi(r)\). The division by `abs_dr` then yields \(\phi(r)\) for the energy sum. This is a common convention for some EAM file formats.)
    4.  **Forces**: Calculates forces on atoms `i` and `j` from the derivatives of the embedding energy (\(dF/d\rho\)), density functions (\(d\rho/dr\)), and pair potentials (\(d\phi/dr\)).
    5.  Accumulates total energy `epot` and virial `wpot`. OpenMP parallelized.
*   `register(this, cfg, m)`: Registers `elements`, `fn` (filename), and `dump` flag with `ptrdict`.

## Important Variables/Constants
*   `MAX_EAM_ELS`: Maximum number of elements supported in an EAM file by this module.
*   The interpretation of the EAM file format is crucial, especially for the pair potential term (\(\phi_{ij}\)). The code scales the pair potential spline by 0.5 upon reading and then in the kernel, adds `spline_value / r` to the energy. This is consistent if the EAM file stores \(r \phi_{ij}(r)\) for the pair term, making the effective pair energy contribution \(0.5 \phi_{ij}(r)\).

## Usage Examples
This module is used for simulations of metallic alloys where EAM is an appropriate model. The EAM parameters are read from an external file.

```fortran
! Conceptual usage:
! USE tabulated_alloy_eam_module
! TYPE(tabulated_alloy_eam_t) :: alloy_pot
!
! CALL alloy_pot%init(fn="my_alloy.eam.fs", elements="Ni,Al")
! ! ... setup particles (p_data), neighbor_list_handler (nl_data) ...
! CALL alloy_pot%bind_to(p_data, nl_data)
!
! ! Calculate energy, forces, virial
! CALL alloy_pot%energy_and_forces(p_data, nl_data, E_total, F_array, V_tensor)
!
! CALL alloy_pot%del()
```

## Dependencies and Interactions
*   **`supplib`**: For file I/O, string manipulation, logging, and physical constants.
*   **`particles`, `neighbors`, `filter`**: For core simulation data structures and utilities.
*   **`simple_spline_t`**: The module relies on a spline library (presumably part of Atomistica or `supplib`) for creating, evaluating (`func`, `f_and_df`, `dfunc`), and managing tabulated functions.
*   The specific EAM file format ("setfl" or "fs" type) dictates how data is read and interpreted.
```
