# src/potentials/eam/tabulated_eam.f90

## Overview

This file defines the `tabulated_eam` module, which implements the Embedded Atom Method (EAM) potential using tabulated values for the embedding function, electron density, and pair potential. It provides the `tabulated_eam_t` type to store potential parameters and routines to calculate energies and forces.

## Key Components

### Modules

*   `tabulated_eam`: Implements the Tabulated EAM potential. It reads potential parameters from a file and provides functions for energy and force calculations.

### Classes/Types

*   `tabulated_eam_t`: Represents a Tabulated EAM potential.
    *   `elements :: character(MAX_EL_STR)`: Specifies the chemical symbols of the elements this potential applies to (e.g., "Al", "Ni,Cu"). Default is "*" (all elements).
    *   `els :: integer`: Internal integer representation of the element filter derived from `elements`.
    *   `fn :: character(100)`: The filename from which the EAM potential parameters (embedding function, density function, pair potential) are read. Defaults to "default.in".
    *   `comment :: character(100)`: A descriptive comment read from the potential file.
    *   `Z :: integer`: Atomic number of the element, read from the potential file.
    *   `mass :: real(DP)`: Atomic mass of the element, read from the potential file.
    *   `a0 :: real(DP)`: Lattice constant of the element in its ground state, read from the potential file.
    *   `lattice :: character(100)`: Ground-state lattice structure (e.g., "FCC"), read from the potential file.
    *   `fF :: type(simple_spline_t)`: A spline object representing the embedding function F(rho).
    *   `fZ :: type(simple_spline_t)`: A spline object representing the effective charge function Z(r) used for the repulsive pair potential part.
    *   `frho :: type(simple_spline_t)`: A spline object representing the electron density function rho(r).
    *   `cutoff :: real(DP)`: The cutoff radius for the interactions.

### Functions/Subroutines

*   `init(this, elements, fn, ierror)`: Constructor for `tabulated_eam_t`. Initializes the potential by reading parameters from the specified file `fn` for the given `elements`.
*   `del(this)`: Destructor for `tabulated_eam_t`. Frees resources associated with the splines.
*   `bind_to(this, p, nl, ierror)`: Binds the potential to a particle set `p` and neighbor list `nl`. It sets up the element filter and requests interaction ranges from the neighbor list handler.
*   `energy_and_forces(this, p, nl, epot, f, wpot, epot_per_at, wpot_per_at, ierror)`: Computes the total potential energy (`epot`), forces on atoms (`f`), and virial stress tensor (`wpot`). Optionally computes per-atom energy and virial. This is a wrapper that calls the kernel.
*   `register(this, cfg, m)`: Registers properties of the `tabulated_eam_t` potential (like `elements` and `fn`) with a configuration dictionary, likely for input parsing.
*   `energy_and_forces_kernel(this, p, nl, epot, f, wpot, maxneb, epot_per_at, wpot_per_at, ierror)`: The core computational routine that calculates EAM energy and forces. It iterates over atoms and their neighbors to compute embedding densities, embedding energies, and pair interactions.

## Important Variables/Constants

*   `MAX_EL_STR`: (Imported from `filter.inc` via `filter` module) Maximum length of the `elements` string.
*   Most important parameters are part of the `tabulated_eam_t` derived type.

## Usage Examples

The following example illustrates how to use the `tabulated_eam_t` potential (extracted from comments in the source file):

```fortran
! type(tabulated_eam_t)  :: pot
! ...
! ! Assuming Cleri_PRB_48_22_Al_Ag_Au is a character variable
! ! holding the path to a valid EAM potential file.
! character(len=200) :: Cleri_PRB_48_22_Al_Ag_Au
! Cleri_PRB_48_22_Al_Ag_Au = "path/to/your/eam_potential_file.eam"
!
! call init(pot, fn = Cleri_PRB_48_22_Al_Ag_Au)
! ! Further setup for particles (p) and neighbor lists (nl) would be needed here
! ! call bind_to(pot, p, nl)
! ! ...
! ! call energy_and_forces(pot, p, nl, epot, f, wpot)
! ! ...
! call del(pot)
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `supplib`: Provides supplementary utilities and file I/O (e.g., `fopen`, `read_line`, `fclose`).
    *   `particles`: Provides the `particles_t` type for storing atom data (positions, types, etc.).
    *   `neighbors`: Provides the `neighbors_t` type for managing neighbor lists required for force calculation.
    *   `filter`: Provides functionalities for filtering atoms by element type (e.g., `filter_from_string`, `IS_EL2`).
    *   `macros.inc`: Contains preprocessor macros, likely for error handling (`INIT_ERROR`, `PASS_ERROR`) and property assignment (`ASSIGN_PROPERTY`).
    *   `filter.inc`: Likely contains definitions related to atom filtering, such as `MAX_EL_STR`.
    *   `spline.inc`: Provides definitions and macros for spline interpolation (`simple_spline_t`, `SPLINE_INLINE`, etc.), crucial for evaluating tabulated functions.
*   **External Libraries:** None explicitly listed beyond standard Fortran. Potential use of MPI via macros if compiled with `_MP` flag, typical for parallel execution.
*   **Interactions:** The `tabulated_eam` module is a key component for performing molecular dynamics or statics simulations using EAM potentials. It interacts with:
    *   The main simulation engine by providing energy and force calculation routines.
    *   Particle data structures to get atomic positions and types.
    *   Neighbor list generation routines to efficiently find interacting atoms.
    *   Configuration systems via the `register` subroutine to allow setting potential parameters from input files.
    Its primary role is to compute the potential energy and forces for a given configuration of atoms according to the EAM formalism, using numerically tabulated functions.
```
