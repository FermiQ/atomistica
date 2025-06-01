# src/potentials/bop/rebo2/bop_kernel_rebo2.f90

## Overview

This file contains a specialized `BOP_KERNEL` subroutine, tailored for the second-generation REBO (Reactive Empirical Bond-Order) potential, primarily based on the work by Brenner et al., J. Phys. Cond. Mat. 14, 783 (2002). This kernel is responsible for calculating the potential energy, atomic forces, and virial stress tensor for systems interacting via the REBO2 potential.

Compared to a more generic BOP kernel, the REBO2 kernel incorporates significantly more complex physics, including:
*   Sophisticated bond order calculations that depend on local coordination numbers and bond angles.
*   Multi-body interactions, including dihedral (torsional) angle terms, which are crucial for describing the structure of hydrocarbons.
*   Optional screening functions to modulate interactions based on the wider atomic environment.
*   Terms dependent on the number and type of neighbors (`#ifdef NUM_NEIGHBORS`).

The kernel's behavior and the terms included in the calculation are heavily influenced by a range of preprocessor flags such as `DIHEDRAL`, `ALT_DIHEDRAL` (alternative dihedral formulation), `SCREENING`, and `NUM_NEIGHBORS`.

## Key Components

### Functions/Subroutines

*   `BOP_KERNEL(this, cell, maxnat, natloc, nat, r, ktyp, nebmax, nebavg, aptr, a2ptr, bptr, ptrmax, dc, epot, f_inout, wpot_inout, ...)`
    *   **Description**: This recursive subroutine is the core computational engine for the REBO2 potential. It iterates over atoms and their pre-calculated neighbor lists to compute interactions.
        1.  **Initialization**: Sets up internal neighbor list buffers if not already allocated or if sizes have changed. These buffers store detailed information about each bond and its local environment.
        2.  **Neighbor List Processing**: It first processes the raw neighbor lists provided by `aptr`, `a2ptr`, `bptr` to build more detailed internal lists (`this%neb`, `this%bndlen`, `this%bndnm`, etc.), applying cutoff functions (`fCin`, `fCar`, `fCbo`, `fCnc`) to determine interaction strengths and store relevant values. Screening logic (`#ifdef SCREENING`) is applied here to determine effective interaction cutoffs and identify screening neighbors.
        3.  **Coordination Dependent Terms** (`#ifdef NUM_NEIGHBORS`): If enabled, it pre-calculates coordination numbers (`this%nn`) and conjugation terms (`nconj`) for each atom, as these influence bond orders and energies.
        4.  **Main Calculation Loop**: Iterates over each atom `i` and its neighbors `j`.
            *   Calculates pair energies (attractive `VA`, repulsive `VR`).
            *   Calculates bond order `b_ij` and `b_ji`. This involves:
                *   Summing angular contributions (`g` function) and length-dependent contributions (`h` function) from neighbors `k` of `i` (excluding `j`) to get `zij`.
                *   Summing similar terms for neighbors `l` of `j` (excluding `i`) to get `zji`.
                *   If `NUM_NEIGHBORS` is active, terms like `Pij` (from `Pcc`, `Pch` tables/functions) are added to `zij` and `zji`.
                *   The final bond order also incorporates `Fij` (from `Fcc`, `Fhh`, `Fch` tables/functions) which depends on coordination and conjugation.
            *   The total energy for the i-j bond includes `VR + b_ave_ij * VA`, where `b_ave_ij` is an average of `bij` and `bji`.
            *   **Dihedral Term** (`#ifdef DIHEDRAL` or `#ifdef ALT_DIHEDRAL`): If enabled, it calculates torsional energy contributions. This involves finding sequences of four atoms (i-j-k-l or similar, depending on the dihedral style) and evaluating a torsional potential based on the angle between the i-j-k plane and j-k-l plane. This term often uses `Tcc` tables/functions.
            *   Forces and virial contributions from all these terms are calculated and accumulated.
        5.  **Screening Force Loop** (`#ifdef SCREENING`): A separate loop to add force contributions arising from the derivatives of screening functions with respect to the positions of screening atoms.
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE)`: The REBO2 potential object (`rebo2_t`) containing parameters and state.
        *   `ktyp(maxnat) :: integer, intent(in)`: Array of element types for each atom.
        *   Other arguments (`cell`, `maxnat`, `natloc`, `nat`, `r`, `nebmax`, `nebavg`, `aptr`, `a2ptr`, `bptr`, `ptrmax`, `dc`, `epot`, `f_inout`, `wpot_inout`, optional per-atom/bond terms, `ierror`) are similar to those in the generic BOP kernel. The LAMMPS version has slightly different arguments (e.g., `tag` instead of `cell` and `dc`).

### Modules Used
*   `tls`: For thread-local storage, essential for OpenMP parallelization.
*   `omp_lib`: If compiled with OpenMP.

### Classes/Types
*   Uses `BOP_TYPE` (which is `rebo2_t` in this context), expected to be defined in `rebo2_type.f90`.

## Important Variables/Constants

*   `typemax :: integer, parameter`: Set to 3, likely for C, H, and potentially a third element type if the parameterization supports it.
*   The kernel uses numerous internal arrays (many allocatable and stored in `this`) to cache neighbor information, bond properties, cutoff values, and intermediate terms for force calculations.
*   Preprocessor flags extensively control the compiled code paths:
    *   `SCREENING`: Enables environmental screening of interactions.
    *   `NUM_NEIGHBORS`: Enables terms dependent on coordination numbers (e.g., `Pij`, `Fij`).
    *   `DIHEDRAL` / `ALT_DIHEDRAL`: Enable different styles of torsional interaction terms.
    *   `LAMMPS` / `PYTHON`: Adapt the interface and some internal calculations for different simulation environments.

## Usage Examples

This is a low-level computational kernel. It is called by higher-level routines within the REBO2 potential module (e.g., by `rebo2_energy_and_forces` which is typically an alias for `COMPUTE_FUNC`).

## Dependencies and Interactions

*   **Internal Dependencies (within Atomistica project):**
    *   `BOP_TYPE` (`rebo2_t`): Relies heavily on the parameters and structure defined for the REBO2 potential.
    *   REBO2-specific functions (defined in `rebo2_func.f90`): `VA`, `VR`, `g`, `h`, `bo`, `fconj`, `eval` (for tables like `Pcc`, `Fcc`, `Tcc`).
    *   Cutoff functions (`fCin`, `fCar`, `fCbo`, `fCnc`) used for smoothly terminating interactions.
*   **External Libraries:** `omp_lib` if OpenMP is used.
*   **Interactions:**
    *   This kernel is the heart of the REBO2 potential calculation.
    *   It takes particle data and neighbor lists, computes interactions based on the complex REBO2 formalism, and returns energy, forces, and virials.
    *   Its internal logic is substantially more involved than a simple pairwise or basic BOP kernel due to the many-body nature and conditional terms of the REBO2 potential.
```
