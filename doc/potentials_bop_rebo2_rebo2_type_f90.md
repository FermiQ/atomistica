# src/potentials/bop/rebo2/rebo2_type.f90

## Overview

This file defines the `BOP_TYPE` derived data type, which is specifically structured for the second-generation REBO (REBO2) potential. When included by `rebo2.f90` or `rebo2_scr.f90`, `BOP_TYPE` is effectively aliased to `rebo2_t` or `rebo2_scr_t`, respectively. This is arguably the most complex potential type definition among the BOPs in Atomistica, reflecting the intricate nature of the REBO2 formalism.

The `BOP_TYPE` for REBO2 serves as a comprehensive container for:
*   A large set of scalar parameters defining pair interactions, angular terms, and bond order modifications for C-C, C-H, and H-H systems. Many of these have default values set directly in the type definition.
*   Parameters and computed coefficients for spline functions, particularly for the C-C angular `g` function.
*   Optionally, spline objects for pair potentials (`VA`, `VR`) and various cutoff functions if `SPLINE_POTENTIAL` or `SPLINE_CUTOFF` macros are defined. If not using splines, it stores objects of `CUTOFF_T` for cutoffs.
*   Objects for 2D and 3D lookup tables (`table2d_t`, `table3d_t`) that store correction terms like `Fcc`, `Pcc`, `Tcc`.
*   Raw input arrays for these tables (e.g., `in_Fcc`) which are populated from `rebo2_default_tables` and then processed into the table objects.
*   Allocatable arrays for internal neighbor lists and bond-specific data required by the `bop_kernel_rebo2.f90`.
*   Flags and counters for managing state (e.g., `tables_allocated`, `neighbor_list_allocated`).

The structure is also conditional on macros like `SCREENING` and `NUM_NEIGHBORS`. This file also declares public interfaces for standard potential operations and for the default table-loading routines.

## Key Components

### Constants
*   `rebo2_C_ = 1`, `rebo2_H_ = 3`: Integer parameters representing internal type indices for Carbon and Hydrogen.
*   `possible_elements(2) = (/ rebo2_C_, rebo2_H_ /)`: Array defining the elements this potential primarily handles.
*   `C_C = 1, C_H = 3, H_H = 6`: Integer parameters for pair type indexing.
*   `H_ = 1, C_ = 6`: Likely atomic numbers or alternative global type indices for Hydrogen and Carbon.

### Derived Data Types
*   `g_coeff_t`: A small helper type holding a `REAL(DP) :: c(6, 3)` array, used for storing spline coefficients for the C-C `g(theta)` function.

*   `BOP_TYPE` (Public)
    *   **Description**: The central data type for REBO2 potentials.
    *   **Selected Fields**:
        *   `elements :: character(MAX_EL_STR)`: String for element filtering (default "C,H").
        *   `els :: integer`: Integer representation of the element filter.
        *   `internal_el(:) :: integer, allocatable`: Array for mapped element types of atoms in the current simulation.
        *   **Scalar Parameters**: Numerous `REAL(DP)` parameters with default values for C-C, C-H, H-H interactions (e.g., `cc_B1, cc_beta1, ch_Q, hh_A, hhh_lambda, cc_re, ch_re, hh_re`) and cutoff radii (e.g., `cc_in_r1, ch_r1`).
        *   **Screening Parameters** (`#ifdef SCREENING`): `Cmin, Cmax`, `screening_threshold`, `dot_threshold`.
        *   **C-C g(theta) Spline Data**: `cc_g_theta(6)` (node points), `cc_g_g1(6), cc_g_dg1(6), cc_g_d2g1(6), cc_g_g2(6)` (target values/derivatives at nodes). `cc_g1_coeff, cc_g2_coeff :: type(g_coeff_t)` store the computed spline coefficients.
        *   **H-H g(theta) Data**: `SPGH(6,3)` (coefficients), `IGH(25)` (indices for polynomial selection).
        *   `with_dihedral :: logical(C_BOOL)`: Flag to enable/disable dihedral term calculations.
        *   **Derived Constants**: Arrays like `conalp, conear(6,6), conpe(3), cut_in_h(10), max_cut_sq(10)`, etc., populated during initialization.
        *   **Table Objects**: `Fcc, Fch, Fhh :: type(table3d_t)`; `Pcc, Pch :: type(table2d_t)`; `Tcc :: type(table3d_t)`.
        *   **Spline Objects** (conditional on `SPLINE_POTENTIAL`, `SPLINE_CUTOFF`): `spl_VA(10), spl_VR(10), spl_fCin(10)`, etc., of `simple_spline_t`. If `SPLINE_CUTOFF` is not defined, `spl_fCin` etc. are of `type(CUTOFF_T)`.
        *   **State Flags**: `neighbor_list_allocated, tables_allocated, it`.
        *   **Neighbor List Arrays**: Allocatable arrays for detailed neighbor information, similar to other BOPs but including specific arrays for REBO2 like `cutfcnnc`, `sfacnc` if `SCREENING` and `NUM_NEIGHBORS` are active.
        *   `zero_tables :: logical(C_BOOL)`: Input flag to initialize tables to zero.
        *   **Input Table Data Arrays**: `in_Fcc(0:4,0:4,0:9)`, `in_dFdi`, etc., used to hold raw table data before processing into table objects.

### Public Interfaces
*   Standard potential operations: `init`, `del`, `bind_to`, `energy_and_forces`, `register`.
*   Also makes public the table loading routines from `rebo2_default_tables`: `rebo2_default_Fcc_table`, `rebo2_default_Fch_table`, `rebo2_default_Fhh_table`, `rebo2_default_Pcc_table`, `rebo2_default_Pch_table`, `rebo2_default_Tcc_table`.

## Important Variables/Constants
*   The type definition itself embeds many default parameter values for the REBO2 potential for C and H.
*   `CUTOFF_T` (defined in `rebo2.f90` or `rebo2_scr.f90`) determines the type of cutoff objects if splined cutoffs are not used.
*   Various preprocessor flags (`SCREENING`, `NUM_NEIGHBORS`, `SPLINE_POTENTIAL`, `SPLINE_CUTOFF`, `LAMMPS`) control the inclusion of specific fields.

## Usage Examples
This file defines the `BOP_TYPE` for REBO2. Instances are managed by procedures in `rebo2_module.f90`.

## Dependencies and Interactions
*   **Internal Dependencies:** Uses `table2d_type`, `table3d_type` (from Atomistica's table utilities), `simple_spline_t` (if splines enabled), and `CUTOFF_T`. Also depends on constants like `MAX_EL_STR`, `MAX_Z`.
*   **External Libraries:** `iso_c_binding`.
*   **Interactions:**
    *   This is the central data structure for the REBO2 potential.
    *   It's populated by `INIT_FUNC` and `BIND_TO_FUNC` (defined in `rebo2_module.f90`), which in turn use routines from `rebo2_db.f90` (which uses `rebo2_default_tables.f90`).
    *   The `bop_kernel_rebo2.f90` reads from and uses the parameters, tables, and splines stored in instances of this type.
    *   The direct inclusion of many default parameters within the type definition is a notable feature.
```
