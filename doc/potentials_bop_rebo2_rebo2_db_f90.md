# src/potentials/bop/rebo2/rebo2_db.f90

## Overview

This file contains subroutines crucial for the initialization and parameter setup of the second-generation REBO (REBO2) potential within the Atomistica framework. Unlike simpler parameter files that just list constants, this file manages the setup of more complex data structures used by REBO2, including 2D and 3D lookup tables for various interaction terms and, optionally, coefficients for spline-interpolated functions.

The main routines are:
*   `rebo2_db_init`: Initializes a REBO2 potential object with default parameters, primarily by fetching table data from `rebo2_default_tables`.
*   `rebo2_db_init_with_parameters`: Performs the detailed setup of the potential object using provided table data and other parameters already set in the `BOP_TYPE` instance. This includes calculating derived constants, setting up cutoff distances, initializing table objects, and preparing splines.
*   `rebo2_db_del`: A destructor that cleans up allocated tables, splines, and internal neighbor lists.
*   Helper routines for generating spline coefficients for specific parts of the potential.

The file also contains functions defining the analytical forms of pair potentials and cutoff functions, which are used if spline interpolation for these terms is enabled (`#ifdef SPLINE_POTENTIAL` or `#ifdef SPLINE_CUTOFF`).

## Key Components

### Initialization Subroutines

*   `rebo2_db_init(this)`
    *   **Description**: Initializes a `BOP_TYPE` object (`this`, which is `rebo2_t`) with default REBO2 parameters.
        1.  Declares local arrays to hold raw data for various 2D and 3D tables (e.g., `in_Fcc`, `in_Pcc`, `in_Tcc`).
        2.  Calls routines from the `rebo2_default_tables` module (e.g., `rebo2_default_Fcc_table`) to populate these local arrays with default numerical values. If the `ZERO_TABLES` macro is defined, these tables are initialized to zero instead.
        3.  Passes the populated local arrays and the `this` object to `rebo2_db_init_with_parameters` for the main setup.
    *   **Arguments**: `this :: type(BOP_TYPE)`: The REBO2 potential object to be initialized.

*   `rebo2_db_init_with_parameters(this, in_Fcc, in_dFdi, in_dFdj, in_dFdk, in_Fch, in_Fhh, in_Pcc, in_Pch, in_Tcc)`
    *   **Description**: This is the core routine for setting up the `rebo2_t` object. It uses parameters already set in `this` (e.g., from a parameter file or defaults in `rebo2_type.f90`) and the provided table data.
        1.  **Screening Setup**: If `#ifdef SCREENING` is active, calculates derived screening parameters like `this%dC` and `this%C_dr_cut`.
        2.  **Logging**: Prints many of the input and derived parameters if `ilog` (logging unit) is active.
        3.  **Bond Order Constants**: Calculates and stores precomputed factors for bond order calculations (e.g., `this%conpe`, `this%conan`).
        4.  **Penalty Function Constants**: Calculates `this%conear` terms used in bond order penalty functions.
        5.  **Cutoff Distance Setup**: Initializes various cutoff distance arrays within `this` (e.g., `this%cut_in_l`, `this%cut_in_h`, `this%cut_in_h2`, `this%cut_in_m`) for C-C, C-H, and H-H interactions based on base parameters like `this%cc_in_r1`, `this%cc_in_r2`. If screening is active, it sets up additional cutoff distance arrays (`cut_ar_`, `cut_bo_`, `cut_nc_`).
        6.  **Max Cutoff**: Determines `this%max_cut_sq` (maximum squared cutoff distance) for each pair type.
        7.  **C-C g-Spline**: Calls `rebo2_db_make_cc_g_spline(this)` to compute and store spline coefficients for the C-C angular dependent `g` function.
        8.  **Table Initialization**: Initializes 2D and 3D table objects within `this` (e.g., `this%Fcc`, `this%Pcc`, `this%Tcc`) using the input arrays `in_Fcc`, etc. These table objects are likely from `table2d_type` and `table3d_type`.
        9.  **Potential/Cutoff Splines**: Calls `rebo2_db_make_splines(this)` to generate spline fits for pair potentials (VA, VR) and cutoff functions if `#ifdef SPLINE_POTENTIAL` or `#ifdef SPLINE_CUTOFF` are defined.
        10. Sets `this%tables_allocated = .true.`.
    *   **Arguments**: `this :: type(BOP_TYPE), intent(inout)`; multiple `intent(in)` real arrays for table data.

### Destruction Subroutine

*   `rebo2_db_del(this)`
    *   **Description**: Deallocates resources held by the `rebo2_t` object.
        1.  Deallocates internal neighbor list arrays if `this%neighbor_list_allocated` is true (similar to `default_del_func.f90`).
        2.  Deallocates the 2D/3D table objects (`this%Fcc`, `this%Fch`, etc.) if `this%tables_allocated` is true.
        3.  Deallocates spline objects for potentials and cutoffs if they were allocated (conditional on `SPLINE_POTENTIAL`, `SPLINE_CUTOFF`).
        4.  Resets `this%tables_allocated` and `this%neighbor_list_allocated` to `.false.`.
    *   **Arguments**: `this :: type(BOP_TYPE), intent(inout)`.

### Spline Generation Subroutines

*   `rebo2_db_make_cc_g_spline(this, error)`
    *   **Description**: Calculates coefficients for a piecewise polynomial spline representation of the C-C angular `g` function. It sets up and solves a system of linear equations (using `gauss1`) based on known values and derivatives of the `g` function at specific angles (`this%cc_g_theta`, `this%cc_g_g1`, etc., which are parameters of the REBO2 potential). The resulting coefficients are stored in `this%cc_g1_coeff` and `this%cc_g2_coeff`.
*   `rebo2_db_make_splines(this)`
    *   **Description**: If `SPLINE_POTENTIAL` or `SPLINE_CUTOFF` are defined, this routine initializes spline objects (e.g., `this%spl_VA(C_C)`, `this%spl_fCin(C_C)`). It uses helper functions (like `cc_VA`, `cutoff_f`, also defined in this file) that provide the analytical form of these functions to generate the necessary data points for creating the spline fits.
    *   The `cutoff_f` helper function itself can be either an exponential or cosine form based on `#ifdef EXP_CUT`.

### Helper Functions for Spline Generation
*   `cc_VA, ch_VA, hh_VA, cc_VR, ch_VR, hh_VR`: These functions return the value of the attractive or repulsive pair potential for C-C, C-H, H-H interactions given a distance `dr` and relevant parameters. They are used by `rebo2_db_make_splines` if `#ifdef SPLINE_POTENTIAL` is active.
*   `cutoff_f(dr, l, h, m)`: Defines a generic cutoff function form (either exponential if `#ifdef EXP_CUT` or cosine-based) used by `rebo2_db_make_splines` if `#ifdef SPLINE_CUTOFF` is active.

## Important Variables/Constants
*   The subroutines extensively modify the `this` object of `BOP_TYPE` (`rebo2_t`), populating its parameter fields, derived constants, cutoff distances, table objects (`table2d_type`, `table3d_type`), and spline objects (`spline_object_t`).
*   Integer constants like `C_C`, `C_H`, `H_H` (likely defined in `rebo2_type.f90` or globally) are used as indices for pair-specific parameters.
*   Conditional compilation (`SCREENING`, `SPLINE_POTENTIAL`, `SPLINE_CUTOFF`, `ZERO_TABLES`, `EXP_CUT`) greatly influences the behavior and initialization paths.

## Usage Examples
These are internal setup routines called by the main `INIT_FUNC` (e.g., `rebo2_init`) of the REBO2 module. They are not directly called by end-users.

## Dependencies and Interactions
*   **Internal Dependencies:**
    *   Relies on the `BOP_TYPE` (`rebo2_t`) structure.
    *   Depends on procedures from `rebo2_default_tables` to get raw numerical data for `Fcc`, `Pcc`, etc.
    *   Uses `table2d_type` and `table3d_type` (and their `init`/`del` methods) from Atomistica's table utilities.
    *   Uses `spline_object_t` (and its `init`/`del` methods) if splines are enabled.
    *   Uses `gauss1` for solving linear equations in spline coefficient calculation.
*   **External Libraries:** None explicitly.
*   **Interactions:**
    *   These routines are responsible for the complex and detailed setup of a REBO2 potential object, transforming raw parameters and table data into a usable state for the `bop_kernel_rebo2.f90`.
    *   The choice of using analytical functions, tabulated values, or splines for different parts of the potential is managed here based on preprocessor flags.
```
