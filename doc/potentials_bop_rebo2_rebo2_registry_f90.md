# src/potentials/bop/rebo2/rebo2_registry.f90

## Overview

This file contains the `REGISTER_FUNC` subroutine, which is typically aliased to `rebo2_register` via preprocessor macros in `rebo2.f90` (or `rebo2_scr.f90`). The primary function of this subroutine is to register the second-generation REBO (REBO2) potential and its extensive set of parameters with the `ptrdict` configuration management system. This registration allows users of the Atomistica simulation package to inspect and modify REBO2 parameters through external configuration files or scripts.

A key feature of this registration process is that it first populates the input table arrays within the `this` (potential object) with default values by calling routines from `rebo2_default_tables`. This ensures that the `ptrdict` system has access to these default values, which can then be selectively overridden by user configurations.

## Key Components

### Functions/Subroutines

*   `REGISTER_FUNC(this, cfg, m)`
    *   **Description**: This subroutine interfaces with the `ptrdict` library to make the parameters of a REBO2 potential instance (`this` of type `rebo2_t`) accessible to the configuration system.
        1.  **Load Default Table Data**: Before registering, it calls:
            *   `rebo2_default_Fcc_table(this%in_Fcc, this%in_dFdi, this%in_dFdj, this%in_dFdk)`
            *   `rebo2_default_Fch_table(this%in_Fch)`
            *   `rebo2_default_Fhh_table(this%in_Fhh)`
            *   `rebo2_default_Pcc_table(this%in_Pcc)`
            *   `rebo2_default_Pch_table(this%in_Pch)`
            *   `rebo2_default_Tcc_table(this%in_Tcc)`
            This populates the `in_...` arrays within the `this` object with standard REBO2 values.
        2.  **Register Section**: It registers a new section in the configuration dictionary (`cfg`) using a name derived from `BOP_STR` (e.g., "Rebo2"). The description notes if it's the screened version.
        3.  **Register Parameters**: It then registers numerous properties:
            *   `elements`: A string for element types this potential applies to.
            *   **Scalar Parameters for Pair Interactions**: Individual parameters for H-H, C-H, and C-C interactions like `HH_Q`, `HH_A`, `HH_alpha`, `HH_B1`, `HH_beta1`, `CH_Q`, ..., `CC_alpha`, `CC_B1`, `CC_B2`, `CC_B3`, `CC_beta1`, `CC_beta2`, `CC_beta3`. These directly modify fields in the `this` object (e.g., `this%hh_Q`, `this%cc_B1`).
            *   **Screening Parameters** (if `#ifdef SCREENING`): `Cmin`, `Cmax`.
            *   **Cutoff Radii**: Parameters like `CC_in_r1`, `CC_in_r2`, and if screening is enabled, `CC_ar_r1`, `CC_ar_r2`, `CC_bo_r1`, `CC_bo_r2`, `CC_nc_r1`, `CC_nc_r2`. Also `CH_r1`, `CH_r2`, `HH_r1`, `HH_r2`.
            *   **Boolean Flags**: `with_dihedral` (to include dihedral terms) and `zero_tables` (to initialize tables to zero, overriding defaults).
            *   **Table Data Arrays**: The raw data for the lookup tables are registered as 2D or 3D arrays: `in_Fcc`, `in_dFdi`, `in_dFdj`, `in_dFdk`, `in_Fch`, `in_Fhh`, `in_Pcc`, `in_Pch`, `in_Tcc`. This allows users to potentially provide entire custom tables via the configuration system.
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE), target`: The REBO2 potential object (`rebo2_t`).
        *   `cfg :: type(c_ptr), intent(in)`: Pointer to the parent configuration dictionary.
        *   `m :: type(c_ptr), intent(out)`: Pointer to the newly created REBO2 section in `cfg`.

## Important Variables/Constants

*   `BOP_STR`: Preprocessor macro for the section name (e.g., "Rebo2").
*   `MAX_EL_STR`: Constant defining the maximum length of the `elements` string.
*   The registered parameters are fields within the `this` object of `BOP_TYPE` (`rebo2_t`), including scalar values and the `in_...` arrays for table data.

## Usage Examples

This `REGISTER_FUNC` is an internal setup routine called during Atomistica's initialization. It's crucial for enabling detailed runtime customization of the REBO2 potential.

```fortran
! Conceptual call during Atomistica initialization:
! TYPE(rebo2_t) :: template_rebo2_potential
! TYPE(c_ptr) :: main_config_dict, rebo2_config_section
!
! ! Initialize main_config_dict via ptrdict
! ! ...
!
! ! Register a REBO2 potential instance. This also loads default table data into template_rebo2_potential.
! CALL rebo2_register(template_rebo2_potential, main_config_dict, rebo2_config_section)
!
! ! Parameters in template_rebo2_potential, including table arrays like in_Fcc,
! ! can now be set or overridden via input files processed by ptrdict.
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on the structure of `BOP_TYPE` (`rebo2_t`).
    *   Critically calls the table-filling subroutines from the `rebo2_default_tables` module (e.g., `rebo2_default_Fcc_table`).
    *   Utilizes `ptrdict` interface functions for registration (e.g., `ptrdict_register_section`, `ptrdict_register_real_property`, `ptrdict_register_array3d_property`).
*   **External Libraries:**
    *   `iso_c_binding` (Standard Fortran module).
    *   `ptrdict`: The dictionary and configuration management library.
*   **Interactions:**
    *   This routine makes nearly all parameters of the REBO2 potential, including the detailed lookup table data, accessible and configurable at runtime via the `ptrdict` system.
    *   By first loading default values into the `this` object's `in_...` table arrays, it ensures that these defaults are available to `ptrdict` and can be modified by the user if needed. These `in_...` arrays are then used by `rebo2_db_init_with_parameters` (called during `BIND_TO_FUNC`) to set up the actual `table2d_type` and `table3d_type` objects used in calculations.
```
