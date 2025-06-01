# src/potentials/bop/juslin/juslin_registry.f90

## Overview

This file contains the `REGISTER_FUNC` subroutine, which is typically aliased to `juslin_register` via preprocessor macros in `juslin.f90`. The purpose of this subroutine is to register the Juslin Bond-Order Potential (BOP) and its comprehensive set of parameters with the `ptrdict` configuration management system. This registration allows users of the Atomistica simulation package to select the Juslin potential and customize its parameters through external configuration files or scripts.

The registration includes standard BOP parameters as well as parameters specific to the Juslin formulation, such as those for the triplet-dependent `h` function (`alpha`, `omega`, `m`).

## Key Components

### Functions/Subroutines

*   `REGISTER_FUNC(this, cfg, m)`
    *   **Description**: This subroutine interfaces with the `ptrdict` library to make the parameters of a Juslin potential instance (`this` of type `juslin_t`) accessible to the configuration system.
        1.  It registers a new section in the configuration dictionary (`cfg`) using a name derived from `BOP_STR` (e.g., "Juslin"). The description of this section notes whether the potential was compiled with screening capabilities (`#ifdef SCREENING`).
        2.  It then registers various properties associated with the Juslin potential, primarily from the `this%db` component (the parameter database). For each property, it provides `ptrdict` with its name, memory location, maximum size, the location of its count variable, and a brief description.
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE), target`: The Juslin potential object (`juslin_t`) whose parameters are to be registered. The `target` attribute allows memory addresses of its components to be taken using `c_loc`.
        *   `cfg :: type(c_ptr), intent(in)`: A C pointer to the parent configuration dictionary.
        *   `m :: type(c_ptr), intent(out)`: A C pointer that receives the address of the newly created section for this Juslin potential within `cfg`.
    *   **Registered Parameters**:
        *   `el`: List of element symbols (max size `JUSLIN_MAX_EL`).
        *   `ref`: The reference string for selecting a predefined parameter set (max size `JUSLIN_MAX_REF`). *Note: The source `c_loc(this%ref(1:1))` might only register the first character of the reference string if `this%ref` is a character array. Assuming `this%ref` is `CHARACTER(JUSLIN_MAX_REF)` as per `juslin_type.f90`, `c_loc(this%ref)` would be more appropriate for the full string.*
        *   **Pair Parameters** (max size `JUSLIN_MAX_PAIRS` which is `JUSLIN_MAX_EL**2`): `D0`, `r0`, `S`, `beta`, `gamma`, `c`, `d`, `h` (angular term parameter), `n`, `r1`, `r2`.
        *   **Triplet Parameters for `h` function** (max size `JUSLIN_MAX_EL**3`):
            *   `alpha :: real(DP)`
            *   `omega :: real(DP)`
            *   `m :: integer` (registered as an integer list)
        *   **Screening Parameters** (conditional on `#ifdef SCREENING`, max size `JUSLIN_MAX_PAIRS`): `or1`, `or2`, `bor1`, `bor2`, `Cmin`, `Cmax`.
        Each parameter is typically registered with the description "See functional form," prompting users to consult documentation on the Juslin potential for details.

## Important Variables/Constants

*   `BOP_STR`: A preprocessor macro defined in `juslin.f90` (e.g., "Juslin"), used as the name for the section in the configuration dictionary.
*   `JUSLIN_MAX_EL`, `JUSLIN_MAX_PAIRS`, `JUSLIN_MAX_REF`: Constants (defined in `juslin_params.f90`) specifying maximum dimensions for parameter arrays and strings. `JUSLIN_MAX_EL**3` is used for triplet parameter arrays.
*   The subroutine uses `iso_c_binding` for Fortran-C interoperability (`c_ptr`, `c_loc`).
*   `CSTR()`: A macro (likely from `macros.inc`) for converting Fortran strings to C strings.

## Usage Examples

The `REGISTER_FUNC` is an internal setup routine called during Atomistica's initialization phase when potentials are made available. It is not directly invoked by users preparing simulation inputs but is essential for enabling parameter customization through those inputs.

```fortran
! Conceptual call during Atomistica initialization:
! TYPE(juslin_t) :: template_juslin_potential
! TYPE(c_ptr) :: main_config_dict, juslin_config_section
!
! ! Initialize main_config_dict via ptrdict
! ! ...
!
! ! Register a Juslin potential instance
! CALL juslin_register(template_juslin_potential, main_config_dict, juslin_config_section)
!
! ! Parameters in template_juslin_potential%db can now be set via input files
! ! processed by ptrdict, using names like "el", "ref", "D0", "alpha", etc.
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on the structure of `BOP_TYPE` (`juslin_t`) and its `db` component (`juslin_db_t`).
    *   Uses constants `JUSLIN_MAX_EL`, `JUSLIN_MAX_PAIRS`, `JUSLIN_MAX_REF` from `juslin_params.f90`.
    *   Utilizes `ptrdict` interface functions for registration.
*   **External Libraries:**
    *   `iso_c_binding` (Standard Fortran module).
    *   `ptrdict`: The dictionary and configuration management library.
*   **Interactions:**
    *   This registration routine makes the Juslin potential's parameters, including its unique triplet terms, discoverable and configurable by the `ptrdict` system.
    *   This allows for flexible runtime customization of the potential without requiring code recompilation.
```
