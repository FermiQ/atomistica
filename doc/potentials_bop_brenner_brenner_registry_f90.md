# src/potentials/bop/brenner/brenner_registry.f90

## Overview

This file contains the `REGISTER_FUNC` subroutine (which is typically aliased to a potential-specific name like `brenner_register` via preprocessor macros in `brenner.f90`). The primary purpose of this subroutine is to register the Brenner Bond-Order Potential (BOP) and its numerous parameters with a configuration management system, referred to as `ptrdict`. This registration process makes the Brenner potential selectable and its parameters tunable through external configuration files or scripts used by the Atomistica simulation package.

## Key Components

### Functions/Subroutines

*   `REGISTER_FUNC(this, cfg, m)`
    *   **Description**: This subroutine interfaces with the `ptrdict` library (via C interoperability provided by `iso_c_binding`) to make the parameters of a Brenner potential instance (`this`) accessible to the configuration system.
        1.  It first registers a new section in the configuration dictionary (`cfg`) using a name derived from `BOP_STR` (e.g., "Brenner"). The description of this section indicates whether the potential was compiled with screening capabilities (`#ifdef SCREENING`).
        2.  It then registers various properties associated with the Brenner potential. For each property, it provides the `ptrdict` system with:
            *   A name for the property (e.g., "el", "ref", "D0", "r0").
            *   The memory location (C address) of the corresponding parameter within the `this%db` structure (the database part of the `brenner_t` object).
            *   The maximum size of the parameter array (e.g., `BRENNER_MAX_PAIRS`).
            *   The memory location of the variable that stores the actual number of elements in that array for the current parameter set (e.g., `this%db%nD0` for the `D0` array).
            *   A brief description of the parameter.
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE), target`: The Brenner potential object (`brenner_t`) whose parameters are to be registered. The `target` attribute is crucial as `c_loc` is used to get memory addresses of its components.
        *   `cfg :: type(c_ptr), intent(in)`: A C pointer to the parent configuration dictionary provided by the `ptrdict` system.
        *   `m :: type(c_ptr), intent(out)`: A C pointer that will receive the address of the newly created section for this Brenner potential within the `cfg` dictionary.
    *   **Registered Parameters**:
        *   `el`: List of element symbols (e.g., "C", "H").
        *   `ref`: The reference string used to select a predefined parameter set from the internal database.
        *   Core Brenner potential parameters: `D0`, `r0`, `S`, `beta`, `gamma`, `c`, `d`, `h`, `mu`, `n`, `m` (integer list), `r1`, `r2`.
        *   If compiled with screening (`#ifdef SCREENING`): `or1`, `or2` (outer cutoff ranges), `bor1`, `bor2` (bond-order specific cutoff ranges), `Cmin`, `Cmax` (screening function parameters).
        Each parameter is registered with a descriptive string, often "See functional form," implying users should refer to Brenner potential documentation for details on each parameter's role.

## Important Variables/Constants

*   `BOP_STR`: A preprocessor macro defined in `brenner.f90` (e.g., "Brenner"), used as the name for the section in the configuration dictionary.
*   `BRENNER_MAX_EL`, `BRENNER_MAX_PAIRS`, `BRENNER_MAX_REF`: Constants (defined in `brenner_params.f90`) that specify the maximum dimensions of parameter arrays and strings. These are used when calling `ptrdict` registration functions.
*   The subroutine uses `iso_c_binding` intrinsic module for Fortran-C interoperability, specifically the `c_ptr` type and `c_loc` function.
*   `CSTR()`: A macro (likely from `macros.inc`) used to convert Fortran character strings to null-terminated C strings required by `ptrdict` functions.

## Usage Examples

The `REGISTER_FUNC` subroutine is part of the setup phase of the Atomistica simulation code. It is typically called once when the program starts or when potentials are being made available to the simulation engine. It is not a function that end-users interact with directly when preparing a simulation input file, but it enables them to customize potential parameters through those input files.

```fortran
! Conceptual call during Atomistica initialization:
! TYPE(brenner_t) :: template_brenner_potential
! TYPE(c_ptr) :: main_config_dictionary, brenner_config_section
!
! ! Initialize main_config_dictionary via ptrdict
! ! ...
!
! ! Register a Brenner potential instance
! CALL brenner_register(template_brenner_potential, main_config_dictionary, brenner_config_section)
!
! ! Now, parameters in template_brenner_potential%db can be set via input files
! ! processed by ptrdict, using the names "el", "ref", "D0", etc.
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on the structure of `BOP_TYPE` (`brenner_t`) and its embedded `db` member (of `BOP_DB_TYPE` / `brenner_db_t`), which holds the parameters.
    *   Uses constants `BRENNER_MAX_EL`, `BRENNER_MAX_PAIRS`, `BRENNER_MAX_REF` from `brenner_params.f90`.
    *   Utilizes `ptrdict` interface functions: `ptrdict_register_section`, `ptrdict_register_string_list_property`, `ptrdict_register_string_property`, `ptrdict_register_list_property`, and `ptrdict_register_integer_list_property`.
*   **External Libraries:**
    *   `iso_c_binding` (Standard Fortran module).
    *   `ptrdict`: An external or internal library providing the dictionary and configuration management services.
*   **Interactions:**
    *   This registration routine makes the Brenner potential's parameters discoverable and configurable by the `ptrdict` system.
    *   When the simulation input is parsed, `ptrdict` can use the registered information to locate and modify the parameters in the `brenner_t` object before the simulation starts. This allows for runtime customization of potentials without recompiling the code.
```
