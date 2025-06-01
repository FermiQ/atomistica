# src/potentials/bop/kumagai/kumagai_registry.f90

## Overview

This file contains the `REGISTER_FUNC` subroutine, which is typically aliased to `kumagai_register` via preprocessor macros in `kumagai.f90`. The primary purpose of this subroutine is to register the Kumagai Bond-Order Potential (BOP) and its parameters with the `ptrdict` configuration management system. This registration allows users of the Atomistica simulation package to inspect and modify Kumagai potential parameters through external configuration files or scripts if such modifications are exposed by the `ptrdict` setup.

## Key Components

### Functions/Subroutines

*   `REGISTER_FUNC(this, cfg, m)`
    *   **Description**: This subroutine interfaces with the `ptrdict` library (using Fortran's `iso_c_binding` for C interoperability) to make the parameters of a Kumagai potential instance (`this` of type `kumagai_t`) accessible to the configuration system.
        1.  It first registers a new section in the configuration dictionary (`cfg`) using a name derived from `BOP_STR` (e.g., "Kumagai"). The description of this section indicates whether the potential was compiled with screening capabilities (`#ifdef SCREENING`).
        2.  It then registers various numerical parameters associated with the Kumagai potential, which are stored within the `this%db` component (the parameter database of the `kumagai_t` object). For each parameter, it provides `ptrdict` with its name, memory location, maximum size (defined by `KUMAGAI_MAX_EL` or `KUMAGAI_MAX_PAIRS`), the location of its count variable (e.g., `this%db%nA`), and a brief description (typically "See functional form").
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE), target`: The Kumagai potential object (`kumagai_t`) whose parameters are to be registered. The `target` attribute allows memory addresses of its components to be obtained using `c_loc`.
        *   `cfg :: type(c_ptr), intent(in)`: A C pointer to the parent configuration dictionary provided by the `ptrdict` system.
        *   `m :: type(c_ptr), intent(out)`: A C pointer that will receive the address of the newly created section for this Kumagai potential within `cfg`.
    *   **Registered Parameters**:
        *   `el`: List of element symbols (max size `KUMAGAI_MAX_EL`).
        *   **Pair Parameters** (max size `KUMAGAI_MAX_PAIRS`): `A`, `B`, `lambda1`, `lambda2`, `alpha`, `r1`, `r2`.
        *   `beta :: integer`: An integer pair parameter (max size `KUMAGAI_MAX_PAIRS`).
        *   **Single-Element Parameters** (max size `KUMAGAI_MAX_EL`): `eta`, `delta`, `c1`, `c2`, `c3`, `c4`, `c5`, `h`. These are used in the bond order and angular `g` functions.
        *   **Screening Parameters** (conditional on `#ifdef SCREENING`, max size `KUMAGAI_MAX_PAIRS`): `or1`, `or2`, `bor1`, `bor2`, `Cmin`, `Cmax`.
    *   *Notable Absence*: Unlike the registry functions for Brenner and Juslin potentials, this `REGISTER_FUNC` for Kumagai does **not** register the `ref` (reference string) parameter. This implies that the specific Kumagai parameter set (e.g., for Si) is likely chosen at compile time or hardcoded as the default, and cannot be switched using a reference string via the `ptrdict` configuration system.

## Important Variables/Constants

*   `BOP_STR`: A preprocessor macro defined in `kumagai.f90` (e.g., "Kumagai"), used as the name for the section in the configuration dictionary.
*   `KUMAGAI_MAX_EL`, `KUMAGAI_MAX_PAIRS`: Constants (defined in `kumagai_params.f90`) specifying the maximum dimensions for parameter arrays.
*   The subroutine uses `iso_c_binding` for Fortran-C interoperability (`c_ptr`, `c_loc`).
*   `CSTR()`: A macro (likely from `macros.inc`) for converting Fortran character strings to C strings.

## Usage Examples

The `REGISTER_FUNC` is an internal setup routine called during Atomistica's initialization phase. It is not directly invoked by users but is essential for enabling parameter inspection and modification through configuration files if the `ptrdict` system is used for this purpose.

```fortran
! Conceptual call during Atomistica initialization:
! TYPE(kumagai_t) :: template_kumagai_potential
! TYPE(c_ptr) :: main_config_dict, kumagai_config_section
!
! ! Initialize main_config_dict via ptrdict
! ! ...
!
! ! Register a Kumagai potential instance
! CALL kumagai_register(template_kumagai_potential, main_config_dict, kumagai_config_section)
!
! ! Parameters in template_kumagai_potential%db can now potentially be
! ! viewed or overridden via ptrdict, using names like "el", "A", "eta", etc.
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on the structure of `BOP_TYPE` (`kumagai_t`) and its embedded `db` component (`kumagai_db_t`).
    *   Uses constants `KUMAGAI_MAX_EL`, `KUMAGAI_MAX_PAIRS` from `kumagai_params.f90`.
    *   Utilizes `ptrdict` interface functions for registration.
*   **External Libraries:**
    *   `iso_c_binding` (Standard Fortran module).
    *   `ptrdict`: The dictionary and configuration management library.
*   **Interactions:**
    *   This registration routine makes the Kumagai potential's numerical parameters accessible via the `ptrdict` system.
    *   The lack of `ref` string registration suggests that parameter set selection for Kumagai potentials might be less flexible at runtime through `ptrdict` compared to other potentials like Brenner or Juslin.
```
