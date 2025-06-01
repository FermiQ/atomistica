# src/potentials/bop/tersoff/tersoff_registry.f90

## Overview

This file contains the `REGISTER_FUNC` subroutine, which is typically aliased to `tersoff_register` via preprocessor macros in `tersoff.f90`. The primary purpose of this subroutine is to register the Tersoff Bond-Order Potential (BOP) and its parameters with the `ptrdict` configuration management system. This registration allows users of the Atomistica simulation package to inspect and modify Tersoff potential parameters through external configuration files or scripts. The parameters registered are explicitly named following the conventions in J. Tersoff, Phys. Rev. B 39, 5566 (1989).

## Key Components

### Functions/Subroutines

*   `REGISTER_FUNC(this, cfg, m)`
    *   **Description**: This subroutine interfaces with the `ptrdict` library (using Fortran's `iso_c_binding` for C interoperability) to make the parameters of a Tersoff potential instance (`this` of type `tersoff_t`) accessible to the configuration system.
        1.  It first registers a new section in the configuration dictionary (`cfg`) using a name derived from `BOP_STR` (e.g., "Tersoff"). The description of this section explicitly mentions the parameter naming convention based on Tersoff's 1989 paper and indicates if the potential was compiled with screening capabilities (`#ifdef SCREENING`).
        2.  It then registers various properties associated with the Tersoff potential, primarily from the `this%db` component (the parameter database of the `tersoff_t` object). For each property, it provides `ptrdict` with its name (e.g., "el", "A", "B", "xi"), memory location, maximum size (defined by `TERSOFF_MAX_EL` or `TERSOFF_MAX_PAIRS`), the location of its count variable (e.g., `this%db%nA`), and a brief description (typically "See functional form").
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE), target`: The Tersoff potential object (`tersoff_t`) whose parameters are to be registered. The `target` attribute allows memory addresses of its components to be obtained using `c_loc`.
        *   `cfg :: type(c_ptr), intent(in)`: A C pointer to the parent configuration dictionary provided by the `ptrdict` system.
        *   `m :: type(c_ptr), intent(out)`: A C pointer that will receive the address of the newly created section for this Tersoff potential within `cfg`.
    *   **Registered Parameters**:
        *   `el`: List of element symbols (max size `TERSOFF_MAX_EL`).
        *   **Pair Parameters** (max size `TERSOFF_MAX_PAIRS`): `A`, `B`, `xi`, `lambda`, `mu`, `omega`, `mubo`, `r1`, `r2`.
        *   `m :: integer`: An integer pair parameter (max size `TERSOFF_MAX_PAIRS`).
        *   **Single-Element Parameters** (max size `TERSOFF_MAX_EL`): `beta`, `n`, `c`, `d`, `h`.
        *   **Screening Parameters** (conditional on `#ifdef SCREENING`, max size `TERSOFF_MAX_PAIRS`): `or1`, `or2`, `bor1`, `bor2`, `Cmin`, `Cmax`.
    *   *Notable Absence*: Similar to the Kumagai potential's registry, this function does not register the `ref` (reference string) parameter. This suggests that selection of specific Tersoff parameter sets (like SiC vs. pure Si) might be handled differently, possibly by using distinct class names or relying on the default set in `tersoff_params.f90`.

## Important Variables/Constants

*   `BOP_STR`: A preprocessor macro defined in `tersoff.f90` (e.g., "Tersoff"), used as the name for the section in the configuration dictionary.
*   `TERSOFF_MAX_EL`, `TERSOFF_MAX_PAIRS`: Constants (defined in `tersoff_params.f90`) specifying the maximum dimensions for parameter arrays.
*   The subroutine uses `iso_c_binding` for Fortran-C interoperability (`c_ptr`, `c_loc`).
*   `CSTR()`: A macro (likely from `macros.inc`) for converting Fortran character strings to C strings.

## Usage Examples

The `REGISTER_FUNC` is an internal setup routine called during Atomistica's initialization phase. It is not directly invoked by users preparing simulation inputs but is essential for enabling parameter customization through those inputs via the `ptrdict` system.

```fortran
! Conceptual call during Atomistica initialization:
! TYPE(tersoff_t) :: template_tersoff_potential
! TYPE(c_ptr) :: main_config_dict, tersoff_config_section
!
! ! Initialize main_config_dict via ptrdict
! ! ...
!
! ! Register a Tersoff potential instance
! CALL tersoff_register(template_tersoff_potential, main_config_dict, tersoff_config_section)
!
! ! Parameters in template_tersoff_potential%db can now potentially be
! ! viewed or overridden via ptrdict, using names like "el", "A", "beta", etc.
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on the structure of `BOP_TYPE` (`tersoff_t`) and its embedded `db` component (`tersoff_db_t`).
    *   Uses constants `TERSOFF_MAX_EL`, `TERSOFF_MAX_PAIRS` from `tersoff_params.f90`.
    *   Utilizes `ptrdict` interface functions for registration.
*   **External Libraries:**
    *   `iso_c_binding` (Standard Fortran module).
    *   `ptrdict`: The dictionary and configuration management library.
*   **Interactions:**
    *   This registration routine makes the Tersoff potential's parameters accessible and configurable by the `ptrdict` system, adhering to the naming convention of Tersoff's 1989 paper.
    *   The absence of `ref` string registration means parameter set selection relies on other mechanisms.
```
