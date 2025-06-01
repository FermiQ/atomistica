# src/potentials/bop/tersoff/tersoff_module.f90

## Overview

This file provides the high-level module procedures that form the public interface for the Tersoff Bond-Order Potential (BOP) within the Atomistica framework. The primary subroutine defined directly in this file is `INIT_FUNC`, which serves as the constructor for initializing the Tersoff potential parameters.

For other standard potential operations, this file follows a common pattern in Atomistica by including default implementations from the parent `bop` directory:
*   `DEL_FUNC` (destructor) is included from `../default_del_func.f90`.
*   `BIND_TO_FUNC` (binding to particle data and neighbor lists) is included from `../default_bind_to_func.f90`.
*   `COMPUTE_FUNC` (main energy/force calculation routine) is included from `../default_compute_func.f90`.

The names `INIT_FUNC`, `DEL_FUNC`, etc., are preprocessor macros that are defined in `tersoff.f90` to point to Tersoff-specific implementations (e.g., `tersoff_init` for `INIT_FUNC`) or the included defaults.

## Key Components

### Functions/Subroutines

*   `INIT_FUNC(this, db)`
    *   **Description**: This is the constructor for the `tersoff_t` (BOP_TYPE) potential object. It initializes the `this` variable with the parameters for the Tersoff potential.
        *   If a `db` argument (of type `tersoff_db_t`) is provided, its parameters are copied into `this%db`. This allows for passing a specific, pre-configured set of Tersoff parameters.
        *   If `db` is not provided, the routine uses the parameters already present in `this%db`. This `this%db` component would have been initialized either by a default value set in `tersoff_type.f90` (if `tersoff_params.f90` defines `tersoff_db_type` and `tersoff_type.f90` initializes `db` with a specific parameter instance from `tersoff_params.f90`) or by a prior mechanism if this function is called on an already partially initialized object.
        *   After ensuring `this%db` is populated, the subroutine logs the values of all key parameters of the Tersoff potential being used. This includes element types (`el`), pair potential parameters (`A`, `B`, `lambda`, `mu`), bond order parameters (`xi`, `omega`, `mubo`, `m`, `beta`, `n`), angular function parameters (`c`, `d`, `h`), and cutoff radii (`r1`, `r2`).
        *   If compiled with screening (`#ifdef SCREENING`), it also logs screening-specific parameters (`or1`, `or2`, `bor1`, `bor2`, `Cmin`, `Cmax`).
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE), intent(inout)`: The Tersoff potential object (`tersoff_t`) to be initialized.
        *   `db :: type(BOP_DB_TYPE), optional, intent(in)`: An optional specific database of Tersoff parameters (`tersoff_db_t`).

### Included Files

*   `../default_del_func.f90`: Provides the standard implementation for `DEL_FUNC` (e.g., `tersoff_del`). This handles the deallocation of internal neighbor list arrays stored within the `BOP_TYPE` object.
*   `../default_bind_to_func.f90`: Provides the standard implementation for `BIND_TO_FUNC` (e.g., `tersoff_bind_to`). This routine prepares the potential for calculations on a specific atomic system by setting up cutoffs (based on `r1`, `r2` from `this%db` and the `CUTOFF_T` type), mapping element types, and requesting appropriate interaction ranges from the neighbor list handler.
*   `../default_compute_func.f90`: Provides the standard implementation for `COMPUTE_FUNC` (e.g., `tersoff_energy_and_forces`). This routine updates neighbor lists and then dispatches the main energy and force calculation to the `BOP_KERNEL` (which is `tersoff_kernel`, an alias for the generic BOP kernel defined in `../bop_kernel.f90`, used with Tersoff-specific functions from `tersoff_func.f90`).

## Important Variables/Constants

*   The `INIT_FUNC` primarily operates on the `this` variable of `BOP_TYPE` (`tersoff_t`) and its `db` component, which holds the actual parameter values defined in `tersoff_params.f90`.
*   `BOP_NAME_STR` (defined in `tersoff.f90` as "tersoff") is used by `INIT_FUNC` for logging purposes.

## Usage Examples

The subroutines in this file are part of the public interface of the Tersoff potential module, typically exposed via aliases defined in `tersoff.f90`.

```fortran
! Conceptual sequence for using the Tersoff potential:
! USE tersoff_module ! (Assuming tersoff.f90 defines this module)
! TYPE(tersoff_t) :: pot_tersoff
! TYPE(tersoff_db_t) :: specific_params ! Optional, can load from file or use defaults
! ! ... potentially load/set parameters in specific_params ...
!
! ! 1. Initialize the Tersoff potential
! CALL tersoff_init(pot_tersoff, db=specific_params) ! Or use default internal db
!
! ! 2. Bind the potential to the system (uses default_bind_to_func)
! CALL tersoff_bind_to(pot_tersoff, atoms_config, neighbor_handler_config)
!
! ! 3. Compute energy and forces (uses default_compute_func)
! CALL tersoff_energy_and_forces(pot_tersoff, atoms_config, neighbor_handler_config, &
!                               total_potential_energy, forces, virial_tensor)
!
! ! 4. Clean up (uses default_del_func)
! CALL tersoff_del(pot_tersoff)
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Types: `BOP_TYPE` (`tersoff_t`), `BOP_DB_TYPE` (`tersoff_db_t`).
    *   Parameter database: `this%db` is populated with parameters defined in `tersoff_params.f90`.
    *   Utility functions: `prlog` (logging), `a2s` (array to string).
    *   The included default function files (`default_del_func.f90`, `default_bind_to_func.f90`, `default_compute_func.f90`) are direct dependencies for providing the full potential interface, which in turn use the generic `bop_kernel.f90` and Tersoff-specific functions from `tersoff_func.f90`.
*   **External Libraries:** None explicitly used in this file.
*   **Interactions:**
    *   `INIT_FUNC` configures the `tersoff_t` object with specific Tersoff parameters.
    *   The default `DEL_FUNC`, `BIND_TO_FUNC`, and `COMPUTE_FUNC` provide standardized mechanisms for memory management, system binding, and computation dispatch, tailored for BOP-style potentials and utilizing the generic BOP kernel.
```
