# src/potentials/bop/kumagai/kumagai_module.f90

## Overview

This file provides the high-level module procedures that form the public interface for the Kumagai Bond-Order Potential (BOP). The primary subroutine defined directly in this file is `INIT_FUNC`, which acts as the constructor for initializing the Kumagai potential parameters.

For other standard potential operations, this file includes default implementations from the parent `bop` directory:
*   `DEL_FUNC` (destructor) is included from `../default_del_func.f90`.
*   `BIND_TO_FUNC` (binding to particle data and neighbor lists) is included from `../default_bind_to_func.f90`.
*   `COMPUTE_FUNC` (main energy/force calculation routine) is included from `../default_compute_func.f90`.

The names `INIT_FUNC`, `DEL_FUNC`, etc., are preprocessor macros that are defined in `kumagai.f90` to point to Kumagai-specific implementations (e.g., `kumagai_init`) or the included defaults.

## Key Components

### Functions/Subroutines

*   `INIT_FUNC(this, db)`
    *   **Description**: This is the constructor for the `kumagai_t` (BOP_TYPE) potential object. It initializes the `this` variable with the parameters for the Kumagai potential.
        *   If a `db` argument (of type `kumagai_db_t`) is provided, its parameters are copied into `this%db`.
        *   If `db` is not provided, the routine uses the parameters already present in `this%db` (which might have been set by a default initialization in `kumagai_type.f90` or by a previous call).
        *   After ensuring `this%db` is populated, the subroutine logs the values of all key parameters of the Kumagai potential being used. This includes element types (`el`), pair potential parameters (`A`, `B`, `lambda1`, `lambda2`), bond order parameters (`eta`, `delta`), length-dependent `h` function parameters (`alpha`, `beta`), angular `g` function parameters (`c1` through `c5`, `h`), and cutoff radii (`r1`, `r2`).
        *   If compiled with screening (`#ifdef SCREENING`), it also logs screening-specific parameters (`or1`, `or2`, `bor1`, `bor2`, `Cmin`, `Cmax`).
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE), intent(inout)`: The Kumagai potential object (`kumagai_t`) to be initialized.
        *   `db :: type(BOP_DB_TYPE), optional, intent(in)`: An optional specific database of Kumagai parameters (`kumagai_db_t`).

### Included Files

*   `../default_del_func.f90`: Provides the standard implementation for `DEL_FUNC` (e.g., `kumagai_del`). This handles deallocation of internal neighbor list arrays.
*   `../default_bind_to_func.f90`: Provides the standard implementation for `BIND_TO_FUNC` (e.g., `kumagai_bind_to`). This routine prepares the potential for calculations on a specific atomic system by setting up cutoffs, mapping element types, and requesting interaction ranges from the neighbor list handler.
*   `../default_compute_func.f90`: Provides the standard implementation for `COMPUTE_FUNC` (e.g., `kumagai_energy_and_forces`). This routine updates neighbor lists and dispatches the main calculation to the `BOP_KERNEL` (which is `kumagai_kernel`).

## Important Variables/Constants

*   The `INIT_FUNC` primarily operates on the `this` variable of `BOP_TYPE` (`kumagai_t`) and its `db` component, which holds the actual parameter values.
*   `BOP_NAME_STR` (defined in `kumagai.f90` as "kumagai") is used by `INIT_FUNC` for logging purposes.

## Usage Examples

The subroutines in this file are part of the public interface of the Kumagai potential module, typically exposed via aliases defined in `kumagai.f90`.

```fortran
! Conceptual sequence for using the Kumagai potential:
! USE kumagai_module ! (Assuming kumagai.f90 defines this module)
! TYPE(kumagai_t) :: pot_kumagai
! TYPE(kumagai_db_t) :: specific_params ! Optional, can load from file or use defaults
! ! ... potentially load/set parameters in specific_params ...
!
! ! 1. Initialize the Kumagai potential
! CALL kumagai_init(pot_kumagai, db=specific_params) ! Or use default internal db
!
! ! ... (Binding and Computation would use the included default handlers) ...
!
! ! 2. Bind the potential to the system (uses default_bind_to_func)
! CALL kumagai_bind_to(pot_kumagai, atoms_config, neighbor_handler_config)
!
! ! 3. Compute energy and forces (uses default_compute_func)
! CALL kumagai_energy_and_forces(pot_kumagai, atoms_config, neighbor_handler_config, &
!                               total_potential_energy, forces, virial_tensor)
!
! ! 4. Clean up (uses default_del_func)
! CALL kumagai_del(pot_kumagai)
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Types: `BOP_TYPE` (`kumagai_t`), `BOP_DB_TYPE` (`kumagai_db_t`).
    *   Parameter database: `this%db` is populated with parameters defined in `kumagai_params.f90`.
    *   Utility functions: `prlog` (logging), `a2s` (array to string).
    *   The included default function files (`default_del_func.f90`, `default_bind_to_func.f90`, `default_compute_func.f90`) are direct dependencies for providing the full potential interface.
*   **External Libraries:** None explicitly used in this file.
*   **Interactions:**
    *   `INIT_FUNC` sets up the `kumagai_t` object with specific Kumagai parameters.
    *   The default `DEL_FUNC`, `BIND_TO_FUNC`, and `COMPUTE_FUNC` provide standardized mechanisms for memory management, system binding, and computation dispatch, respectively. These defaults call the appropriate `BOP_KERNEL` (i.e., `kumagai_kernel`) for the actual physics.
```
