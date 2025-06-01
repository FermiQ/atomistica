# src/potentials/bop/brenner/brenner_module.f90

## Overview

This file provides key high-level module procedures that form the public interface for the Brenner Bond-Order Potential (BOP). These subroutines manage the lifecycle and setup of the Brenner potential object (`BOP_TYPE`, aliased to `brenner_t`). Specifically, it includes:
*   An initialization routine (`INIT_FUNC`) to set up parameters for the Brenner potential, either from a provided database structure or by looking up a named parameter set.
*   A binding routine (`BIND_TO_FUNC`) to connect the potential to the particle data and neighbor list handler, performing necessary pre-calculations and requesting interaction ranges.

This file also includes default implementations for:
*   The destructor (`DEL_FUNC`) via `../default_del_func.f90`.
*   The primary computation dispatch routine (`COMPUTE_FUNC`) via `../default_compute_func.f90`, which typically calls the specialized `BOP_KERNEL`.

The names `INIT_FUNC`, `BIND_TO_FUNC`, etc., are preprocessor macros that are defined in `brenner.f90` to point to Brenner-specific implementations like `brenner_init`, `brenner_bind_to`.

## Key Components

### Functions/Subroutines

*   `INIT_FUNC(this, db, el, D0, r0, S, beta, gamma, c, d, h, mu, n, m, r1, r2, ..., ierror)`
    *   **Description**: This is the constructor for the `brenner_t` (BOP_TYPE) potential object. It populates the `this` variable with the necessary parameters for the Brenner potential.
        *   It can accept a fully populated parameter database (`db` of type `brenner_db_t`).
        *   Alternatively, if `this%ref` (a reference string) is set, it attempts to find a matching parameter set from a predefined global list of Brenner databases (`BOP_DB`).
        *   Individual parameters (e.g., `D0`, `r0`, `S`, `beta`, etc., corresponding to terms in the Brenner potential formula) can be passed as optional arguments to override values from the selected database.
        *   The routine logs the parameters being used and checks for consistency.
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE), intent(inout)`: The Brenner potential object to be initialized.
        *   `db :: type(BOP_DB_TYPE), intent(in), optional`: A specific database of Brenner parameters.
        *   `el :: character(2), intent(in), optional`: Array of element symbols (e.g., "C", "H").
        *   `D0, r0, S, beta, ... :: real(DP), intent(in), optional`: Arrays of Brenner potential parameters. The list is extensive and includes all terms needed for the Brenner formulation (pair interactions, angular terms, bond order parameters, cutoff ranges). If `SCREENING` is enabled, additional parameters (`or1, or2, Cmin, Cmax`, etc.) are included.
        *   `ierror :: integer, intent(inout), optional`: Error status flag.

*   `BIND_TO_FUNC(this, p, nl, ierror)`
    *   **Description**: This subroutine "binds" an initialized Brenner potential object (`this`) to a specific physical system described by `particles_t` (`p`) and a neighbor list generator (`neighbors_t`, `nl`). This is a crucial setup step before any energy or force calculations can be performed. Its tasks include:
        1.  Performing consistency checks on the dimensions of parameter arrays stored in `this%db`.
        2.  Creating a mapping from atomic numbers (Z values from `p%el2Z`) to the internal indices used by the Brenner parameter database (`this%Z2db`).
        3.  Pre-calculating several derived parameters from the base set in `this%db` to optimize later computations. These include terms like `bo_exp` (bond order exponent), `expR` and `expA` (Morse potential exponents), `c_sq`, `d_sq` (related to angular terms), `VR_f`, `VA_f` (pair potential prefactors).
        4.  Initializing cutoff function objects (`this%cut_in`, `this%cut_out`, `this%cut_bo`) and storing plain cutoff distances (`cut_in_l`, `cut_in_h`, `max_cut_sq`, etc.) based on the parameters `r1`, `r2`, `or1`, `or2`, `bor1`, `bor2`.
        5.  Informing the neighbor list handler (`nl`) about the maximum interaction distance required for each pair of element types. This is done by calling `request_interaction_range`. If screening is enabled (`#ifdef SCREENING`), the requested range might be larger to include potential screening atoms, adjusted by `this%C_dr_cut`.
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE), intent(inout)`: The initialized Brenner potential object.
        *   `p :: type(particles_t), intent(inout)`: The structure containing particle information (types, positions, etc.).
        *   `nl :: type(neighbors_t), intent(inout)`: The neighbor list handler object.
        *   `ierror :: integer, optional, intent(out)`: Error status flag.

### Included Files

*   `../default_del_func.f90`: This included file provides the implementation for `DEL_FUNC` (e.g., `brenner_del`). This is the destructor for the `brenner_t` object, responsible for deallocating any dynamically allocated memory within the object, particularly for internal neighbor list buffers if they were managed by the potential itself.
*   `../default_compute_func.f90`: This included file provides the implementation for `COMPUTE_FUNC` (e.g., `brenner_energy_and_forces`). This subroutine is the primary entry point for calculating the potential energy and forces. It typically prepares arguments and then calls the main computational kernel (`BOP_KERNEL`, which is `brenner_kernel` for this potential).

## Important Variables/Constants

These subroutines primarily manipulate the `this` variable (of `BOP_TYPE` / `brenner_t`), which stores:
*   `this%db`: A pointer or allocatable component of `BOP_DB_TYPE` (`brenner_db_t`) holding the raw parameters.
*   `this%ref`: A reference string for selecting a parameter set.
*   Arrays for pre-calculated parameters (e.g., `this%bo_exp`, `this%expR`, `this%expA`).
*   Cutoff function objects and distances (e.g., `this%cut_in`, `this%cut_in_l`, `this%max_cut_sq`).
*   Mapping tables like `this%Z2db`.

## Usage Examples

The subroutines in this file are part of the internal mechanics of using the Brenner potential within the Atomistica framework. Users would typically interact with them via the aliased names (`brenner_init`, `brenner_bind_to`, `brenner_energy_and_forces`).

```fortran
! Conceptual sequence:
! USE brenner_module ! (Assuming brenner.f90 defines brenner_module which makes these available)
! TYPE(brenner_t) :: pot
! TYPE(particles_t) :: atoms
! TYPE(neighbors_t) :: neighbor_handler
! REAL(DP) :: total_energy, forces_on_atoms(3,num_atoms), virial(3,3)
!
! ! 1. Initialize the potential
! CALL brenner_init(pot, ref_string="Brenner_CH_1990_original") ! Or pass db/parameters directly
!
! ! 2. Setup atom data (positions, types) in 'atoms'
! ! ...
!
! ! 3. Initialize neighbor handler
! ! ...
!
! ! 4. Bind potential to particles and neighbor handler
! CALL brenner_bind_to(pot, atoms, neighbor_handler)
!
! ! 5. Compute energy and forces
! CALL brenner_energy_and_forces(pot, atoms, neighbor_handler, total_energy, forces_on_atoms, virial)
!
! ! 6. Clean up
! CALL brenner_del(pot)
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Types: `BOP_TYPE` (`brenner_t`), `BOP_DB_TYPE` (`brenner_db_t`), `particles_t`, `neighbors_t`.
    *   Parameter database: `BOP_DB` (global array of `brenner_db_t`).
    *   Utility functions: `prlog` (logging), `a2s` (array/string conversion), `atomic_number` (string to Z), `Z2pair` (element types to pair index), `request_interaction_range`.
    *   Macros: `ASSIGN_STRING_ARRAY_PROPERTY`, `ASSIGN_ARRAY_PROPERTY`, `RAISE_ERROR`, `INIT_ERROR`.
    *   The included `default_del_func.f90` and `default_compute_func.f90` are direct dependencies for providing the full potential interface.
*   **External Libraries:** None explicitly used in this file.
*   **Interactions:**
    *   `INIT_FUNC` populates the `brenner_t` object with parameters.
    *   `BIND_TO_FUNC` prepares this object for calculations on a specific atomic system by setting up derived parameters and cutoffs, and by configuring the neighbor list handler.
    *   The `COMPUTE_FUNC` (from the include) will subsequently use this prepared `brenner_t` object and call the `BOP_KERNEL` (specialized as `brenner_kernel`) to perform the actual physics calculations.
```
