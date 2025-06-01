# src/potentials/bop/juslin/juslin_module.f90

## Overview

This file provides the core module procedures for the Juslin W-C-H (Tungsten-Carbon-Hydrogen) Bond-Order Potential (BOP). These subroutines define the public Application Programming Interface (API) for creating, configuring, and utilizing the Juslin potential within the Atomistica simulation framework. The file includes implementations for:
*   `INIT_FUNC`: The constructor, responsible for initializing the Juslin potential parameters.
*   `DEL_FUNC`: The destructor, responsible for deallocating any dynamically allocated memory associated with the potential object.
*   `BIND_TO_FUNC`: A crucial setup routine that connects the potential definition to the physical system (particle data and neighbor lists) and performs necessary pre-calculations.
*   `COMPUTE_FUNC`: The main routine that triggers the calculation of potential energy, atomic forces, and virial stress.

The names `INIT_FUNC`, `DEL_FUNC`, etc., are preprocessor macros that are defined in `juslin.f90` to point to Juslin-specific implementations (e.g., `juslin_init`, `juslin_bind_to`).

## Key Components

### Functions/Subroutines

*   `INIT_FUNC(this, db, el, D0, r0, S, beta, gamma, c, d, h, n, m, alpha, omega, r1, r2, ..., ierror)`
    *   **Description**: This is the constructor for the `juslin_t` (BOP_TYPE) potential object. It populates the `this` variable with parameters for the Juslin potential.
        *   It can accept a fully populated parameter database (`db` of type `juslin_db_t`).
        *   If `db` is not provided, it can look up a parameter set from a global list (`BOP_DB` defined in `juslin_params.f90`) using a reference string stored in `this%ref`.
        *   Individual parameters (e.g., `D0`, `r0`, `S`, `beta`, `gamma`, `c`, `d`, `h` (for angular term), `n` (bond order exponent), `alpha`, `omega`, `m` (for Juslin's length-dependent `h` function), `r1`, `r2` (cutoff radii)) can be passed as optional arguments to override values from the selected database.
        *   Parameters specific to the Juslin potential's `h` function, `alpha` and `omega`, are included in the argument list.
        *   Screening-related parameters (`or1, or2, bor1, bor2, Cmin, Cmax`) are included if compiled with `#ifdef SCREENING`.
        *   The routine logs the parameters being used.
    *   **Key Arguments**: `this`, `db`, `el`, and numerous optional arrays for Juslin potential parameters.

*   `DEL_FUNC(this)`
    *   **Description**: This is the destructor for the `juslin_t` object. It deallocates internal, dynamically allocated arrays used for neighbor lists and bond data, provided `this%neighbor_list_allocated` is true. This implementation is identical to the one found in `default_del_func.f90`.
    *   **Key Arguments**: `this :: type(BOP_TYPE), intent(inout)`.

*   `BIND_TO_FUNC(this, p, nl, ierror)`
    *   **Description**: This subroutine "binds" an initialized Juslin potential object (`this`) to particle data (`p`) and a neighbor list handler (`nl`). Key setup tasks include:
        1.  **Parameter Validation**: Checks for consistent numbers of entries for parameter arrays in `this%db`.
        2.  **Parameter Symmetrization**: Iterates through pairs of element types (i,j) and if `this%db%r0(PAIR_INDEX_NS(i,j,...))` is negative (used as a flag), it copies all parameters from pair (j,i) to (i,j). This ensures that parameters for A-B interactions are defined even if only B-A is given in the input file. `PAIR_INDEX_NS` (non-symmetric) is used for indexing.
        3.  **Screening Setup**: If `#ifdef SCREENING` is active, initializes `this%Cmin`, `this%Cmax`, `this%dC`, and `this%C_dr_cut` from `this%db`.
        4.  **Element Mapping**: Maps atomic numbers (Z) to internal database indices (`this%Z2db`).
        5.  **Pre-calculation of Derived Parameters**: Computes terms like `bo_exp`, `expR`, `expA`, `c_sq`, `d_sq`, `c_d`, `VR_f`, `VA_f` for efficiency.
        6.  **Cutoff Initialization**: Sets up `cut_in_l`, `cut_in_h`, `cut_in_h2`. Critically, for the cosine-based cutoffs used in `juslin_func.f90`, it pre-calculates factors `this%cut_in_fca = PI / (r2 - r1)` and `this%cut_in_fc = -0.5 * this%cut_in_fca`. Similar setup for `cut_out_` and `cut_bo_` parameters if screening is enabled.
        7.  **Request Interaction Ranges**: Informs the neighbor list handler `nl` of the required interaction cutoff distances for each pair of element types present in the system `p`.
    *   **Key Arguments**: `this`, `p`, `nl`, `ierror`.

*   `COMPUTE_FUNC(this, p, nl, epot, f, wpot, ..., ierror)`
    *   **Description**: This is the main routine for calculating energy, forces, and virial. It first calls `update(nl, p, ...)` to ensure neighbor lists are current. Then, it maps element types to internal indices, determines `nebmax` and `nebavg` (maximum and average number of neighbors), and finally calls `BOP_KERNEL` (which is `juslin_kernel`) to perform the detailed calculations. This implementation is identical to the one in `default_compute_func.f90`.
    *   **Key Arguments**: `this`, `p`, `nl`, `epot`, `f`, `wpot`, optional per-atom/per-bond terms, `ierror`.

## Important Variables/Constants

*   These subroutines primarily operate on the `this` variable (of `BOP_TYPE` / `juslin_t`), which stores the potential's parameters (in `this%db`) and state.
*   The `PAIR_INDEX_NS` macro is used in `BIND_TO_FUNC` for parameter symmetrization.
*   The `INIT_FUNC` specifically handles Juslin parameters like `alpha` and `omega`.
*   The `BIND_TO_FUNC` pre-calculates `_fca` and `_fc` factors for the cosine cutoff functions defined in `juslin_func.f90`.

## Usage Examples

The subroutines in this file are part of the public interface of the Juslin potential module.

```fortran
! Conceptual sequence for using the Juslin potential:
! USE juslin_module ! (Assuming juslin.f90 defines this module)
! TYPE(juslin_t) :: pot_juslin
! TYPE(particles_t) :: atoms_config
! TYPE(neighbors_t) :: neighbor_handler_config
! REAL(DP) :: total_potential_energy, forces(3,num_atoms), virial_tensor(3,3)
!
! ! 1. Initialize the Juslin potential
! CALL juslin_init(pot_juslin, ref_string="Juslin_WCH_2005")
!
! ! 2. Setup atomic configuration in 'atoms_config'
! ! ...
!
! ! 3. Initialize neighbor list handler 'neighbor_handler_config'
! ! ...
!
! ! 4. Bind the potential to the system
! CALL juslin_bind_to(pot_juslin, atoms_config, neighbor_handler_config)
!
! ! 5. Compute energy and forces
! CALL juslin_energy_and_forces(pot_juslin, atoms_config, neighbor_handler_config, &
!                               total_potential_energy, forces, virial_tensor)
!
! ! 6. Clean up
! CALL juslin_del(pot_juslin)
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Types: `BOP_TYPE` (`juslin_t`), `BOP_DB_TYPE` (`juslin_db_t`), `particles_t`, `neighbors_t`.
    *   Parameter database: `BOP_DB` (global array of `juslin_db_t` from `juslin_params.f90`).
    *   Utility functions: `prlog`, `a2s`, `atomic_number`, `Z2pair`, `request_interaction_range`.
    *   Macros: `ASSIGN_STRING_ARRAY_PROPERTY`, `ASSIGN_ARRAY_PROPERTY`, `RAISE_ERROR`, `INIT_ERROR`, `PASS_ERROR`.
    *   Relies on `juslin_kernel` (via `BOP_KERNEL` macro) for actual computations.
*   **External Libraries:** None explicitly used in this file.
*   **Interactions:**
    *   `INIT_FUNC` populates the `juslin_t` object with parameters, including those specific to Juslin's formulation (e.g., `alpha`, `omega`).
    *   `DEL_FUNC` provides memory management for internal arrays.
    *   `BIND_TO_FUNC` prepares the `juslin_t` object for calculations on a specific atomic system, performing parameter symmetrization and pre-calculating factors for its specific cosine cutoff functions.
    *   `COMPUTE_FUNC` orchestrates the call to the `juslin_kernel`.
```
