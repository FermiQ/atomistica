# src/potentials/pair_potentials/born_mayer.f90

## Overview

The `born_mayer` module implements the Born-Mayer potential, a simple pair potential typically used to model short-range repulsive interactions between ions in ionic crystals, or as a general-purpose repulsive core in other contexts. The functional form of the Born-Mayer potential is:

\[ V(r_{ij}) = A \exp\left(-\frac{r_{ij}}{\rho}\right) \]

where \(r_{ij}\) is the distance between atoms \(i\) and \(j\), \(A\) is an energy prefactor, and \(\rho\) is a characteristic length scale determining the steepness of the repulsion.

This implementation applies the potential between two specified element types and includes a cutoff mechanism where the potential is shifted to be zero at the `cutoff` distance.

## Key Components

### Modules

*   `born_mayer`
    *   **Uses**: `libAtoms_module`, `ptrdict`, `logging`, `timer`, `particles`, `neighbors`, `filter`.

### Data Types

*   `born_mayer_t` (Public)
    *   **Description**: Holds the parameters for the Born-Mayer potential.
    *   **Fields**:
        *   `A :: real(DP)`: The energy prefactor for the exponential repulsion. Default: `1.0`.
        *   `rho :: real(DP)`: The characteristic length scale in the exponent. Default: `1.0`.
        *   `cutoff :: real(DP)`: The cutoff distance for the interaction. Default: `1.0`.
        *   `element1 :: character(MAX_EL_STR)`: The chemical symbol of the first element type involved in this specific Born-Mayer interaction. Default: "C".
        *   `element2 :: character(MAX_EL_STR)`: The chemical symbol of the second element type. Default: "C".
        *   `el1 :: integer`: Internal filter integer corresponding to `element1`.
        *   `el2 :: integer`: Internal filter integer corresponding to `element2`.
        *   `shift :: real(DP)`: The value of the potential at the `cutoff` distance (\(A \exp(-\text{cutoff}/\rho)\)), used to shift the calculated potential to zero at and beyond the cutoff.

### Public Subroutines & Interfaces

*   `init(this, element1, element2, A, rho, cutoff)` (maps to `born_mayer_init`)
    *   **Description**: Constructor. Initializes the `born_mayer_t` object with the specified parameters: `element1`, `element2`, `A`, `rho`, and `cutoff`.
*   `del(this)` (maps to `born_mayer_del`)
    *   **Description**: Destructor. This is currently an empty subroutine as `born_mayer_t` does not dynamically allocate memory that it owns.
*   `bind_to(this, p, nl, ierror)` (maps to `born_mayer_bind_to`)
    *   **Description**: Binds the potential to a particle set (`p`) and a neighbor list handler (`nl`).
        1.  Logs the potential parameters (`A`, `rho`, `cutoff`, element types).
        2.  Converts the `element1` and `element2` character strings into integer filter IDs (`this%el1`, `this%el2`) using `filter_from_string`.
        3.  Calls `request_interaction_range(nl, this%cutoff)` to inform the neighbor list system of the interaction range.
        4.  Calculates `this%shift = this%A * exp(-this%cutoff / this%rho)`.
*   `energy_and_forces(this, p, nl, epot, for, wpot, ...)` (maps to `born_mayer_energy_and_forces`)
    *   **Description**: Calculates the potential energy and forces due to the Born-Mayer interaction between specified element pairs.
        1.  Iterates through each local atom `i`.
        2.  Checks if atom `i` matches `this%el1` or `this%el2`.
        3.  For each neighbor `j` of atom `i`:
            *   Checks if the pair `(i,j)` is of type (`el1,el2`) or (`el2,el1`).
            *   If the distance `abs_dr` is less than `this%cutoff`:
                *   Calculates the energy contribution: `e_contrib = this%A * exp(-abs_dr/this%rho) - this%shift`.
                *   Calculates the force contribution: `f_contrib_vec = (this%A / this%rho) * exp(-abs_dr/this%rho) * (dr_vec / abs_dr)`.
                *   Adds `e_contrib` to the total potential energy `epot`.
                *   Adds `f_contrib_vec` to `for(i)` and subtracts it from `for(j)`.
    *   **Note**: This routine currently **does not calculate the virial tensor** (`wpot`). `wpot` is an argument but is not modified.
*   `register(this, cfg, m)` (maps to `born_mayer_register`)
    *   **Description**: Registers the parameters `A`, `rho`, `cutoff`, `element1`, and `element2` with the `ptrdict` configuration system.

## Important Variables/Constants

*   `A`: Energy prefactor of the potential.
*   `rho`: Characteristic length scale of the repulsion.
*   `cutoff`: Distance at which the interaction is truncated (after shifting).
*   `element1`, `element2`: Define the specific pair of element types this potential instance applies to.

## Usage Examples

The Born-Mayer potential is typically used for repulsive interactions in ionic models or as a general short-range repulsion.

```fortran
! Conceptual usage:
! USE born_mayer_module
! TYPE(born_mayer_t) :: bm_potential
!
! ! Define parameters for a C-O repulsion
! CALL bm_potential%init(element1="C", element2="O", A=1000.0_DP, rho=0.2_DP, cutoff=5.0_DP)
!
! ! ... setup particles (p_data), neighbor_list_handler (nl_data) ...
! CALL bm_potential%bind_to(p_data, nl_data)
!
! ! Calculate energy and forces
! CALL bm_potential%energy_and_forces(p_data, nl_data, E_total, F_array, V_tensor)
! ! Note: V_tensor will not be updated by this specific potential
!
! CALL bm_potential%del()
```

## Dependencies and Interactions

*   Uses standard Atomistica modules: `libAtoms_module`, `ptrdict`, `logging`, `timer`, `particles`, `neighbors`, `filter`.
*   The interaction is pairwise and only applies between atoms matching the specified `element1` and `element2`.
*   The potential is shifted to be zero at the cutoff, which ensures continuity of energy but not necessarily of the force if the force at cutoff is non-zero before shifting (though for a purely repulsive potential like Born-Mayer, the force would be small if the potential is near zero).
```
