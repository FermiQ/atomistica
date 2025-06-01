# src/potentials/bop/default_cutoff.f90

## Overview

This file provides default implementations for various cutoff function calls used within the Bond-Order Potential (BOP) framework in Atomistica. Cutoff functions are essential for smoothly truncating interatomic interactions at a defined range, ensuring that the potential energy and forces go to zero at the cutoff distance. This helps in managing computational cost, especially for systems with many atoms.

The file defines wrapper subroutines for different types of cutoffs:
*   `fCin`: For the inner cutoff, always present.
*   `fCar`: For the attractive/repulsive range cutoff (typically the outer cutoff when screening is active), compiled only if `#ifdef SCREENING` is defined.
*   `fCbo`: For a cutoff specifically applied to bond-order terms, compiled only if `#ifdef SCREENING` is defined.

These subroutines delegate the actual calculation of the cutoff function value and its derivative to a generic `fc` method of a specific cutoff object (of type `CUTOFF_T`) stored within the main `BOP_TYPE` potential object.

## Key Components

### Functions/Subroutines

*   `fCin(this, ijpot, dr, val, dval)`
    *   **Description**: Calculates the value (`val`) and derivative (`dval`) of the **inner cutoff function**. This cutoff is typically applied to all interactions to define their primary range.
    *   **Arguments**:
        *   `this :: type(BOP_TYPE), intent(in)`: The main BOP object, which contains the specific cutoff function object for the inner range (`this%cut_in`).
        *   `ijpot :: integer, intent(in)`: An index representing the type of atom pair (e.g., C-C, C-H), used to select the correct set of cutoff parameters.
        *   `dr :: real(DP), intent(in)`: The interatomic distance at which the cutoff function is to be evaluated.
        *   `val :: real(DP), intent(out)`: The calculated value of the inner cutoff function. This value typically ranges from 1 (well within the interaction range) to 0 (at or beyond the cutoff distance).
        *   `dval :: real(DP), intent(out)`: The derivative of the inner cutoff function with respect to `dr`.
    *   **Implementation**: This subroutine calls `call fc(this%cut_in(ijpot), dr, val, dval)`. The actual calculation is performed by the `fc` type-bound procedure of the `this%cut_in(ijpot)` object, whose specific type is determined by the `CUTOFF_T` macro (e.g., `trig_off_t`, `exp_cutoff_t`).

*   `fCar(this, ijpot, dr, val, dval)` (Compiled only if `#ifdef SCREENING` is defined)
    *   **Description**: Calculates the value (`val`) and derivative (`dval`) of the cutoff function typically associated with the **attractive/repulsive (pair) part** of the potential, especially when screening mechanisms are active. This often corresponds to an "outer" cutoff if different ranges are used for pair and bond-order terms under screening.
    *   **Arguments**: Same as `fCin`.
    *   **Implementation**: Calls `call fc(this%cut_out(ijpot), dr, val, dval)`, using the `this%cut_out(ijpot)` object.

*   `fCbo(this, ijpot, dr, val, dval)` (Compiled only if `#ifdef SCREENING` is defined)
    *   **Description**: Calculates the value (`val`) and derivative (`dval`) of a cutoff function applied specifically to the **bond-order terms** of the potential. This is relevant when screening is enabled and bond-order interactions might have a different cutoff range or shape than the pair terms.
    *   **Arguments**: Same as `fCin`.
    *   **Implementation**: Calls `call fc(this%cut_bo(ijpot), dr, val, dval)`, using the `this%cut_bo(ijpot)` object.

## Important Variables/Constants

*   The functionality of these subroutines depends critically on the `BOP_TYPE` instance (`this`) and the cutoff objects stored within it:
    *   `this%cut_in(:)`: Array of inner cutoff function objects.
    *   `this%cut_out(:)`: Array of outer/pair-term cutoff function objects (if screening is active).
    *   `this%cut_bo(:)`: Array of bond-order term cutoff function objects (if screening is active).
*   The actual mathematical form of the cutoff (e.g., polynomial, exponential) is determined by the specific type of these cutoff objects (defined by the `CUTOFF_T` macro at compile time of the main potential module like `brenner.f90`).

## Usage Examples

These cutoff subroutines are low-level utilities. They are not typically called directly by end-users but are used internally by the `BOP_KERNEL` or by potential-specific function files (like `brenner_func.f90`) during the calculation of interaction energies and forces. The kernel needs these functions to ensure that interactions and their effects diminish smoothly to zero at the defined cutoff radii.

```fortran
! Conceptual usage within BOP_KERNEL or brenner_func.f90:
! REAL(DP) :: distance, cutoff_value, cutoff_derivative
! INTEGER :: pair_type_index
! TYPE(BOP_TYPE) :: potential_object
!
! ! ... obtain distance, pair_type_index, and potential_object ...
!
! CALL fCin(potential_object, pair_type_index, distance, cutoff_value, cutoff_derivative)
!
! ! Use cutoff_value and cutoff_derivative to modulate interaction energy and forces
! energy_term = energy_term * cutoff_value
! force_term  = force_term * cutoff_value + energy_term_raw * cutoff_derivative
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   These subroutines are tightly coupled with the `BOP_TYPE` data structure, as they access its cutoff function objects (`cut_in`, `cut_out`, `cut_bo`).
    *   The behavior depends on the type-bound `fc` procedure of the `CUTOFF_T` objects. Different cutoff types (e.g., `trig_off_t` from `cutoff_types/trig_cutoff.f90` or `exp_cutoff_t` from `cutoff_types/exp_cutoff.f90`) will have different `fc` implementations.
*   **External Libraries:** None.
*   **Interactions:**
    *   These functions provide a consistent interface for applying cutoff functions within the BOP calculations.
    *   The `#ifdef SCREENING` preprocessor directive controls the availability of `fCar` and `fCbo`, tailoring the cutoff mechanism to whether screening effects are being modeled.
    *   The `BOP_KERNEL` relies on these (or similar direct calls to `cutoff_object%fc`) to ensure that all calculated energies and forces are properly zeroed out beyond the interaction range.
```
