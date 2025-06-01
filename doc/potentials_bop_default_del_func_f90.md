# src/potentials/bop/default_del_func.f90

## Overview

This file provides a default implementation for the `DEL_FUNC` subroutine, a standardized interface function used by various Bond-Order Potentials (BOPs) within the Atomistica simulation package. When a specific BOP module (e.g., `brenner.f90`, `tersoff.f90`) is compiled, its `DEL_FUNC` macro (which defines the public destructor like `brenner_del`) typically points to this default implementation.

The primary responsibility of this subroutine is to act as a destructor for `BOP_TYPE` objects. It handles the deallocation of dynamically allocated arrays that were used to store internal neighbor lists and related bond data during computations. This is crucial for proper memory management and preventing memory leaks.

## Key Components

### Functions/Subroutines

*   `DEL_FUNC(this)`
    *   **Description**: This subroutine checks if the internal neighbor list arrays within the `BOP_TYPE` object `this` have been allocated. This is determined by the boolean flag `this%neighbor_list_allocated`. If the flag is true, the subroutine proceeds to deallocate all known allocatable arrays associated with these lists. After deallocation, the `neighbor_list_allocated` flag is set back to `.false.`.
    *   **Deallocated Arrays**:
        *   Core neighbor and bond indexing: `this%neb`, `this%nbb`.
        *   Displacement vectors for periodic images: `this%dcell` (deallocation is skipped if compiled for LAMMPS via `#ifndef LAMMPS`).
        *   Properties of each bond: `this%bndtyp` (bond type), `this%bndlen` (bond length), `this%bndnm` (normalized bond vector).
        *   Stored values of cutoff functions and their derivatives: `this%cutfcnar`, `this%cutdrvar`.
    *   **Conditionally Deallocated Arrays** (if `#ifdef SCREENING` is defined):
        *   Cutoff values for bond-order terms: `this%cutfcnbo`, `this%cutdrvbo`.
        *   Data structures for screening neighbors: `this%sneb_seed`, `this%sneb_last` (indices for screening neighbor lists), `this%sneb` (screening atom indices), `this%sbnd` (screening bond indices).
        *   Screening-related force components: `this%sfacbo`.
        *   Derivatives of cutoff functions with respect to screening atom positions: `this%cutdrarik`, `this%cutdrarjk`, `this%cutdrboik`, `this%cutdrbojk`.
    *   **Key Arguments**:
        *   `this :: type(BOP_TYPE), intent(inout)`: The Bond-Order Potential object whose dynamically allocated internal arrays are to be deallocated.

## Important Variables/Constants

*   The subroutine primarily interacts with the `this` object of `BOP_TYPE`.
*   `this%neighbor_list_allocated :: logical`: A flag that tracks whether the internal arrays are currently allocated. This prevents errors from attempting to deallocate already deallocated memory or unallocated arrays.
*   The preprocessor directives `#ifndef LAMMPS` and `#ifdef SCREENING` control which arrays are deallocated, ensuring the destructor matches the compiled features and structure of `BOP_TYPE`.

## Usage Examples

This `DEL_FUNC` is not typically called directly by an end-user. It is the subroutine that gets executed when a user calls the public destructor of a specific BOP module (e.g., `call brenner_del(my_potential)`).

```fortran
! Conceptual usage within a specific potential module (e.g., brenner.f90)
! MODULE brenner_module
!   ! ... other definitions ...
!   PUBLIC :: brenner_del
!   INTERFACE brenner_del
!     MODULE PROCEDURE DEL_FUNC ! DEL_FUNC is this default implementation
!   END INTERFACE
!   ! ...
! END MODULE brenner_module

! User code:
! USE brenner_module
! TYPE(brenner_t) :: my_brenner_potential
! ! ... initialize and use my_brenner_potential ...
! CALL brenner_del(my_brenner_potential) ! This call executes the default DEL_FUNC
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on the structure of `BOP_TYPE` as defined in files like `brenner_type.f90`, particularly its allocatable array members and the `neighbor_list_allocated` flag.
*   **External Libraries:** None explicitly used.
*   **Interactions:**
    *   This subroutine is essential for cleaning up memory used by a `BOP_TYPE` instance, typically after it's no longer needed or before it's re-initialized.
    *   The arrays deallocated here are usually allocated within the `BOP_KERNEL` during its first execution or potentially within `BIND_TO_FUNC` if some pre-sizing is done.
    *   Properly calling this destructor prevents memory leaks in long simulations or when frequently changing potentials.
```
