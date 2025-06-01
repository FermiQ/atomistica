# src/potentials/dispersion/dispdftd3.f90

## Overview

The `dispdftd3` module provides an interface to Grimme's DFT-D3 method for calculating dispersion corrections to energy and forces. This method is widely used to add van der Waals interactions, which are often missing or poorly described in standard Density Functional Theory (DFT) functionals or some classical potentials.

This module acts as a wrapper around an external library, `dftd3-lib` (expected to be available if the `HAVE_DFTD3` preprocessor macro is defined). It handles parameter settings, unit conversions, and communication with the DFT-D3 library to obtain the dispersion energy, atomic forces, and stress tensor contributions.

The implementation supports:
*   Becke-Johnson (BJ) damping or zero-damping schemes.
*   Optional inclusion of three-body interaction terms.
*   Customizable parameters for the D3 model (e.g., `s6`, `s8`, `a1`, `a2`).

References:
*   Grimme et al., J. Chem. Phys. 132, 154104 (2010)
*   Grimme et al., J. Comp. Chem. 32, 1456 (2011)

## Key Components

### Modules

*   `dispdftd3`
    *   **Uses**: `libAtoms_module`, `ptrdict`, `logging`, `timer`, `particles`, `neighbors`, `filter`.
    *   **Uses (conditional)**: `dftd3_api` (if `HAVE_DFTD3` is defined). This is the Fortran API for the external `dftd3-lib`.

### Data Types

*   `dispdftd3_t` (Public)
    *   **Description**: Holds the parameters and state for the DFT-D3 dispersion correction.
    *   **Fields**:
        *   `a1, a2 :: real(DP)`: Parameters for Becke-Johnson damping. Defaults: `0.5719`, `3.6017`.
        *   `s6, s8 :: real(DP)`: Global scaling factors for \(C_6\) and \(C_8\) terms. Defaults: `1.0000`, `0.5883`.
        *   `sr6, sr8 :: real(DP)`: Damping parameters for the zero-damping scheme (used if `BeckeJohnson` is false). Defaults: `0.7461`, `1.0000`.
        *   `alpha6 :: real(DP)`: Damping parameter for the zero-damping scheme. Default: `14.000`.
        *   `cutoff :: real(DP)`: Cutoff for two-body dispersion terms. Default: `80.0` (units depend on Atomistica's system, converted to Bohr for the library).
        *   `cutoffCN :: real(DP)`: Cutoff for coordination numbers in the three-body term. Default: `40.0` (units converted to Bohr).
        *   `BeckeJohnson :: logical`: If true, use Becke-Johnson damping; otherwise, use zero-damping. Default: `.true.`.
        *   `threebody :: logical`: If true, include the three-body dispersion term. Default: `.false.`.
        *   `calculator :: type(dftd3_calc), allocatable`: (If `HAVE_DFTD3`) An object from the `dftd3_api` module that encapsulates the DFT-D3 calculation state and parameters for the external library.

### Public Subroutines & Interfaces

*   `init(this)` (maps to `dispdftd3_init`)
    *   **Description**: Constructor. If `HAVE_DFTD3` is defined, it allocates `this%calculator`. It then prepares a `dftd3_input` structure (setting `threebody`, `numgrad=.false.`, cutoffs converted to Bohr units) and calls `dftd3_init(this%calculator, input)` from the DFT-D3 library. Finally, it calls `dftd3_set_params` to pass the chosen damping parameters (`s6,a1,s8,a2` for BJ or `s6,sr6,s8,sr8,alpha6` for zero-damping) to the library.
*   `del(this)` (maps to `dispdftd3_del`)
    *   **Description**: Destructor. If `HAVE_DFTD3` and `this%calculator` is allocated, it deallocates `this%calculator` (which should internally call the library's cleanup for that calculator instance).
*   `bind_to(this, p, nl, ierror)` (maps to `dispdftd3_bind_to`)
    *   **Description**: Binds the potential to particle and neighbor list data. If `HAVE_DFTD3` is not defined, this routine raises an error. Otherwise, it primarily logs the current DFT-D3 parameters. It does not call `request_interaction_range` as the DFT-D3 library handles its own interaction ranges based on its internal cutoff parameters.
*   `energy_and_forces(this, p, nl, epot, f, wpot, mask, epot_per_at, wpot_per_at, ierror)` (maps to `dispdftd3_energy_and_forces`)
    *   **Description**: Calculates the DFT-D3 dispersion energy, forces, and virial contribution.
        1.  Checks if `HAVE_DFTD3` is defined.
        2.  Allocates local arrays `coords` (for positions) and `grads` (for gradients/forces).
        3.  Copies particle atomic numbers (`p%Z`) to `elem` and particle positions (`POS3(p,...)`) to `coords`, converting positions to Bohr radii. Cell vectors `p%Abox` are also converted to Bohr.
        4.  Determines unit conversion factors (`unit_energy`, `unit_forces`, `unit_virial`) to convert results from atomic units (used by `dftd3-lib`) to the current Atomistica system of units.
        5.  Calls `dftd3_pbc_dispersion(this%calculator, coords, elem, latvecs, edisp, grads, stress)` from the DFT-D3 library. This is the core call that performs the dispersion calculation.
        6.  Adds the returned dispersion energy `edisp` (after unit conversion) to `epot`.
        7.  Subtracts the returned gradients `grads` (after unit conversion, as forces are \(-\nabla E\)) from `f`.
        8.  Subtracts the returned stress tensor `stress` (after unit conversion and scaling by volume) from `wpot`.
        9.  Deallocates local `coords` and `grads`.
*   `register(this, cfg, m)` (maps to `dispdftd3_register`)
    *   **Description**: Registers all tunable DFT-D3 parameters (`BeckeJohnson`, `threebody`, `a1`, `a2`, `s6`, `s8`, `sr6`, `sr8`, `alpha6`, `cutoff`, `cutoffCN`) with the `ptrdict` configuration system.

## Important Variables/Constants
*   The various parameters controlling the D3 model: `s6, s8` (global scaling), `a1, a2` (BJ damping), `sr6, sr8, alpha6` (zero-damping).
*   `BeckeJohnson`, `threebody`: Logical flags to switch model features.
*   `cutoff`, `cutoffCN`: Cutoff radii for two-body and three-body terms.
*   Unit conversion factors: `length_to_Bohr`, `Hartree`.

## Usage Examples
This module is typically added to a primary energy model (like a DFT calculation or a classical potential that lacks dispersion) to include van der Waals interactions.

```fortran
! Conceptual usage with another potential:
! USE dispdftd3_module
! TYPE(dispdftd3_t) :: d3_correction
! ! ... other potential setup (e.g., DFT) ...
!
! CALL d3_correction%init()
! ! Optionally, modify parameters from defaults if not using ptrdict:
! ! d3_correction%s6 = 1.1_DP
! ! ...
! CALL d3_correction%bind_to(particles, neighbors)
!
! ! In the main energy/force call:
! ! CALL dft_potential%energy_and_forces(..., E_dft, F_dft, V_dft)
! CALL d3_correction%energy_and_forces(particles, neighbors, E_d3, F_d3, V_d3)
!
! ! E_total = E_dft + E_d3
! ! F_total = F_dft + F_d3
! ! V_total = V_dft + V_d3
!
! CALL d3_correction%del()
```

## Dependencies and Interactions
*   **External Library**: Critically depends on the `dftd3-lib` being available and Atomistica being compiled with `HAVE_DFTD3` defined. If the library is not available, using this module will result in an error at the `bind_to` stage.
*   **Atomistica Core Modules**: `libAtoms_module` (for `system_of_units`, unit constants), `ptrdict`, `logging`, `timer`, `particles`, `neighbors`, `filter`.
*   **Unit Conversion**: The module handles conversion of input coordinates and cell parameters to Bohr (for `dftd3-lib`) and conversion of output energy, forces, and stress back to the active unit system in Atomistica.
*   It acts as an additive correction; the primary interactions are expected to be calculated by another potential module.
```
