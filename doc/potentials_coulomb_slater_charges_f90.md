# src/potentials/coulomb/slater_charges.f90

## Overview

The `slater_charges` module implements corrections to electrostatic interactions arising from the use of Slater-type (exponential) atomic charge distributions instead of simple point charges. The charge density around an atom `i` is modeled as a point-like effective nuclear charge \(Z_i\) plus an electronic part \((q_i-Z_i)\) distributed according to a Slater orbital \(f_i(\vec{r}-\vec{r_i})\).

This module calculates the **correction terms** to the energy, potential, and forces due to these finite-sized, overlapping charge distributions. It does **not** calculate the standard \(1/r\) point-charge Coulomb interaction itself but provides the necessary adjustments. Thus, it's designed to be used in conjunction with a primary Coulomb solver (like Ewald sum or direct summation for the \(1/r\) part).

The module is essential for methods like Density Functional Tight Binding (DFTB), particularly:
*   **DFTB2-like models**: Can include damping of interactions involving hydrogen (`damp_gamma` flag) via functions from `damp_short_gamma`.
*   **DFTB3-like models**: Can include terms related to the derivative of Hubbard U parameters with respect to charge (`dftb3` flag and `dU` parameters) using functions from `coulomb_short_gamma`.

## Key Components

### Modules

*   `slater_charges`
    *   **Uses**: `supplib`, `particles`, `neighbors`, `filter`.
    *   **Uses (for DFTB2/3 features)**: `coulomb_short_gamma`, `damp_short_gamma`.

### Data Types

*   `slater_charges_db_t` (Internal to the module)
    *   **Description**: Stores database parameters before processing for specific particle types.
    *   **Fields**:
        *   `nel, nU, nZ :: integer`: Counts of elements, U values, and Z values.
        *   `el(2, SLATER_CHARGES_MAX_EL) :: character`: Element symbols.
        *   `U(SLATER_CHARGES_MAX_EL) :: real(DP)`: Hubbard U-like parameters (related to Slater orbital exponents \(\zeta_{slater}\)).
        *   `Z(SLATER_CHARGES_MAX_EL) :: real(DP)`: Effective nuclear charges.
        *   `dU(0:116) :: real(DP)`: Derivative of Hubbard U with respect to atomic charge, indexed by atomic number (for DFTB3). Default: 0.0.

*   `slater_charges_t` (Public)
    *   **Description**: Holds parameters and state for the Slater charges potential component.
    *   **Fields**:
        *   `elements :: character(MAX_EL_STR)`: Element filter string.
        *   `els :: integer`: Compiled element filter.
        *   `cutoff, cutoff_sq :: real(DP)`: Cutoff radius for these correction terms. Default: 5.0.
        *   `dftb3 :: logical`: Flag to enable DFTB3 specific terms. Default: `.false.`.
        *   `damp_gamma :: logical`: Flag to enable XH damping (DFTB2-like). Default: `.false.`.
        *   `zeta :: real(DP)`: Exponent for the damping function from `damp_short_gamma`. Default: 4.0.
        *   `U(:) :: real(DP), allocatable`: Processed Slater exponents \(\zeta_{slater}\) for each element type in the simulation. Derived from `db%U` after unit conversion and scaling (`db%U / (Hartree*Bohr) * 16.0/5.0`).
        *   `Z(:) :: real(DP), allocatable`: Effective nuclear charges for each element type.
        *   `dU(:) :: real(DP), allocatable`: Processed \(dU/dq\) for each element type, in atomic units.
        *   `db :: type(slater_charges_db_t)`: Internal storage for database parameters.

### Public Subroutines & Interfaces

*   `init(this, p, U, elements, cutoff, error)`: Constructor.
*   `del(this)`: Destructor (deallocates `U`, `Z`, `dU`).
*   `set_Hubbard_U(this, p, U, Z, error)`: Sets `db%U` and `db%Z` (optional, defaults to 0). `U` and `Z` are arrays per element type in `p`.
*   `bind_to(this, p, nl, ierror)`: Binds to particle/neighbor data. Processes filters, validates DB counts, maps `db%U, db%Z, db%dU` to `this%U, this%Z, this%dU` including unit conversions and the \(U_{Hubbard} \rightarrow \zeta_{slater}\) scaling for `this%U`. Requests interaction range.
*   `potential(this, p, nl, q, phi, ierror)`: Calculates corrections to the electrostatic potential `phi`.
    *   Includes nuclear-electron attraction terms based on Slater overlaps (e.g., `~ -Z_j * (0.5*U_i + 1/r_ij) * exp(-U_i*r_ij)`).
    *   Includes electron-electron repulsion terms between Slater charge clouds \((q_i-Z_i)\) and \((q_j-Z_j)\), with functional forms depending on whether `U_i` and `U_j` are similar.
    *   Applies `hij` damping from `damp_short_gamma` to e-e terms if `damp_gamma` is true and H is involved.
    *   Adds DFTB3 corrections involving `capital_short_gamma` and `dU` if `dftb3` is true.
    *   Includes an on-site term `5/16*q_i*U_i` (using the scaled `this%U` which is \(\zeta_{slater}\)). If `dftb3`, also adds `-0.5*(q_i-Z_i)^2*dU_i`.
*   `energy_and_forces(this, p, nl, q, epot, f, wpot, error)`: Calculates corrections to energy, forces, and virial, following similar logic as `potential` for the different interaction components (nuclear-electron, electron-electron, DFTB2 damping, DFTB3 terms, on-site).
*   `register(this, cfg, m)`: Registers parameters including `elements`, `cutoff`, `dftb3`, `damp_gamma`, `zeta`, `db%el`, `db%U`, `db%Z`, and the array `db%dU` (with per-element entries).

## Important Variables/Constants

*   `U` parameter in `slater_charges_t`: Represents the Slater orbital exponent \(\zeta_{slater}\). It's derived from input "Hubbard U" values, typically via a relation like \(U_{Hubbard} = \frac{5}{16} \zeta_{slater}\) (in atomic units). The code applies this scaling: `this%U = this%db%U / (Hartree*Bohr) * 16.0/5.0`.
*   `Z`: Effective nuclear charge for an atom type.
*   `dU`: Derivative \(dA/dq_A\) (where A is related to Hubbard U) used in DFTB3.
*   `dftb3`, `damp_gamma`, `zeta`: Flags and parameter controlling DFTB3/DFTB2 specific behaviors.

## Usage Examples

This module is intended for use in DFTB or similar variable charge/tight-binding models.

```fortran
! Conceptual usage within a DFTB simulation:
! USE slater_charges_module
! TYPE(slater_charges_t) :: slater_corr
! ! ... setup particles, neighbors, charges q ...
! ! ... load U, Z, dU parameters from a database into slater_corr%db ...
!
! CALL slater_corr%init(cutoff=..., elements="C,H")
! CALL slater_corr%bind_to(particles, neighbors) ! This processes db params into this%U, etc.
!
! ! Calculate corrections
! CALL slater_corr%energy_and_forces(particles, neighbors, q, E_slater, F_slater, V_slater)
!
! ! E_total = E_point_charge_coulomb + E_slater + E_short_range_repulsion + ...
```

## Dependencies and Interactions

*   **`supplib`**: For constants (`Hartree`, `Bohr`, `PI`) and utilities.
*   **`particles`, `neighbors`, `filter`**: For core simulation infrastructure.
*   **`coulomb_short_gamma`**: Used if `dftb3 = .true.` for specific DFTB3 energy/potential terms.
*   **`damp_short_gamma`**: Used if `damp_gamma = .true.` to apply damping to interactions involving hydrogen.
*   This module calculates **corrections** to a primary point-charge Coulomb interaction. The main \(1/r\) term should be handled by another module (e.g., `pme`, `direct_coulomb`).
*   The interaction terms are complex, involving exponential functions and error functions (implicitly via `coulomb_short_gamma` which uses them, although this module's direct electron-electron terms are simpler exponentials).
```
