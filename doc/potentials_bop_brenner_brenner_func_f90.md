# src/potentials/bop/brenner/brenner_func.f90

## Overview

This file provides the specific mathematical functions that define the interactions for the Brenner Bond-Order Potential (BOP). These functions are designed to be called by the generic `BOP_KERNEL` when the Brenner potential is active. They calculate various components of the potential, such as pair interaction energies, bond order terms, and angular corrections, along with their derivatives. The file also includes `../default_cutoff.f90`, which provides standard cutoff function implementations used to smoothly truncate interactions at a certain distance.

The functions in this file typically take `this` (an instance of `BOP_TYPE`, specifically `brenner_t`) as an argument, which contains the necessary parameters for the Brenner potential, often stored in a database-like structure `this%db`.

## Key Components

### Functions/Subroutines

*   `VA(this, ijpot, dr, val, dval)`
    *   **Type**: Elemental Subroutine
    *   **Description**: Calculates the attractive part of the Morse-style pair potential (`val`) and its derivative with respect to distance (`dval`).
    *   **Arguments**:
        *   `this :: type(BOP_TYPE), intent(in)`: The Brenner potential object containing parameters.
        *   `ijpot :: integer, intent(in)`: Index representing the type of atom pair (e.g., C-C, C-H).
        *   `dr :: real(DP), intent(in)`: The distance between the two atoms.
        *   `val :: real(DP), intent(out)`: The calculated attractive potential energy.
        *   `dval :: real(DP), intent(out)`: The calculated derivative of the attractive potential energy with respect to `dr`.
    *   **Implementation**: `val = -VA_f * exp(-expA * (dr - r0))`, where `VA_f`, `expA`, and `r0` are parameters for the given `ijpot`.

*   `VR(this, ijpot, dr, val, dval)`
    *   **Type**: Elemental Subroutine
    *   **Description**: Calculates the repulsive part of the Morse-style pair potential (`val`) and its derivative (`dval`).
    *   **Arguments**: Similar to `VA`.
    *   **Implementation**: `val = VR_f * exp(-expR * (dr - r0))`, where `VR_f`, `expR`, and `r0` are parameters for the `ijpot`.

*   `g(this, ktypj, ktypi, ktypk, ijpot, ikpot, costh, val, dval_dcosth)`
    *   **Type**: Elemental Subroutine
    *   **Description**: Computes the angular contribution (`val`) to the bond order term for an angle formed by atoms j-i-k. It also returns the derivative of this contribution with respect to `costh` (`dval_dcosth`).
    *   **Arguments**:
        *   `ktypj, ktypi, ktypk :: integer, intent(in)`: Element types of atoms j, i, and k.
        *   `ijpot, ikpot :: integer, intent(in)`: Pair type indices for bonds i-j and i-k.
        *   `costh :: real(DP), intent(in)`: Cosine of the angle j-i-k.
        *   `val :: real(DP), intent(out)`: The calculated angular contribution.
        *   `dval_dcosth :: real(DP), intent(out)`: Derivative of `val` with respect to `costh`.
    *   **Implementation**: `val = gamma * (1 + c_d - c_sq / (d_sq + (h + costh)^2))`, where `gamma`, `c_d`, `c_sq`, `d_sq`, and `h` are parameters for the `ikpot` type.

*   `bo(this, ktypi, ijpot, zij, fcij, faij, bij, dfbij)`
    *   **Type**: Subroutine
    *   **Description**: Calculates the bond order parameter (`bij`) for a bond of type `ijpot` given the summed angular and length contributions from its environment (`zij`). It also computes the derivative of the bond order (`dfbij`).
    *   **Arguments**:
        *   `ktypi :: integer, intent(in)`: Element type of atom i.
        *   `ijpot :: integer, intent(in)`: Pair type index for the bond.
        *   `zij :: real(DP), intent(in)`: Sum of angular (`g`) and length (`h`) dependent terms from neighboring atoms.
        *   `fcij, faij :: real(DP), intent(in)`: Values of cutoff/switching functions for the current bond.
        *   `bij :: real(DP), intent(out)`: The calculated bond order parameter.
        *   `dfbij :: real(DP), intent(out)`: The derivative of `bij` (used in force calculations).
    *   **Implementation**: `bij = (1 + zij^n)^bo_exp` or `(1 + zij)^bo_exp` depending on `n`. `dfbij` involves derivatives of this expression. Parameters `n`, `bo_exp`, `bo_fac`, `bo_exp1` are used.

*   `h(this, ktypj, ktypi, ktypk, ijpot, ikpot, dr, val, dval)`
    *   **Type**: Elemental Subroutine
    *   **Description**: Computes the length-dependent contribution (`val`) to the bond order term. This function typically depends on the difference in bond lengths or a scaled bond length (`dr`). It also returns its derivative (`dval`).
    *   **Arguments**:
        *   `dr :: real(DP), intent(in)`: A measure related to bond length or difference in bond lengths.
        *   Other arguments similar to `g`.
    *   **Implementation**: `val = exp((2*mu*dr)^m)` or similar exponential forms, where `mu` and `m` are parameters for the `ikpot` type. If `mu` is zero, `val` is 1.

*   `Z2pair(this, i, j)`
    *   **Type**: Elemental Function
    *   **Description**: Generates a unique integer index for a given pair of element types `i` and `j`. This index is then used to look up pair-specific parameters.
    *   **Arguments**:
        *   `this :: type(BOP_TYPE), intent(in)`: The Brenner potential object.
        *   `i, j :: integer, intent(in)`: Element types (as integers).
    *   **Returns**: `integer` - the pair index.
    *   **Implementation**: Uses the `PAIR_INDEX(i, j, num_elements)` macro.

### Includes

*   `../default_cutoff.f90`: This file provides standard implementations for cutoff functions (e.g., `fCin` for inner cutoff, `fCar` for attractive/repulsive range, `fCbo` for bond order range). These functions ensure that interactions smoothly go to zero at the cutoff distance.

## Important Variables/Constants

The functions in this file are parameterized by values stored in the `this` argument, specifically within `this%db` (a structure of type `BOP_DB_TYPE`, e.g., `brenner_db_t`). Key parameters include:
*   For pair potentials (`VA`, `VR`): `expA`, `expR` (exponential decay factors), `r0` (equilibrium distance part), `VA_f`, `VR_f` (energy prefactors).
*   For angular function (`g`): `gamma`, `c_d`, `c_sq`, `d_sq`, `h`.
*   For bond order function (`bo`): `n` (exponent for `zij`), `bo_exp` (bond order exponent), `bo_fac` (prefactor for derivative), `bo_exp1` (exponent for derivative).
*   For length contribution (`h`): `mu`, `m` (parameters for exponential term).
*   `this%db%nel`: Number of element types considered in the parameter database.

## Usage Examples

These functions are not meant to be called directly by end-users. They are called by the `BOP_KERNEL` subroutine during the computation of energies and forces. The `BOP_KERNEL` iterates over atoms and their neighbors, calling these functions with appropriate arguments (distances, angles, element types) to evaluate the different parts of the Brenner potential.

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on the `BOP_TYPE` (specifically `brenner_t`) structure for parameters.
    *   Uses the `PAIR_INDEX` macro (likely defined in `macros.inc` or a similar core include file).
    *   Depends on the functions provided by the included `../default_cutoff.f90`.
*   **External Libraries:** None.
*   **Interactions:**
    *   These functions are the building blocks of the Brenner potential energy expression.
    *   They are orchestrated by `BOP_KERNEL` to calculate the total energy and forces.
    *   The `elemental` attribute for many of these routines allows them to be called with scalar or array arguments, which can be useful for vectorized calculations if the compiler supports it well.
```
