# src/potentials/bop/kumagai/kumagai_func.f90

## Overview

This file provides the specific mathematical functions that define the interactions for the Kumagai Bond-Order Potential (BOP), as described in Kumagai et al., Comp. Mater. Sci. 39, 457 (2007). These functions are designed to be called by the generic `BOP_KERNEL` when the Kumagai potential is active. The file includes implementations for attractive and repulsive pair potentials (`VA`, `VR`), the angular contribution to the bond order (`g`), the bond order function (`bo`), and the length-dependent contribution to bond order (`h`). Default cutoff function implementations are included via `../default_cutoff.f90`.

## Key Components

### Included Files
*   `../default_cutoff.f90`: This file is included to provide standard cutoff function implementations (like `fCin`, and conditionally `fCar`, `fCbo` if screening is enabled for a Kumagai variant). These functions ensure smooth truncation of interactions.

### Pair Potentials

*   `VA(this, ijpot, dr, val, dval)`
    *   **Type**: Subroutine
    *   **Description**: Calculates the attractive part (`val`) of the pair potential and its derivative (`dval`) with respect to distance `dr`.
    *   **Implementation**: `val = -B(ijpot) * exp(-lambda2(ijpot) * dr)`. Parameters `B` and `lambda2` are specific to the atom pair type `ijpot` and are retrieved from `this%db`.
*   `VR(this, ijpot, dr, val, dval)`
    *   **Type**: Subroutine
    *   **Description**: Calculates the repulsive part (`val`) of the pair potential and its derivative (`dval`).
    *   **Implementation**: `val = A(ijpot) * exp(-lambda1(ijpot) * dr)`. Parameters `A` and `lambda1` are specific to the atom pair type `ijpot` and are retrieved from `this%db`.

### Angular Contribution to Bond Order

*   `g(this, ktypj, ktypi, ktypk, ijpot, ikpot, costh, val, dval_dcosth)`
    *   **Type**: Subroutine
    *   **Description**: Computes the angular contribution (`val`) to the bond order term for an angle formed by atoms j-i-k, where `i` is the central atom. It also returns the derivative of this contribution with respect to `costh` (`dval_dcosth`).
    *   **Implementation**: The function uses a complex form: `val = c1 + h_cos * go * (1 + ga1)`, where `go = c2 * tmp`, `tmp = h_cos / (c3 + h_cos_sq)`, `ga1 = c4 * exp(-c5 * h_cos_sq)`, and `h_cos = h - costh`. The parameters `c1, c2, c3, c4, c5, h` are specific to the central atom type `ktypi` and are retrieved from `this%db`.

### Bond Order Function

*   `bo(this, ktypi, ijpot, zij, fcij, faij, bij, dfbij)`
    *   **Type**: Subroutine
    *   **Description**: Calculates the bond order parameter (`bij`) for a bond of type `ijpot` given the summed angular and length contributions from its environment (`zij`). It also computes the derivative of the bond order (`dfbij`). `fcij` and `faij` are cutoff function values.
    *   **Implementation**: `bij = (1 + zij^eta)^delta_neg`, where `eta` and `delta_neg` (negative of `delta` parameter) are specific to the central atom type `ktypi` and retrieved from `this%db`. If `zij <= 0`, `bij = 1` and `dfbij = 0`.

### Length-Dependent Contribution to Bond Order

*   `h(this, ktypj, ktypi, ktypk, ijpot, ikpot, dr, val, dval)`
    *   **Type**: Subroutine
    *   **Description**: Computes the length-dependent contribution (`val`) to the bond order term for the bond i-k (where `i` is the central atom of the j-i-k angle). It also returns its derivative (`dval`).
    *   **Implementation**: The function takes the form `val = exp(alpha * dr^beta_int)` or variants depending on the integer parameter `beta_int`. If `alpha` is zero, `val = 1`. Parameters `alpha` (real) and `beta` (integer) are specific to the i-k pair type (`ikpot`) and retrieved from `this%db`.

### Pair Indexing

*   `Z2pair(this, ktypi, ktypj)`
    *   **Type**: Function
    *   **Description**: Generates a unique integer index for a given pair of (internal) element types `ktypi` and `ktypj`.
    *   **Implementation**: Uses the `PAIR_INDEX(ktypi, ktypj, this%db%nel)` macro. `this%db%nel` provides the number of elements.

## Important Variables/Constants

These functions are parameterized by values stored in the `this` argument (an instance of `kumagai_t`, via the `BOP_TYPE` macro), specifically within `this%db`:
*   Pair parameters (`A, B, lambda1, lambda2`) for `VA` and `VR` are per `ijpot` pair type.
*   Angular parameters (`c1, c2, c3, c4, c5, h`) for `g` are per central atom type `ktypi`.
*   Bond order parameters (`eta, delta`) for `bo` are per central atom type `ktypi`.
*   Length contribution parameters (`alpha, beta`) for `h` are per `ikpot` pair type.

## Usage Examples

These functions are low-level routines called by the `BOP_KERNEL` (specifically `kumagai_kernel`) during the computation of energies and forces. They are not intended for direct use by end-users.

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Relies on the `BOP_TYPE` (`kumagai_t`) structure for parameters.
    *   Uses the `PAIR_INDEX` macro.
    *   Includes `../default_cutoff.f90` for standard cutoff function implementations (`fCin`, etc.), which are then used by the `BOP_KERNEL`.
*   **External Libraries:** None.
*   **Interactions:**
    *   These functions collectively define the energy landscape according to the Kumagai potential.
    *   The `BOP_KERNEL` orchestrates calls to these functions based on atomic geometries and neighbor lists.
```
