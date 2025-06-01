# src/potentials/coulomb/gaussian_charges.f90

## Overview

The `gaussian_charges` module implements the calculation of corrections to Coulomb interactions arising from the finite size of charge distributions, modeled as Gaussian functions. Instead of point charges, each atom `i` can be assigned a Gaussian charge distribution \(\rho_i(\vec{r}) = q_i f_i(\vec{r}-\vec{r_i})\), where \(f_i\) is a Gaussian function whose width is related to an effective Hubbard U-like parameter.

This module calculates the **difference** in electrostatic energy and potential due to these Gaussian shapes compared to a point-charge model. It does **not** compute the standard long-range \(1/r\) Coulomb interaction itself; rather, it provides the short-range correction terms. Therefore, it's typically used in conjunction with another Coulomb solver (like Ewald sum or direct summation for point charges).

This functionality is noted as being required for Tight-Binding (DFTB) and Variable Charge model implementations. For pure tight-binding, the nuclear charges \(Z_i\) might be considered zero in this context, and only the electronic charge distributions contribute.

## Key Components

### Modules

*   `gaussian_charges`
    *   **Uses**: `supplib` (for constants like `Hartree`, `Bohr`, and utilities), `particles`, `neighbors`, `filter`. Implicitly uses `ptrdict` for parameter registration.

### Data Types

*   `gaussian_charges_db_t` (Internal to the module)
    *   **Description**: Stores the database of parameters before they are processed for specific particle types.
    *   **Fields**:
        *   `nel, nU :: integer`: Number of element types and U values defined in the database.
        *   `el(2, GAUSSIAN_CHARGES_MAX_EL) :: character`: Element symbols (e.g., "Si").
        *   `U(GAUSSIAN_CHARGES_MAX_EL) :: real(DP)`: Hubbard U-like parameters for each element type in the database. These are related to the inverse width of the Gaussian.

*   `gaussian_charges_t` (Public)
    *   **Description**: Holds the parameters and state for the Gaussian charges potential component.
    *   **Fields**:
        *   `elements :: character(MAX_EL_STR)`: A string to filter which elements this potential component applies to (default "*").
        *   `els :: integer`: The compiled element filter from `elements`.
        *   `cutoff, cutoff_sq :: real(DP)`: Cutoff radius for the Gaussian correction terms (and its square). Default: 5.0.
        *   `U(:) :: real(DP), allocatable`: Array storing the processed Hubbard U-like parameters for each actual element type present in the simulation. These are derived from `db%U` and converted to atomic units (inverse Bohr radii, by dividing by `Hartree*Bohr`). A value of `U(type) > 0` indicates atom type `type` has a Gaussian distribution.
        *   `db :: type(gaussian_charges_db_t)`: Internal storage for database parameters.

### Public Subroutines & Interfaces

*   `init(this, p, U, elements, cutoff, error)`: Constructor. Initializes `elements` filter string, `cutoff`. Optionally, if `p` (particles) and `U` (array of U values per element type in `p`) are provided, it calls `set_Hubbard_U`.
*   `del(this)`: Destructor. Deallocates `this%U`.
*   `set_Hubbard_U(this, p, U, error)`: Sets the Hubbard U values in `this%db%U`. `U` is an array corresponding to the element types in `p`. This must be called before `bind_to`.
*   `bind_to(this, p, nl, ierror)`: Binds the potential to particle data `p` and neighbor list `nl`.
    1.  Processes the `elements` filter string.
    2.  Validates `db%nel` and `db%nU`.
    3.  Allocates `this%U` based on the number of element types in `p`.
    4.  Maps parameters from `this%db%el` and `this%db%U` to `this%U(j)` for each particle type `j` in `p`, converting `U` values by dividing by `(Hartree*Bohr)`.
    5.  Sets `this%cutoff_sq`.
    6.  Calls `request_interaction_range(nl, this%cutoff)`.
*   `potential(this, p, nl, q, phi, ierror)`: Calculates the correction to the electrostatic potential `phi(i)` at each atom `i` due to the Gaussian shape of its own charge `q(i)` and surrounding charges `q(j)`.
    *   For a pair i-j within `cutoff`:
        *   If both i and j are Gaussian: adds `q(j) * -erfc(sqrt(PI/2*Ui^2*Uj^2/(Ui^2+Uj^2)) * rij) / rij` to `phi(i)` (and symmetrically for `phi(j)`).
        *   If i is Gaussian, j is point: adds `q(j) * -erfc(sqrt(PI/2)*Ui*rij) / rij` to `phi(i)` (and symmetrically).
    *   A self-term `q(i)*this%U(p%el(i))` is added to `phi(i)`. (This is related to `erf(0)/0` limit and self-energy).
    *   Uses OpenMP for parallelization.
*   `energy_and_forces(this, p, nl, q, epot, f, wpot, error)`: Calculates the energy correction `epot` and corresponding forces `f` and virial `wpot`.
    *   For each pair i-j within `cutoff` where at least one atom is Gaussian, it adds a term `q(i)*q(j)*hlp` to `epot`, where `hlp` is `-erfc(effective_width * rij) / rij`. The derivative of this term gives the force correction.
    *   A self-energy term `0.5*q(i)*q(i)*this%U(p%el(i))` is added to `epot` for each Gaussian atom `i`.
*   `register(this, cfg, m)`: Registers `elements`, `cutoff`, and the database parameters `db%el` and `db%U` with `ptrdict`.

## Important Variables/Constants

*   `U`: The Hubbard U-like parameter, inversely related to the Gaussian width. A positive `U` value for an atom type signifies a Gaussian distribution; zero or negative might imply a point charge for this correction model. It undergoes unit conversion to `1/(Hartree*Bohr)`.
*   `cutoff`: Radius at which the Gaussian correction terms are truncated.
*   `erfc`: The complementary error function, which arises from the electrostatic interaction of two Gaussian charge distributions.
*   `Hartree`, `Bohr`: Physical constants from `supplib` used for unit conversion of `U`.

## Usage Examples

This module is typically used in DFTB or variable charge simulations where accounting for the spatial extent of charge distributions is important. It provides corrections that are added to a primary point-charge Coulomb solver.

```fortran
! Conceptual usage within a DFTB code:
! USE gaussian_charges_module
! TYPE(gaussian_charges_t) :: gc_correction
! TYPE(particles_t) :: atoms
! TYPE(neighbors_t) :: nl_handler
! REAL(DP), ALLOCATABLE :: charges(:), potential_correction(:), force_correction(:,:)
! REAL(DP) :: energy_correction, virial_correction(3,3)
! REAL(DP), ALLOCATABLE :: U_params_for_types(:)
!
! ! ... setup atoms, nl_handler, charges, U_params_for_types ...
!
! CALL gc_correction%init(cutoff=10.0_DP)
! CALL gc_correction%set_Hubbard_U(atoms, U_params_for_types)
! CALL gc_correction%bind_to(atoms, nl_handler)
!
! ! Calculate potential correction (phi_gc)
! potential_correction = 0.0_DP
! CALL gc_correction%potential(atoms, nl_handler, charges, potential_correction)
!
! ! Calculate energy and force corrections (E_gc, F_gc)
! energy_correction = 0.0_DP
! force_correction = 0.0_DP
! virial_correction = 0.0_DP
! CALL gc_correction%energy_and_forces(atoms, nl_handler, charges, energy_correction, force_correction, virial_correction)
!
! ! Total potential = phi_point_charge + potential_correction
! ! Total energy    = E_point_charge   + energy_correction
! ! Total forces    = F_point_charge   + force_correction
```

## Dependencies and Interactions

*   **`supplib`**: For constants and utility functions.
*   **`particles`, `neighbors`, `filter`**: For core simulation data structures and element filtering.
*   **`ptrdict`**: For parameter registration.
*   This module calculates **corrections**. The main \(1/r\) Coulomb interaction must be handled by another module (e.g., `direct_coulomb`, `pme`).
*   The interpretation of `U` as related to Gaussian width is key. Larger `U` implies a more point-like (narrower) Gaussian, and smaller `U` a more diffuse (wider) one. The self-energy term `0.5*q^2*U` is characteristic of such models.
```
