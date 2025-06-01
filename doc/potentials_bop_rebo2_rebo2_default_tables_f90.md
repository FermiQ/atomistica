# src/potentials/bop/rebo2/rebo2_default_tables.f90

## Overview

The `rebo2_default_tables` module is a crucial component for the second-generation REBO (REBO2) potential in Atomistica. Its primary purpose is to provide the default numerical data for various multi-dimensional lookup tables that form integral parts of the REBO2 potential energy expression. These tables typically store correction factors or complex functional dependencies that are not easily represented by simple analytical functions, such as bond energy corrections based on local coordination (`Fcc`, `Fch`, `Fhh`), pi-conjugation terms (`Pcc`, `Pch`), and dihedral interaction terms (`Tcc`).

The numerical values are generally taken from the original REBO2 publication by Brenner et al., J. Phys.: Condens. Matter 14, 783 (2002). The module contains separate subroutines for populating each specific table. An option to initialize these tables with zeros (using the `#ifdef ZERO_TABLES` preprocessor macro) is available, likely for testing or custom parameterization workflows.

## Key Components

### Modules

*   `rebo2_default_tables`
    *   **Uses**: `supplib` (presumably for utility functions, although not directly evident in the table-filling routines themselves).
    *   **Contains**: A set of subroutines, each dedicated to providing data for a specific REBO2 table.

### Subroutines

*   `rebo2_default_Fcc_table(F, dFdi, dFdj, dFdk)`
    *   **Description**: Populates the 3D array `F` with values for the `Fcc` table, which represents a correction term to C-C bond energies based on the coordination numbers of the two carbon atoms involved in the bond and the conjugation number of the bond. It also populates arrays `dFdi`, `dFdj`, `dFdk` with analytical or finite-difference derivatives of this term.
    *   **Data Source**: Values are primarily from Table IV of the Brenner et al. (2002) paper. Some values are explicitly set, while others appear to be interpolated or adjusted (as indicated by comments like "Refit for proper graphene elastic constants").
    *   **Symmetrization**: The routine explicitly symmetrizes the table values and their derivatives.
    *   **Arguments**:
        *   `F(0:4, 0:4, 0:9) :: real(DP), intent(out)`: The `Fcc` table.
        *   `dFdi(0:4, 0:4, 0:9), dFdj(0:4, 0:4, 0:9), dFdk(0:4, 0:4, 0:9) :: real(DP), intent(out)`: Derivatives of `Fcc`.

*   `rebo2_default_Fch_table(F)`
    *   **Description**: Populates the 3D array `F` with values for the `Fch` table, a correction term for C-H bond energies.
    *   **Data Source**: Values from Table IX of Brenner et al. (2002).
    *   **Symmetrization**: Includes a step to symmetrize table values.
    *   **Arguments**: `F(0:4, 0:4, 0:9) :: real(DP), intent(out)`.

*   `rebo2_default_Fhh_table(F)`
    *   **Description**: Populates the 3D array `F` with values for the `Fhh` table, a correction term for H-H bond energies.
    *   **Data Source**: Values from Table VI of Brenner et al. (2002).
    *   **Arguments**: `F(0:4, 0:4, 0:9) :: real(DP), intent(out)`.

*   `rebo2_default_Pcc_table(P)`
    *   **Description**: Populates the 2D array `P` with values for the `Pcc` table, which is part of the bond-order correction term `Pij` accounting for pi-conjugation effects in C-C bonds.
    *   **Data Source**: Values from Table VIII of Brenner et al. (2002).
    *   **Arguments**: `P(0:5, 0:5) :: real(DP), intent(out)`.

*   `rebo2_default_Pch_table(P)`
    *   **Description**: Populates the 2D array `P` with values for the `Pch` table, part of the `Pij` correction for C-H bonds.
    *   **Data Source**: Values from Table VIII of Brenner et al. (2002).
    *   **Arguments**: `P(0:5, 0:5) :: real(DP), intent(out)`.

*   `rebo2_default_Tcc_table(T)`
    *   **Description**: Populates the 3D array `T` with values for the `Tcc` table, which defines the dihedral (torsional) interaction term for C-C-C-C linkages as a function of conjugation numbers.
    *   **Data Source**: Values from Table V of Brenner et al. (2002).
    *   **Arguments**: `T(0:4, 0:4, 0:9) :: real(DP), intent(out)`.

**Note on `#ifdef ZERO_TABLES`**: If this preprocessor macro is defined at compile time, all the above subroutines will fill their respective output arrays with zeros instead of the default REBO2 values.

## Important Variables/Constants

The primary content of this module is the numerical data itself, hardcoded into the assignments within each subroutine. These represent the default parameterization of the tabulated components of the REBO2 potential.

## Usage Examples

These subroutines are not intended for direct use by end-users. They are called by `rebo2_db_init` (located in `rebo2_db.f90`) during the initialization phase of a REBO2 potential object.

```fortran
! Conceptual call within rebo2_db_init:
! REAL(DP) :: Fcc_data(0:4, 0:4, 0:9)
! REAL(DP) :: Pcc_data(0:5, 0:5)
! ! ... and so on for other tables
!
! CALL rebo2_default_Fcc_table(Fcc_data, Fcc_deriv_i, Fcc_deriv_j, Fcc_deriv_k)
! CALL rebo2_default_Pcc_table(Pcc_data)
! ! ...
! ! Then, this data is used to initialize table_2d_type or table_3d_type objects
! ! within the rebo2_t potential object.
! CALL init(this_rebo2_potential%Fcc, 4, 4, 9, Fcc_data, Fcc_deriv_i, Fcc_deriv_j, Fcc_deriv_k)
! CALL init(this_rebo2_potential%Pcc, 5, 5, Pcc_data)
```

## Dependencies and Interactions

*   **Internal Dependencies:** None beyond standard Fortran.
*   **External Libraries/Modules:** Uses `supplib` (though its direct use is not apparent in these specific routines, it might be used by other parts of the module or for common constants/macros not shown).
*   **Interactions:**
    *   This module serves as the primary source of default numerical data for the tabulated parts of the REBO2 potential.
    *   The subroutine `rebo2_db_init` in `rebo2_db.f90` calls the routines in this module to get the data, which is then used to initialize `table2d_type` and `table3d_type` objects within the main `rebo2_t` potential object.
    *   The data provided here is essential for the correct evaluation of REBO2 energies and forces by `bop_kernel_rebo2.f90`.
```
