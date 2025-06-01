# src/potentials/bop/brenner/brenner_params.f90

## Overview

This file defines the parameter sets for various Brenner-type Bond-Order Potentials (BOPs) used within the Atomistica simulation package. It introduces the `BOP_DB_TYPE` (aliased as `brenner_db_t` within the Brenner module context) derived data type, which encapsulates all necessary parameters for a specific Brenner potential parameterization. The file then declares several constant instances of `BOP_DB_TYPE`, each corresponding to a published parameter set for different materials (e.g., SiC, PtC, FeC, AlN). These instances are collected into a global private array `BOP_DB`, allowing the simulation to select and load these predefined parameter sets by a reference string.

## Key Components

### Constants

*   `BRENNER_MAX_REF :: integer, parameter`: Defines the maximum length of the reference string (e.g., a citation) for a parameter set. Value: 1000.
*   `BRENNER_MAX_EL :: integer, parameter`: Maximum number of distinct chemical elements that a single parameter set can describe (e.g., for a binary or ternary system). Value: 3.
*   `BRENNER_MAX_PAIRS :: integer, parameter`: Maximum number of unique chemical element pairs that can be formed from `BRENNER_MAX_EL` elements. Calculated using the `PAIR_INDEX` macro. This determines the size of parameter arrays like `D0`, `r0`, etc.

### Data Types

*   `BOP_DB_TYPE` (Public, typically aliased to `brenner_db_t`)
    *   **Description**: A Fortran derived data type designed to store a complete set of parameters defining a specific Brenner-like interatomic potential.
    *   **Fields**:
        *   `nel :: integer`: The number of chemical elements defined in this parameter set (e.g., 2 for SiC).
        *   `nD0, nr0, nS, ... :: integer`: Integer counts for each parameter array, used to store the actual number of defined values for validation (e.g., `nD0` would be 3 for a 2-element system: A-A, A-B, B-B).
        *   `ref :: character(BRENNER_MAX_REF)`: A character string containing the bibliographic reference for the parameter set.
        *   `el(2, BRENNER_MAX_EL) :: character`: An array storing the chemical symbols of the elements (e.g., `el(:,1) = "C ", el(:,2) = "Si"`).
        *   **Pair Interaction Parameters** (arrays dimensioned to `BRENNER_MAX_PAIRS`):
            *   `D0 :: real(DP)`: Equilibrium binding energy of the dimer.
            *   `r0 :: real(DP)`: Equilibrium bond distance of the dimer.
            *   `S :: real(DP)`: Parameter related to the slope of the Pauling plot, affecting the range of interaction.
            *   `beta :: real(DP)`: Parameter related to the stiffness of the dimer bond (vibrational frequency).
        *   **Bond Order Parameters**:
            *   `gamma :: real(DP)`: Scaling factor for the bond-order contribution.
            *   `c, d, h :: real(DP)`: Parameters defining the angular dependence (g_ijk term) of the bond order.
            *   `mu :: real(DP)`: Parameter for the exponential length-dependent contribution to the bond order term.
            *   `n :: real(DP)`: Exponent in the bond order function, typically `(1 + zij^n)`.
            *   `m :: integer`: Exponent for the distance-dependent part of the bond order, `exp((2*mu*dr)**m)`.
        *   **Cutoff Parameters**:
            *   `r1, r2 :: real(DP)`: Start and end distances for the primary interaction cutoff (or inner cutoff if screening is active).
        *   **Screening Parameters** (compiled if `#ifdef SCREENING` is active):
            *   `or1, or2 :: real(DP)`: Start and end distances for the outer cutoff function (used in screening).
            *   `bor1, bor2 :: real(DP)`: Start and end distances for the bond-order specific cutoff function (used in screening).
            *   `Cmin, Cmax :: real(DP)`: Parameters defining the range of the screening function `S(C)`.

### Predefined Parameter Sets

The file defines the following `BOP_DB_TYPE` instances as `parameter` constants, with values taken from literature:

*   `Erhart_PRB_71_035211_SiC`: Parameters for Silicon Carbide (SiC) from Erhart, P., & Albe, K. (2005). Phys. Rev. B, 71(3), 035211.
*   `Albe_PRB_65_195124_PtC`: Parameters for Platinum Carbide (PtC) from Albe, K., Nordlund, K., & Averback, R. S. (2002). Phys. Rev. B, 65(19), 195124.
*   `Henriksson_PRB_79_144107_FeC`: Parameters for Iron Carbide (FeC) from Henriksson, K. O. E., & Nordlund, K. (2009). Phys. Rev. B, 79(14), 144107.
*   `Kioseoglou_PSSb_245_1118_AlN`: Parameters for Aluminium Nitride (AlN) from Kioseoglou, J., Komninou, P., & Karakostas, T. (2008). Phys. Stat. Sol. (b), 245(6), 1118–1126.

These parameter sets are initialized with specific values for all fields of the `BOP_DB_TYPE` structure, including element types, reference strings, and all numerical parameters for pair interactions, bond order terms, and cutoffs. The initialization uses `FILL1`, `FILL3`, `FILL3i` macros for padding unused array elements when `BRENNER_MAX_EL` is larger than the number of elements in a specific parameter set.

### Global Parameter Database

*   `BOP_DB(4) :: type(BOP_DB_TYPE), parameter, private`: A private array that aggregates all the predefined parameter sets (`Erhart_PRB_71_035211_SiC`, etc.). This array serves as an internal database from which the `INIT_FUNC` (in `brenner_module.f90`) can retrieve a specific parameter set based on a reference string.

## Important Variables/Constants

This file is primarily dedicated to the definition of the `BOP_DB_TYPE` and the constant parameter data. The `FILL*` macros are helper macros for array initialization.

## Usage Examples

This file itself does not contain executable code but provides the data that is consumed by other parts of the Brenner potential module. The `INIT_FUNC` in `brenner_module.f90` uses the `BOP_DB` array to look up and load the parameters into a `brenner_t` object.

```fortran
! Conceptual lookup in INIT_FUNC (brenner_module.f90):
! IF (this%ref_string == Erhart_PRB_71_035211_SiC%ref) THEN
!   this%db_parameters = Erhart_PRB_71_035211_SiC
! END IF
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   Uses the `PAIR_INDEX` macro (likely defined in `macros.inc` or a similar core include file) to calculate array sizes.
*   **External Libraries:** None.
*   **Interactions:**
    *   This file is central to defining the actual behavior of the Brenner potential, as it holds the numerical values for all its terms.
    *   The `brenner_module.f90::INIT_FUNC` reads from the `BOP_DB` array or uses a directly passed `BOP_DB_TYPE` structure to configure a `brenner_t` potential instance.
    *   The structure of `BOP_DB_TYPE` must be consistent with how `brenner_func.f90` and `brenner_kernel.f90` access and use these parameters.
    *   The `#ifdef SCREENING` directive conditionally includes parameters relevant to screening functionalities, ensuring the data structure adapts to compiled features.
```
