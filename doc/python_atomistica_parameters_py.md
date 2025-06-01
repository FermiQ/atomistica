# src/python/atomistica/parameters.py

## Overview

This Python module, `parameters.py`, serves as a database of pre-defined parameter sets for various interatomic potentials commonly used within the Atomistica framework. These potentials include several variants of Tersoff, Brenner-type (Albe, Erhart, Henriksson, Kioseoglou), Juslin, and Kumagai potentials.

The parameters are stored as Python dictionaries, with keys corresponding to the parameter names used by the respective potential models (e.g., `A`, `B`, `lambda`, `mu` for Tersoff). For multi-element systems, pair and triplet interaction parameters are typically stored in flat lists, and this module provides helper functions (`pair_index`, `triplet_index`) to map element indices to the correct list index, consistent with Atomistica's Fortran backend. It also includes utilities for applying standard mixing rules (arithmetic and geometric) to generate cross-interaction parameters.

## Key Components

### Indexing Functions

These functions are crucial for correctly ordering parameters in flat lists for multi-element systems, matching the indexing expected by the Fortran backend of Atomistica.

*   `pair_index(i, j, maxval)`
    *   **Description**: Calculates a unique flat list index for an atom pair `(i, j)`. `i` and `j` are 0-based indices of the element types, and `maxval` is the total number of element types in the parameter set. The indexing ensures that `(i,j)` and `(j,i)` map to the same location if the potential's parameters are symmetric.
    *   Corresponds to the `PAIR_INDEX` macro in Atomistica's Fortran code.
*   `triplet_index(i, j, k, maxval)`
    *   **Description**: Calculates a unique flat list index for an atom triplet `(i, j, k)`. `i`, `j`, `k` are 0-based element type indices, and `maxval` is the total number of element types.
    *   Corresponds to the `TRIPLET_INDEX_NS` (non-symmetric) macro in Atomistica's Fortran code, implying the order of i, j, k matters for parameters like Juslin's `alpha` and `omega`.

### Mixing Rule Functions

These functions are used to generate parameters for interactions between different element types when only same-element parameters are explicitly provided, or to fill in missing cross-terms in a multi-element parameter set.

*   `mix(p, key, rule)`
    *   **Description**: A generic function that applies a given binary `rule` (a lambda function) to generate parameters for mixed pairs `(i,j)` from homo-elemental parameters `(i,i)` and `(j,j)`.
    *   **Arguments**:
        *   `p :: dict`: The parameter dictionary.
        *   `key :: str`: The parameter name (e.g., "A", "sigma") whose mixed terms are to be calculated.
        *   `rule :: function`: A lambda function `lambda x, y: ...` that defines how to combine parameters `x` (for type `i`) and `y` (for type `j`).
*   `mix_arithmetic(p, key)`
    *   **Description**: Applies the arithmetic mean mixing rule: \(P_{ij} = (P_{ii} + P_{jj}) / 2\).
*   `mix_geometric(p, key)`
    *   **Description**: Applies the geometric mean mixing rule: \(P_{ij} = \sqrt{P_{ii} \cdot P_{jj}}\).

### Parameter Dictionaries

The module defines numerous Python dictionaries, each representing a specific published or commonly used parameter set. Examples include:

*   **Tersoff Potentials**:
    *   `Tersoff_PRB_39_5566_Si_C`: For Si-C systems (Tersoff 1989).
    *   `Tersoff_PRB_39_5566_Si_C__Scr`: A screened version.
    *   `Goumri_Said_ChemPhys_302_135_Al_N`: For Al-N systems.
    *   `Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N`: For B-C-N systems, demonstrating the use of mixing rules to fill out cross-terms.
    *   `Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr`: Screened version.

*   **Brenner-type Potentials (for use with `Brenner` ASE calculator)**:
    *   `Erhart_PRB_71_035211_SiC`: Erhart & Albe's parameters for SiC.
    *   `Erhart_PRB_71_035211_SiC__Scr`: Screened version.
    *   `Albe_PRB_65_195124_PtC`: Albe et al.'s parameters for PtC.
    *   `Henriksson_PRB_79_114107_FeC`: Henriksson & Nordlund's parameters for FeC.
    *   `Kioseoglou_PSSb_245_1118_AlN`: Kioseoglou et al.'s parameters for AlN.
    *   `Brenner_PRB_42_9458_C_I` & `Brenner_PRB_42_9458_C_II`: Brenner's original 1990 potentials for Carbon.

*   **Juslin Potentials (for use with `Juslin` ASE calculator)**:
    *   `Juslin_JAP_98_123520_WCH`: Juslin et al.'s W-C-H parameters.
    *   `Juslin_JAP_98_123520_WCH__Scr`: Screened version.
    *   `Kuopanportti_CMS_111_525_FeCH`: Kuopanportti et al.'s Fe-C-H (Juslin-like form).
    *   `Kuopanportti_CMS_111_525_FeCH__Scr`: Screened version.

*   **Kumagai Potential (for use with `Kumagai` ASE calculator)**:
    *   `Kumagai_CompMaterSci_39_457_Si`: For Silicon.
    *   `Kumagai_CompMaterSci_39_457_Si__Scr`: Screened version.

Each parameter dictionary typically contains:
*   `"__ref__" :: str`: A bibliographic reference for the parameter set.
*   `"el" :: list[str]`: A list of element symbols covered by the parameter set (e.g., `["Si", "C"]`).
*   Keys for each potential-specific parameter (e.g., `A`, `B`, `lambda`, `mu`, `r1`, `r2` for Tersoff; `D0`, `r0`, `S`, `beta` for Brenner-type). Parameter values are often lists, ordered according to `pair_index` or `triplet_index`.
*   Screened versions (often denoted with `__Scr`) include additional or modified parameters relevant to screening functions (e.g., `or1`, `or2`, `Cmin`, `Cmax`).

## Important Variables/Constants
The dictionaries containing the parameter sets are the primary "constants" provided by this module.

## Usage Examples

These parameter dictionaries are designed to be used with the `convpar` function (from `atomistica.aseinterface`) and then passed to the constructor of the corresponding Atomistica ASE calculator class.

```python
from atomistica.parameters import Tersoff_PRB_39_5566_Si_C
from atomistica.aseinterface import convpar, Tersoff # Assuming Tersoff is dynamically created

# Get the parameter dictionary
tersoff_params_SiC = Tersoff_PRB_39_5566_Si_C

# Convert for use with native code (if needed by the calculator's __init__)
# The Atomistica ASE calculator usually handles this internally via convpar.
# native_params = convpar(tersoff_params_SiC) # Not always directly needed by user

# Initialize the ASE calculator with the parameters
# The calculator's __init__ typically expects the Python dictionary directly.
try:
    calculator = Tersoff(**tersoff_params_SiC)
    # For calculators that take a 'param_file' or 'ref_string' argument,
    # one might use the "__ref__" or a name derived from the dictionary key.
    # Alternatively, the calculator might look up these dictionaries internally.
except Exception as e:
    print(f"Error initializing calculator: {e}")
    print("Note: Specific calculator constructors might vary.")

```

## Dependencies and Interactions

*   **`copy.deepcopy`**: Used for creating modified copies of parameter sets (e.g., for screened versions).
*   **`math.log`, `math.sqrt`**: Used for calculating some parameter values directly within this file (e.g., mixed terms, or derived `mubo` in `Tersoff_PRB_39_5566_Si_C__Scr`).
*   **`atomistica.aseinterface.convpar`**: This function (defined elsewhere) is the intended consumer of these parameter dictionaries to format them for the native Atomistica backend.
*   The structure and parameter names within these dictionaries must align with the expectations of the specific Fortran potential kernels and the `convpar` function.
```
