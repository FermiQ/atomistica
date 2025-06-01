# src/python/atomistica/snippets.py

## Overview

This Python module, `snippets.py`, is intended to house a collection of small, reusable utility functions ("code snippets") that simplify common tasks within the Atomistica framework or related atomistic simulations. Currently, it provides one key function: `mic` for applying the Minimum Image Convention.

## Key Components

### Functions

*   `mic(dr, cell, pbc=None)`
    *   **Description**: Applies the Minimum Image Convention (MIC) to an array of distance vectors, considering periodic boundary conditions. Given a raw distance vector (or an array of such vectors) \(\vec{dr}\) and the simulation cell vectors, this function returns the shortest vector \(\vec{dr}_{\text{MIC}}\) that connects two points (or pairs of points) in a periodic system.
    *   **Algorithm**:
        1.  Calculates the inverse of the simulation cell matrix: `rec = cell_inverse`.
        2.  If `pbc` flags (a boolean array indicating periodicity along each cell vector) are provided, it effectively nullifies non-periodic components of `rec` to ensure wrapping only occurs along periodic directions.
        3.  Projects the input distance vector(s) `dr` into the basis of the cell vectors (scaled coordinates): `dr_scaled = dr @ rec`.
        4.  Rounds the scaled coordinates to the nearest integer: `dri = round(dr_scaled)`. These integers represent the number of cell image shifts required.
        5.  Calculates the MIC vector: `dr_mic = dr - dri @ cell`.
    *   **Arguments**:
        *   `dr :: numpy.ndarray`: A distance vector (shape (3,)) or an array of distance vectors (shape (N, 3)).
        *   `cell :: numpy.ndarray`: A 3x3 matrix where rows (or columns, depending on convention, but ASE typically uses rows as cell vectors `cell[i,:]`) are the simulation cell vectors. The function uses `np.linalg.inv(cell)` suggesting `cell` should be invertible and its columns are basis vectors if `dr` is a row vector, or rows are basis vectors if `dr` is a column vector or `dot(dr,rec)` means `dr @ rec.T`. Given `np.dot(dri, cell)` at the end, `cell` rows are likely the cell vectors.
        *   `pbc :: array_like (bool), optional`: A 3-element boolean array indicating periodicity along each cell vector direction (e.g., `[True, True, False]`). If `None`, assumes periodicity in all directions where cell vectors are non-zero (though current implementation applies `rec` scaling regardless of `pbc` before rounding, then `pbc` only affects `rec` for the rounding step if `rec` was modified by `pbc` - the provided code multiplies `rec` by `pbc` shaped as `(3,1)` which means `rec` columns are scaled by `pbc` values. This effectively makes `dr_scaled` components for non-periodic directions zero if `pbc` for that direction is `False` *before* rounding, leading to `dri` being zero for that component).
    *   **Returns**: `numpy.ndarray`: The distance vector(s) according to the Minimum Image Convention, with the same shape as the input `dr`.

### Classes/Types
This module does not define any classes.

## Important Variables/Constants
This module does not define any public module-level constants.

## Usage Examples

```python
import numpy as np
from atomistica.snippets import mic

# Example 1: Orthorhombic cell
cell_ortho = np.array([[10.0, 0.0, 0.0],
                       [0.0, 12.0, 0.0],
                       [0.0, 0.0, 15.0]])
pbc_true = [True, True, True]

# Distance vector pointing outside the primary cell
dr_original = np.array([11.0, -7.0, 16.0])
dr_mic_ortho = mic(dr_original, cell_ortho, pbc=pbc_true)
# Expected: dr_mic_ortho = [1.0, 5.0, 1.0] (or [-1.0, 5.0, 1.0] depending on rounding for exact 0.5)
# Correction: The rounding logic np.round(np.dot(dr, rec)) and then dr - np.dot(dri, cell)
# For dr_original[0] = 11.0, scaled = 1.1, round(1.1)=1.0. dr_mic[0] = 11.0 - 1.0*10.0 = 1.0.
# For dr_original[1] = -7.0, scaled = -7/12 = -0.5833, round(-0.5833)=-1.0. dr_mic[1] = -7.0 - (-1.0*12.0) = 5.0.
# For dr_original[2] = 16.0, scaled = 16/15 = 1.066, round(1.066)=1.0. dr_mic[2] = 16.0 - 1.0*15.0 = 1.0.
print(f"Original: {dr_original}, MIC (ortho): {dr_mic_ortho}")


# Example 2: Array of distance vectors
dr_array = np.array([[11.0, -7.0, 16.0],
                     [-1.0, 13.0, -0.5]])
dr_mic_array = mic(dr_array, cell_ortho, pbc=pbc_true)
print(f"Original array:\n{dr_array}\nMIC array (ortho):\n{dr_mic_array}")


# Example 3: Triclinic cell
cell_triclinic = np.array([[10.0, 0.0, 0.0],
                           [2.0, 12.0, 0.0],
                           [1.0, 1.5, 15.0]])
dr_triclinic = np.array([11.0, 1.0, 1.0]) # A point just outside along first cell vector
# scaled = dr_triclinic @ inv(cell_triclinic)
# scaled approx = [1.1, -0.1833, -0.0055] -> round -> [1., 0., 0.]
# dr_mic = dr_triclinic - [1.,0.,0.] @ cell_triclinic = [11,1,1] - [10,0,0] = [1,1,1]
dr_mic_triclinic = mic(dr_triclinic, cell_triclinic, pbc=pbc_true)
print(f"Original: {dr_triclinic}, MIC (triclinic): {dr_mic_triclinic}")

# Example 4: Non-periodic dimension
pbc_mixed = [True, True, False]
dr_non_periodic = np.array([11.0, 7.0, 16.0]) # z-component should not wrap
# rec will have its 3rd column (or row, depending on inv) zeroed out by pbc scaling.
# So dri for z will be 0.
# dr_mic_mixed should be [1.0, -5.0, 16.0] (assuming 7.0 -> -5.0 for y)
# For y: 7.0 / 12.0 = 0.5833 -> round(0.5833) = 1.0.  7.0 - 1.0*12.0 = -5.0
dr_mic_mixed = mic(dr_non_periodic, cell_ortho, pbc=pbc_mixed)
print(f"Original: {dr_non_periodic}, MIC (mixed pbc): {dr_mic_mixed}")

```

## Dependencies and Interactions

*   **`numpy`**: This module heavily relies on NumPy for all numerical computations, including:
    *   `numpy.linalg.inv` for calculating the inverse of the cell matrix.
    *   `numpy.dot` for matrix-vector and matrix-matrix products.
    *   `numpy.round` for rounding scaled coordinates to the nearest integer.
    *   Basic array creation and manipulation.
*   The function is essential for correctly calculating distances and displacement vectors in simulations employing periodic boundary conditions, which is a common requirement in many analysis tasks (e.g., radial distribution function, mean squared displacement, strain calculation).
```
