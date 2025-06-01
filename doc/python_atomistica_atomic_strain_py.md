# src/python/atomistica/atomic_strain.py

## Overview

This Python module, `atomic_strain.py`, provides tools to compute local atomic strain measures, specifically the per-atom deformation gradient tensor and the \(D^2_{\text{min}}\) measure of non-affine displacements. The \(D^2_{\text{min}}\) measure quantifies the deviation of local atomic displacements from a best-fit affine transformation, providing insight into plastic deformation, defects, and other non-linear material responses. The methodology is based on the work by Falk and Langer, Phys. Rev. B 57, 7192 (1998).

## Key Components

### Functions

*   `get_XIJ(nat, i_now, dr_now, dr_old)`
    *   **Description**: Calculates the matrix \( \mathbf{X}_k \) for each atom \(k\). This matrix is defined as \( \mathbf{X}_k = \sum_{j \in \text{neigh}(k)} \vec{dr}_{kj}^{\text{now}} \otimes \vec{dr}_{kj}^{\text{old}} \), where \(\vec{dr}^{\text{now}}\) are current bond vectors and \(\vec{dr}^{\text{old}}\) are reference bond vectors to neighbors.
    *   **Arguments**:
        *   `nat :: int`: Total number of atoms.
        *   `i_now :: array_like`: 1D array of central atom indices for each bond.
        *   `dr_now :: array_like`: Array of current bond vectors (N_bonds x 3).
        *   `dr_old :: array_like`: Array of reference bond vectors (N_bonds x 3).
    *   **Returns**: `numpy.ndarray`: Array of shape (nat, 3, 3) containing the \( \mathbf{X}_k \) matrix for each atom.

*   `get_YIJ(nat, i_now, dr_old)`
    *   **Description**: Calculates the matrix \( \mathbf{Y}_k \) for each atom \(k\). This matrix is defined as \( \mathbf{Y}_k = \sum_{j \in \text{neigh}(k)} \vec{dr}_{kj}^{\text{old}} \otimes \vec{dr}_{kj}^{\text{old}} \).
    *   **Arguments**:
        *   `nat :: int`: Total number of atoms.
        *   `i_now :: array_like`: 1D array of central atom indices for each bond.
        *   `dr_old :: array_like`: Array of reference bond vectors (N_bonds x 3).
    *   **Returns**: `numpy.ndarray`: Array of shape (nat, 3, 3) containing the \( \mathbf{Y}_k \) matrix for each atom.

*   `array_inverse(A)`
    *   **Description**: Computes the inverse for each 3x3 matrix in a list (or an array of shape (N, 3, 3)). It uses `numpy.linalg.lapack_lite.dgesv` for efficient inversion of multiple small matrices.
    *   **Arguments**:
        *   `A :: numpy.ndarray`: Array of shape (N, 3, 3) containing N matrices to be inverted.
    *   **Returns**: `numpy.ndarray`: Array of shape (N, 3, 3) containing the inverted matrices. Raises `numpy.linalg.LinAlgError` if any matrix is singular.

*   `get_delta_plus_epsilon(nat, i_now, dr_now, dr_old)`
    *   **Description**: Calculates the local deformation gradient tensor \( \mathbf{F}_k = \mathbf{I} + \mathbf{\epsilon}_k \) for each atom \(k\). This tensor describes the best-fit affine transformation that maps the reference neighborhood of atom \(k\) to its current neighborhood. It is calculated as \( \mathbf{F}_k = \mathbf{X}_k (\mathbf{Y}_k)^{-1} \).
    *   **Arguments**: Same as `get_XIJ`.
    *   **Returns**: `numpy.ndarray`: Array of shape (nat, 3, 3) containing the deformation gradient \( \mathbf{F}_k \) for each atom.

*   `get_D_square_min(atoms_now, atoms_old, i_now, j_now, delta_plus_epsilon=None)`
    *   **Description**: Calculates the \(D^2_{\text{min}}\) norm as defined by Falk and Langer for each atom. This scalar value quantifies the degree of non-affine displacement in the neighborhood of an atom.
        1.  It computes current bond vectors \(\vec{dr}^{\text{now}}\) and reference bond vectors \(\vec{dr}^{\text{old}}\) using the Minimum Image Convention (`mic`) for periodic systems.
        2.  If `delta_plus_epsilon` (the per-atom deformation gradient \( \mathbf{F}_k \)) is not provided, it calculates it using `get_delta_plus_epsilon`.
        3.  For each atom \(k\), it computes \(D^2_k = \sum_{j \in \text{neigh}(k)} |\vec{dr}_{kj}^{\text{now}} - \mathbf{F}_k \vec{dr}_{kj}^{\text{old}}|^2 \).
    *   **Arguments**:
        *   `atoms_now`: ASE `Atoms` object for the current configuration.
        *   `atoms_old`: ASE `Atoms` object for the reference configuration.
        *   `i_now :: array_like`: 1D array of central atom indices for each bond in the neighbor list.
        *   `j_now :: array_like`: 1D array of neighbor atom indices for each bond.
        *   `delta_plus_epsilon :: numpy.ndarray, optional`: Pre-calculated per-atom deformation gradients.
    *   **Returns**: A tuple `(delta_plus_epsilon, d_sq)`, where `delta_plus_epsilon` is the (N_atoms, 3, 3) array of deformation gradients and `d_sq` is a 1D array (N_atoms) of \(D^2_{\text{min}}\) values.

*   `atomic_strain(atoms_now, atoms_old, cutoff=None, i_now=None, j_now=None)`
    *   **Description**: This is the main user-facing function to calculate the per-atom deformation gradient tensor (\(\mathbf{F}_k\)) and the \(D^2_{\text{min}}\) measure of non-affine displacements.
    *   **Arguments**:
        *   `atoms_now`: ASE `Atoms` object for the current (deformed) configuration.
        *   `atoms_old`: ASE `Atoms` object for the reference (undeformed) configuration.
        *   `cutoff :: float, optional`: If neighbor lists `i_now, j_now` are not provided, this cutoff distance is used to compute them using Atomistica's native neighbor list routines.
        *   `i_now, j_now :: array_like, optional`: Pre-computed neighbor lists (indices of central and neighbor atoms for each bond). If `None`, they are computed using `cutoff`.
    *   **Returns**: A tuple `(delta_plus_epsilon, d_sq)`:
        *   `delta_plus_epsilon`: NumPy array of shape (N_atoms, 3, 3) containing the deformation gradient \( \mathbf{F}_k = \mathbf{I} + \mathbf{\epsilon}_k \) for each atom.
        *   `d_sq`: NumPy array of shape (N_atoms,) containing the \(D^2_{\text{min}}\) value for each atom.
    *   **Raises**: `ValueError` if neither `cutoff` nor neighbor lists (`i_now`, `j_now`) are provided.

### Classes/Types
This module does not define any classes.

## Important Variables/Constants
This module does not define any public module-level constants.

## Usage Examples

```python
import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from atomistica.analysis import atomic_strain
from atomistica.deformation import F_matrix_to_positions # For applying deformation
from atomistica.native import Neighbors # For neighbor list if needed

# 1. Create a reference configuration
atoms_old = bulk('Cu', 'fcc', a=3.6, cubic=True)
atoms_old.set_calculator(EMT()) # Calculator needed for some operations if not using native NL

# 2. Create a deformed configuration (example: simple shear)
F_def = np.eye(3)
F_def[0,1] = 0.05 # Apply a small shear
atoms_now = atoms_old.copy()
F_matrix_to_positions(atoms_now, F_def, scale_positions=True, scale_cell=True)
atoms_now.set_calculator(EMT())


# 3. Calculate atomic strain and D^2_min
# Option A: Provide a cutoff for neighbor list calculation
try:
    # Note: A realistic cutoff would depend on the material.
    # For atomistica.native.Neighbors, a calculator might not be strictly needed
    # on atoms_now if only positions and cell are used.
    deformation_gradient, d_square_min = atomic_strain(atoms_now, atoms_old, cutoff=3.0)

    print("Deformation gradient for first atom:\n", deformation_gradient[0])
    print("D^2_min for first atom:", d_square_min[0])

except ImportError:
    print("Atomistica native modules not available for neighbor list calculation.")
except Exception as e:
    print(f"An error occurred during atomic_strain: {e}")

# Option B: Pre-calculate neighbor list (e.g., using ASE or Atomistica's tools)
# This example assumes you have a way to get i_now, j_now.
# For instance, using atomistica.native.Neighbors:
# p_now = native.from_atoms(atoms_now)
# nl = native.Neighbors(avgn=12) # Ensure avgn is reasonable
# nl.request_interaction_range(3.0)
# i_list, j_list, _ = nl.get_neighbors(p_now)
# deformation_gradient, d_square_min = atomic_strain(atoms_now, atoms_old, i_now=i_list, j_now=j_list)

```

## Dependencies and Interactions

*   **`numpy`**: Extensively used for numerical computations, array storage, and linear algebra operations (including LAPACK interface for matrix inversion).
*   **`atomistica.native`**: The `atomic_strain` function can use `atomistica.native.from_atoms` and `atomistica.native.Neighbors` to compute neighbor lists if they are not explicitly provided. This implies a dependency on the compiled core of Atomistica.
*   **`atomistica.snippets.mic`**: Used by `get_D_square_min` to correctly calculate interatomic displacement vectors respecting periodic boundary conditions (Minimum Image Convention).
*   **ASE (Atomistic Simulation Environment)**: The primary input data structures are ASE `Atoms` objects.

The module provides a robust way to quantify local, non-affine deformations in atomic systems, which is essential for analyzing plastic events, defect structures, and mechanical responses beyond the linear elastic regime.
```
