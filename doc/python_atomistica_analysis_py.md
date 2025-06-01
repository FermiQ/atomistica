# src/python/atomistica/analysis.py

## Overview

This Python module, `analysis.py`, provides a collection of tools for analyzing atomic configurations within the Atomistica framework. Key functionalities include performing Voronoi tessellation analysis using the external `voro++` program and calculating stress tensor invariants.

## Key Components

### Functions

*   `get_enclosing_orthorhombic_box(cell)`
    *   **Description**: Calculates and returns the lower and upper bounds of an orthorhombic bounding box that perfectly encloses a given simulation cell (parallelepiped).
    *   **Arguments**:
        *   `cell`: A 3x3 matrix or list of lists representing the three cell vectors.
    *   **Returns**: A tuple `(lower, upper)`, where `lower` and `upper` are 3-element NumPy arrays representing the minimum and maximum XYZ coordinates of the bounding box.

*   `voropp_for_non_orthorhombic_cells(_a, q='%v', voropp_path=VOROPP_PATH, fast=False, dump=None)`
    *   **Description**: Performs Voronoi analysis using the `voro++` program, specifically handling non-orthorhombic cells or cells with Lees-Edwards shear conditions. It achieves this by constructing a supercell (2x2x2 if `fast=True`, otherwise 3x3x3) to ensure periodic images are correctly handled by `voro++` when operating within an orthorhombic bounding box of this supercell. Results corresponding to the original central cell are then extracted.
    *   **Arguments**:
        *   `_a`: An ASE (Atomistic Simulation Environment) `Atoms` object.
        *   `q :: str or list`: Specifies the output quantities desired from `voro++` (e.g., '%v' for cell volume, '%n' for face list). See `voro++ -hc` for options. Default: `'%v'`.
        *   `voropp_path :: str`: Path to the `voro++` executable. Default: Value of `VOROPP_PATH` constant.
        *   `fast :: bool`: If true, uses a smaller 2x2x2 supercell for non-orthorhombic cases. Default: `False`.
        *   `dump :: str`: If a filename is provided, dumps the constructed supercell (before Voronoi analysis) to this file using `atomistica.io.write`. Default: `None`.
    *   **Returns**: A NumPy array. If `q` requests a single quantity, it's a 1D array of that quantity for each original atom. If `q` requests multiple quantities, it's a 2D array (num_quantities x num_atoms).

*   `voropp(a, q='%v', voropp_path=VOROPP_PATH, fast=False, dump=None)`
    *   **Description**: The main interface for Voronoi analysis. It checks if the input cell `a.get_cell()` is orthorhombic and has no shear (from `a.info['shear_dx']`). If so, it calls `voro++` directly on the system. Otherwise, it delegates to `voropp_for_non_orthorhombic_cells` to handle the geometric complexities. The function writes atomic positions to a temporary file (`tmp.voronoi`), executes `voro++` via `os.system`, reads the output from `tmp.voronoi.vol`, and then removes these temporary files.
    *   **Arguments**: Same as `voropp_for_non_orthorhombic_cells`.
    *   **Returns**: Same as `voropp_for_non_orthorhombic_cells`.

*   `stress_invariants(s)`
    *   **Description**: Calculates and returns the three invariants of a given stress tensor(s). It can process a single stress tensor (Voigt 6-component vector or 3x3 matrix) or an array of stress tensors.
    *   **Arguments**:
        *   `s`: A NumPy array representing the stress tensor(s). Can be shape (6,), (3,3), (N,6), or (N,3,3).
    *   **Returns**: A tuple `(hydrostatic_pressure, octahedral_shear_stress, J3)`:
        *   `hydrostatic_pressure = -I1/3`
        *   `octahedral_shear_stress = sqrt(2*J2/3)`
        *   `J3` (the third deviatoric stress invariant)
        where I1, I2, I3 are the principal invariants of the stress tensor, and J2, J3 are invariants of the deviatoric stress tensor.

### Classes/Types
This module does not define any classes.

## Important Variables/Constants

*   `VOROPP_PATH :: str`: Module-level constant defining the default command/path for the `voro++` executable. Value: `'voro++'`.

## Usage Examples

```python
import numpy as np
from ase import Atoms
from atomistica.analysis import voropp, stress_invariants

# Example for voropp (requires voro++ installed)
# Create a simple ASE Atoms object
atoms = Atoms('Ar2', positions=[[0,0,0], [3,0,0]], cell=[10,10,10], pbc=True)
try:
    volumes = voropp(atoms, q='%v') # Get Voronoi volumes
    print("Voronoi volumes:", volumes)
    # To get volumes and number of faces:
    # data = voropp(atoms, q=['v', 'n'])
    # print("Volumes:", data[0])
    # print("Face orders:", data[1]) # This would need further parsing if 'n' gives list of face orders
except OSError:
    print("voro++ command not found or failed to run.")


# Example for stress_invariants
# Voigt notation: [s_xx, s_yy, s_zz, s_xy, s_yz, s_zx]
stress_voigt = np.array([10, 12, 11, 1, 0.5, 0.2])
p, tau_oct, j3 = stress_invariants(stress_voigt)
print(f"Hydrostatic Pressure: {p}")
print(f"Octahedral Shear Stress: {tau_oct}")
print(f"J3: {j3}")

# Matrix notation
stress_matrix = np.array([[10, 1, 0.2],
                          [1, 12, 0.5],
                          [0.2, 0.5, 11]])
p_m, tau_oct_m, j3_m = stress_invariants(stress_matrix)
print(f"Hydrostatic Pressure (from matrix): {p_m}")
```

## Dependencies and Interactions

*   **External Program**:
    *   `voro++`: The `voropp` and `voropp_for_non_orthorhombic_cells` functions require the `voro++` command-line tool to be installed and accessible in the system's PATH or at the path specified by `VOROPP_PATH`.
*   **Python Libraries**:
    *   `os`: Used for executing the `voro++` command via `os.system` and for file operations (creating `tmp.voronoi`, removing temporary files).
    *   `numpy` (as `np`): Extensively used for array manipulations, mathematical operations (e.g., `np.min`, `np.max`, `np.linalg.norm`), and loading data from text files (`np.loadtxt`).
*   **Internal Atomistica Dependencies**:
    *   `atomistica.io.write`: Used by `voropp_for_non_orthorhombic_cells` if the `dump` argument is provided, for writing ASE `Atoms` objects to a file.
*   The `voropp` functions interact with the file system by creating and deleting temporary files (`tmp.voronoi`, `tmp.voronoi.vol`) in the current working directory during their execution.
```
