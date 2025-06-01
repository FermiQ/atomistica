# src/python/atomistica/native.py

## Overview

This Python module, `native.py`, serves as a direct interface to some of Atomistica's core data structures implemented in the native C/Fortran extension module (`_atomistica`). It provides utility functions to:
1.  Convert ASE (Atomistic Simulation Environment) `Atoms` objects into Atomistica's native `Particles` objects.
2.  Facilitate the creation and population of Atomistica's native `Neighbors` list objects.

These functions are primarily used by other Python components within Atomistica (such as the ASE calculator interface or analysis scripts) when they need to interact with the underlying compiled Atomistica routines that operate on these native types.

## Key Components

### Functions

*   `from_atoms(atoms)`
    *   **Description**: Converts an ASE `Atoms` object into an instance of Atomistica's native `_atomistica.Particles` class. This involves transferring information such as the number of atoms, cell dimensions, periodic boundary conditions, atomic numbers (mapping symbols to Z), and atomic positions.
    *   **Steps**:
        1.  Creates an `_atomistica.Particles()` object.
        2.  Allocates memory within the `Particles` object for the number of atoms in the input `atoms` object.
        3.  Sets the simulation cell and periodic boundary conditions using `particles.set_cell()`.
        4.  Populates the atomic numbers (`particles.Z`) by converting element symbols from the ASE `Atoms` object to integers using `ase.data.atomic_numbers`.
        5.  Copies the atomic positions from `atoms.get_positions()` to `particles.coordinates`.
        6.  Calls `particles.I_changed_positions()` to signal that positions have been set and internal data structures (like scaled coordinates) might need an update.
        7.  Calls `particles.update_elements()` to finalize element-specific setup within the `Particles` object.
    *   **Arguments**:
        *   `atoms :: ase.Atoms`: The input ASE `Atoms` object.
    *   **Returns**: `_atomistica.Particles`: An instance of the native Atomistica `Particles` class, populated with data from the ASE `Atoms` object.

*   `neighbor_list(particles, cutoff, avgn=100)`
    *   **Description**: Creates, builds, and returns an Atomistica native neighbor list object (`_atomistica.Neighbors`).
    *   **Steps**:
        1.  Instantiates an `_atomistica.Neighbors(avgn)` object. The `avgn` parameter typically provides an estimate for the average number of neighbors per atom, which can help in pre-allocating memory.
        2.  Calls `neighbors.request_interaction_range(cutoff)` to set the cutoff distance for pair searching.
        3.  Calls `neighbors.update(particles)` to build the neighbor list based on the atomic positions in the provided `particles` object and the set cutoff.
    *   **Arguments**:
        *   `particles :: _atomistica.Particles`: An instance of Atomistica's native `Particles` class.
        *   `cutoff :: float`: The cutoff distance for including pairs in the neighbor list.
        *   `avgn :: int, optional`: An estimate of the average number of neighbors per atom. Default: `100`.
    *   **Returns**: `_atomistica.Neighbors`: An instance of the native Atomistica `Neighbors` class, populated with the neighbor list for the given `particles` and `cutoff`.

### Classes/Types
This module does not define any new classes; it uses classes imported from `_atomistica`.

## Important Variables/Constants
This module does not define any public module-level constants.

## Usage Examples

These functions are typically used as internal utilities within Atomistica's Python scripts when interfacing with the compiled core.

```python
import numpy as np
from ase.build import bulk
from atomistica.native import from_atoms, neighbor_list # Assuming _atomistica is compiled and accessible

# 1. Create an ASE Atoms object
ase_cu = bulk('Cu', 'fcc', a=3.6)

# 2. Convert ASE Atoms to Atomistica Particles
native_particles = from_atoms(ase_cu)
print(f"Created Atomistica Particles object with {len(native_particles.Z)} atoms.")

# 3. Create a neighbor list for the Atomistica Particles
#    (Example cutoff, adjust based on material and needs)
cutoff_distance = 3.0
try:
    nl = neighbor_list(native_particles, cutoff_distance)
    # The 'nl' object can now be passed to other native Atomistica routines
    # that require a pre-built neighbor list.
    # Example: i_indices, j_indices, _ = nl.get_neighbors(native_particles)
    print(f"Neighbor list created with cutoff {cutoff_distance} Å.")
except Exception as e:
    # This might happen if _atomistica.Neighbors or its methods are not found
    # or if there's an issue during neighbor list construction.
    print(f"Error creating neighbor list: {e}")

```

## Dependencies and Interactions

*   **`_atomistica` (C-extension/compiled library)**: This is the primary dependency. The module directly imports and uses the `Particles` and `Neighbors` classes from `_atomistica` (e.g., `from _atomistica import *`). Without a correctly compiled and accessible `_atomistica` module, `native.py` cannot function.
*   **`numpy`**: Used for array operations, such as handling periodic boundary condition flags (`np.array(atoms.get_pbc())`).
*   **`ase.data.atomic_numbers`**: Used to convert chemical symbols from ASE `Atoms` objects into integer atomic numbers for the `Particles` object.
*   This module acts as a low-level bridge, enabling Python scripts that primarily use ASE data structures to also leverage Atomistica's compiled routines that expect data in the native `Particles` and `Neighbors` format.
```
