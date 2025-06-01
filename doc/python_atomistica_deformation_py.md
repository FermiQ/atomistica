# src/python/atomistica/deformation.py

## Overview

This Python module, `deformation.py`, provides tools for handling and analyzing homogeneously deformed simulation volumes, with a particular focus on simple shear deformations. Such deformations are common in simulations studying mechanical properties, but output from simulation codes (like LAMMPS) can sometimes involve cell "flipping" or remapping when shear strains become large, making direct analysis of strain history complex. This module offers utilities to address these issues.

## Key Components

### Functions

*   `get_shear_distance(a)`
    *   **Description**: This function calculates the accumulated shear displacement in a simulation cell. It is designed to work with cells undergoing simple shear, typically in the xy-plane where the x and y components of the third cell vector (`a.cell[2]`) represent the shear. It can also account for Lees-Edwards boundary conditions if the shear information is stored in `a.info['shear_dx']`.
    *   The function asserts that other off-diagonal components of the cell matrix are zero, consistent with simple xy-shear.
    *   **Arguments**:
        *   `a`: An ASE (Atomistic Simulation Environment) `Atoms` object.
    *   **Returns**: A tuple `(dx, dy)` representing the accumulated shear distances in the x and y directions.

### Classes

*   `RemoveSimpleShearDeformation`
    *   **Description**: This class wraps an iterable ASE trajectory (e.g., a list of `Atoms` objects or an ASE trajectory reader object). Its primary purpose is to "unwrap" the shear deformation history from trajectories where the simulation cell might have been remapped or "flipped" when shear strains exceeded values like +/- 0.5 of the cell dimension. This is a common behavior in simulation codes like LAMMPS to keep the cell reasonably "close" to orthorhombic.
    *   The class reconstructs a continuous history of shear deformation and allows access to `Atoms` objects where atomic positions have been remapped to an unsheared reference frame, while storing the true (unwrapped) sheared cell information.
    *   **`__init__(self, traj)`**:
        *   The constructor takes an ASE trajectory object (`traj`).
        *   Initializes internal lists: `self.last_d` (to store the last unwrapped (dx, dy) shear distances), `self.sheared_cells` (to store the reconstructed true sheared cell vectors for each frame), and `self.unsheared_cells` (to store the diagonal/orthorhombic cell vectors corresponding to each frame).
    *   **`_fill_cell_info_upto(self, i)`**:
        *   A private helper method that processes the trajectory up to frame `i`.
        *   It iteratively calls `get_shear_distance` on each frame.
        *   It "unwraps" the raw shear distances `(cur_dx, cur_dy)` obtained from each frame by comparing them to the `last_dx, last_dy` from the previously processed frame. If a jump larger than half the cell dimension in x or y is detected (e.g., from +0.45 strain to -0.45 strain), it adds/subtracts the cell dimension (sx or sy) to reconstruct the continuous, unwrapped shear (`dx`, `dy`).
        *   Stores the unwrapped `(dx, dy)` in `self.last_d`.
        *   Stores the reconstructed true sheared cell (e.g., `[[sx,0,0],[0,sy,0],[dx,dy,sz]]`) in `self.sheared_cells` and the corresponding diagonal cell `[sx,sy,sz]` in `self.unsheared_cells`.
    *   **`__getitem__(self, i)`**:
        *   Allows accessing a processed `Atoms` object for frame `i` (supports negative indexing).
        *   Ensures trajectory information is processed up to frame `i` by calling `_fill_cell_info_upto(i)`.
        *   Retrieves the original `Atoms` object for frame `i` from the input trajectory.
        *   Sets the cell of this `Atoms` object to the true, unwrapped sheared cell (`self.sheared_cells[i]`) without scaling atoms.
        *   Then, sets the cell to the diagonal, unsheared cell (`self.unsheared_cells[i]`) with `scale_atoms=True`. This operation remaps atomic coordinates from the sheared frame to the unsheared orthorhombic reference frame.
        *   Wraps atomic positions into the unit cell using `a.set_scaled_positions(a.get_scaled_positions()%1.0)`.
        *   Stores the true sheared cell in `a.info['true_cell']` for reference.
        *   Returns the modified `Atoms` object.
    *   **`__len__(self)`**: Returns the total number of frames in the trajectory.

## Important Variables/Constants
This module does not define any public module-level constants.

## Usage Examples

```python
from ase.io import read
from atomistica.deformation import RemoveSimpleShearDeformation
# Assuming 'shear_trajectory.xyz' is a trajectory file from a simulation
# with simple shear that might involve cell flipping.

original_traj = read('shear_trajectory.xyz', index=':') # Read all frames
unwrapped_traj = RemoveSimpleShearDeformation(original_traj)

# Access a specific frame (e.g., the 10th frame)
# atoms_frame_10_unsheared will have atomic positions remapped to an
# orthorhombic cell, and atoms_frame_10_unsheared.info['true_cell']
# will contain the actual, continuously deformed cell for that frame.
if len(unwrapped_traj) > 10:
    atoms_frame_10_unsheared = unwrapped_traj[10]
    print("Original cell from file for frame 10:", original_traj[10].cell)
    print("Unwrapped true sheared cell for frame 10:", atoms_frame_10_unsheared.info['true_cell'])
    print("Cell of returned Atoms object (unsheared):", atoms_frame_10_unsheared.cell)

# Iterate through the whole trajectory
# for frame_index in range(len(unwrapped_traj)):
#     atoms_unsheared = unwrapped_traj[frame_index]
#     # Perform analysis on atoms_unsheared, using atoms_unsheared.info['true_cell']
#     # if the actual deformed cell is needed.
```

## Dependencies and Interactions

*   **`numpy`**: Used for numerical operations, particularly for handling cell vectors and arrays.
*   **ASE (Atomistic Simulation Environment)**: The module operates on ASE `Atoms` objects and trajectory objects (anything iterable that yields `Atoms` objects). It uses `a.cell`, `a.get_cell()`, `a.set_cell()`, `a.get_scaled_positions()`, `a.set_scaled_positions()`, `a.info`.
*   The `RemoveSimpleShearDeformation` class is particularly useful for post-processing simulation trajectories where simple shear has been applied, ensuring that analysis (like strain calculation or defect tracking) is performed with a consistent reference frame that accounts for the total applied shear, irrespective of cell remapping artifacts.
```
