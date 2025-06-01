# src/python/atomistica/io.py

## Overview

This Python module, `io.py`, provides convenience functions for file input/output (I/O) operations related to atomic configurations within the Atomistica framework. It offers unified `read` and `write` functions that attempt to automatically detect the file format based on the filename extension. Depending on the detected format, these functions delegate the actual I/O operation to either Atomistica's native routines (from `atomistica.mdcore_io`) or to the versatile I/O capabilities of the Atomistic Simulation Environment (ASE).

The module includes specific handling for:
*   Atomistica's native binary formats (`.out`, `.dat`).
*   NetCDF trajectory files (`.nc`), including reading specific frames.
*   LAMMPS data files (`.lammps`).
*   Other formats supported by ASE.

## Key Components

### Functions

*   `read(fn, **kwargs)`
    *   **Description**: Reads an atomic configuration from the file specified by `fn`. It determines the reading method based on the file extension.
    *   **Arguments**:
        *   `fn :: str`: The path to the file to be read. For NetCDF files, a specific frame can be requested by appending `@<frame_number>` (e.g., `"trajectory.nc@5"`).
        *   `**kwargs`: Additional keyword arguments to be passed to the underlying ASE I/O function (`ase.io.read` or `NetCDFTrajectory`).
    *   **Behavior**:
        *   If `fn` ends with `.out` or `.dat`:
            *   It uses `atomistica.mdcore_io.read_atoms(fn, ...)`.
            *   It checks for the existence of a `cyc.dat` file in the same directory as `fn`. If found, it's passed as the `cycfn` argument to `read_atoms`, which typically contains cycle/time information for Atomistica snapshots.
        *   If `fn` ends with `.nc`:
            *   It uses `ase.io.NetCDFTrajectory`.
            *   If a frame number is specified (e.g., `traj.nc@5`), that specific frame is read and returned.
            *   Otherwise (e.g., `traj.nc`), the last frame in the trajectory is read and returned.
        *   For any other file extension:
            *   It delegates to `ase.io.read(fn, **kwargs)`.
    *   **Returns**: An ASE `Atoms` object or a list of `Atoms` objects, depending on the file format and ASE's behavior.

*   `write(fn, a, **kwargs)`
    *   **Description**: Writes an ASE `Atoms` object `a` to the file specified by `fn`. The writing method is chosen based on the file extension.
    *   **Arguments**:
        *   `fn :: str`: The path to the file to be written.
        *   `a :: ase.Atoms`: The ASE `Atoms` object (or a list of `Atoms` objects for trajectory formats) to be written.
        *   `**kwargs`: Additional keyword arguments to be passed to the underlying ASE I/O function (`ase.io.write` or specific writers like `write_lammps_data`).
    *   **Behavior**:
        *   If `fn` ends with `.out` or `.dat`:
            *   It uses `atomistica.mdcore_io.write_atoms(fn, a)`.
        *   If `fn` ends with `.lammps`:
            *   It uses `write_lammps_data(fn, a, velocities=True, **kwargs)`. The function attempts to import `write_lammps_data` first from `ase.calculators.lammps` and falls back to `ase.calculators.lammpsrun` for compatibility with older ASE versions.
        *   If `fn` ends with `.nc`:
            *   It uses `ase.io.NetCDFTrajectory(fn, 'w').write(a)` to write or append to a NetCDF trajectory.
        *   For any other file extension:
            *   It delegates to `ase.io.write(fn, a, **kwargs)`.

### Classes/Types
This module does not define any classes.

## Important Variables/Constants
This module does not define any public module-level constants.

## Usage Examples

```python
from ase.build import bulk
from atomistica.io import read, write

# Create an example Atoms object
si_atoms = bulk('Si', 'diamond', a=5.43)

# Write to Atomistica native format
try:
    write('si_atomistica.out', si_atoms)
    re_read_atoms = read('si_atomistica.out')
    print(f"Read back {len(re_read_atoms)} atoms from si_atomistica.out")
except Exception as e:
    print(f"Error with Atomistica native I/O: {e}")

# Write to XYZ format (delegates to ASE)
write('si_ase.xyz', si_atoms)
re_read_xyz = read('si_ase.xyz')
print(f"Read back {len(re_read_xyz)} atoms from si_ase.xyz")

# Write to LAMMPS data format
try:
    write('si.lammps', si_atoms, atom_style='atomic') # atom_style is a kwarg for write_lammps_data
    print("Wrote si.lammps data file.")
except Exception as e:
    print(f"Error writing LAMMPS data: {e}")

# Example for NetCDF (if ASE has NetCDF support)
try:
    traj_writer = NetCDFTrajectory('test_traj.nc', 'w')
    traj_writer.write(si_atoms)
    si_atoms.positions[0,0] += 0.1
    traj_writer.write(si_atoms)
    traj_writer.close()

    # Read last frame
    last_frame = read('test_traj.nc')
    print(f"Last frame from NetCDF has {len(last_frame)} atoms.")
    # Read specific frame
    first_frame = read('test_traj.nc@0')
    print(f"First frame from NetCDF has {len(first_frame)} atoms.")
except NameError: # NetCDFTrajectory might not be imported if ASE lacks it
    print("NetCDFTrajectory not available in ASE.")
except Exception as e:
    print(f"Error with NetCDF I/O: {e}")

```

## Dependencies and Interactions

*   **`os`**: Used for path manipulations (`os.path.dirname`, `os.path.exists`) to locate auxiliary files like `cyc.dat`.
*   **`ase.io`**: A major dependency for reading and writing a wide variety of common atomic simulation file formats (e.g., XYZ, PDB, CIF) and for `NetCDFTrajectory`.
*   **`ase.calculators.lammps.write_lammps_data`** (or `ase.calculators.lammpsrun.write_lammps_data`): Specifically used for writing LAMMPS data files. The module includes a try-except block to handle different import paths for this function across ASE versions.
*   **`atomistica.mdcore_io`**: Provides `read_atoms` and `write_atoms` for handling Atomistica's native binary file formats (`.out`, `.dat`).
*   The module attempts to import `NetCDFTrajectory` in a try-except block, allowing it to function (for other formats) even if NetCDF support is not available in the installed ASE version.
```
