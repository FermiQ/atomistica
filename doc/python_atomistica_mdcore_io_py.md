# src/python/atomistica/mdcore_io.py

## Overview

This Python module, `mdcore_io.py`, provides input and output functionalities for a specific, somewhat deprecated, file format associated with the MDCORE or tblmd simulation codes, previously used at Fraunhofer IWM. The primary files handled are typically named `atoms.dat` (for atomic coordinates, types, and other per-atom data) and `cyc.dat` (for simulation cell parameters).

These functions are often used as backends by the more general `atomistica.io.read` and `atomistica.io.write` convenience functions when `.out` or `.dat` file extensions are detected.

## Key Components

### Functions

*   `read_atoms(fn, cycfn=None, pos_only=False, conv=1.0)`
    *   **Description**: Reads atomic configuration data from an `atoms.dat`-like file.
    *   **Arguments**:
        *   `fn :: str`: Path to the `atoms.dat` file.
        *   `cycfn :: str, optional`: Path to an associated `cyc.dat` file containing cell information. If provided, `read_cyc` is called.
        *   `pos_only :: bool, optional`: If `True`, only atomic positions and types are read. If `False` (default), the function attempts to read subsequent sections like "VELOCITIES", "FORCES", "CHARGES", "CELL", and other custom arrays.
        *   `conv :: float, optional`: A conversion factor applied to the coordinates read from the file. Default: `1.0`.
    *   **File Format Details Handled**:
        *   Skips initial header lines (starting with '#', '<').
        *   Reads the number of atoms.
        *   Reads atom lines: `type x y z group gamma T`. `type` can be an integer (atomic number Z) or a chemical symbol.
        *   Parses optional sections marked by headers like `<--- VELOCITIES --->`.
        *   Stores read data into an ASE `Atoms` object, including standard arrays (`positions`, `momenta`, `forces`, `charges`) and any custom arrays found.
    *   **Returns**: An `ase.Atoms` object.

*   `read_cyc(this, fn, conv=1.0)`
    *   **Description**: Reads simulation cell (lattice vectors) information from a `cyc.dat` file and applies it to an existing ASE `Atoms` object `this`.
    *   **Arguments**:
        *   `this :: ase.Atoms`: The ASE `Atoms` object to which the cell information will be applied.
        *   `fn :: str`: Path to the `cyc.dat` file.
        *   `conv :: float, optional`: A conversion factor applied to the cell vector components. Default: `1.0`.
    *   **File Format Details Handled**: Reads the 9 components of the three cell vectors, skipping header lines and placeholder lines for barostat settings. Sets periodic boundary conditions (`pbc=True`) on the `Atoms` object.

*   `write_atoms(fn, this, cycfn=None, conv=1.0, symbols=True)`
    *   **Description**: Writes an ASE `Atoms` object `this` to the `atoms.dat` format.
    *   **Arguments**:
        *   `fn :: str`: Path to the output `atoms.dat` file.
        *   `this :: ase.Atoms`: The ASE `Atoms` object to write.
        *   `cycfn :: str, optional`: If a path is provided, `write_cyc` is called to write a corresponding `cyc.dat` file.
        *   `conv :: float, optional`: A conversion factor applied to coordinates before writing. Default: `1.0`.
        *   `symbols :: bool, optional`: If `True` (default), writes element types as chemical symbols. If `False`, writes them as atomic numbers (Z).
    *   **File Format Details Generated**:
        *   Writes standard header lines and the number of atoms.
        *   Writes atom data: `type mass x y z group gamma T`. `group`, `gamma`, and `T` values are taken from `this.get_array("groups")` etc., if present, otherwise defaults (1, 0.0, 0.0) are used. Mass is taken from `ase.data.atomic_masses`.
        *   Writes a "VELOCITIES" section (calculating velocities from momenta and masses).
        *   Writes a "FORCES" section (currently writes "0 0 0" for all forces, meaning it's a placeholder or expects forces to be added by other means if needed in this format).
        *   Writes a "CELL" section with the cell vectors.
        *   Writes any other per-atom arrays found in `this.arrays` (excluding those in `atoms_default_fields`) under their respective names as headers.

*   `write_cyc(fn, this, conv=1.0)`
    *   **Description**: Writes the cell information from an ASE `Atoms` object `this` to the `cyc.dat` format.
    *   **Arguments**:
        *   `fn :: str`: Path to the output `cyc.dat` file.
        *   `this :: ase.Atoms`: The ASE `Atoms` object supplying the cell vectors.
        *   `conv :: float, optional`: A conversion factor applied to cell vector components before writing. Default: `1.0`.
    *   **File Format Details Generated**: Writes the specific formatted structure of a `cyc.dat` file, including placeholders for barostat status, initial and final box vectors (written as the same), and initial/final stress tensors (written as zeros).

### Module-Level Constants

*   `atoms_default_fields :: numpy.ndarray`: An array of strings listing ASE array names that are handled in a specific way by `write_atoms` (e.g., "positions", "momenta", "numbers") and are not to be written out as generic auxiliary arrays.

## Usage Examples

These functions are primarily intended to be used by the main `atomistica.io.read()` and `atomistica.io.write()` wrapper functions when files with `.out` or `.dat` extensions are encountered.

```python
from ase.build import bulk
from atomistica.mdcore_io import read_atoms, write_atoms, read_cyc, write_cyc

# Create an example Atoms object
atoms = bulk('Si', 'diamond', a=5.43, cubic=True)
atoms.set_array("groups", np.ones(len(atoms), dtype=int)) # Example auxiliary data

# Write to Atomistica native format
write_atoms('atoms.dat', atoms, cycfn='cyc.dat')

# Read back
read_back_atoms = read_atoms('atoms.dat', cycfn='cyc.dat')

print(f"Original cell:\n{atoms.cell}")
print(f"Read back cell:\n{read_back_atoms.cell}")
print(f"Original number of atoms: {len(atoms)}")
print(f"Read back number of atoms: {len(read_back_atoms)}")
if read_back_atoms.has("groups"):
    print("Read back 'groups' array successfully.")
```

## Dependencies and Interactions

*   **`os`**: Used for path manipulation.
*   **`numpy`**: Used for array creation and data handling.
*   **`ase`**: Core ASE modules (`ase.Atoms`, `ase.Atom`) are used for creating the atoms objects.
*   **`ase.data`**: Used for accessing `atomic_masses` and `chemical_symbols`.
*   **`ase.parallel.paropen`**: Used for file I/O, ensuring that in an MPI parallel run, only the master process performs the actual file operations.
*   The format is described as "deprecated," suggesting that while these I/O routines provide compatibility, newer or more standard formats might be preferred for general use.
```
