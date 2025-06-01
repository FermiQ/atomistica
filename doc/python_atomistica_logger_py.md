# src/python/atomistica/logger.py

## Overview

This Python module, `logger.py`, provides the `MDLogger` class, designed for logging various physical quantities during molecular dynamics (MD) simulations performed with the Atomistic Simulation Environment (ASE). It allows for flexible configuration of what data is logged and how it is formatted. The logger can be attached to an ASE dynamics object and will write data at specified intervals to a file or standard output.

## Key Components

### Classes

*   `MDLogger`
    *   **Description**: This class facilitates the logging of MD simulation data. It can track and record time, total energy, potential energy, kinetic energy, temperature, stress tensor, cell parameters, and volume. Logging can be configured to be per-atom for energies.
    *   **`__init__(self, dyn, atoms, logfile, header=True, stress=False, cell=False, volume=False, peratom=False, hiprec=False, mode="a")`**:
        *   **Arguments**:
            *   `dyn`: The ASE dynamics object (e.g., `VelocityVerlet`, `Langevin`). A weak reference is kept to avoid circular dependencies. This can be `None` if time logging from a dynamics object is not needed.
            *   `atoms :: ase.Atoms`: The `Atoms` object representing the system being simulated.
            *   `logfile :: str or file-object`: The destination for the log data. Can be a filename, an already open file object, or "-" to write to standard output (`sys.stdout`).
            *   `header :: bool, optional`: If `True` (default), a header line describing the columns is written at the beginning of the log.
            *   `stress :: bool, optional`: If `True`, the six components of the stress tensor (in GPa) are logged. Default: `False`.
            *   `cell :: bool, optional`: If `True`, the cell vectors (xx, yy, zz, xy, yz, zx components) are logged. Default: `False`.
            *   `volume :: bool, optional`: If `True`, the cell volume (in Å³) is logged. Default: `False`.
            *   `peratom :: bool, optional`: If `True`, energies (total, potential, kinetic) are logged as per-atom values. Default: `False`.
            *   `hiprec :: bool, optional`: If `True`, floating-point numbers are written with higher precision (`%20.12e`). Otherwise, a standard precision (`%12.4f`) is used. Default: `False`.
            *   `mode :: str, optional`: File opening mode if `logfile` is a filename (e.g., "a" for append, "w" for write). Default: "a".
        *   **Functionality**:
            *   Stores references to `dyn` and `atoms`.
            *   Opens the `logfile` if a filename is given. If `ase.parallel.rank > 0` (i.e., running with MPI and not the master process), `logfile` is redirected to `/dev/null`.
            *   Dynamically constructs a header string (`self.hdr`) and a Python format string (`self.fmt`) based on the logging options selected (time, energies, stress, cell, etc.).
    *   **`__del__(self)`**:
        *   Calls `self.close()` to ensure the log file is closed if this instance opened it.
    *   **`close(self)`**:
        *   Closes `self.logfile` if `self.ownlogfile` is true (i.e., if the `MDLogger` instance was responsible for opening the file).
    *   **`__call__(self)`**:
        *   This method is designed to be called by ASE's dynamics loop (e.g., via `AttachLogger` or as an observer function).
        *   It retrieves the current simulation data:
            *   Time: from `self.dyn.get_time()` (converted to ps) if `self.dyn` is available.
            *   Potential energy: from `self.atoms.get_potential_energy()`.
            *   Kinetic energy: from `self.atoms.get_kinetic_energy()`.
            *   Temperature: calculated from kinetic energy.
            *   Stress: from `self.atoms.get_stress()` (converted to GPa) if `self.stress` is true.
            *   Cell: from `self.atoms.get_cell()` if `self.cell` is true (logs diagonal and off-diagonal components).
            *   Volume: from `self.atoms.get_volume()` if `self.volume` is true.
        *   If `self.peratom` is true, energies are divided by the number of atoms.
        *   Formats the collected data using `self.fmt` and writes it to `self.logfile`.
        *   Flushes the `self.logfile` to ensure data is written immediately.

## Important Variables/Constants
This module does not define any public module-level constants.

## Usage Examples

```python
from ase import Atoms
from ase.md.verlet import VelocityVerlet
from ase.units import fs, kB
from atomistica.logger import MDLogger
from ase.calculators.emt import EMT # Example calculator

# Setup a simple system and dynamics
atoms = Atoms('Ar10', positions=[(i,i,i) for i in range(10)], pbc=True)
atoms.set_cell([10,10,10])
atoms.set_calculator(EMT())
atoms.set_momenta([(0.01, 0, 0) for _ in range(10)]) # Give some initial momenta

dyn = VelocityVerlet(atoms, timestep=1.0*fs)

# Create an MDLogger instance
# Log to standard output, include stress, high precision
logger = MDLogger(dyn, atoms, '-', header=True, stress=True, hiprec=True)

# Attach the logger to the dynamics (ASE typically calls attached loggers)
# For direct use as an observer, one might do:
# dyn.attach(logger, interval=10)
# Or, if the dynamics object doesn't have AttachLogger, one might call it manually in a loop:

print("# Manual logging example every step:")
for i in range(5):
    dyn.run(1) # Run 1 step
    logger()   # Call the logger instance to log data

logger.close() # Close if it was a file
```

## Dependencies and Interactions

*   **`weakref`**: Used to store a `weakref.proxy` to the dynamics object (`dyn`) to prevent circular references that might delay garbage collection.
*   **`sys`**: Used for `sys.stdout` when logging to the standard output.
*   **`ase.units`**: Used for physical constants like `kB` (Boltzmann constant) and `fs` (femtosecond) for unit conversions (e.g., time from ASE internal units to picoseconds).
*   **`ase.parallel`**: Used to check the MPI rank (`ase.parallel.rank`) so that logging is typically performed only by the master process in an MPI parallel simulation.
*   **ASE `Atoms` object**: The logger relies on the `atoms` object passed to it to provide methods like `get_potential_energy()`, `get_kinetic_energy()`, `get_stress()`, `get_cell()`, `get_volume()`, and `get_number_of_atoms()`.
*   **ASE dynamics object**: If a dynamics object (`dyn`) is provided and has a `get_time()` method, the logger will use it to record simulation time.
```
