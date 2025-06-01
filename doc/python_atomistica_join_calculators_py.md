# src/python/atomistica/join_calculators.py

## Overview

This Python module, `join_calculators.py`, provides utilities for combining the effects of multiple ASE (Atomistic Simulation Environment) calculators or for applying simple, analytically defined external forces to an atomic system. It contains two main classes:
*   `JoinCalculators`: Allows multiple individual ASE calculators to be treated as a single composite calculator, where total energies, forces, and stresses are the sum of contributions from each component calculator.
*   `LinearPotential`: Implements a potential corresponding to a constant force field, useful for simulating the effect of a uniform external field or applying a constant bias force.

## Key Components

### Classes

*   `JoinCalculators`
    *   **Description**: This class acts as a composite ASE-compatible calculator. It takes a list of individual calculator instances and aggregates their results. For example, one could combine a DFT calculator with a classical dispersion correction potential (like DFT-D3). The potential energies, forces, and stress tensors obtained from each calculator in the list are summed to produce the total quantities.
    *   **`__init__(self, calcs)`**:
        *   **Arguments**:
            *   `calcs :: list`: A list of ASE-compatible calculator objects.
        *   Initializes the `JoinCalculators` instance by storing the provided list of calculators in `self.calcs`.
    *   **`get_forces(self, a)`**:
        *   Calculates the total atomic forces for the given ASE `Atoms` object `a`.
        *   It iterates through each calculator in `self.calcs`, calls its `get_forces(a)` method, and sums the resulting force arrays.
        *   **Returns**: `numpy.ndarray` - The total forces on atoms.
    *   **`get_potential_energy(self, a)`**:
        *   Calculates the total potential energy for the `Atoms` object `a`.
        *   It iterates through `self.calcs`, calls `get_potential_energy(a)` for each, and sums the energies.
        *   **Returns**: `float` - The total potential energy.
    *   **`get_stress(self, a)`**:
        *   Calculates the total stress tensor for the `Atoms` object `a`.
        *   It iterates through `self.calcs`, calls `get_stress(a)` for each, and sums the stress tensors (expected in Voigt notation, 6 components).
        *   **Returns**: `numpy.ndarray` - The total stress tensor in Voigt form.
    *   **`set_atoms(self, a)`**:
        *   Assigns the `Atoms` object `a` to each underlying calculator in `self.calcs` that has a `set_atoms` method.
    *   **NotImplemented Methods**:
        *   `get_magnetic_moments(self, a)`
        *   `get_potential_energies(self, a)` (per-atom energies)
        *   `get_spin_polarized(self)`
        *   `get_stresses(self, a)` (per-atom stresses)
        These methods will raise `NotImplementedError` if called.

*   `LinearPotential`
    *   **Description**: This class implements a potential that results in a constant force applied to atoms. This can be used to simulate, for example, the effect of a uniform electric field on charged particles or a constant pulling/pushing force in a specific direction.
    *   **`__init__(self, force, mask=None)`**:
        *   **Arguments**:
            *   `force :: array_like`: A 3-element vector (e.g., `[fx, fy, fz]`) representing the constant force to be applied.
            *   `mask :: array_like, optional`: A boolean or integer array. If provided, the force is applied only to atoms where `mask` is true (or non-zero). If `None`, the force is applied to all atoms.
    *   **`get_forces(self, a=None)`**:
        *   Calculates the forces on atoms. If `a` (ASE `Atoms` object) is not provided, it uses `self.a` (set by `set_atoms`).
        *   Returns a zero force array if no atoms are selected by the mask or if `self.force` is zero. Otherwise, applies `self.force` to the selected atoms.
        *   **Returns**: `numpy.ndarray` - The array of forces.
    *   **`get_potential_energy(self, a=None)`**:
        *   Calculates the potential energy associated with the constant force field: \(E = - \sum_i \vec{r}_i \cdot \vec{F}_i\), where the sum is over atoms selected by the mask.
        *   **Returns**: `float` - The potential energy.
    *   **`get_stress(self, a=None)`**:
        *   Returns a zero stress tensor (`numpy.zeros(6, dtype=float)`), as a uniform constant force field does not inherently contribute to stress in the standard virial definition.
    *   **`set_atoms(self, a)`**:
        *   Assigns the provided ASE `Atoms` object `a` to `self.a` for later use by `get_forces` and `get_potential_energy`.
    *   **NotImplemented Methods**:
        *   `get_potential_energies(self, a)` (per-atom energies).

## Important Variables/Constants
This module does not define any public module-level constants.

## Usage Examples

**Joining two calculators (e.g., a primary calculator and a dispersion correction):**
```python
from ase.build import bulk
from ase.calculators.lj import LennardJones # Example primary calculator
from atomistica.join_calculators import JoinCalculators
# Assuming DispCorrCalc is another ASE-compatible calculator for dispersion
# from some_dispersion_module import DispCorrCalc

atoms = bulk('Ar', 'fcc', a=5.2)

calc1 = LennardJones(epsilon=0.0103, sigma=3.40)
# calc_disp = DispCorrCalc(parameters_for_Ar) # Fictitious dispersion calculator

# combined_calc = JoinCalculators([calc1, calc_disp])
combined_calc = JoinCalculators([calc1]) # Using only one for a simpler run

atoms.set_calculator(combined_calc)
total_energy = atoms.get_potential_energy()
total_forces = atoms.get_forces()

print(f"Total Energy: {total_energy}")
```

**Applying a constant force (e.g., simulating an electric field on one atom):**
```python
import numpy as np
from ase import Atoms
from atomistica.join_calculators import LinearPotential

atoms = Atoms('Na', positions=[[0,0,0]])
# Apply a force of [0, 0, 1.0] eV/A to the first atom (index 0)
mask = np.array([True]) # Mask for the first atom
const_force_potential = LinearPotential(force=np.array([0,0,1.0]), mask=mask)

# This LinearPotential can be used alone or within JoinCalculators
atoms.set_calculator(const_force_potential)
forces = atoms.get_forces()
energy = atoms.get_potential_energy()

print("Applied force:", forces) # Will be [[0,0,1.0]]
print("Energy from constant force:", energy) # Will be 0 if atom is at origin
```

## Dependencies and Interactions

*   **`numpy`**: Used for array initializations (e.g., `np.zeros`) and numerical operations (e.g., `np.sum`, `np.dot`).
*   **ASE (Atomistic Simulation Environment)**: Both classes are designed to be compatible with the ASE calculator interface. They take ASE `Atoms` objects as arguments to their methods (e.g., `get_forces(self, a)`).
*   The `JoinCalculators` class relies on the component calculators in `self.calcs` to correctly implement the standard ASE calculator methods (`get_forces`, `get_potential_energy`, `get_stress`).
```
