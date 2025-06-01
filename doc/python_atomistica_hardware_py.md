# src/python/atomistica/hardware.py

## Overview

This Python module, `hardware.py`, provides utilities for interacting with high-performance computing (HPC) cluster environments. Its main features include:
1.  Auto-detection of the cluster architecture based on hostname or system files.
2.  A database (`_hardware_info`) containing configurations for known clusters (e.g., "bwUniCluster", "jureca", "nemo", "justus"), including scheduler type (MOAB, SLURM), cores per node, etc.
3.  Generation of job submission scripts tailored to the detected or specified cluster architecture and its scheduler.
4.  Helper functions for time conversion (seconds to DHMS/HMS format).

The goal is to simplify the process of submitting Atomistica (or other Python-based) computations to different HPC systems.

## Key Components

### Module-Level Dictionaries

*   `moab :: dict`:
    *   A dictionary defining templates and commands for the MOAB Torque/PBS job scheduler. Includes keys like `cmdstr` (e.g., `#MSUB `), `jobid` (e.g., `$MOAB_JOBID`), `mpirun`, `name` (job name flag), `nodes` (nodes request flag), `ppn` (processors-per-node flag), `walltime`.

*   `_hardware_info :: dict`:
    *   A database where keys are strings representing known HPC cluster architectures (e.g., "bwUniCluster", "jureca").
    *   Each value is a dictionary containing:
        *   `"cores_per_node" :: int`: The number of physical CPU cores per compute node.
        *   `"loginnodes" :: list[str]`: A list of regular expression patterns used to identify the cluster from system hostnames.
        *   `"modules" :: list[str]` (optional): A list of environment modules that might be relevant for that cluster.
        *   `"scheduler" :: dict`: A dictionary describing the job scheduler, similar in structure to the `moab` dictionary, but specific to that cluster (e.g., it might contain SLURM directives for "jureca").

### Time Conversion Functions

*   `dhms(secs)`:
    *   Converts a duration in seconds into a list of four integers: `[days, hours, minutes, seconds]`.
*   `hms(secs)`:
    *   Converts a duration in seconds into a list of three integers: `[hours, minutes, seconds]`.
*   `hms_string(secs)`:
    *   Converts a duration in seconds into a formatted string "HH:MM:SS" (e.g., "02:00:45").

### Classes

*   `ComputeCluster`
    *   **Description**: A class to represent and interact with a specific compute cluster architecture.
    *   **`__init__(self, architecture=None)`**:
        *   Constructor. If `architecture` (a string like "jureca") is provided, it attempts to load the configuration from `_hardware_info`.
        *   If `architecture` is `None`, it tries to auto-detect the current cluster by:
            1.  Checking for a Juelich-specific system name file (`/etc/FZJ/systemname`).
            2.  Reading the `HOSTNAME` environment variable.
            3.  Using `socket.gethostname()`.
            4.  Executing the `hostname -s` command.
        *   The obtained hostname is matched against the `loginnodes` regex patterns in `_hardware_info` to identify the cluster.
        *   If the architecture cannot be determined or is not found in `_hardware_info`, a `KeyError` is raised.
    *   **`list_architectures(self)`**:
        *   Returns a string listing the names of all known architectures defined in `_hardware_info`.
    *   **`write(self, filename=None, **set)`**:
        *   **Description**: Generates a job submission script for the cluster represented by the `ComputeCluster` instance.
        *   **Arguments**:
            *   `filename :: str, optional`: The name of the job script file to be written. Defaults to `run.<self.arch>`.
            *   `**set`: Keyword arguments specifying job parameters:
                *   `name :: str`: Job name.
                *   `cores :: int`: Total number of CPU cores requested for the job.
                *   `smt :: bool`: If true, indicates that Simultaneous Multi-Threading (SMT) or Hyper-Threading is to be used, effectively doubling the `cores_per_node` for node calculation.
                *   `time :: int`: Requested walltime in seconds.
                *   `mail :: str, optional`: User's email address for job notifications.
                *   `wd :: str`: The working directory where the job should run.
                *   `script :: str`: The Python script to be executed by the job.
                *   `parameters :: str, optional`: Command-line parameters to be passed to the `script`.
                *   `out :: str`: Basename for the standard output file (job ID will be appended).
                *   `err :: str`: Basename for the standard error file (job ID will be appended).
        *   **Functionality**:
            1.  Calculates the number of nodes (`nodes`) and processes per node (`ppn`) based on `cores` and the cluster's `cores_per_node` (adjusted for `smt`).
            2.  Writes the job script line by line, starting with `#!/bin/bash -x`.
            3.  Uses the `cmdstr` (e.g., `#MSUB ` or `#SBATCH `) from the cluster's scheduler configuration to write directives for job name, nodes, ppn (or tasks per node for SLURM), walltime (formatted using `hms_string`), and mail notifications.
            4.  Sets the working directory (`cd <wd>`).
            5.  Sets `OMP_NUM_THREADS` environment variable (Note: sets it to `$PBS_NP`, which is specific to PBS/MOAB; for SLURM, a different variable like `$SLURM_CPUS_PER_TASK` or a calculated value might be more appropriate).
            6.  Exports the current `MODULEPATH` and reloads `LOADEDMODULES` to replicate the user's environment.
            7.  Constructs the command line to execute the Python `script` with its `parameters`, redirecting stdout and stderr to files named with `out`/`err` basenames and the scheduler's job ID variable.
        *   **Returns**: The `filename` of the generated script.

## Usage Examples

```python
from atomistica.hardware import ComputeCluster, hms_string

# Assuming the script is run on a known cluster, or specify architecture:
# cluster = ComputeCluster(architecture="jureca")
try:
    cluster = ComputeCluster() # Auto-detect
    print(f"Detected cluster: {cluster.arch}")

    job_settings = {
        'name': 'MyAtomisticaJob',
        'cores': 48,
        'smt': False,
        'time': 3600 * 4,  # 4 hours in seconds
        'mail': 'user@example.com',
        'wd': '/path/to/my/simulation/directory',
        'script': 'run_simulation.py',
        'parameters': '--input params.in',
        'out': 'sim_output',
        'err': 'sim_error'
    }
    script_file = cluster.write(**job_settings)
    print(f"Job script written to: {script_file}")

except KeyError as e:
    print(e)
    # print("Available architectures:")
    # print(ComputeCluster().list_architectures()) # Need instance if __init__ fails
```

## Dependencies and Interactions

*   **`os`**: Used for accessing environment variables (`HOSTNAME`, `MODULEPATH`, `LOADEDMODULES`) and for `os.path.isfile`, `os.popen4` (fallback hostname).
*   **`re`**: Used for regular expression matching of hostnames against patterns in `_hardware_info`.
*   **`socket`** (optional, try-except import): Used as one method to get the system hostname.
*   The module's main functionality, job script generation, is highly dependent on the structure and content of the `_hardware_info` dictionary and the `moab` template for scheduler command formatting.
*   The environment variable for setting `OMP_NUM_THREADS` might need adjustments for SLURM systems if `$PBS_NP` is not available or appropriate.
```
