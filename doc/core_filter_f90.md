# src/core/filter.f90

## Overview

This module, `filter`, is designed to filter atoms based on their type. It allows other parts of the Atomistica simulation code to operate on specific subsets of atoms.

## Key Components

### Functions/Subroutines

*   `filter_from_string`: Convert a string to a filter
*   `filter_count`: Count how many atoms that match this filter we have locally
*   `filter_sum`: Sum field x
*   `filter_average`: Sum field x
*   `filter_mask`: Sum field x
*   `filter_pack`: Pack property into an array
*   `filter_unpack`: Pack property into an array
*   `filter_prlog`: Dump filter information to log file

### Modules

*   `filter`: Filter a type of atom

### Classes/Types

*   N/A

## Important Variables/Constants

*   `MAX_EL_STR`: Parameter defining the maximum length for strings holding element information.

## Usage Examples

Specific examples of how to use the `filter` module and its routines can be found in the test suites or higher-level modules within the Atomistica project that utilize this filtering capability.

```fortran
! Example:
! integer :: my_filter
! type(particles_t) :: my_particles
!
! ! Initialize my_particles appropriately
!
! ! Create a filter for Carbon atoms
! my_filter = filter_from_string("C", my_particles)
!
! ! Count Carbon atoms
! print *, "Number of Carbon atoms: ", filter_count(my_filter, my_particles)
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `supplib`: Provides supplementary utilities.
    *   `logging`: Used for logging messages.
    *   `data`: Provides core data structures, likely including particle data types.
    *   `particles`: Manages particle information.
*   **External Libraries:** None explicitly listed beyond standard Fortran capabilities.
*   **Interactions:** The `filter` module interacts with other components of the Atomistica system by providing a mechanism to select atoms based on their chemical element type. This allows for targeted operations, analysis, or modifications on specific atom types within a simulation. For example, it can be used to apply different potentials to different elements or to analyze properties of a particular subset of atoms.
