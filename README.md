# Wave Equation Solver (MPI + deal.II)

This repository contains a fully implemented numerical solver for the wave equation, developed as part of a PDE / Numerical Analysis project.

The solver combines finite element spatial discretization (via deal.II) with Newmark time integration, and is designed to run in parallel using MPI.

A key design choice of this project is the strict separation between code and configuration: all physical and numerical parameters are specified in an external parameter file, allowing the simulation setup to be modified without recompiling the code.

## Overview

- **Equation:** Wave equation
- **Spatial discretization:** Finite Elements (deal.II)
- **Time integration:** Newmark method (second-order accurate)
- **Parallelization:** MPI
- **Execution environment:** Apptainer container

The code is designed to be modular, extensible, and reproducible, making it suitable for both experimentation and further numerical analysis.

## Key Features

- Second-order accurate Newmark time integration
- Fully MPI-parallel solver
- Clean separation between numerical code and runtime configuration
- No recompilation required when changing physical or numerical parameters
- Output compatible with ParaView for visualization
- Designed for execution inside a controlled Apptainer container

## Requirements

The code must be compiled and executed inside an Apptainer container with the following modules available:

- `gcc-glibc`
- `dealii`
- MPI (provided inside the container)

Required build and runtime tools:

- `cmake`
- `make`
- `mpirun`

## Build Instructions

From the root directory of the repository, create a build directory and compile the project using CMake:

```bash
mkdir build
cd build
cmake ..
make

```

This will generate the executable: `wave_solver`

The `build/` directory is intentionally excluded from version control.

## Running the Solver (MPI)

The solver is executed using `mpirun` and reads all simulation parameters at runtime from an external configuration file.

### Execution syntax

```bash
mpirun -n <N_PROCS> ./wave_solver ../include/parameters.prm

```

* `<N_PROCS>`: number of MPI processes
* `parameters.prm`: runtime configuration file

### Example

```bash
mpirun -n 4 ./wave_solver ../include/parameters.prm

```

## Runtime Configuration (parameters.prm)

All modifiable parameters are defined in the following file:

```text
include/parameters.prm

```

The parameter file contains a structured and documented configuration that allows the user to modify:

* Time step size
* Final simulation time
* Mesh resolution
* Physical parameters
* Output and solver settings

No recompilation is required when modifying these values. This design enables rapid testing and experimentation, while keeping the numerical code clean, modular, and independent of the simulation setup.

## Repository Structure

The repository is organized to clearly separate source code, configuration files, and documentation:

```text
.
├── CMakeLists.txt          # Build configuration
├── cmake-common.cmake      # Shared CMake settings
├── main.cpp                # Entry point
├── plot_energy.py          # Post-processing and energy plots
├── README.md               # Documentation
├── include/                # Configuration and helpers
│   ├── input_parser.cpp    # Parameter parsing implementation
│   ├── input_parser.hpp    # Parameter parsing header
│   ├── parameters.prm      # Runtime configuration file
│   └── test_cases.hpp      # Test cases definitions
├── wave_equation/          # Core solver logic
│   ├── wave_equation.cpp   # Numerical scheme implementation
│   └── wave_equation.hpp   # Solver class declaration
└── repos/                  # Reports and references
    └── report_finale_pde_parziale.pdf

─ build/                  # Build directory (not versioned)
```

## Output and Visualization

The solver produces output files compatible with ParaView, enabling efficient post-processing and visualization of the solution.

### Visualization steps

1. Open ParaView
2. Navigate to the output directory generated inside `build/`
3. Load the `.pvtu` file
4. Click Apply

### 2D visualization tips

* Apply "Warp By Scalar"
* Scalar: displacement
* Normal: 0 0 1
* Disable 2D mode to freely rotate the view

## Mathematical Background

The project includes a detailed mathematical formulation of the problem, covering:

* Strong and weak forms of the wave equation
* Spatial discretization using the finite element method
* Semi-discrete Newmark time integration scheme

For readability reasons, the full mathematical derivation is not reproduced directly in this README. Refer instead to:

* `repos/`: detailed mathematical derivations and notes
* `wave_equation.cpp`: implementation details of the numerical scheme

## Notes

This repository represents a complete and functional implementation of a numerical solver for the wave equation.

The focus of the project is on:

* correctness of the numerical scheme,
* clarity of the software structure,
* reproducibility of results,
* separation between implementation and configuration.

The code is intended primarily for educational and research purposes, and can serve as a solid starting point for extensions such as alternative time integration schemes, different boundary conditions, or more advanced parallel strategies.

## Scope and Intended Use

This project was developed as part of a PDE / Numerical Analysis coursework, with the following goals:

* Implement a physically consistent wave equation solver
* Combine finite element discretization with second-order time integration
* Exploit MPI parallelism in a clean and scalable way
* Enable fast experimentation through runtime configuration files

While the solver is not optimized for large-scale production runs, it is designed to be robust, readable, and extensible, making it suitable for academic use and further development.

## Final Remarks

* All numerical and physical parameters are controlled via `parameters.prm`
* No recompilation is required to change the simulation setup
* The build system is fully based on CMake
* Visualization is handled externally via ParaView

This repository provides a self-contained and well-documented reference implementation of a parallel wave equation solver using deal.II and MPI

