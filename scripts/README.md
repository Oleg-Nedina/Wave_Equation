# Scripts and Post-Processing Tools

This directory contains a suite of Python tools designed to automate simulation execution and validate the numerical and physical properties of the implemented wave equation solver.

## Directory Structure

### 1. Performance Analysis (Scaling)

These scripts evaluate the MPI parallel efficiency of the C++ solver.

* **`strong_scaling.py`**: Measures the reduction in execution time for a fixed global problem size () as the number of processors increases from 1 to 4. It calculates the **Speedup** () and generates a plot comparing measured performance against ideal linear speedup.
* **`weak_scaling.py`**: Evaluates parallel efficiency by increasing both the problem size and the computational resources proportionally (e.g., 1 core at Refinement 7 vs. 4 cores at Refinement 8). It calculates **Parallel Efficiency** and generates a corresponding visualization.

### 2. Physical and Numerical Validation

These scripts use the solver's output to verify energy properties and numerical behavior.

* **`plot_energy_conservation.py`**: Verifies that the total energy of the system remains constant over time using the **Implicit Newmark** scheme with conservative parameters ().
* **`plot_energy_dissipation.py`**: Demonstrates **algorithmic numerical dissipation** by setting  (specifically ). It illustrates how total energy artificially decreases over time due to the integration scheme.
* **`plot_dispersion_comparison.py`**: Analyzes **numerical dispersion** by comparing phase errors on a sparse grid () between the **Explicit FE(1)** and **Implicit FE(2)** solvers.
* **`plot_energy_work_overlap.py`**: Designed for **Scenario 6 (Pumping Wall)**, it compares the total computed energy  with the theoretical work  performed by the moving boundary, validating physical consistency in forced systems.

## Requirements and Usage

### Python Dependencies

The scripts require the following packages:

* `pandas`
* `matplotlib`
* `numpy`

### Execution Instructions

The scaling scripts are designed to be executed from the `scripts/` directory to properly locate the compiled solver and parameter files.

1. **Generate Scaling Data**:

```bash
python3 strong_scaling.py
python3 weak_scaling.py

```

1. **Generate Validation Plots**:
Ensure the solver has been compiled in the `../build/` directory before running:

```bash
python3 plot_energy_conservation.py
python3 plot_energy_dissipation.py

```

## Theoretical Background

The validation scripts specifically target properties of the **Newmark-Beta** method and the **Finite Element** discretization.

The algorithmic dissipation is controlled by the  parameter; while  is second-order accurate and energy-conserving, values  introduce high-frequency damping often used to stabilize numerical oscillations at the cost of energy loss. Similarly, the dispersion plot highlights how higher-order elements (FE(2)) significantly reduce the phase lag common in low-order wave simulations.
