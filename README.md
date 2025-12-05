# 🌊 NMPDE Project 2: Parallel Wave Equation Solver

**Course:** Numerical Methods for Partial Differential Equations (A.Y. 2025/2026)  
**Topic:** Scalable Finite Element Solver for the 2D/3D Wave Equation using deal.II and MPI.

## 📌 Project Overview

This project implements a parallel Finite Element Method (FEM) solver for the scalar wave equation in 2D and 3D domains. The code is built upon the **deal.II** library and leverages **MPI** for distributed memory parallelization, making it suitable for High-Performance Computing (HPC) clusters.

The solver currently focuses on the **Implicit Newmark-Beta** time integration scheme, ensuring unconditional stability for stiff problems.

### 📐 Mathematical Formulation
The problem solves the strong form:
$$
\begin{cases}
\frac{\partial^2 u}{\partial t^2} - \Delta u = f & \text{in } \Omega \times (0, T] \\
u = g(t) & \text{on } \partial\Omega \\
u(\mathbf{x},0) = u_0, \quad \dot{u}(\mathbf{x},0) = v_0 & \text{in } \Omega
\end{cases}
$$

The semi-discrete algebraic system solved at each time step is:
$$
\left( \frac{1}{\beta \Delta t^2} \mathbf{M} + \mathbf{K} \right) \mathbf{U}_{n+1} = \mathbf{F}_{n+1} + \mathbf{M} \left[ \frac{1}{\beta \Delta t^2} \mathbf{U}_n + \frac{1}{\beta \Delta t} \mathbf{V}_n + \frac{1-2\beta}{2\beta} \mathbf{A}_n \right]
$$

---

## 🚀 Implemented Features

* **MPI Parallelism:** Fully distributed mesh handling (`parallel::distributed::Triangulation`) and linear algebra using **Trilinos** wrappers.
* **Dimension Independence:** The same codebase compiles and runs for **2D** and **3D** problems via C++ templates.
* **Implicit Time Stepping:** Implementation of the A-Stable Newmark-$\beta$ scheme ($\beta=0.25, \gamma=0.5$).
* **Consistent Initialization:** Solves $M \mathbf{a}_0 = \mathbf{F}(0) - \mathbf{K} \mathbf{u}_0$ to compute the correct initial acceleration consistent with the FEM formulation.
* **Scenario Strategy Pattern:** Easy switching between different physical test cases without recompiling the core solver logic.
* **Parallel Output:** Efficient `.pvtu` / `.vtu` output generation in a dedicated directory.

---

## 🧪 Test Scenarios

The solver includes 5 built-in scenarios to verify correctness and demonstrate physics. You can switch scenarios by changing the `scenario_id` in `wave_equation.cpp` (setup section).

| ID | Name | Description | Purpose |
| :--- | :--- | :--- | :--- |
| **1** | **Manufactured Solution** | $u = \sin(\pi t)\prod \sin(\pi x_i)$ | **Verification**. Checks $L^2$ error convergence order. |
| **2** | **Gaussian Pulse** | Single pulse starting from center. | Visualizing propagation and reflection. |
| **3** | **Collision** | Two pulses starting from corners. | Testing wave superposition/interference. |
| **4** | **Raindrop Symphony** | Multiple drops of varying amplitude. | Complex interference patterns (Chaos). |
| **5** | **Ring Implosion** | A ring collapsing into a singularity. | Stress-testing the solver stability. |

---

## 🛠️ Build and Run Instructions

### Prerequisites
* **deal.II** (v9.3 or newer) with MPI and Trilinos support.
* **MPI** (OpenMPI or MPICH).
* **CMake** and a C++17 compatible compiler.

### Compilation
```bash
mkdir build
cd build
cmake ..
make -j4  # Compiles using a suitable number (just let -j for the defoult)
```


### Running the Solver
The executable `wave_solver` accepts command-line arguments to select the dimension and the parameter file.

**Syntax:**
```bash
mpirun -n <N_PROCS> ./wave_solver [DIMENSION] [PARAM_FILE]
```

### Examples:
**Run in 2D (Default):**
```bash
mpirun -n 4 ./wave_solver
```

**Run in 3D:**
```bash
mpirun -n 4 ./wave_solver 3
```

**Run in 3D with specific parameters:**
```bash
mpirun -n 8 ./wave_solver 3 parameters.prm
```

---

## 📊 Visualization (ParaView)
1. Open ParaView.  
2. Navigate to the `build/output_results` folder.  
3. Open the group file `solution_..pvtu` (ParaView should group them automatically).  
4. Click **Apply**.

### For 2D visualization:
- Add filter: **Warp By Scalar**.  
- Scalars: `displacement`.  
- Normal: `0 0 1`.  
- Uncheck **"2D" view mode** to rotate the camera.

---

## 📝 To-Do List (Future Development)
The following features are planned or currently under development:

- [ ] **Explicit Time Stepping:** Implement the Explicit Newmark scheme (β = 0) with Mass Lumping to avoid linear solvers and increase efficiency for wave propagation.  
- [ ] **Parameter Parsing:** Fully connect the `parameters.prm` file to the C++ code (currently parameters are hardcoded in the constructor for rapid testing).  
- [ ] **Time-Dependent BCs:** Implement non-homogeneous Dirichlet boundary conditions (g(t) ≠ 0).  
- [ ] **CFL Condition Check:** Automate time-step calculation based on mesh size (h_min) for the explicit solver.

---

## 📁 Repository Structure
```
├── CMakeLists.txt         # Build configuration
├── cmake-common.cmake     # Shared CMake settings
├── main.cpp               # Entry point (MPI init & Arg parsing)
├── wave_equation.hpp      # Class declaration & Documentation
├── wave_equation.cpp      # Implementation (Logic & Math)
└── parameters.prm         # Runtime parameters configuration
```

