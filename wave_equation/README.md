## 1. File Structure: `wave_equation.hpp`

The header file defines the `WaveEquation` class template, which encapsulates the entire simulation state and logic.

### Key Components

* **Enums and Types**: Defines `TimeType` (`IMPLICIT` or `EXPLICIT`) to determine the integration scheme at runtime.
* **Class Declaration**: `template <int dim> class WaveEquation` supports both 2D and 3D simulations.
* **Parallel Infrastructure**:
* `mpi_communicator`: Manages MPI ranks.
* `pcout`: A conditional stream that only allows Rank 0 to print, preventing terminal clutter.
* `triangulation`: A `parallel::distributed::Triangulation` for domain partitioning across nodes.

* **Finite Element & DoF Management**:
* `fe`: Defines the element type (e.g.,  or ).
* `dof_handler`: Manages the global enumeration of degrees of freedom.
* `IndexSet`: Distinguishes between **Locally Owned** and **Locally Relevant** (ghost) DoFs.

* **Linear Algebra**:
* `TrilinosWrappers`: Parallel sparse matrices (, , ) and vectors (, , ).

* **Function Pointers (Strategy Pattern)**:
* `std::shared_ptr<Function<dim>>` pointers for initial conditions, forcing terms, and boundary conditions, allowing the solver to switch scenarios without code changes.

---

## 2. File Structure: `wave_equation.cpp`

The implementation file contains the logic for the simulation lifecycle, from grid generation to the final time-stepping loop.

### Method Workflow

#### A. Setup and Assembly

1. **Constructor**: Parses parameters via `InputParser` and initializes the chosen physical scenario.
2. **`make_grid()`**: Generates a hypercube and performs global refinement.
3. **`setup_system()`**: Distributes DoFs and initializes parallel matrices and vectors based on the sparsity pattern.
4. **`assemble_matrices()`**: Computes the static Mass () and Stiffness () matrices using numerical quadrature.

#### B. Initial Conditions

* **`solve_initial_conditions()`**: Computes the acceleration at  by solving . This ensures the simulation starts without numerical shocks.

#### C. Solvers (The Time Loop)

* **`solve_IMPLICIT()`**: Executes the Newmark-beta scheme. It includes:
* Predictor step.
* Solving the linear system  using **CG + AMG**.
* Corrector step for velocity and acceleration.

* **`solve_EXPLICIT()`**: Executes the Velocity Verlet scheme. It features:
* **Mass Lumping**: Diagonalizing  for  efficiency.
* **CFL Check**: Automatically calculating a stable .
* Two-stage updates for velocity and displacement.

#### D. Post-Processing and I/O

* **`output_results()`**: Synchronizes ghost vectors and writes parallel `.pvtu` files for ParaView.
* **`compute_energy()`**: Calculates the sum of kinetic () and potential () energy at each step for conservation analysis.
* **`compute_error_L2()`**: For verification, compares the numerical solution against the exact analytical function.

---

## Summary of Execution Flow

1. **Run**:

* `make_grid()`
* `setup_system()`
* `assemble_matrices()`
* `solve_initial_conditions()`
* **Loop**: `solve_EXPLICIT()` or `solve_IMPLICIT()`
* `output_results()` (at defined frequencies)
