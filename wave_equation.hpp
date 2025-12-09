#ifndef WAVE_EQUATION_HPP
#define WAVE_EQUATION_HPP

// ============================================================================
// INCLUDE FILES
// ============================================================================

// --- Standard deal.II Utilities ---
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/parameter_handler.h> 

// --- Grid and DoF Handling ---
#include <deal.II/distributed/tria.h>      // Distributed mesh for MPI
#include <deal.II/grid/grid_generator.h>   // For hyper_cube
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

// --- Finite Elements ---
#include <deal.II/fe/fe_q.h>       // Lagrangian Qp elements
#include <deal.II/fe/fe_values.h>  // Evaluation of basis functions

// --- Linear Algebra (Trilinos Wrappers) ---
// We use Trilinos for parallel sparse matrices and vectors.
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/affine_constraints.h> // Handles hanging nodes and BCs

// --- Numerics and Output ---
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

// --- C++ Standards ---
#include <memory>
#include <string>
#include <fstream>
#include <iostream>

namespace WaveEquationProject {

    using namespace dealii;

    /**
     * Enumerator to select the time integration scheme.
     * - IMPLICIT: Unconditionally stable, requires linear solver (Newmark beta > 0).
     * - EXPLICIT: Conditionally stable (CFL), fast diagonal solve (Newmark beta = 0).
     */
    enum TimeType { IMPLICIT, EXPLICIT };

    /**
     * Main class for the Parallel Wave Equation Solver.
     * Solves the scalar wave equation d2u/dt2 - Laplacian(u) = f
     * using the Finite Element Method (FEM) in space and Newmark-beta in time.
     *
     * @tparam dim Spatial dimension (2 or 3).
     */
    template <int dim>
    class WaveEquation {
    public:
        /**
         * Constructor.
         * @param param_file Name of the parameter file (currently unused, parameters are hardcoded for testing).
         */
        WaveEquation(const std::string &param_file);

        /**
         * Main driver function.
         * Orchestrates the simulation: grid generation, setup, assembly,
         * initial conditions, time loop, and output.
         */
        void run();

    private:
        // SETUP & MESH GENERATION
        
        /**
         * Initializes the linear algebra structures (matrices, vectors)
         * and distributes Degrees of Freedom (DoFs) across MPI processes.
         */
        void setup_system();

        /**
         * Generates the mesh (HyperCube) and performs global refinement.
         */
        void make_grid();

        // ====================================================================
        // ASSEMBLY ROUTINES
        // ====================================================================

        /**
         * Assembles the static matrices M (Mass) and K (Stiffness/Laplace).
         * Since the mesh is static, these are computed only once at t=0.
         */
        void assemble_matrices();         

        /**
         * Builds the specific system matrix for the Implicit Newmark scheme.
         * A_sys = M + beta * dt^2 * K
         */
        void assemble_matrix_IMPLICITO(); 

        /**
         * (Placeholder) Builds the lumped mass matrix for the Explicit scheme.
         */
        void assemble_matrix_ESPLICITO();

        // SOLVERS & ALGORITHMS

        /**
         * Computes consistent initial acceleration (A0).
         * Solves M * A0 = F(0) - K * U0.
         * Necessary because standard initial data only provides U0 and V0.
         */
        void solve_initial_conditions();

        /**
         * Computes the L2 error norm against the analytical solution.
         * Used for convergence analysis (Scenario 1).
         * @return The global L2 error across all processors.
         */
        double compute_error_L2();

        /**
         * Executes the Implicit Time Stepping Loop (Newmark-beta).
         * Includes prediction, linear solve (CG + AMG), and correction.
         */
        void solve_IMPLICITO(); 

        /**
         * (Placeholder) Executes the Explicit Time Stepping Loop.
         * Uses Mass Lumping for fast inversion.
         */
        void solve_ESPLICITO(); 

        /**
         * Outputs the solution in .vtu/.pvtu format for ParaView.
         * @param time_step The current time step number (used for filename generation).
         */
        void output_results(const unsigned int time_step);

        void compute_lumped_mass_matrix();
        void auto_check_cfl_condition();
        void solve_EXPLICIT();

        // MEMBER VARIABLES

        // --- MPI & Parallelism ---
        MPI_Comm mpi_communicator;
        ConditionalOStream pcout;       // stream that prints only on process 0
        TimerOutput computing_timer;    // Profiling tool

        // --- FEM & Grid ---
        parallel::distributed::Triangulation<dim> triangulation; // Distributed mesh
        FE_Q<dim> fe;                                            // Q1 elements (linear)
        DoFHandler<dim> dof_handler;

        // --- Index Sets (Parallel partitioning) ---
        IndexSet locally_owned_dofs;    // DoFs strictly owned by this processor
        IndexSet locally_relevant_dofs; // Owned + Ghost DoFs needed for visualization/gradients

        // --- Constraints ---
        AffineConstraints<double> constraints; // Handles hanging nodes and Dirichlet BCs

        // --- Vectors (Solution state) ---
        // U: Displacement, V: Velocity, A: Acceleration
        TrilinosWrappers::MPI::Vector U;      
        TrilinosWrappers::MPI::Vector U_old; // Previous step displacement
        TrilinosWrappers::MPI::Vector V;      
        TrilinosWrappers::MPI::Vector A;     
        
        // --- Vectors (Linear Algebra) ---
        TrilinosWrappers::MPI::Vector system_rhs; 
        TrilinosWrappers::MPI::Vector lumped_mass_matrix; // Diagonal approximation of M

        // --- Matrices ---
        // system_matrix: The Jacobian matrix to invert at each step
        TrilinosWrappers::SparseMatrix system_matrix;
        // mass_matrix: Standard mass matrix (M)
        TrilinosWrappers::SparseMatrix mass_matrix;   
        // laplace_matrix: Stiffness matrix (K)
        TrilinosWrappers::SparseMatrix laplace_matrix;

        // --- Physical Functions (Strategy Pattern) ---
        // These pointers allow switching between different scenarios (Gaussian, Rain, Ring, etc.)
        // without changing the core solver logic.
        std::shared_ptr<Function<dim>> initial_forcing_term; // f(x,t)
        std::shared_ptr<Function<dim>> initial_u0;           // u(x,0)
        std::shared_ptr<Function<dim>> initial_v0;           // v(x,0)
        std::shared_ptr<Function<dim>> boundary_function;    // g(t) - Dirichlet values
        std::shared_ptr<Function<dim>> exact_solution;       // u_ex (for verification)

        // --- Parameters ---
        unsigned int scenario_id; 
        double time;
        double time_step;
        double final_time;
        unsigned int output_time_step;

        // Newmark-beta parameters
        TimeType method_type;
        double beta;   
        double gamma; 
    };
}
#endif
