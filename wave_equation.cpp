#include "./wave_equation.hpp"
#include "test_cases.hpp" 

// Necessary includes for grid generation, output, and file system operations
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <filesystem>

namespace WaveEquationProject {



    // WAVE EQUATION CLASS IMPLEMENTATION

    /**
     * @brief Constructor.
     * Initializes MPI, Timers, Triangulation, and Finite Elements.
     * Selects the test scenario and sets up initial/boundary conditions.
     */
    template <int dim>
    WaveEquation<dim>::WaveEquation(const std::string &param_file)
        : mpi_communicator(MPI_COMM_WORLD),
        
        // pcout: Wraps std::cout so that only process 0 prints to the terminal.
        // This prevents 100 processes from printing the same line simultaneously.
          pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
         

        // computing_timer: Measures wall-clock time for each section (Assemble, Solve, Output).
          computing_timer(mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times),
        // triangulation: Distributed parallel mesh.
        // MeshSmoothing ensures that when refining, adjacent cells don't differ too much in size.
          triangulation(mpi_communicator,
                        typename Triangulation<dim>::MeshSmoothing(
                            Triangulation<dim>::smoothing_on_refinement |
                            Triangulation<dim>::smoothing_on_coarsening)),
         
        // fe: Lagrangian Q1 elements (linear basis functions).
          fe(1), 
        // dof_handler: Manages the enumeration of Degrees of Freedom on the mesh.
          dof_handler(triangulation) 
    {
        // --- PARAMETER HANDLING ---
        // Ideally, we would read parameters from the .prm file here.
        // WE WILL IMPLEMENT THIS LATER. FOR NOW, WE SET DEFAULTS DIRECTLY.
        ParameterHandler prm;
        prm.enter_subsection("Global parameters");

        // Suppress "unused parameter" warning (we will use this later).
        (void)param_file;


        // --- SIMULATION SETTINGS ---
        // dt: Time step size. Critical for stability (Explicit) and accuracy (Implicit).
        time_step = 0.01;


        // T: Final simulation time.
        final_time = 1.0;

        // Output freq: 1 = save every step (smooth video), 10 = save every 10th step (save disk space).
        output_time_step = 1;
        prm.leave_subsection();
        

        // --- SCENARIO SELECTION ---
        // Change this ID to switch between tests (we will later read from prm file) 
        // 1: Verification (Sinusoidal)
        // 2: Gaussian Pulse
        // 3: Collision
        // 4: Raindrops
        // 5: Ring Implosion
       // 6: Time-Dependent BCs (Pumping Wall)
        scenario_id = 1; //it will provide a convergence test



        // Strategy Pattern: We assign the correct Function object to the shared pointers.
        // The rest of the code is agnostic to which specific function is being used.
        if (scenario_id == 1) {
            pcout << "Scenario 1: Manufactured Solution (Convergence Test)" << std::endl;
            initial_u0 = std::make_shared<InitialValuesU_Case1<dim>>();
            initial_v0 = std::make_shared<InitialValuesV_Case1<dim>>();
            initial_forcing_term = std::make_shared<RightHandSide_Case1<dim>>();
            boundary_function = std::make_shared<Functions::ZeroFunction<dim>>(); // Homogeneous Dirichlet
            exact_solution = std::make_shared<ExactSolution_Case1<dim>>();
        } 
        else if (scenario_id == 2) {
            pcout << "Scenario 2: Gaussian Pulse Propagation" << std::endl;
            initial_u0 = std::make_shared<InitialValuesU_Case2<dim>>();
            initial_v0 = std::make_shared<InitialValuesV_Case2<dim>>();
            initial_forcing_term = std::make_shared<RightHandSide_Case2<dim>>(); 
            boundary_function = std::make_shared<Functions::ZeroFunction<dim>>(); 
        }
        else if (scenario_id == 3) {
            pcout << "Scenario 3: Double Gaussian Collision" << std::endl;
            initial_u0 = std::make_shared<InitialValuesU_Case3<dim>>();
            initial_v0 = std::make_shared<InitialValuesV_Case2<dim>>(); //reuse V=0
            initial_forcing_term = std::make_shared<RightHandSide_Case2<dim>>();//reuse F=0
            boundary_function = std::make_shared<Functions::ZeroFunction<dim>>(); 
        }
        else if (scenario_id == 4) {
            pcout << "Scenario 4: The Raindrop Symphony" << std::endl;
            initial_u0 = std::make_shared<InitialValuesU_Raindrops<dim>>();
            
            initial_v0 = std::make_shared<InitialValuesV_Case2<dim>>(); 
            initial_forcing_term = std::make_shared<RightHandSide_Case2<dim>>();
            
            boundary_function = std::make_shared<Functions::ZeroFunction<dim>>(); 
        }
        else if (scenario_id == 5) {
            pcout << "Scenario 5: Ring Implosion" << std::endl;
            initial_u0 = std::make_shared<InitialValuesU_Ring<dim>>();
            initial_v0 = std::make_shared<InitialValuesV_Case2<dim>>(); // V=0
            initial_forcing_term = std::make_shared<RightHandSide_Case2<dim>>(); // F=0
            boundary_function = std::make_shared<Functions::ZeroFunction<dim>>(); 
        }
        else if (scenario_id == 6) {
            pcout << "Scenario 6: Time-Dependent BCs (Pumping Wall)" << std::endl;
            
            initial_u0 = std::make_shared<Functions::ZeroFunction<dim>>();
            initial_v0 = std::make_shared<Functions::ZeroFunction<dim>>(); 
            initial_forcing_term = std::make_shared<Functions::ZeroFunction<dim>>();
            
            boundary_function = std::make_shared<BoundaryFunction_MovingWall<dim>>(); 
        }
        else {
            throw std::runtime_error("Invalid scenario ID specified.");
    }

        //--- NEWMARK SCHEME PARAMETERS ---
        time = 0.0;
        method_type = EXPLICIT;

        // Beta = 0.25, Gamma = 0.5 corresponds to the "Constant Average Acceleration" method.
        // It is unconditionally stable (A-stable) and second-order accurate.
        beta = 0.25;
        gamma = 0.5;
    }


    // GRID GENERATION

    /**
     * @brief Generates the computational domain.
     * We use a simple HyperCube (square in 2D) and refine it globally.
     * In a real application, this would read a .msh file.
     */
    template <int dim>
    void WaveEquation<dim>::make_grid() {

        TimerOutput::Scope t(computing_timer, "Make grid");
        pcout << "Generating Grid..." << std::endl;

        // Create a unit square [0,1]^d
        GridGenerator::hyper_cube(triangulation, 0, 1);



        int tmp_int_ref;
        if(dim==2){
            tmp_int_ref=6;
        }
        else{
            tmp_int_ref=4;
        }
        // Refine the mesh globally 'n' times. 
        // 6 refinements = 2^6 = 4096 cells in 2D. (Adjust as needed for convergence tests the numer of refinements)
        triangulation.refine_global(tmp_int_ref); 

        pcout << "  Number of active cells: " << triangulation.n_global_active_cells() << std::endl;
    }


    // ERROR ANALYSIS

    /**
     * @brief Computes the L2 error norm against the exact solution.
     * E = || u_h - u_ex ||_L2
     * Used to verify the order of convergence (should be 2 for Q1 elements).
     */
    template <int dim>
     double WaveEquation<dim>::compute_error_L2() {
        pcout << "  Computing L2 Error against Exact Solution..." << std::endl;

        // Synchronize the exact solution function to the current simulation time
        exact_solution->set_time(time);

        // Vector to store the error norm squared for each cell
        Vector<float> difference_per_cell(triangulation.n_active_cells());

        // Numerical integration of the difference (u_h - u_ex)^2
        VectorTools::integrate_difference(dof_handler,
                                          U,
                                          *exact_solution,
                                          difference_per_cell,
                                          QGauss<dim>(fe.degree + 2), // Higher order quadrature for error
                                          VectorTools::L2_norm);

        // Sum the local errors squared
        const double local_error_sq = difference_per_cell.norm_sqr();


        // MPI Reduction: Sum errors from all processors to get the global error            
        const double global_error_sq = Utilities::MPI::sum(local_error_sq, mpi_communicator);

        return std::sqrt(global_error_sq);
    }


    // SYSTEM SETUP

    /**
     * @brief Initializes DoFs, Matrices, and Vectors.
     * Sets up the parallel sparsity patterns and resizes all linear algebra objects.
     */
    template <int dim>
    void WaveEquation<dim>::setup_system() {
        TimerOutput::Scope t(computing_timer, "Setup system");
        pcout << "Setting up system..." << std::endl;

        // Distribute Degrees of Freedom based on the FE selection
        dof_handler.distribute_dofs(fe);

        // --- MPI PARTITIONING ---
        // locally_owned: DoFs this processor is responsible for computing.
        locally_owned_dofs = dof_handler.locally_owned_dofs();

        // locally_relevant: Owned DoFs + Ghost DoFs on adjacent cells.
        // Required for computing gradients and writing output.
        locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

        pcout << "  Number of DoFs: " << dof_handler.n_dofs() << std::endl;

        // --- CONSTRAINTS ---
        // Handling hanging nodes (if any) and periodic boundaries.
        // Note: Dirichlet BCs are handled dynamically in the time loop or via AffineConstraints here.
        constraints.clear();
        constraints.reinit(locally_relevant_dofs);
        DoFTools::make_hanging_node_constraints(dof_handler, constraints);
        constraints.close();

        // --- MATRIX INITIALIZATION ---
        // Build the sparsity pattern using Trilinos (parallel efficient).
        TrilinosWrappers::SparsityPattern dsp(locally_owned_dofs, mpi_communicator);
        DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
        dsp.compress();

        // Initialize system matrices with the sparsity pattern
        mass_matrix.reinit(dsp);
        laplace_matrix.reinit(dsp);
        system_matrix.reinit(dsp); // Used for the implicit solver (Jacobian)

        // --- VECTOR INITIALIZATION ---
        // Solution vectors need 'locally_relevant_dofs' to store ghosts.
        U.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
        U_old.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
        V.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
        A.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

        // RHS and temporary vectors only need 'locally_owned_dofs' for algebra.
        system_rhs.reinit(locally_owned_dofs, mpi_communicator);
        lumped_mass_matrix.reinit(locally_owned_dofs, mpi_communicator);
    }



    // MATRIX ASSEMBLY

    /**
     * @brief Assembles the time-independent Mass (M) and Stiffness (K) matrices.
     * M_ij = \int \phi_i \phi_j dx
     * K_ij = \int \nabla\phi_i \cdot \nabla\phi_j dx
     */
    template <int dim>
    void WaveEquation<dim>::assemble_matrices() {

        TimerOutput::Scope t(computing_timer, "Assemble matrices");
        pcout << "Assembling Mass and Laplace matrices..." << std::endl;

        mass_matrix = 0;
        laplace_matrix = 0;

        // Quadrature rule: Gauss formula exact for polynomial degree 2*p + 1
        const QGauss<dim> quadrature_formula(fe.degree + 1);

        // FEValues extracts shape function values and gradients at quadrature point
        FEValues<dim> fe_values(fe, quadrature_formula,
                                update_values | update_gradients | 
                                update_JxW_values | update_quadrature_points);

        const unsigned int dofs_per_cell = fe.dofs_per_cell;
        const unsigned int n_q_points = quadrature_formula.size();

        FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
        FullMatrix<double> cell_laplace_matrix(dofs_per_cell, dofs_per_cell);
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

        // Loop over all cells
        for (const auto &cell : dof_handler.active_cell_iterators()) {
           
            // MPI Check: Assemble only cells owned by this processor
            if (cell->is_locally_owned()) { 
                fe_values.reinit(cell);
                cell_mass_matrix = 0;
                cell_laplace_matrix = 0;

                // Loop over quadrature points
                for (unsigned int q = 0; q < n_q_points; ++q) {

                    // Loop over test functions (i) and trial functions (j)
                    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                        for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                            
                            // Mass Matrix: (phi_i, phi_j)
                            cell_mass_matrix(i, j) += 
                                (fe_values.shape_value(i, q) *
                                 fe_values.shape_value(j, q) *
                                 fe_values.JxW(q));

                            // Laplace (Stiffness) Matrix: (grad_phi_i, grad_phi_j)
                            // Note: No minus sign here, it is added in the weak form equation
                            cell_laplace_matrix(i, j) += 
                                (fe_values.shape_grad(i, q) *
                                 fe_values.shape_grad(j, q) *
                                 fe_values.JxW(q));
                        }
                    }
                }

                // Map local DoF indices to global indices
                cell->get_dof_indices(local_dof_indices);
                
                // Add local contributions to the global matrices (thread-safe for MPI/Constraints)
                constraints.distribute_local_to_global(cell_mass_matrix, local_dof_indices, mass_matrix);
                constraints.distribute_local_to_global(cell_laplace_matrix, local_dof_indices, laplace_matrix);
            }
        }

        // Finalize assembly: communicate data between processors
        mass_matrix.compress(VectorOperation::add);
        laplace_matrix.compress(VectorOperation::add);
    }


    // INITIAL CONDITIONS (CONSISTENT ACCELERATION)

    /**
     * @brief Computes the initial acceleration U''(0) consistent with the strong form.
     * * In second-order PDE solvers, we are given u(0) and u'(0). 
     * However, the Newmark scheme also requires the initial acceleration u''(0).
     * We cannot just assume u''(0)=0 (unless the system is at rest with no forces).
     * * Algorithm 1 (from Math PDF):
     * 1. Interpolate initial displacement U0 and velocity V0.
     * 2. Construct the RHS at t=0: b = F(0) - K * U0.
     * 3. Solve the linear system: M * A0 = b.
     */
template <int dim>

    void WaveEquation<dim>::solve_initial_conditions() {
        TimerOutput::Scope t(computing_timer, "Initial Conditions");
        pcout << "Computing Initial Conditions..." << std::endl;

        // --- 1. Interpolate U0 and V0 ---
        // Project analytical functions onto the FEM space
        VectorTools::interpolate(dof_handler, *initial_u0, U);
        VectorTools::interpolate(dof_handler, *initial_v0, V);
        
       // Exchange ghost data to ensure consistency across processors
        U.compress(VectorOperation::insert);
        V.compress(VectorOperation::insert);
        
        // Store U0 for the first time step logic
        U_old = U; 
 
        // --- 2. Construct RHS: b = F(0) - K * U0 ---
        
        // Prepare a vector for the RHS (non-ghosted, fully distributed)
        TrilinosWrappers::MPI::Vector rhs_initial(locally_owned_dofs, mpi_communicator);
        
        // Compute F(0): The forcing term vector
        initial_forcing_term->set_time(0.0);
        VectorTools::create_right_hand_side(dof_handler, 
                                            QGauss<dim>(fe.degree + 1), 
                                            *initial_forcing_term, 
                                            rhs_initial);
        
        // Compute K * U0
        TrilinosWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator);
        
        // Create a non-ghosted copy of U for matrix-vector multiplication validity
        TrilinosWrappers::MPI::Vector U_owned(locally_owned_dofs, mpi_communicator);
        U_owned = U; 
        
        // Perform sparse matrix-vector multiplication: tmp = K * U
        laplace_matrix.vmult(tmp, U_owned);

        // Subtract K*U from F: rhs_initial = F - K*U
        rhs_initial.add(-1.0, tmp); 

        // --- 3. Solve M * A0 = RHS ---
        
        // Setup Solver Control: Max 1000 iterations, tolerance relative to RHS norm
        SolverControl solver_control(1000, 1e-12 * rhs_initial.l2_norm());
        SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
        
        // Identity preconditioner is sufficient because Mass Matrix is well-conditioned
        TrilinosWrappers::PreconditionIdentity preconditioner; 
        
        pcout << "  Solving for initial acceleration..." << std::endl;

        // Use a non-ghosted vector for the solve operation
        TrilinosWrappers::MPI::Vector A_owned(locally_owned_dofs, mpi_communicator);
        A_owned = A; 
        
        // Solve the system
        solver.solve(mass_matrix, A_owned, rhs_initial, preconditioner);
        
        // Update the ghosted acceleration vector
        A = A_owned; 
        
        pcout << "  Initial acceleration computed. Iterations: " << solver_control.last_step() << std::endl;
    }


    // IMPLICIT SOLVER HELPERS

    /**
     * @brief Assembles the system matrix for the Implicit Newmark-beta scheme.
     * * The discrete equation to solve for U_{n+1} is:
     * (M + beta * dt^2 * K) * U_{n+1} = RHS
     * * This function constructs the matrix on the left hand side:
     * A_sys = scaling * M + K
     * where scaling = 1 / (beta * dt^2)
     * * We modify the equation slightly by dividing everything by (beta * dt^2) 
     * to improve conditioning or match standard implementations:
     * A_sys_actual = (1/(beta*dt^2)) * M + K
     */
    template <int dim>
    void WaveEquation<dim>::assemble_matrix_IMPLICITO() {
        pcout << "  Building Implicit System Matrix..." << std::endl;
        
        // A_sys = scaling * M + K
        const double scaling_factor = 1.0 / (beta * time_step * time_step);

        // Copy Mass matrix structure and values to System matrix
        system_matrix.copy_from(mass_matrix);

        // Scale Mass matrix: M <- M * (1/beta*dt^2)
        system_matrix *= scaling_factor;
       
        // Add Stiffness matrix: M_scaled + K
        system_matrix.add(1.0, laplace_matrix);
    }




    // IMPLICIT SOLVER (NEWMARK-BETA)

    /**
     * @brief Executes the time loop using the Implicit Newmark scheme.
     * * Solves the system: (M + beta*dt^2*K) * U_{n+1} = RHS
     * where RHS includes F_{n+1} and inertial terms from the previous step.
     * * This method is unconditionally stable (A-stable) for beta=0.25, gamma=0.5.
     */
    template <int dim>
    void WaveEquation<dim>::solve_IMPLICITO() {
        pcout << "Starting Implicit Time Loop (Newmark-Beta)..." << std::endl;

        // 1. Build the system matrix (constant if dt is constant)
        // Matrix structure: A_sys = (1/(beta*dt^2)) * M + K_ij
        assemble_matrix_IMPLICITO();

        // --- SOLVER CONFIGURATION ---
        // We use a Conjugate Gradient (CG) solver because the system matrix is 
        // Symmetric Positive Definite (SPD).
        SolverControl solver_control(5000 /*maxiter*/ , 1e-10 /*use 1e-6 if the problem is the 5-ith one becouse we have an "explosion" and if possible -10 */);
        SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);

        // --- PRECONDITIONER ---
        // Algebraic Multigrid (AMG) is highly effective for elliptic-like operators
        // such as the one we are solving here (Mass + Stiffness).
        TrilinosWrappers::PreconditionAMG preconditioner;
        preconditioner.initialize(system_matrix);

        // --- TEMPORARY VECTORS ---
        // "Owned" vectors do not store ghost elements. They are used for 
        // algebraic operations (Solver, Matrix-Vector products) to ensure MPI consistency.
        TrilinosWrappers::MPI::Vector U_owned(locally_owned_dofs, mpi_communicator);
        TrilinosWrappers::MPI::Vector RHS_owned(locally_owned_dofs, mpi_communicator);
        TrilinosWrappers::MPI::Vector TMP_owned(locally_owned_dofs, mpi_communicator);

        // Temporary vector for acceleration update (needs ghosts for later use in output)
        TrilinosWrappers::MPI::Vector A_new(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

        unsigned int time_step_number = 0;


        // TIME LOOP
        while (time < final_time) {
            time += time_step;
            time_step_number++;


            //debug print 1:
           // pcout << "  Step " << time_step_number << " (t=" << time << ")" << std::flush;

            // 1. CONSTRUCT RHS
            // RHS = F_{n+1} + M * (Predictor Terms)
            
            // a. Compute Forcing Term F_{n+1}
            initial_forcing_term->set_time(time);
            VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(fe.degree + 1), 
                                                *initial_forcing_term, RHS_owned);

            // b. Compute Newmark Predictor Coefficients
            // Derived from rearranging the Newmark displacement formula to solve for A_{n+1}.
            // See Eq (4.4) in Project Documentation[cite: 317].
            double c1 = 1.0 / (beta * time_step * time_step);
            double c2 = 1.0 / (beta * time_step);
            double c3 = (1.0 - 2.0 * beta) / (2.0 * beta);

            // c. Accumulate Predictors in TMP_owned
            // TMP = c1*U_n + c2*V_n + c3*A_n
            U_owned = U; // Copy previous state (U_n)
            TMP_owned = U_owned; 
            TMP_owned *= c1;
            
            U_owned = V; 
            TMP_owned.add(c2, U_owned);

            U_owned = A; 
            TMP_owned.add(c3, U_owned);

            // d. Apply Mass Matrix: RHS += M * TMP
            mass_matrix.vmult(U_owned, TMP_owned);
            RHS_owned.add(1.0, U_owned);


            // 2. APPLY BOUNDARY CONDITIONS
            
            // We apply Dirichlet BCs strongly to the linear system.
            // Note: In Trilinos, modifying the matrix after assembly can be expensive,
            // but MatrixTools::apply_boundary_values handles the local rows correctly.
            std::map<types::global_dof_index, double> boundary_values;
            boundary_function->set_time(time);
            VectorTools::interpolate_boundary_values(dof_handler, 0, *boundary_function, boundary_values);
            
            // Apply to System Matrix and RHS:
            // Rows corresponding to boundary DoFs are cleared and set to identity diagonal.
            // RHS is set to the boundary value.
            MatrixTools::apply_boundary_values(boundary_values, system_matrix, U_owned, RHS_owned, false);

            // 3. SOLVE LINEAR SYSTEM

            // Solve: A_sys * U_{n+1} = RHS
            solver.solve(system_matrix, U_owned, RHS_owned, preconditioner);
            
            // Distribute results to ghost layers
            constraints.distribute(U_owned);
            
            // Save U_n into U_old before overwriting U
            U_old = U; 
            U = U_owned; // U now holds U_{n+1}

      
            // 4. CORRECTOR STEP (UPDATE V and A)
            
            // Calculate A_{n+1} using the displacement difference.
            // Formula: A_{n+1} = c1*(U_{n+1} - U_n) - c2*V_n - c3*A_n
            A_new = U;       // U_{n+1}
            A_new -= U_old;  // U_{n+1} - U_n
            A_new *= c1;     
            
            A_new.add(-c1 * time_step, V); // Subtract velocity term
            A_new.add(-(1.0 - 2.0 * beta) / (2.0 * beta), A); // Subtract acceleration term

            // Calculate V_{n+1}
            // Formula: V_{n+1} = V_n + dt * [ (1-gamma)*A_n + gamma*A_{n+1} ]
            V.add(time_step * (1.0 - gamma), A);       
            V.add(time_step * gamma, A_new);           

            // Update Acceleration for next step
            A = A_new;
            
            //debug print 2:
            //pcout << " -> CG Iter: " << solver_control.last_step() << std::endl;

            // Output results at specified intervals
            if (time_step_number % output_time_step == 0) {

            //debug print 3:
            /*pcout << "Step " << time_step_number 
                      << " (t=" << time << ")"
                      << " -> CG Iter: " << solver_control.last_step() 
                      << " [Output Saved]" << std::endl; */

                output_results(time_step_number);
            }
        } // End of Time Loop


        pcout << "Simulation loop finished." << std::endl;

    
        // POST-PROCESSING (VERIFICATION)
        // --------------------------------------------------------------------
        // If running the manufactured solution scenario, check the L2 error
        // to verify the implementation correctness.
        if (scenario_id == 1) {
            double error = compute_error_L2();
            pcout << "========================================" << std::endl;
            pcout << "CONVERGENCE CHECK (t=" << time << ")" << std::endl;
            pcout << "L2 Error: " << error << std::endl;
            pcout << "========================================" << std::endl;
        }
    }

    /**
     * @brief Computes the Lumped Mass Matrix (Diagonal) and its Inverse.
     * Technique: Row-Summing the consistent mass matrix.
     * M_lumped_ii = Sum_j (M_consistent_ij)
     */
    template <int dim>
    void WaveEquation<dim>::compute_lumped_mass_matrix() {
        pcout << "  Computing Lumped Mass Matrix..." << std::endl;

        // 1. Create ones vector
        TrilinosWrappers::MPI::Vector ones(locally_owned_dofs, mpi_communicator);
        ones = 1.0;

        // 2. Perform Matrix-Vector multiplication: M_lumped = M_consistent * ones
        // This effectively sums rows.
        mass_matrix.vmult(lumped_mass_matrix, ones);

        // 3. Compute Inverse Mass Vector (for explicit division)
        // Reusing ones vector or create new one to store 1/M.
        // Note: Must iterate manually because Trilinos vectors don't support element-wise division easily.

        // Store inverse in 'lumped_mass_matrix' directly to save memory
        // For safety, invert in place.

        for (unsigned int i = 0; i < lumped_mass_matrix.local_size(); ++i) {
            if (lumped_mass_matrix[i] > 1e-15) {
                lumped_mass_matrix[i] = 1.0 / lumped_mass_matrix[i];
            } else {
                // Handle potential zero mass (should not happen with standard Q1 elements)
                lumped_mass_matrix[i] = 0.0;
            }
        }

        pcout << "  Mass Lumping complete (stored as Inverse Mass)." << std::endl;
    }

    /**
     * @brief Automatically selects a stable time step based on the mesh size.
     * Sets time_step = 0.9 * h_min / (c * sqrt(dim))
     */
    template <int dim>
    void WaveEquation<dim>::auto_check_cfl_condition() {
        pcout << "  Auto_Checking CFL Condition..." << std::endl;

        // Find local minimum h
        double min_h_local = std::numeric_limits<double>::max();

        for (const auto &cell : dof_handler.active_cell_iterators()) {
            if (cell->is_locally_owned()) {
                min_h_local = std::min(min_h_local, cell->diameter());
            }
        }

        // Find global minimum h
        const double min_h_global = Utilities::MPI::min(min_h_local, mpi_communicator);

        // Define Wave Speed
        const double wave_speed = 1.0;

        // Calculate Critical Time Step
        // CFL Limit: dt <= h / (c * sqrt(d))
        const double dt_critical = min_h_global / (wave_speed * std::sqrt(dim));

        // Apply Safety Factor (0.5) to ensure stability
        this->time_step = dt_critical * 0.5;

        // Report
        pcout << "  ------------------------------------------" << std::endl;
        pcout << "  | AUTOMATIC TIME STEP SELECTION          |" << std::endl;
        pcout << "  ------------------------------------------" << std::endl;
        pcout << "  | Min Element Size (h): " << min_h_global << std::endl;
        pcout << "  | Critical Limit:       " << dt_critical << std::endl;
        pcout << "  | -> NEW Time Step:     " << this->time_step << " (Safety: 90%)" << std::endl;
        pcout << "  ------------------------------------------" << std::endl;
    }

    // MAIN DRIVER

    /**
     * @brief High-level orchestration of the simulation.
     * 1. Creates output directories (handling MPI race conditions).
     * 2. Calls grid generation and system setup.
     * 3. Assembles static matrices.
     * 4. Computes initial conditions.
     * 5. Selects and runs the appropriate time-stepping solver.
     */
   template <int dim>
    void WaveEquation<dim>::run() {


    pcout << "===========================================" << std::endl;
        pcout << "   WAVE EQUATION SOLVER " << std::endl;
        pcout << "===========================================" << std::endl;

        // --- OUTPUT DIRECTORY SETUP ---
        const std::string output_folder = "output_results";



        // Filesystem operations must be done by a single process to avoid race conditions.
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0) {

           // Check if directory exists before creating it
            if (!std::filesystem::exists(output_folder)) {
                std::filesystem::create_directories(output_folder);
            }
        }

        // CRITICAL: All other processes must wait here until Rank 0 has finished
        // creating the directory. If we don't wait, Rank 1 might try to write
        // a file to a non-existent folder, causing a crash.
        MPI_Barrier(mpi_communicator);


        // --- INITIALIZATION PHASE ---

        make_grid(); // Generate geometry
        setup_system(); // Distribute DoFs and init matrices
        assemble_matrices(); // Compute M and K (static)
        if (method_type != IMPLICIT) {auto_check_cfl_condition();} //Automate time-step selection based on mesh size
        solve_initial_conditions(); // Compute A0 consistent with U0, V0

        pcout << "Initial conditions computed. Starting simulation..." << std::endl;

        // Save the state at t=0
        output_results(0);

        // --- SOLVER SELECTION ---
        if (method_type == IMPLICIT) {
            solve_IMPLICITO();
        } else if (method_type == EXPLICIT) {
            solve_EXPLICIT();
        } else {
            // Placeholder for Phase 3 , we will implement the explicit solver later soon (nick be fast please)
            pcout << "Please select IMPLICIT or EXPLICIT as option" << std::endl;
        }

        pcout << "===========================================" << std::endl;
        pcout << "   SIMULATION COMPLETED " << std::endl;
        pcout << "===========================================" << std::endl;

   }

    // OUTPUT ROUTINES

    /**
     * @brief Executes the Explicit Newmark scheme (Central Difference).
     * Parameters: Beta = 0.0, Gamma = 0.5
     * Stability: Conditionally stable (CFL condition required).
     * Cost: Very low per step (No linear solver).
     */
    template <int dim>
    void WaveEquation<dim>::solve_EXPLICIT() {
        pcout << "Starting Explicit Time Loop (Velocity Verlet)..." << std::endl;

        // --- SETUP ---
        compute_lumped_mass_matrix();

        // Initialize Owned Vectors
        TrilinosWrappers::MPI::Vector U_owned(locally_owned_dofs, mpi_communicator);
        TrilinosWrappers::MPI::Vector V_owned(locally_owned_dofs, mpi_communicator);
        TrilinosWrappers::MPI::Vector A_owned(locally_owned_dofs, mpi_communicator);
        TrilinosWrappers::MPI::Vector F_ext(locally_owned_dofs, mpi_communicator);
        TrilinosWrappers::MPI::Vector F_int(locally_owned_dofs, mpi_communicator);

        // Sync initial state (U0 and V0)
        U_owned = U;
        V_owned = V;

        // --- WARM START: RECOMPUTE A0 consistent with LUMPED MASS ---
        // The previous 'A' was computed using the Consistent Mass Matrix.
        // We must recompute it using M_lumped to avoid an initial kick.

        // 1. Compute Forces at t=0
        initial_forcing_term->set_time(time);
        VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(fe.degree + 1),
                                            *initial_forcing_term, F_ext);
        laplace_matrix.vmult(F_int, U); // K * U0

        // 2. Solve M_lumped * A0 = F_ext - F_int
        A_owned = F_ext;
        A_owned.add(-1.0, F_int);
        for (unsigned int i = 0; i < A_owned.local_size(); ++i) {
            A_owned[i] *= lumped_mass_matrix[i];
        }

        // Save consistent A0
        A = A_owned;
        constraints.distribute(A);

        unsigned int time_step_number = 0;

        // --- TIME LOOP ---
        while (time < final_time) {
            time += time_step;
            time_step_number++;

            // 1. UPDATE VELOCITY (First Half-Step)
            // v_{n+0.5} = v_n + 0.5 * dt * a_n
            V_owned.add(0.5 * time_step, A_owned);

            // 2. UPDATE DISPLACEMENT
            // u_{n+1} = u_n + dt * v_{n+0.5}
            U_owned.add(time_step, V_owned);

            // Apply Boundary Conditions to U
            std::map<types::global_dof_index, double> boundary_values;
            boundary_function->set_time(time);
            VectorTools::interpolate_boundary_values(dof_handler, 0, *boundary_function, boundary_values);

            for (auto const& [dof_index, value] : boundary_values) {
                if (locally_owned_dofs.is_element(dof_index)) {
                    U_owned(dof_index) = value;
                }
            }

            // Distribute U to Ghosted Vector
            U = U_owned;
            constraints.distribute(U);

            // 3. COMPUTE ACCELERATION (a_{n+1})
            // M * a = F_ext - K * u

            // External Force
            initial_forcing_term->set_time(time);
            VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(fe.degree + 1),
                                                *initial_forcing_term, F_ext);

            // Internal Force
            laplace_matrix.vmult(F_int, U);

            // Residual
            A_owned = F_ext;
            A_owned.add(-1.0, F_int);

            // Solve Diagonal System
            for (unsigned int i = 0; i < A_owned.local_size(); ++i) {
                A_owned[i] *= lumped_mass_matrix[i];
            }

            // Apply BCs to Acceleration (0 on fixed boundaries)
            for (auto const& [dof_index, value] : boundary_values) {
                if (locally_owned_dofs.is_element(dof_index)) {
                    A_owned(dof_index) = 0.0;
                }
            }
            A = A_owned;
            constraints.distribute(A);

            // 4. UPDATE VELOCITY (Second Half-Step)
            // v_{n+1} = v_{n+0.5} + 0.5 * dt * a_{n+1}
            V_owned.add(0.5 * time_step, A_owned);

            V = V_owned;
            constraints.distribute(V);

            // --- OUTPUT ---
            if (time_step_number % output_time_step == 0) {
                output_results(time_step_number);
            }
        }

        pcout << "Explicit simulation finished." << std::endl;

        // POST-PROCESSING (VERIFICATION)
        if (scenario_id == 1) {
            double error = compute_error_L2();
            pcout << "========================================" << std::endl;
            pcout << "CONVERGENCE CHECK (t=" << time << ")" << std::endl;
            pcout << "L2 Error: " << error << std::endl;
            pcout << "========================================" << std::endl;
        }
    }

    /**
     * @brief Writes the current solution state to disk.
     * Uses the VTK format (Paraview compatible).
     * In parallel, this generates:
     * - One master file: solution_X.pvtu
     * - N piece files:   solution_X_p.vtu (where p is the processor rank)
     * * @param step The current time step number (used for filename numbering).
     */
    template <int dim>
    void WaveEquation<dim>::output_results(const unsigned int step) {
        TimerOutput::Scope t(computing_timer, "Output");

        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler);

        // Attach the solution vectors.
        // Note: These must be ghosted vectors (locally_relevant) to allow
        // DataOut to access values on the boundary of the processor's subdomain.
        data_out.add_data_vector(U, "displacement");
        data_out.add_data_vector(V, "velocity");
        data_out.add_data_vector(A, "acceleration");

        // Attach Subdomain ID (useful to visualize MPI partitioning in Paraview)
        Vector<float> subdomain(triangulation.n_active_cells());
        for (unsigned int i = 0; i < subdomain.size(); ++i)
            subdomain(i) = triangulation.locally_owned_subdomain();
        data_out.add_data_vector(subdomain, "subdomain");

        // Prepare the patches for writing
        data_out.build_patches();

        // Write the parallel file structure.
        // 1st arg: Directory
        // 2nd arg: Base filename (deal.II appends _step.pvtu)
        data_out.write_vtu_with_pvtu_record("output_results/", "solution", step, mpi_communicator);

    }

    // --- TEMPLATE INSTANTIATION ---
    // Forces the compiler to generate code for the 2D case.
    template class WaveEquation<3>;
    template class WaveEquation<2>;

} // namespace WaveEquationProject
