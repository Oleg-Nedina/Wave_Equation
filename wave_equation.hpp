#define WAVE_EQUATION_HPP

// --- Include Standard ---
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

// --- Include Grid e DoF ---
#include <deal.II/distributed/tria.h> 
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h> // STANDARD:FE_Q per quadrilateri/esaedri
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

namespace WaveEquationProject {

    using namespace dealii;

    enum TimeType { IMPLICITO, ESPLICITO };

    template <int dim>
    class WaveEquation {
    public:
        WaveEquation(const std::string &param_file);
        void run();

    private:
        // Setup 
        void setup_system();
        void make_grid();

        // Assemblaggio
        void assemble_matrices();         
        void assemble_rhs_IMPLICITO(); 
        void assemble_rhs_ESPLICITO();

        // Solver
        void solve_IMPLICITO(); 
        void solve_ESPLICITO(); 

        // Output
        void output_results(const unsigned int time_step);

        

        MPI_Comm mpi_communicator;
        
        ConditionalOStream pcout; 
        TimerOutput computing_timer;

        parallel::distributed::Triangulation<dim> triangulation;
        FE_Q<dim> fe; 
        DoFHandler<dim> dof_handler;

        IndexSet locally_owned_dofs;
        IndexSet locally_relevant_dofs;

        AffineConstraints<double> constraints;

        TrilinosWrappers::MPI::Vector U;      
        TrilinosWrappers::MPI::Vector U_old; 
        TrilinosWrappers::MPI::Vector V;      
        TrilinosWrappers::MPI::Vector A;     
        
        TrilinosWrappers::MPI::Vector system_rhs; 
        TrilinosWrappers::MPI::Vector lumped_mass_matrix;

        TrilinosWrappers::SparseMatrix system_matrix;
        TrilinosWrappers::SparseMatrix mass_matrix;   
        TrilinosWrappers::SparseMatrix laplace_matrix;

        double time;
        double time_step;
        double final_time;
        unsigned int output_time_step;

        TimeType method_type;
        double beta;
        double gamma;
    };
}
#endif 
