
#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/parameter_handler.h>

#include "wave_equation.hpp"

using namespace dealii;

int main(int argc, char *argv[])
{
    // ========================================================================
    // MPI INITIALIZATION
    // ========================================================================
    // Utilities::MPI::MPI_InitFinalize automatically handles MPI_Init and 
    // MPI_Finalize via RAII. The third argument (1) restricts the number of 
    // Threading Building Blocks (TBB) threads to minimize conflicts with MPI.
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    try
    {
        // ====================================================================
        // ARGUMENT PARSING
        // ====================================================================
        // Usage: ./wave_solver [dimension] [parameter_file]
        // Default: dimension = 2, file = "parameters.prm"
        
        unsigned int dimension = 2; // Default dimension
        std::string parameter_filename = "parameters.prm"; // Default filename

        // 1. Read Dimension (if provided)
        if (argc > 1)
        {
            // Convert string argument to integer
            dimension = std::stoi(argv[1]);
        }

        // 2. Read Parameter File (if provided)
        if (argc > 2)
        {
            parameter_filename = argv[2];
        }

        // ====================================================================
        // PROBLEM INSTANTIATION
        // ====================================================================
        // Since templates are compile-time constructs, we must switch on the 
        // runtime value of 'dimension' to instantiate the correct class.

        if (dimension == 2)
        {
            // Instantiate and run the 2D solver
            WaveEquationProject::WaveEquation<2> wave_problem(parameter_filename);
            wave_problem.run();
        }
        else if (dimension == 3)
        {
            // Instantiate and run the 3D solver
            WaveEquationProject::WaveEquation<3> wave_problem(parameter_filename);
            wave_problem.run();
        }
        else
        {
            // Error handling for invalid dimensions
            std::cerr << "Error: Dimension " << dimension 
                      << " is not supported. Use 2 or 3." << std::endl;
            return 1;
        }
    }
    // ========================================================================
    // EXCEPTION HANDLING
    // ========================================================================
    catch (std::exception &exc)
    {
        // Catch standard C++ exceptions and deal.II exceptions
        std::cerr << std::endl 
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Exception on processing: " << std::endl
                  << exc.what() << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    }
    catch (...)
    {
        // Catch unknown exceptions
        std::cerr << std::endl 
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Unknown exception!" << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    }

    return 0;
}
