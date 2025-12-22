#include "wave_equation.hpp"

// Main function
int main(int argc, char *argv[])
{
    try
    {
        using namespace dealii;
        using namespace WaveEquationProject;

        // Inizializza MPI
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

        // Controllo argomenti:
        // Se l'utente passa un file, lo usa. Altrimenti cerca "parameters.prm" di default.
        std::string parameter_file;
        if (argc > 1) {
            parameter_file = argv[1];
        } else {
            parameter_file = "parameters.prm";
        }

        // Lancia la simulazione in 2D passando il nome del file
        WaveEquation<2> wave_equation_solver(parameter_file);
        wave_equation_solver.run();
    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl << std::endl
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
        std::cerr << std::endl << std::endl
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
