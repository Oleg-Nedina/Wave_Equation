 //esempio di main.cpp


#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/parameter_handler.h>

#include "wave_equation.hpp"

using namespace dealii;

int main(int argc, char *argv[])
{
    // Inizializzazione MPI
    // Utilities::MPI::MPI_InitFinalize gestisce automaticamente init e finalize , sto comando in teoria imposta i tbb al minimo per evitare conflitti
       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    try
    {
        //bisogna definire un file di parametri in modo tale da non dover ricompilare ogni volta , limportante è il formato con cui decidiamo di passare i file  , in modo tale che il programma sappia leggerli (altra opzione è passare i parametri da riga di comando ma è scomodo) 
         //ho riportato una linea di codice che permette di gestire enrambi i casi
        std::string parameter_filename = "parameters.prm";
        if (argc > 1)
        {
            parameter_filename = argv[1];
        }

        // usiamo direttamente il 2 perchè il problema è 2D come richiesto , possiamo pero fare anche con 3d dato che la gestione dovrebbe divenrarea automatica grazie a dealii (o quasi)
        WaveEquationProject::WaveEquation<2> wave_problem(parameter_filename);
        
        wave_problem.run();
    }
     //stile java nic , catch per la gestione degli errori durante lesecuzione del programma
    catch (std::exception &exc)
    {
        // Gestione errori standard C++ (ho cercato i comandi per fare un log "carino")
        std::cerr << "Exception: " << std::endl
                  << exc.what() << std::endl
                  << "Esploso!" << std::endl
                  << std::endl;
        return 1;
    }
// Gestione errori non standard
    catch (...)
    {
        std::cerr << "Unknown exception!" << std::endl;
        return 1;
    }
    return 0;
}
