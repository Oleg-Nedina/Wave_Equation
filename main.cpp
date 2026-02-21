#include "wave_equation.hpp"

// Main function
int main(int argc, char *argv[]) {
  try {
    using namespace dealii;
    using namespace WaveEquationProject;

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    std::string parameter_file;
    if (argc > 1) {
      parameter_file = argv[1];
    } else {
      parameter_file = "parameters.prm";
    }

    WaveEquation<2> wave_equation_solver(parameter_file);
    wave_equation_solver.run();
  } catch (std::exception &exc) {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  } catch (...) {
    std::cerr << std::endl
              << std::endl
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
