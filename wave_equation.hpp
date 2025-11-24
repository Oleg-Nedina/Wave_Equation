//ho preso ispirazione dal laboratorio 3 , quello con il calcolo parallelo 
#ifndef WAVE_EQUATION_HPP
#define WAVE_EQUATION_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <filesystem>
#include <fstream>
#include <iostream>


namespace WaveEquationProject {

    using namespace dealii;
 
    enum TimeType { IMPLICITO, ESPLICITO};

    template <int dim>
    class WaveEquation {
           public:

    // Costruttore base
                WaveEquation(const std::string &param_file);
    //metodo per la run del programma
                void run();

           private:
    //metodi per il setup e per la cstruzione della griglia
                void setup_system();
                void make_grid();
    //metodi per l'assemblaggio di Brhs a ogni step
                void assemble_system_IMPLICITO();
                void assemble_system_ESPLICITO();
    //metodi per la risoluzione del sistema lineare
                void solve_IMPLICITO();
                void solve_ESPLICITO();
    //metodo per l'output dei risultati
                void output_results(const unsigned int time_step);


    //vado a definire degli oggetti MPI e degli oggetti per la gestione dei data
    //ho preso ispirazione dai laboratori e dal resto

    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    DoFHandler<dim> dof_handler;

    //dettaglio tecnico , ho utilizzato LinearAlgebra per la gestione dei vettori perche piu rubusto , dovendo invece usare un precodinzionatore scalabile (AMG multi grid) per la risoluzione del sistema lineare ho dovuto utilizzare le librerie corrette per dealii
    //
    LinearAlgebra::distributed::Vector<double> U; //spostamento a tempo t_n
    LinearAlgebra::distributed::Vector<double> U_t; //velocita a tempo t_n
    LinearAlgebra::distributed::Vector<double> U_tt; //accelerazione a tempo t_n
    LinearAlgebra::distributed::Vector<double> B_rhs; //spostamento a tempo
    
    //Definisco le matrici adesso !

    TrilinosWrappers::SparseMatrix A_sys; //matrice di sistema
    TrilinosWrappers::SparseMatrix M; //matrice di massa
    TrilinosWrappers::SparseMatrix K; //matrice di laplace
    //
    //per la matrice di massa devo definire la lumpded mass matrix nel xaso ESPLICITO

    LinearAlgebra::distributed::Vector<double> M_lumped;

    //dat di simulazione
    //
    double time;
    double time_step;
    double final_time;
    unsigned int output_time_step; //devo salvarlo ogni tot step se no mi ritrovo il disco del pc pieno
    
    //parametri di Newmark
    
    TimeType method_type;
    double beta;
    double gamma;

   };
}
#endif


