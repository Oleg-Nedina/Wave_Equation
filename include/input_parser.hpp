#ifndef INPUT_PARSER_HPP
#define INPUT_PARSER_HPP

#include <deal.II/base/parameter_handler.h>
#include <string>
#include <iostream>
#include <fstream> 

namespace WaveEquationProject {
    using namespace dealii;

    class InputParser {
    public:
        InputParser();

        // Funzione principale di parsing
        void parse_parameters(const std::string &filename);

        // Getters
        double get_time_step() const;
        double get_final_time() const;
        unsigned int get_output_frequency() const;
        unsigned int get_scenario_id() const;
        
        std::string get_solver_type() const; 

        double get_beta() const;
        double get_gamma() const;

    private:
        // 'mutable' permette ai metodi const di modificare lo stato interno
        mutable ParameterHandler prm;

        // Dichiarazione della funzione privata (deve stare qui, PRIMA della chiusura della classe)
        void declare_parameters();
    };
}

#endif
