#include "input_parser.hpp"

namespace WaveEquationProject {

    InputParser::InputParser() {
        declare_parameters();
    }

    void InputParser::declare_parameters() {
        
        prm.enter_subsection("Time Settings");
        {
            prm.declare_entry("Time step", "0.01",
                              Patterns::Double(0.0),
                              "The time step size (dt).");

            prm.declare_entry("Final time", "1.0",
                              Patterns::Double(0.0),
                              "The end time of the simulation.");
            
            prm.declare_entry("Output frequency", "1",
                              Patterns::Integer(1),
                              "Save results every N time steps.");
        }
        prm.leave_subsection();

    
        prm.enter_subsection("Scenario");
        {
            prm.declare_entry("Scenario ID", "1",
                              Patterns::Integer(1, 6),
                              "Select test case: 1=Conv, 2=Gauss, 3=Collision, 4=Rain, 5=Ring, 6=Wall");
        }
        prm.leave_subsection();

        
        prm.enter_subsection("Solver Settings");
        {
            prm.declare_entry("Method", "IMPLICIT",
                              Patterns::Selection("IMPLICIT|EXPLICIT"),
                              "Time integration scheme: IMPLICIT or EXPLICIT");

            prm.declare_entry("Newmark Beta", "0.25",
                              Patterns::Double(0.0),
                              "Beta parameter for Newmark scheme.");

            prm.declare_entry("Newmark Gamma", "0.5",
                              Patterns::Double(0.0),
                              "Gamma parameter for Newmark scheme.");
        }
        prm.leave_subsection();
    }

    void InputParser::parse_parameters(const std::string &filename) {
        if (filename.empty()) {
             std::cout << "No parameter file given, using defaults." << std::endl;
             return;
        }
        
        try {
            prm.parse_input(filename);
        } catch (std::exception &e) {
            std::cerr << "Error parsing parameter file: " << e.what() << std::endl;
            std::cerr << "Generating a default 'parameters.prm' file for you." << std::endl;
            std::ofstream out("parameters.prm");
            prm.print_parameters(out, ParameterHandler::Text);
            throw; 
        }
    }

    
    double InputParser::get_time_step() const {
        prm.enter_subsection("Time Settings");
        double val = prm.get_double("Time step");
        prm.leave_subsection();
        return val;
    }

    double InputParser::get_final_time() const {
        prm.enter_subsection("Time Settings");
        double val = prm.get_double("Final time");
        prm.leave_subsection();
        return val;
    }

    unsigned int InputParser::get_output_frequency() const {
        prm.enter_subsection("Time Settings");
        int val = prm.get_integer("Output frequency");
        prm.leave_subsection();
        return static_cast<unsigned int>(val);
    }

    unsigned int InputParser::get_scenario_id() const {
        prm.enter_subsection("Scenario");
        int val = prm.get_integer("Scenario ID");
        prm.leave_subsection();
        return static_cast<unsigned int>(val);
    }

    std::string InputParser::get_solver_type() const {
        prm.enter_subsection("Solver Settings");
        std::string val = prm.get("Method");
        prm.leave_subsection();
        return val;
    }

    double InputParser::get_beta() const {
        prm.enter_subsection("Solver Settings");
        double val = prm.get_double("Newmark Beta");
        prm.leave_subsection();
        return val;
    }

    double InputParser::get_gamma() const {
        prm.enter_subsection("Solver Settings");
        double val = prm.get_double("Newmark Gamma");
        prm.leave_subsection();
        return val;
    }
}
