cmake_minimum_required(VERSION 3.12.0)
project(wave_equation_solver LANGUAGES CXX C)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake-common.cmake)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/wave_equation
)

add_executable(wave_solver
    main.cpp                          
    wave_equation/wave_equation.cpp  
    include/input_parser.cpp       
)

deal_ii_setup_target(wave_solver)

target_link_libraries(wave_solver MPI::MPI_CXX)
