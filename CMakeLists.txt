cmake_minimum_required(VERSION 3.12.0)
project(wave_equation_solver LANGUAGES CXX C)

# Includiamo le configurazioni comuni
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake-common.cmake)

# Definiamo l'eseguibile con i nostri file sorgente
add_executable(wave_solver
    main.cpp
    wave_equation.cpp
)

# Configurazione automatica di deal.II (include path, flags, link libs)
deal_ii_setup_target(wave_solver)

# --- FIX: Rimosso "PRIVATE" ---
# Dobbiamo usare lo stile "plain" per compatibilità con deal_ii_setup_target
target_link_libraries(wave_solver MPI::MPI_CXX)
