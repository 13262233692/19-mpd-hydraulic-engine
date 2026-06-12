cmake_minimum_required(VERSION 3.18)
project(hydraulic_engine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/O2 /fp:fast /openmp)
else()
    add_compile_options(-O3 -ffast-math -fopenmp -march=native)
endif()

set(HYDRAULIC_ENGINE_SOURCES
    src/cpp/grid.cpp
    src/cpp/variables.cpp
    src/cpp/physics.cpp
    src/cpp/fvm.cpp
    src/cpp/pentadiagonal_solver.cpp
    src/cpp/time_integrator.cpp
    src/cpp/simulator.cpp
    src/cpp/choke_valve.cpp
)

add_library(hydraulic_engine_core STATIC ${HYDRAULIC_ENGINE_SOURCES})
target_include_directories(hydraulic_engine_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(Python3 COMPONENTS Interpreter Development)
if(Python3_FOUND)
    find_package(pybind11 CONFIG)
    if(pybind11_FOUND)
        pybind11_add_module(hydraulic_engine src/python/bindings.cpp)
        target_link_libraries(hydraulic_engine PRIVATE hydraulic_engine_core)
        set_target_properties(hydraulic_engine PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/hydraulic_engine
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/hydraulic_engine
        )
    else()
        message(STATUS "pybind11 not found, skipping Python bindings")
    endif()
else()
    message(STATUS "Python3 not found, skipping Python bindings")
endif()

enable_testing()
add_executable(test_solver tests/test_pentadiagonal.cpp)
target_link_libraries(test_solver PRIVATE hydraulic_engine_core)
add_test(NAME test_pentadiagonal COMMAND test_solver)

add_executable(test_fvm tests/test_fvm.cpp)
target_link_libraries(test_fvm PRIVATE hydraulic_engine_core)
add_test(NAME test_fvm COMMAND test_fvm)

add_executable(test_simulator tests/test_simulator.cpp)
target_link_libraries(test_simulator PRIVATE hydraulic_engine_core)
add_test(NAME test_simulator COMMAND test_simulator)

add_executable(test_kick_stability tests/test_kick_stability.cpp)
target_link_libraries(test_kick_stability PRIVATE hydraulic_engine_core)
add_test(NAME test_kick_stability COMMAND test_kick_stability)

add_executable(test_choke_control tests/test_choke_control.cpp)
target_link_libraries(test_choke_control PRIVATE hydraulic_engine_core)
add_test(NAME test_choke_control COMMAND test_choke_control)
