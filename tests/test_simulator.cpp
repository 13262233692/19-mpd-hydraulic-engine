#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include "hydraulic_engine/simulator.hpp"

using namespace hydraulic_engine;

int main() {
    std::cout << "Testing Hydraulic Engine Simulator..." << std::endl;
    std::cout << "=========================================" << std::endl;
    
    SimulationConfig config;
    config.total_depth = 3000.0;
    config.num_cells = 100;
    config.wellbore_diameter = 0.216;
    config.initial_pressure = 1e5;
    config.initial_void_fraction = 0.0;
    config.initial_temperature = 300.0;
    config.liquid_density_ref = 1200.0;
    config.gas_density_ref = 150.0;
    config.liquid_viscosity = 1e-3;
    config.gas_viscosity = 1.8e-5;
    config.surface_tension = 0.072;
    config.molar_mass_gas = 0.016;
    config.distribution_coefficient = 1.2;
    config.drift_velocity_coeff = 0.35;
    config.gravity = 9.81;
    config.dt_init = 0.001;
    config.dt_min = 1e-5;
    config.dt_max = 0.01;
    config.cfl = 0.05;
    config.max_newton_iter = 30;
    config.newton_tol = 1.0;
    config.adaptive_time_stepping = false;
    config.flux_limiter = FluxLimiter::FIRST_ORDER_UPWIND;
    
    std::cout << "Creating simulator..." << std::endl;
    Simulator simulator(config);
    
    std::cout << "Grid: " << config.num_cells << " cells, " 
              << config.total_depth << " m depth" << std::endl;
    
    const auto& output = simulator.get_output();
    std::cout << "\nInitial conditions:" << std::endl;
    std::cout << "  Top pressure: " << output.pressure[0] / 1e6 << " MPa" << std::endl;
    std::cout << "  Bottom pressure: " << output.pressure[config.num_cells - 1] / 1e6 << " MPa" << std::endl;
    std::cout << "  Max void fraction: " << 
        *std::max_element(output.void_fraction.begin(), output.void_fraction.end()) << std::endl;
    
    std::cout << "\nRunning 10s without kick to verify steady state..." << std::endl;
    simulator.run(10.0, 100000);
    
    auto after_steady = simulator.get_output();
    std::cout << "After steady state:" << std::endl;
    std::cout << "  Top pressure: " << after_steady.pressure[0] / 1e6 << " MPa" << std::endl;
    std::cout << "  Bottom pressure: " << after_steady.pressure[config.num_cells - 1] / 1e6 << " MPa" << std::endl;
    
    std::cout << "\nEnabling kick at 2000m depth..." << std::endl;
    simulator.enable_kick(5.0, 0.8, 2000.0, 10.0, 30.0);
    
    double total_simulation_time = 15.0;
    std::cout << "\nRunning simulation for total " << total_simulation_time << " seconds..." << std::endl;
    std::cout.flush();
    
    auto start = std::chrono::high_resolution_clock::now();
    
    bool success = simulator.run(total_simulation_time, 20);
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    if (!success) {
        std::cout << "ERROR: Simulation failed!" << std::endl;
        return 1;
    }
    
    std::cout << "\nSimulation completed!" << std::endl;
    std::cout << "  Total steps: " << simulator.get_total_steps() << std::endl;
    std::cout << "  Final time: " << simulator.get_current_time() << " s" << std::endl;
    std::cout << "  Elapsed time: " << elapsed.count() << " s" << std::endl;
    std::cout << "  Average time per step: " 
              << elapsed.count() / simulator.get_total_steps() * 1000 << " ms" << std::endl;
    
    const auto& final_output = simulator.get_output();
    
    std::cout << "\nFinal results:" << std::endl;
    std::cout << "  Max void fraction: " << 
        *std::max_element(final_output.void_fraction.begin(), 
                         final_output.void_fraction.end()) << std::endl;
    
    double max_void = 0.0;
    int max_void_idx = 0;
    for (int i = 0; i < config.num_cells; ++i) {
        if (final_output.void_fraction[i] > max_void) {
            max_void = final_output.void_fraction[i];
            max_void_idx = i;
        }
    }
    std::cout << "  Max void fraction location: " 
              << final_output.depth[max_void_idx] << " m" << std::endl;
    
    double top_void = final_output.void_fraction[0];
    double bottom_void = final_output.void_fraction[config.num_cells - 1];
    std::cout << "  Top void fraction: " << top_void << std::endl;
    std::cout << "  Bottom void fraction: " << bottom_void << std::endl;
    
    std::cout << "\nPressure profile (selected points):" << std::endl;
    for (int i = 0; i < config.num_cells; i += config.num_cells / 5) {
        std::cout << "  Depth " << final_output.depth[i] << " m: " 
                  << final_output.pressure[i] / 1e6 << " MPa, "
                  << "void: " << final_output.void_fraction[i] << std::endl;
    }
    
    std::cout << "\nVelocity profile (selected points):" << std::endl;
    for (int i = 0; i < config.num_cells; i += config.num_cells / 5) {
        std::cout << "  Depth " << final_output.depth[i] << " m: "
                  << "vm=" << final_output.mixture_velocity[i] << " m/s, "
                  << "vg=" << final_output.gas_velocity[i] << " m/s, "
                  << "vl=" << final_output.liquid_velocity[i] << " m/s" << std::endl;
    }
    
    double max_velocity = 0.0;
    for (int i = 0; i < config.num_cells; ++i) {
        max_velocity = std::max(max_velocity, std::abs(final_output.mixture_velocity[i]));
    }
    std::cout << "\nMax mixture velocity: " << max_velocity << " m/s" << std::endl;
    
    std::cout << "\nLast Newton iteration: " << simulator.get_last_newton_iterations() 
              << " iterations, residual: " << simulator.get_last_newton_residual() << std::endl;
    
    if (max_void < 0.01) {
        std::cout << "WARNING: Kick may not have propagated properly" << std::endl;
    }
    
    std::cout << "\n=========================================" << std::endl;
    std::cout << "All simulator tests completed successfully!" << std::endl;
    
    return 0;
}
