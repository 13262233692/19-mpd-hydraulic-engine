#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include "hydraulic_engine/simulator.hpp"

using namespace hydraulic_engine;

int main() {
    std::cout << "=========================================" << std::endl;
    std::cout << "Testing Gas Kick Stability (Sharp Front)" << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << std::endl;
    
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
    
    std::vector<FluxLimiter::Type> limiters = {
        FluxLimiter::FIRST_ORDER_UPWIND,
        FluxLimiter::MINMOD,
        FluxLimiter::SUPERBEE,
        FluxLimiter::VAN_LEER,
        FluxLimiter::VAN_ALBADA
    };
    
    std::vector<std::string> limiter_names = {
        "FIRST_ORDER_UPWIND",
        "MINMOD",
        "SUPERBEE",
        "VAN_LEER",
        "VAN_ALBADA"
    };
    
    bool all_passed = true;
    
    for (size_t l = 0; l < limiters.size(); ++l) {
        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "Testing limiter: " << limiter_names[l] << std::endl;
        std::cout << "-----------------------------------------" << std::endl;
        
        config.flux_limiter = limiters[l];
        
        try {
            Simulator simulator(config);
            
            const auto& output_initial = simulator.get_output();
            std::cout << "  Initial conditions:" << std::endl;
            std::cout << "    Top pressure: " << output_initial.pressure[0] / 1e6 
                      << " MPa" << std::endl;
            std::cout << "    Bottom pressure: " 
                      << output_initial.pressure[config.num_cells - 1] / 1e6 
                      << " MPa" << std::endl;
            std::cout << "    Max void fraction: " 
                      << *std::max_element(output_initial.void_fraction.begin(), 
                                           output_initial.void_fraction.end()) 
                      << std::endl;
            
            std::cout << "\n  Enabling extreme kick at 2000m depth..." << std::endl;
            std::cout << "  Kick flow rate: 10.0 m^3/s (extreme)" << std::endl;
            std::cout << "  Kick void fraction: 0.95 (near pure gas)" << std::endl;
            
            simulator.enable_kick(10.0, 0.95, 2000.0, 0.0, 5.0);
            
            std::cout << "\n  Running simulation for 3 seconds..." << std::endl;
            std::cout.flush();
            
            auto start = std::chrono::high_resolution_clock::now();
            
            bool success = simulator.run(3.0, 100);
            
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            
            const auto& output = simulator.get_output();
            
            std::cout << "\n  Simulation completed in " << elapsed.count() 
                      << " seconds" << std::endl;
            std::cout << "  Total steps: " << simulator.get_total_steps() 
                      << std::endl;
            std::cout << "  Final time: " << simulator.get_current_time() 
                      << " s" << std::endl;
            std::cout << "  Avg time/step: " 
                      << elapsed.count() / simulator.get_total_steps() * 1000 
                      << " ms" << std::endl;
            
            bool has_nan = false;
            bool has_negative_density = false;
            bool has_negative_pressure = false;
            bool has_invalid_void_fraction = false;
            
            double max_void = 0.0;
            int max_void_idx = 0;
            double min_density = 1e9;
            double min_pressure = 1e9;
            
            for (int i = 0; i < config.num_cells; ++i) {
                if (!std::isfinite(output.pressure[i])) has_nan = true;
                if (!std::isfinite(output.void_fraction[i])) has_nan = true;
                if (!std::isfinite(output.mixture_density[i])) has_nan = true;
                if (!std::isfinite(output.mixture_velocity[i])) has_nan = true;
                
                if (output.pressure[i] <= 0) has_negative_pressure = true;
                if (output.mixture_density[i] <= 0) has_negative_density = true;
                if (output.void_fraction[i] < -1e-6 || output.void_fraction[i] > 1.0 + 1e-6) 
                    has_invalid_void_fraction = true;
                
                max_void = std::max(max_void, output.void_fraction[i]);
                if (output.void_fraction[i] > max_void) max_void_idx = i;
                min_density = std::min(min_density, output.mixture_density[i]);
                min_pressure = std::min(min_pressure, output.pressure[i]);
            }
            
            std::cout << "\n  Results:" << std::endl;
            std::cout << "    Success: " << (success ? "YES" : "NO") << std::endl;
            std::cout << "    Has NaN: " << (has_nan ? "YES (CRITICAL!)" : "NO") << std::endl;
            std::cout << "    Negative density: " 
                      << (has_negative_density ? "YES (CRITICAL!)" : "NO") << std::endl;
            std::cout << "    Negative pressure: " 
                      << (has_negative_pressure ? "YES (CRITICAL!)" : "NO") << std::endl;
            std::cout << "    Invalid void fraction: " 
                      << (has_invalid_void_fraction ? "YES (CRITICAL!)" : "NO") << std::endl;
            std::cout << "    Max void fraction: " << max_void << std::endl;
            std::cout << "    Max void location: " 
                      << config.total_depth - max_void_idx * (config.total_depth / config.num_cells) 
                      << " m" << std::endl;
            std::cout << "    Min mixture density: " << min_density << " kg/m^3" << std::endl;
            std::cout << "    Min pressure: " << min_pressure / 1e6 << " MPa" << std::endl;
            
            bool test_passed = success && !has_nan && !has_negative_density && 
                              !has_negative_pressure && !has_invalid_void_fraction &&
                              max_void > 0.1;
            
            std::cout << "\n  Test " << (test_passed ? "PASSED" : "FAILED") << std::endl;
            
            if (!test_passed) all_passed = false;
            
        } catch (const std::exception& e) {
            std::cout << "\n  EXCEPTION: " << e.what() << std::endl;
            std::cout << "  Test FAILED" << std::endl;
            all_passed = false;
        }
        
        std::cout << std::endl;
    }
    
    std::cout << "=========================================" << std::endl;
    std::cout << "Overall Result: " << (all_passed ? "ALL TESTS PASSED" : "SOME TESTS FAILED") << std::endl;
    std::cout << "=========================================" << std::endl;
    
    return all_passed ? 0 : 1;
}
