#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <fstream>
#include "hydraulic_engine/simulator.hpp"
#include "hydraulic_engine/choke_valve.hpp"

using namespace hydraulic_engine;

int main() {
    std::cout << "============================================================" << std::endl;
    std::cout << "  CHOKE VALVE CLOSED-LOOP PRESSURE CONTROL TEST" << std::endl;
    std::cout << "============================================================" << std::endl;
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
    config.flux_limiter = FluxLimiter::MINMOD;
    
    double hydrostatic_bottom_pressure = config.initial_pressure + 
        config.liquid_density_ref * config.gravity * config.total_depth;
    
    std::cout << "  Well configuration:" << std::endl;
    std::cout << "    Total depth: " << config.total_depth << " m" << std::endl;
    std::cout << "    Wellbore diameter: " << config.wellbore_diameter << " m" << std::endl;
    std::cout << "    Hydrostatic BHP: " << hydrostatic_bottom_pressure / 1e6 
              << " MPa" << std::endl;
    std::cout << std::endl;
    
    ChokeValveConfig choke_config;
    choke_config.discharge_coefficient = 0.85;
    choke_config.max_area = 0.01;
    choke_config.min_area = 1e-5;
    choke_config.initial_area = 0.005;
    choke_config.max_area_change_rate = 0.5;
    choke_config.target_bhp = hydrostatic_bottom_pressure * 1.05;
    choke_config.bhp_tolerance = 1e5;
    choke_config.kp = 1e-7;
    choke_config.ki = 2e-6;
    choke_config.kd = 5e-9;
    choke_config.integral_min = -5.0;
    choke_config.integral_max = 5.0;
    choke_config.derivative_filter_tau = 0.05;
    
    std::cout << "  Choke valve configuration:" << std::endl;
    std::cout << "    Discharge coefficient: " << choke_config.discharge_coefficient << std::endl;
    std::cout << "    Max area: " << choke_config.max_area << " m^2" << std::endl;
    std::cout << "    Min area: " << choke_config.min_area << " m^2" << std::endl;
    std::cout << "    Initial opening: " 
              << (choke_config.initial_area - choke_config.min_area) / 
                 (choke_config.max_area - choke_config.min_area) * 100.0 
              << " %" << std::endl;
    std::cout << "    Target BHP: " << choke_config.target_bhp / 1e6 << " MPa" << std::endl;
    std::cout << "    BHP tolerance: " << choke_config.bhp_tolerance / 1e6 << " MPa" << std::endl;
    std::cout << "    Kp: " << choke_config.kp << std::endl;
    std::cout << "    Ki: " << choke_config.ki << std::endl;
    std::cout << "    Kd: " << choke_config.kd << std::endl;
    std::cout << std::endl;
    
    std::cout << "  Test scenario: Gas kick with closed-loop choke control" << std::endl;
    std::cout << "    Kick flow rate: 2.0 m^3/s" << std::endl;
    std::cout << "    Kick void fraction: 0.9" << std::endl;
    std::cout << "    Kick depth: 2500 m" << std::endl;
    std::cout << "    Simulation time: 5 seconds" << std::endl;
    std::cout << std::endl;
    
    Simulator simulator(config);
    
    std::cout << "  Enabling choke valve control..." << std::endl;
    simulator.enable_choke_valve_control(choke_config);
    
    std::cout << "  Enabling gas kick at t=1.0s..." << std::endl;
    simulator.enable_kick(2.0, 0.9, 2500.0, 1.0, 10.0);
    
    std::cout << std::endl;
    std::cout << "  Running simulation..." << std::endl;
    std::cout << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    bool success = simulator.run(5.0, 10);
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    std::cout << std::endl;
    std::cout << "  Simulation completed in " << elapsed.count() 
              << " seconds" << std::endl;
    std::cout << "  Total steps: " << simulator.get_total_steps() << std::endl;
    std::cout << "  Final time: " << simulator.get_current_time() << " s" << std::endl;
    std::cout << std::endl;
    
    const auto& output = simulator.get_output();
    
    bool has_nan = false;
    bool has_negative_density = false;
    bool has_negative_pressure = false;
    
    for (int i = 0; i < config.num_cells; ++i) {
        if (!std::isfinite(output.pressure[i])) has_nan = true;
        if (!std::isfinite(output.void_fraction[i])) has_nan = true;
        if (!std::isfinite(output.mixture_density[i])) has_nan = true;
        
        if (output.pressure[i] <= 0) has_negative_pressure = true;
        if (output.mixture_density[i] <= 0) has_negative_density = true;
    }
    
    std::cout << "  Numerical stability check:" << std::endl;
    std::cout << "    Has NaN: " << (has_nan ? "YES (FAILED)" : "NO (PASSED)") << std::endl;
    std::cout << "    Negative density: " 
              << (has_negative_density ? "YES (FAILED)" : "NO (PASSED)") << std::endl;
    std::cout << "    Negative pressure: " 
              << (has_negative_pressure ? "YES (FAILED)" : "NO (PASSED)") << std::endl;
    std::cout << std::endl;
    
    double final_bhp = output.pressure[config.num_cells - 1];
    double final_bhp_error = final_bhp - choke_config.target_bhp;
    double max_void = *std::max_element(output.void_fraction.begin(), 
                                         output.void_fraction.end());
    
    std::cout << "  Final state:" << std::endl;
    std::cout << "    Final BHP: " << final_bhp / 1e6 << " MPa" << std::endl;
    std::cout << "    Target BHP: " << choke_config.target_bhp / 1e6 << " MPa" << std::endl;
    std::cout << "    BHP error: " << final_bhp_error / 1e6 << " MPa" << std::endl;
    std::cout << "    Max void fraction: " << max_void << std::endl;
    std::cout << "    Wellhead pressure: " << output.pressure[0] / 1e6 << " MPa" << std::endl;
    std::cout << std::endl;
    
    simulator.print_choke_control_report();
    
    bool control_passed = false;
    if (simulator.get_choke_controller()) {
        const auto& log = simulator.get_choke_controller()->get_log();
        if (!log.empty()) {
            double max_error = 0.0;
            double avg_error = 0.0;
            int within_tolerance = 0;
            
            size_t start_idx = log.size() / 2;
            for (size_t i = start_idx; i < log.size(); ++i) {
                double err = std::abs(log[i].bhp_error);
                max_error = std::max(max_error, err);
                avg_error += err;
                if (err <= choke_config.bhp_tolerance) {
                    within_tolerance++;
                }
            }
            avg_error /= (log.size() - start_idx);
            
            double pct_within = 100.0 * within_tolerance / (log.size() - start_idx);
            
            std::cout << "  Control performance (second half of simulation):" << std::endl;
            std::cout << "    Max absolute error: " << max_error / 1e6 << " MPa" << std::endl;
            std::cout << "    Average absolute error: " << avg_error / 1e6 << " MPa" << std::endl;
            std::cout << "    Time within tolerance: " << pct_within << " %" << std::endl;
            std::cout << std::endl;
            
            control_passed = pct_within > 50.0 && max_error < 5e6;
        }
    }
    
    bool overall_passed = success && !has_nan && !has_negative_density && 
                          !has_negative_pressure && control_passed;
    
    std::cout << "============================================================" << std::endl;
    std::cout << "  OVERALL RESULT: " << (overall_passed ? "PASSED" : "FAILED") << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;
    
    return overall_passed ? 0 : 1;
}
