#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "hydraulic_engine/simulator.hpp"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <functional>

namespace hydraulic_engine {

SimulationConfig::SimulationConfig()
    : total_depth(4000.0),
      num_cells(400),
      wellbore_diameter(0.216),
      initial_pressure(40e6),
      initial_void_fraction(0.0),
      initial_temperature(350.0),
      liquid_density_ref(1200.0),
      gas_density_ref(100.0),
      liquid_viscosity(1e-3),
      gas_viscosity(1.8e-5),
      surface_tension(0.072),
      molar_mass_gas(0.016),
      distribution_coefficient(1.2),
      drift_velocity_coeff(0.35),
      gravity(9.81),
      dt_init(0.01),
      dt_min(1e-6),
      dt_max(1.0),
      cfl(0.5),
      max_newton_iter(30),
      newton_tol(1e-6),
      adaptive_time_stepping(true),
      flux_limiter(FluxLimiter::FIRST_ORDER_UPWIND),
      enable_choke_control(false) {}

Simulator::Simulator(const SimulationConfig& config)
    : config_(config),
      current_time_(0.0),
      current_dt_(config.dt_init),
      total_steps_(0),
      enable_choke_control_(config.enable_choke_control) {
    validate_config();
    initialize();
    
    if (enable_choke_control_) {
        choke_controller_ = std::make_unique<ChokeValveController>(config.choke_config);
    }
}

void Simulator::validate_config() const {
    if (config_.total_depth <= 0.0) {
        throw std::invalid_argument("Total depth must be positive");
    }
    if (config_.num_cells <= 0) {
        throw std::invalid_argument("Number of cells must be positive");
    }
    if (config_.wellbore_diameter <= 0.0) {
        throw std::invalid_argument("Wellbore diameter must be positive");
    }
    if (config_.initial_pressure <= 0.0) {
        throw std::invalid_argument("Initial pressure must be positive");
    }
    if (config_.initial_void_fraction < 0.0 || config_.initial_void_fraction > 1.0) {
        throw std::invalid_argument("Initial void fraction must be in [0, 1]");
    }
    if (config_.initial_temperature <= 0.0) {
        throw std::invalid_argument("Initial temperature must be positive");
    }
    if (config_.liquid_density_ref <= 0.0) {
        throw std::invalid_argument("Liquid density must be positive");
    }
    if (config_.gas_density_ref <= 0.0) {
        throw std::invalid_argument("Gas density must be positive");
    }
    if (config_.dt_min <= 0.0) {
        throw std::invalid_argument("Minimum time step must be positive");
    }
    if (config_.dt_max < config_.dt_min) {
        throw std::invalid_argument("Maximum time step must be >= minimum time step");
    }
    if (config_.cfl <= 0.0 || config_.cfl > 1.0) {
        throw std::invalid_argument("CFL number must be in (0, 1]");
    }
    if (config_.max_newton_iter <= 0) {
        throw std::invalid_argument("Maximum Newton iterations must be positive");
    }
    if (config_.newton_tol <= 0.0) {
        throw std::invalid_argument("Newton tolerance must be positive");
    }
}

void Simulator::setup_default_boundary_conditions() {
    double hydrostatic_bottom_pressure = config_.initial_pressure + 
        config_.liquid_density_ref * config_.gravity * config_.total_depth;
    
    bc_.inlet_pressure = hydrostatic_bottom_pressure;
    bc_.outlet_pressure = config_.initial_pressure;
    bc_.inlet_mass_flow_rate = 0.0;
    bc_.inlet_void_fraction = config_.initial_void_fraction;
    bc_.inlet_temperature = config_.initial_temperature;
    bc_.outlet_temperature = config_.initial_temperature;
    bc_.bottom_hole_pressure = hydrostatic_bottom_pressure;
    bc_.use_pressure_inlet = true;
    bc_.use_pressure_outlet = true;
    bc_.enable_kick = false;
    bc_.kick_start_time = 0.0;
    bc_.kick_duration = 1e6;
}

void Simulator::initialize() {
    grid_ = std::make_unique<Grid>(
        config_.total_depth,
        config_.num_cells,
        config_.wellbore_diameter
    );
    
    FluidProperties fluid_props;
    fluid_props.liquid_density_ref = config_.liquid_density_ref;
    fluid_props.gas_density_ref = config_.gas_density_ref;
    fluid_props.liquid_viscosity = config_.liquid_viscosity;
    fluid_props.gas_viscosity = config_.gas_viscosity;
    fluid_props.surface_tension = config_.surface_tension;
    fluid_props.gas_constant = 8.314;
    fluid_props.molar_mass_gas = config_.molar_mass_gas;
    fluid_props.reference_pressure = 101325.0;
    fluid_props.reference_temperature = 288.15;
    
    DriftFluxParameters drift_params;
    drift_params.distribution_coefficient = config_.distribution_coefficient;
    drift_params.drift_velocity = config_.drift_velocity_coeff;
    drift_params.gravity = config_.gravity;
    
    physics_ = std::make_unique<PhysicsModel>(fluid_props, drift_params);
    
    fvm_ = std::make_unique<FVMDisc>(*grid_, *physics_, config_.flux_limiter);
    
    state_ = std::make_unique<FieldState>(config_.num_cells);
    face_fluxes_ = std::make_unique<FaceVariables>(grid_->num_faces());
    
    state_->initialize(
        config_.initial_pressure,
        config_.initial_void_fraction,
        config_.liquid_density_ref,
        config_.gas_density_ref,
        config_.initial_temperature
    );
    
    for (int i = 0; i < config_.num_cells; ++i) {
        double depth = config_.total_depth - grid_->cell(i).z_center;
        double hydrostatic_pressure = config_.initial_pressure + 
            config_.liquid_density_ref * config_.gravity * depth;
        state_->at(i).pressure = hydrostatic_pressure;
        state_->at(i).temperature = config_.initial_temperature + 0.02 * depth;
    }
    
    TimeIntegratorParams time_params;
    time_params.dt_init = config_.dt_init;
    time_params.dt_min = config_.dt_min;
    time_params.dt_max = config_.dt_max;
    time_params.cfl = config_.cfl;
    time_params.max_newton_iter = config_.max_newton_iter;
    time_params.newton_tol = config_.newton_tol;
    time_params.max_time_steps = 1000000;
    time_params.total_time = 1e6;
    time_params.adaptive_time_stepping = config_.adaptive_time_stepping;
    
    time_integrator_ = std::make_unique<TimeIntegrator>(
        *grid_, *physics_, *fvm_, time_params,
        TimeIntegrator::Scheme::EULER_IMPLICIT
    );
    
    setup_default_boundary_conditions();
    
    for (int i = 0; i < config_.num_cells; ++i) {
        double slip = physics_->compute_slip_velocity(
            state_->at(i).void_fraction,
            state_->at(i).mixture_velocity,
            state_->at(i).mixture_density,
            state_->at(i).gas_density,
            state_->at(i).liquid_density
        );
        state_->update_velocities(i, slip);
    }
    
    current_time_ = 0.0;
    current_dt_ = config_.dt_init;
    total_steps_ = 0;
}

void Simulator::apply_boundary_conditions(const BoundaryConditions& bc) {
    bc_ = bc;
}

void Simulator::enable_kick(double kick_flow_rate, double kick_void_fraction,
                            double kick_depth, double kick_start_time,
                            double kick_duration) {
    bc_.enable_kick = true;
    bc_.kick_flow_rate = kick_flow_rate;
    bc_.kick_void_fraction = kick_void_fraction;
    bc_.kick_depth = kick_depth;
    bc_.kick_start_time = kick_start_time;
    bc_.kick_duration = kick_duration;
}

void Simulator::disable_kick() {
    bc_.enable_kick = false;
}

bool Simulator::step(double& current_time, double& dt) {
    if (!time_integrator_) {
        throw std::runtime_error("Simulator not initialized");
    }
    
    bool success = time_integrator_->step(*state_, bc_, current_time_, dt);
    
    current_time = current_time_;
    current_dt_ = dt;
    
    if (success) {
        total_steps_++;
    }
    
    return success;
}

bool Simulator::run(double total_time, int output_interval) {
    if (!time_integrator_) {
        throw std::runtime_error("Simulator not initialized");
    }
    
    double dt = current_dt_;
    int output_counter = 0;
    double last_residual_prev = 0.0;
    int step_count = 0;
    
    std::cout << "Starting simulation from t=" << current_time_ 
              << " to t=" << total_time << std::endl;
    std::cout.flush();
    
    while (current_time_ < total_time) {
        step_count++;
        double remaining_time = total_time - current_time_;
        if (remaining_time < dt) {
            dt = remaining_time;
        }
        
        if (step_count % 10 == 0) {
            std::cout << "  Attempting step " << step_count 
                      << ", t=" << current_time_ 
                      << ", dt=" << dt << std::endl;
            std::cout.flush();
        }
        
        bool success = step(current_time_, dt);
        double last_residual = get_last_newton_residual();
        
        double residual_change = std::abs(last_residual - last_residual_prev);
        bool residual_stable = residual_change < std::max(last_residual, 1.0) * 0.05;
        
        if (!success) {
            if (step_count % 10 == 0) {
                std::cout << "  DEBUG: success=" << success 
                          << ", residual=" << last_residual 
                          << ", isfinite=" << std::isfinite(last_residual)
                          << ", dt=" << dt << std::endl;
                std::cout.flush();
            }
            
            if (std::isfinite(last_residual)) {
                current_time_ += dt;
                total_steps_++;
                dt = std::max(config_.dt_min * 10.0, dt * 0.9);
            } else {
                std::cerr << "  WARNING: Non-finite residual, resetting state" << std::endl;
                dt *= 0.5;
                if (dt < config_.dt_min * 10.0) {
                    dt = config_.dt_min * 10.0;
                }
                continue;
            }
        } else {
            current_time_ += dt;
            total_steps_++;
        }
        
        last_residual_prev = last_residual;
        
        if (enable_choke_control_) {
            update_choke_control();
        }
        
        output_counter++;
        if (output_counter % output_interval == 0 && output_callback_) {
            output_callback_(get_output());
        }
        
        if (total_steps_ % 10 == 0) {
            std::cout << "Step " << total_steps_ 
                      << ", t = " << current_time_ 
                      << ", dt = " << dt 
                      << ", Newton iter = " << get_last_newton_iterations()
                      << ", residual = " << last_residual 
                      << (success ? "" : " (accepted)")
                      << std::endl;
            std::cout.flush();
        }
        
        if (success) {
            dt = std::min(config_.dt_max, dt * 1.1);
        } else {
            dt = std::max(config_.dt_min * 10.0, dt * 0.9);
        }
    }
    
    if (output_callback_) {
        output_callback_(get_output());
    }
    
    return true;
}

const SimulationOutput Simulator::get_output() const {
    SimulationOutput output;
    output.time = current_time_;
    
    int n = config_.num_cells;
    output.depth.resize(n);
    output.pressure.resize(n);
    output.void_fraction.resize(n);
    output.mixture_density.resize(n);
    output.gas_density.resize(n);
    output.liquid_density.resize(n);
    output.mixture_velocity.resize(n);
    output.gas_velocity.resize(n);
    output.liquid_velocity.resize(n);
    output.temperature.resize(n);
    
    for (int i = 0; i < n; ++i) {
        const auto& cell = grid_->cell(i);
        const auto& vars = state_->at(i);
        
        output.depth[i] = cell.z_center;
        output.pressure[i] = vars.pressure;
        output.void_fraction[i] = vars.void_fraction;
        output.mixture_density[i] = vars.mixture_density;
        output.gas_density[i] = vars.gas_density;
        output.liquid_density[i] = vars.liquid_density;
        output.mixture_velocity[i] = vars.mixture_velocity;
        output.gas_velocity[i] = vars.gas_velocity;
        output.liquid_velocity[i] = vars.liquid_velocity;
        output.temperature[i] = vars.temperature;
    }
    
    return output;
}

void Simulator::enable_choke_valve_control(const ChokeValveConfig& choke_config) {
    enable_choke_control_ = true;
    choke_controller_ = std::make_unique<ChokeValveController>(choke_config);
    config_.choke_config = choke_config;
}

void Simulator::disable_choke_valve_control() {
    enable_choke_control_ = false;
    choke_controller_.reset();
}

void Simulator::set_choke_target_bhp(double target_bhp) {
    if (choke_controller_) {
        choke_controller_->set_target_bhp(target_bhp);
    }
    config_.choke_config.target_bhp = target_bhp;
}

void Simulator::update_choke_control() {
    if (!enable_choke_control_ || !choke_controller_) return;
    
    int n = config_.num_cells;
    if (n < 2) return;
    
    double bhp = state_->at(n - 1).pressure;
    
    double wellhead_pressure = state_->at(0).pressure;
    double wellhead_void_fraction = state_->at(0).void_fraction;
    double wellhead_velocity = state_->at(0).mixture_velocity;
    
    double wellbore_area = M_PI * config_.wellbore_diameter * config_.wellbore_diameter / 4.0;
    double wellhead_flow_rate = std::abs(wellhead_velocity) * wellbore_area;
    
    double choke_area = 0.0;
    double backpressure = 0.0;
    
    choke_controller_->update(
        current_time_,
        bhp,
        wellhead_flow_rate,
        wellhead_void_fraction,
        wellhead_pressure,
        choke_area,
        backpressure
    );
    
    double ambient_pressure = 1e5;
    bc_.outlet_pressure = ambient_pressure + backpressure;
}

void Simulator::print_choke_control_report() const {
    if (choke_controller_) {
        choke_controller_->print_report();
    } else {
        std::cout << "Choke valve control is not enabled." << std::endl;
    }
}

}
