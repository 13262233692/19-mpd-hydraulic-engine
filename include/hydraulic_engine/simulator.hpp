#pragma once

#include <memory>
#include <vector>
#include <string>
#include <functional>
#include "grid.hpp"
#include "variables.hpp"
#include "physics.hpp"
#include "fvm.hpp"
#include "time_integrator.hpp"
#include "pentadiagonal_solver.hpp"

namespace hydraulic_engine {

struct SimulationConfig {
    double total_depth;
    int num_cells;
    double wellbore_diameter;
    
    double initial_pressure;
    double initial_void_fraction;
    double initial_temperature;
    
    double liquid_density_ref;
    double gas_density_ref;
    double liquid_viscosity;
    double gas_viscosity;
    double surface_tension;
    double molar_mass_gas;
    
    double distribution_coefficient;
    double drift_velocity_coeff;
    double gravity;
    
    double dt_init;
    double dt_min;
    double dt_max;
    double cfl;
    int max_newton_iter;
    double newton_tol;
    bool adaptive_time_stepping;
    
    FluxLimiter::Type flux_limiter;
    
    SimulationConfig();
};

struct SimulationOutput {
    double time;
    std::vector<double> depth;
    std::vector<double> pressure;
    std::vector<double> void_fraction;
    std::vector<double> mixture_density;
    std::vector<double> gas_density;
    std::vector<double> liquid_density;
    std::vector<double> mixture_velocity;
    std::vector<double> gas_velocity;
    std::vector<double> liquid_velocity;
    std::vector<double> temperature;
};

class Simulator {
public:
    explicit Simulator(const SimulationConfig& config);
    
    void initialize();
    
    bool run(double total_time, int output_interval = 10);
    
    bool step(double& current_time, double& dt);
    
    void apply_boundary_conditions(const BoundaryConditions& bc);
    
    void enable_kick(double kick_flow_rate, double kick_void_fraction,
                     double kick_depth, double kick_start_time = 0.0,
                     double kick_duration = 1e6);
    
    void disable_kick();
    
    const SimulationOutput get_output() const;
    
    double get_current_time() const { return current_time_; }
    double get_current_dt() const { return current_dt_; }
    int get_total_steps() const { return total_steps_; }
    
    const Grid& grid() const { return *grid_; }
    const FieldState& state() const { return *state_; }
    const BoundaryConditions& boundary_conditions() const { return bc_; }
    
    int get_last_newton_iterations() const { 
        return time_integrator_ ? time_integrator_->last_newton_iterations() : 0; 
    }
    
    double get_last_newton_residual() const { 
        return time_integrator_ ? time_integrator_->last_newton_residual() : 0.0; 
    }
    
    void set_output_callback(std::function<void(const SimulationOutput&)> callback) {
        output_callback_ = callback;
    }

private:
    SimulationConfig config_;
    
    std::unique_ptr<Grid> grid_;
    std::unique_ptr<PhysicsModel> physics_;
    std::unique_ptr<FVMDisc> fvm_;
    std::unique_ptr<TimeIntegrator> time_integrator_;
    std::unique_ptr<FieldState> state_;
    std::unique_ptr<FaceVariables> face_fluxes_;
    
    BoundaryConditions bc_;
    
    double current_time_;
    double current_dt_;
    int total_steps_;
    
    std::vector<SimulationOutput> output_history_;
    std::function<void(const SimulationOutput&)> output_callback_;
    
    void setup_default_boundary_conditions();
    void validate_config() const;
};

}
