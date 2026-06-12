#include "hydraulic_engine/time_integrator.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace hydraulic_engine {

TimeIntegrator::TimeIntegrator(const Grid& grid,
                               const PhysicsModel& physics,
                               const FVMDisc& fvm,
                               const TimeIntegratorParams& params,
                               Scheme scheme)
    : grid_(grid),
      physics_(physics),
      fvm_(fvm),
      params_(params),
      scheme_(scheme),
      solver_(grid.num_cells()),
      state_old_(grid.num_cells()),
      state_prev_(grid.num_cells()),
      residuals_(grid.num_cells()),
      face_fluxes_(grid.num_faces()),
      last_newton_iter_(0),
      last_newton_residual_(0.0),
      total_time_steps_(0) {
    int n = grid.num_cells();
    a_.resize(n);
    b_.resize(n);
    c_.resize(n);
    d_.resize(n);
    e_.resize(n);
    rhs_.resize(n);
    delta_rho_.resize(n);
    delta_alpha_.resize(n);
    delta_p_.resize(n);
}

double TimeIntegrator::compute_suggested_dt(const FieldState& state) const {
    double max_velocity = 0.0;
    double max_speed_of_sound = 1500.0;
    
    for (int i = 0; i < grid_.num_cells(); ++i) {
        double vm = std::abs(state.at(i).mixture_velocity);
        double vg = std::abs(state.at(i).gas_velocity);
        double vl = std::abs(state.at(i).liquid_velocity);
        max_velocity = std::max(max_velocity, std::max({vm, vg, vl}));
    }
    
    double min_dz = grid_.cell(0).dz;
    for (int i = 1; i < grid_.num_cells(); ++i) {
        min_dz = std::min(min_dz, grid_.cell(i).dz);
    }
    
    double dt_cfl = params_.cfl * min_dz / (max_velocity + max_speed_of_sound + eps_);
    dt_cfl = std::max(params_.dt_min, std::min(params_.dt_max, dt_cfl));
    
    return dt_cfl;
}

void TimeIntegrator::update_physical_properties(FieldState& state) {
    for (int i = 0; i < grid_.num_cells(); ++i) {
        auto& vars = state.at(i);
        
        vars.gas_density = physics_.compute_gas_density(vars.pressure, vars.temperature);
        vars.liquid_density = physics_.compute_liquid_density(vars.pressure, vars.temperature);
        
        vars.void_fraction = std::max(0.0, std::min(0.999, vars.void_fraction));
        
        vars.mixture_density = physics_.compute_mixture_density(
            vars.void_fraction, vars.gas_density, vars.liquid_density
        );
        
        double slip = physics_.compute_slip_velocity(
            vars.void_fraction, vars.mixture_velocity,
            vars.mixture_density, vars.gas_density, vars.liquid_density
        );
        
        physics_.update_phase_velocities(vars, slip);
    }
}

double TimeIntegrator::compute_residual_norm(const FieldState& residuals) const {
    double norm_mass = 0.0;
    double norm_mass_gas = 0.0;
    double norm_mom = 0.0;
    
    for (int i = 0; i < grid_.num_cells(); ++i) {
        norm_mass += residuals.conserved(i).mass_mixture * residuals.conserved(i).mass_mixture;
        norm_mass_gas += residuals.conserved(i).mass_gas * residuals.conserved(i).mass_gas;
        norm_mom += residuals.conserved(i).momentum_mixture * residuals.conserved(i).momentum_mixture;
    }
    
    int n = grid_.num_cells();
    return std::sqrt((norm_mass + norm_mass_gas + norm_mom) / (3.0 * n));
}

void TimeIntegrator::compute_newton_update(const FieldState& current_state,
                                           FieldState& updated_state,
                                           const std::vector<double>& delta_rho,
                                           const std::vector<double>& delta_alpha,
                                           const std::vector<double>& delta_p) {
    for (int i = 0; i < grid_.num_cells(); ++i) {
        auto& vars = updated_state.at(i);
        const auto& curr = current_state.at(i);
        
        vars.mixture_density = curr.mixture_density + delta_rho[i];
        vars.void_fraction = curr.void_fraction + delta_alpha[i];
        vars.pressure = curr.pressure + delta_p[i];
        
        vars.void_fraction = std::max(0.0, std::min(0.999, vars.void_fraction));
        vars.pressure = std::max(1e4, vars.pressure);
        
        vars.gas_density = physics_.compute_gas_density(vars.pressure, vars.temperature);
        vars.liquid_density = physics_.compute_liquid_density(vars.pressure, vars.temperature);
        
        double alpha = vars.void_fraction;
        double rho_g = vars.gas_density;
        double rho_l = vars.liquid_density;
        double rho_m = vars.mixture_density;
        
        if (alpha < 1.0 - 1e-12) {
            double computed_rho_m = alpha * rho_g + (1.0 - alpha) * rho_l;
            if (std::abs(computed_rho_m - rho_m) > 1.0) {
                vars.mixture_density = computed_rho_m;
            }
        }
    }
}

void TimeIntegrator::apply_under_relaxation(FieldState& new_state,
                                            const FieldState& old_state,
                                            double relaxation_factor) {
    for (int i = 0; i < grid_.num_cells(); ++i) {
        auto& new_var = new_state.at(i);
        const auto& old_var = old_state.at(i);
        
        new_var.mixture_density = relaxation_factor * new_var.mixture_density + 
                                  (1.0 - relaxation_factor) * old_var.mixture_density;
        new_var.void_fraction = relaxation_factor * new_var.void_fraction + 
                                (1.0 - relaxation_factor) * old_var.void_fraction;
        new_var.pressure = relaxation_factor * new_var.pressure + 
                           (1.0 - relaxation_factor) * old_var.pressure;
        new_var.mixture_velocity = relaxation_factor * new_var.mixture_velocity + 
                                   (1.0 - relaxation_factor) * old_var.mixture_velocity;
    }
}

bool TimeIntegrator::solve_linear_system(FieldState& state,
                                         const BoundaryConditions& bc,
                                         double current_time,
                                         double dt) {
    int n = grid_.num_cells();
    solver_.resize(n);
    
    fvm_.compute_jacobian_vectors(state, bc, a_, b_, c_, d_, e_, rhs_, dt, current_time, 0);
    if (!solver_.solve_modified_thomas(a_, b_, c_, d_, e_, rhs_, delta_rho_)) {
        return false;
    }
    
    fvm_.compute_jacobian_vectors(state, bc, a_, b_, c_, d_, e_, rhs_, dt, current_time, 1);
    if (!solver_.solve_modified_thomas(a_, b_, c_, d_, e_, rhs_, delta_alpha_)) {
        return false;
    }
    
    fvm_.compute_jacobian_vectors(state, bc, a_, b_, c_, d_, e_, rhs_, dt, current_time, 2);
    if (!solver_.solve_modified_thomas(a_, b_, c_, d_, e_, rhs_, delta_p_)) {
        return false;
    }
    
    return true;
}

bool TimeIntegrator::solve_nonlinear(FieldState& state,
                                     const BoundaryConditions& bc,
                                     double current_time,
                                     double dt) {
    FieldState state_guess = state;
    
    last_newton_iter_ = 0;
    last_newton_residual_ = 0.0;
    
    for (int i = 0; i < grid_.num_cells(); ++i) {
        state_old_.primitive_to_conservative(i, grid_.cell(i).volume);
    }
    
    for (int iter = 0; iter < params_.max_newton_iter; ++iter) {
        last_newton_iter_ = iter + 1;
        
        update_physical_properties(state_guess);
        fvm_.apply_boundary_conditions(state_guess, bc, current_time);
        
        for (int i = 0; i < grid_.num_cells(); ++i) {
            state_guess.primitive_to_conservative(i, grid_.cell(i).volume);
        }
        
        fvm_.compute_face_fluxes(state_guess, face_fluxes_, dt);
        fvm_.compute_residuals(state_guess, face_fluxes_, bc, residuals_, current_time);
        
        for (int i = 0; i < grid_.num_cells(); ++i) {
            double vol = grid_.cell(i).volume;
            residuals_.conserved(i).mass_mixture = 
                (state_guess.conserved(i).mass_mixture - state_old_.conserved(i).mass_mixture) / dt
                - residuals_.conserved(i).mass_mixture;
            
            residuals_.conserved(i).mass_gas = 
                (state_guess.conserved(i).mass_gas - state_old_.conserved(i).mass_gas) / dt
                - residuals_.conserved(i).mass_gas;
            
            residuals_.conserved(i).momentum_mixture = 
                (state_guess.conserved(i).momentum_mixture - state_old_.conserved(i).momentum_mixture) / dt
                - residuals_.conserved(i).momentum_mixture;
        }
        
        last_newton_residual_ = compute_residual_norm(residuals_);
        
        if (last_newton_residual_ < params_.newton_tol) {
            state = state_guess;
            return true;
        }
        
        if (iter == params_.max_newton_iter - 1) {
            break;
        }
        
        fvm_.compute_face_fluxes(state_guess, face_fluxes_, dt);
        fvm_.compute_residuals(state_guess, face_fluxes_, bc, residuals_, current_time);
        
        double relax = 0.8;
        double dt_eff = dt * relax;
        
        for (int i = 0; i < grid_.num_cells(); ++i) {
            double vol = grid_.cell(i).volume;
            
            double dmass_m = state_guess.conserved(i).mass_mixture + dt_eff * residuals_.conserved(i).mass_mixture;
            double dmass_g = state_guess.conserved(i).mass_gas + dt_eff * residuals_.conserved(i).mass_gas;
            double dmom = state_guess.conserved(i).momentum_mixture + dt_eff * residuals_.conserved(i).momentum_mixture;
            
            double new_rho_m = dmass_m / vol;
            double new_alpha_g = dmass_g / (state_guess.at(i).gas_density * vol);
            new_alpha_g = std::max(0.0, std::min(0.999, new_alpha_g));
            
            double new_vm = 0.0;
            if (std::abs(new_rho_m) > 1e-12) {
                new_vm = dmom / (new_rho_m * vol);
            }
            
            double new_rho_l = state_guess.at(i).liquid_density;
            if (new_alpha_g < 0.999) {
                new_rho_l = (new_rho_m - new_alpha_g * state_guess.at(i).gas_density) / (1.0 - new_alpha_g);
            }
            
            double final_relax = 0.5;
            state_guess.at(i).mixture_density = final_relax * new_rho_m + (1.0 - final_relax) * state_guess.at(i).mixture_density;
            state_guess.at(i).void_fraction = final_relax * new_alpha_g + (1.0 - final_relax) * state_guess.at(i).void_fraction;
            state_guess.at(i).mixture_velocity = final_relax * new_vm + (1.0 - final_relax) * state_guess.at(i).mixture_velocity;
            state_guess.at(i).liquid_density = new_rho_l;
        }
        
        update_physical_properties(state_guess);
        fvm_.apply_boundary_conditions(state_guess, bc, current_time);
    }
    
    state = state_guess;
    return last_newton_residual_ < params_.newton_tol * 1000.0;
}

bool TimeIntegrator::step(FieldState& state,
                          const BoundaryConditions& bc,
                          double& current_time,
                          double& dt) {
    int n = grid_.num_cells();
    
    if (params_.adaptive_time_stepping) {
        double suggested_dt = compute_suggested_dt(state);
        dt = std::max(params_.dt_min, std::min(params_.dt_max, suggested_dt));
    }
    
    state_old_ = state;
    
    for (int i = 0; i < n; ++i) {
        state_old_.primitive_to_conservative(i, grid_.cell(i).volume);
    }
    
    bool success = solve_nonlinear(state, bc, current_time, dt);
    
    if (success) {
        current_time += dt;
        total_time_steps_++;
        
        update_physical_properties(state);
        fvm_.apply_boundary_conditions(state, bc, current_time);
    }
    
    return success;
}

}
