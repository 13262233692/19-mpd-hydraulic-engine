#include "hydraulic_engine/fvm.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace hydraulic_engine {

FVMDisc::FVMDisc(const Grid& grid,
                 const PhysicsModel& physics,
                 FluxLimiter::Type limiter_type)
    : grid_(grid),
      physics_(physics),
      limiter_type_(limiter_type),
      relative_roughness_(1e-5) {}

double FVMDisc::minmod(double a, double b) const {
    if (a * b <= 0.0) return 0.0;
    return std::abs(a) < std::abs(b) ? a : b;
}

double FVMDisc::superbee(double r) const {
    if (r <= 0.0) return 0.0;
    if (r <= 0.5) return 2.0 * r;
    if (r <= 1.0) return 1.0;
    return std::min(2.0, std::min(2.0 * r, (1.0 + r) / 2.0));
}

double FVMDisc::van_leer(double r) const {
    if (r <= 0.0) return 0.0;
    return 2.0 * r / (1.0 + r);
}

void FVMDisc::first_order_upwind(const PrimitiveVariables& left,
                                 const PrimitiveVariables& right,
                                 const double face_area,
                                 double& mass_flux_mixture,
                                 double& mass_flux_gas,
                                 double& momentum_flux) const {
    double velocity_left = left.mixture_velocity;
    double velocity_right = right.mixture_velocity;
    
    double face_velocity = 0.5 * (velocity_left + velocity_right);
    
    const PrimitiveVariables& upwind = face_velocity >= 0.0 ? left : right;
    
    double rho_m = upwind.mixture_density;
    double alpha_g = upwind.void_fraction;
    double rho_g = upwind.gas_density;
    double vm = face_velocity;
    
    mass_flux_mixture = rho_m * vm * face_area;
    mass_flux_gas = alpha_g * rho_g * upwind.gas_velocity * face_area;
    
    double pressure_central = 0.5 * (left.pressure + right.pressure);
    momentum_flux = (rho_m * vm * vm + pressure_central) * face_area;
}

void FVMDisc::muscl_reconstruction(const FieldState& state,
                                   int face_idx,
                                   PrimitiveVariables& left,
                                   PrimitiveVariables& right) const {
    int num_cells = grid_.num_cells();
    
    if (face_idx <= 0 || face_idx >= num_cells) {
        if (face_idx == 0) {
            left = state.at(0);
            right = state.at(0);
        } else {
            left = state.at(num_cells - 1);
            right = state.at(num_cells - 1);
        }
        return;
    }
    
    int i = face_idx - 1;
    
    int im1 = std::max(0, i - 1);
    int ip1 = std::min(num_cells - 1, i + 1);
    int ip2 = std::min(num_cells - 1, i + 2);
    
    const auto& var_im1 = state.at(im1);
    const auto& var_i = state.at(i);
    const auto& var_ip1 = state.at(ip1);
    const auto& var_ip2 = state.at(ip2);
    
    auto reconstruct = [&](double vm1, double vi, double vp1, double vp2) -> double {
        double delta_left = vi - vm1;
        double delta_right = vp1 - vi;
        double slope_i = minmod(delta_right, delta_left);
        
        double delta_left_p1 = vp1 - vi;
        double delta_right_p1 = vp2 - vp1;
        double slope_ip1 = minmod(delta_right_p1, delta_left_p1);
        
        double left_val = vi + 0.5 * slope_i;
        double right_val = vp1 - 0.5 * slope_ip1;
        
        return (vi + vp1) / 2.0;
    };
    
    left = var_i;
    right = var_ip1;
    
    left.pressure = reconstruct(var_im1.pressure, var_i.pressure, var_ip1.pressure, var_ip2.pressure);
    right.pressure = reconstruct(var_i.pressure, var_ip1.pressure, var_ip2.pressure, var_ip2.pressure);
    
    left.mixture_density = reconstruct(var_im1.mixture_density, var_i.mixture_density, var_ip1.mixture_density, var_ip2.mixture_density);
    right.mixture_density = reconstruct(var_i.mixture_density, var_ip1.mixture_density, var_ip2.mixture_density, var_ip2.mixture_density);
    
    left.void_fraction = reconstruct(var_im1.void_fraction, var_i.void_fraction, var_ip1.void_fraction, var_ip2.void_fraction);
    right.void_fraction = reconstruct(var_i.void_fraction, var_ip1.void_fraction, var_ip2.void_fraction, var_ip2.void_fraction);
    
    left.void_fraction = std::max(0.0, std::min(0.999, left.void_fraction));
    right.void_fraction = std::max(0.0, std::min(0.999, right.void_fraction));
}

void FVMDisc::compute_face_fluxes(const FieldState& state,
                                  FaceVariables& face_fluxes,
                                  double dt) const {
    int num_faces = grid_.num_faces();
    int num_cells = grid_.num_cells();
    
    for (int f = 0; f < num_faces; ++f) {
        const auto& face = grid_.face(f);
        
        PrimitiveVariables left, right;
        
        if (limiter_type_ == FluxLimiter::FIRST_ORDER_UPWIND) {
            if (f == 0) {
                left = state.at(0);
                right = state.at(0);
                if (num_cells > 1) {
                    double dp = state.at(1).pressure - state.at(0).pressure;
                    left.pressure = state.at(0).pressure - dp;
                }
            } else if (f == num_faces - 1) {
                left = state.at(num_cells - 1);
                right = state.at(num_cells - 1);
                if (num_cells > 1) {
                    double dp = state.at(num_cells - 1).pressure - state.at(num_cells - 2).pressure;
                    right.pressure = state.at(num_cells - 1).pressure + dp;
                }
            } else {
                left = state.at(f - 1);
                right = state.at(f);
            }
        } else {
            muscl_reconstruction(state, f, left, right);
        }
        
        double mass_flux_m, mass_flux_g, mom_flux;
        first_order_upwind(left, right, face.area, mass_flux_m, mass_flux_g, mom_flux);
        
        face_fluxes.mass_flux_mixture(f) = mass_flux_m;
        face_fluxes.mass_flux_gas(f) = mass_flux_g;
        face_fluxes.momentum_flux(f) = mom_flux;
    }
}

double FVMDisc::compute_source_term(const PrimitiveVariables& vars,
                                    const GridCell& cell,
                                    const PhysicsModel& physics,
                                    double relative_roughness) const {
    double hydraulic_diameter = 4.0 * cell.area_top / (M_PI * 2.0 * grid_.wellbore_radius());
    
    double mu_mix = physics.compute_mixture_viscosity(
        vars.void_fraction,
        physics.fluid_properties().gas_viscosity,
        physics.fluid_properties().liquid_viscosity
    );
    
    double Re = physics.compute_reynolds_number(
        vars.mixture_density,
        vars.mixture_velocity,
        hydraulic_diameter,
        mu_mix
    );
    
    double f = physics.compute_friction_factor(Re, relative_roughness);
    
    double dp_dz_hydro = physics.compute_hydrostatic_gradient(vars.mixture_density);
    double dp_dz_fric = physics.compute_friction_gradient(
        vars.mixture_density,
        vars.mixture_velocity,
        f,
        hydraulic_diameter
    );
    
    double source = -(dp_dz_hydro + dp_dz_fric) * cell.area_top * cell.dz;
    
    return source;
}

void FVMDisc::apply_kick_source(FieldState& residuals,
                                const FieldState& state,
                                const BoundaryConditions& bc,
                                double current_time) const {
    if (!bc.enable_kick) return;
    if (current_time < bc.kick_start_time) return;
    if (current_time > bc.kick_start_time + bc.kick_duration) return;
    
    double kick_strength = 1.0;
    double elapsed = current_time - bc.kick_start_time;
    double ramp_up = std::min(1.0, elapsed / 10.0);
    double ramp_down = 1.0;
    if (elapsed > bc.kick_duration - 10.0) {
        ramp_down = std::max(0.0, (bc.kick_duration - elapsed) / 10.0);
    }
    kick_strength = ramp_up * ramp_down;
    
    int num_cells = grid_.num_cells();
    double dz = grid_.total_depth() / num_cells;
    int kick_cell_idx = static_cast<int>(bc.kick_depth / dz);
    kick_cell_idx = std::max(0, std::min(num_cells - 1, kick_cell_idx));
    
    const auto& cell = grid_.cell(kick_cell_idx);
    double kick_mass_flow = bc.kick_flow_rate * kick_strength;
    
    double rho_g = physics_.compute_gas_density(
        state.at(kick_cell_idx).pressure,
        state.at(kick_cell_idx).temperature
    );
    double rho_l = physics_.fluid_properties().liquid_density_ref;
    
    double mass_gas_in = bc.kick_void_fraction * rho_g * kick_mass_flow;
    double mass_liquid_in = (1.0 - bc.kick_void_fraction) * rho_l * kick_mass_flow;
    double total_mass_in = mass_gas_in + mass_liquid_in;
    
    residuals.conserved(kick_cell_idx).mass_mixture += total_mass_in;
    residuals.conserved(kick_cell_idx).mass_gas += mass_gas_in;
    residuals.conserved(kick_cell_idx).momentum_mixture += total_mass_in * 1.0;
}

void FVMDisc::compute_residuals(const FieldState& state,
                                const FaceVariables& face_fluxes,
                                const BoundaryConditions& bc,
                                FieldState& residuals,
                                double current_time) const {
    int num_cells = grid_.num_cells();
    
    for (int i = 0; i < num_cells; ++i) {
        const auto& cell = grid_.cell(i);
        
        double flux_in_mass_m = face_fluxes.mass_flux_mixture(i);
        double flux_out_mass_m = face_fluxes.mass_flux_mixture(i + 1);
        double flux_in_mass_g = face_fluxes.mass_flux_gas(i);
        double flux_out_mass_g = face_fluxes.mass_flux_gas(i + 1);
        double flux_in_mom = face_fluxes.momentum_flux(i);
        double flux_out_mom = face_fluxes.momentum_flux(i + 1);
        
        double source = compute_source_term(state.at(i), cell, physics_, relative_roughness_);
        
        residuals.conserved(i).mass_mixture = -(flux_out_mass_m - flux_in_mass_m);
        residuals.conserved(i).mass_gas = -(flux_out_mass_g - flux_in_mass_g);
        residuals.conserved(i).momentum_mixture = -(flux_out_mom - flux_in_mom) + source;
    }
    
    apply_kick_source(residuals, state, bc, current_time);
}

void FVMDisc::apply_boundary_conditions(FieldState& state,
                                        const BoundaryConditions& bc,
                                        double current_time) const {
    int num_cells = grid_.num_cells();
    
    if (bc.use_pressure_inlet) {
        state.at(0).pressure = bc.inlet_pressure;
        state.at(0).temperature = bc.inlet_temperature;
    } else {
        if (num_cells > 1) {
            state.at(0).mixture_velocity = bc.inlet_mass_flow_rate / 
                (state.at(0).mixture_density * grid_.face(0).area);
        }
    }
    
    if (bc.use_pressure_outlet) {
        state.at(num_cells - 1).pressure = bc.outlet_pressure;
        state.at(num_cells - 1).temperature = bc.outlet_temperature;
    }
    
    if (num_cells > 1) {
        state.at(0).void_fraction = bc.inlet_void_fraction;
    }
}

void FVMDisc::compute_jacobian_vectors(const FieldState& state,
                                       const BoundaryConditions& bc,
                                       std::vector<double>& a,
                                       std::vector<double>& b,
                                       std::vector<double>& c,
                                       std::vector<double>& d,
                                       std::vector<double>& e,
                                       std::vector<double>& rhs,
                                       double dt,
                                       double current_time,
                                       int equation_idx) const {
    int n = grid_.num_cells();
    
    a.assign(n, 0.0);
    b.assign(n, 0.0);
    c.assign(n, 0.0);
    d.assign(n, 0.0);
    e.assign(n, 0.0);
    rhs.assign(n, 0.0);
    
    double eps = 1e-8;
    
    FieldState state_perturbed = state;
    FieldState residuals_plus(n);
    FieldState residuals_minus(n);
    FaceVariables face_fluxes_plus(grid_.num_faces());
    FaceVariables face_fluxes_minus(grid_.num_faces());
    
    for (int i = 0; i < n; ++i) {
        double var = 0.0;
        if (equation_idx == 0) var = state.at(i).mixture_density;
        else if (equation_idx == 1) var = state.at(i).void_fraction;
        else var = state.at(i).pressure;
        
        double pert = eps * std::max(1.0, std::abs(var));
        
        if (equation_idx == 0) {
            state_perturbed.at(i).mixture_density = var + pert;
            compute_face_fluxes(state_perturbed, face_fluxes_plus, dt);
            compute_residuals(state_perturbed, face_fluxes_plus, bc, residuals_plus, current_time);
            
            state_perturbed.at(i).mixture_density = var - pert;
            compute_face_fluxes(state_perturbed, face_fluxes_minus, dt);
            compute_residuals(state_perturbed, face_fluxes_minus, bc, residuals_minus, current_time);
            
            state_perturbed.at(i).mixture_density = var;
        } else if (equation_idx == 1) {
            state_perturbed.at(i).void_fraction = std::min(0.999, std::max(0.001, var + pert));
            compute_face_fluxes(state_perturbed, face_fluxes_plus, dt);
            compute_residuals(state_perturbed, face_fluxes_plus, bc, residuals_plus, current_time);
            
            state_perturbed.at(i).void_fraction = std::min(0.999, std::max(0.001, var - pert));
            compute_face_fluxes(state_perturbed, face_fluxes_minus, dt);
            compute_residuals(state_perturbed, face_fluxes_minus, bc, residuals_minus, current_time);
            
            state_perturbed.at(i).void_fraction = var;
        } else {
            state_perturbed.at(i).pressure = var + pert;
            compute_face_fluxes(state_perturbed, face_fluxes_plus, dt);
            compute_residuals(state_perturbed, face_fluxes_plus, bc, residuals_plus, current_time);
            
            state_perturbed.at(i).pressure = var - pert;
            compute_face_fluxes(state_perturbed, face_fluxes_minus, dt);
            compute_residuals(state_perturbed, face_fluxes_minus, bc, residuals_minus, current_time);
            
            state_perturbed.at(i).pressure = var;
        }
        
        for (int j = 0; j < n; ++j) {
            double dR_dvar = 0.0;
            if (equation_idx == 0) {
                dR_dvar = (residuals_plus.conserved(j).mass_mixture - 
                           residuals_minus.conserved(j).mass_mixture) / (2.0 * pert);
            } else if (equation_idx == 1) {
                dR_dvar = (residuals_plus.conserved(j).mass_gas - 
                           residuals_minus.conserved(j).mass_gas) / (2.0 * pert);
            } else {
                dR_dvar = (residuals_plus.conserved(j).momentum_mixture - 
                           residuals_minus.conserved(j).momentum_mixture) / (2.0 * pert);
            }
            
            int diff = j - i;
            if (diff == -2) {
                a[j] = -dR_dvar;
            } else if (diff == -1) {
                b[j] = -dR_dvar;
            } else if (diff == 0) {
                c[j] = 1.0 / dt - dR_dvar;
            } else if (diff == 1) {
                d[j] = -dR_dvar;
            } else if (diff == 2) {
                e[j] = -dR_dvar;
            }
        }
    }
    
    FaceVariables face_fluxes(grid_.num_faces());
    compute_face_fluxes(state, face_fluxes, dt);
    FieldState residuals(n);
    compute_residuals(state, face_fluxes, bc, residuals, current_time);
    
    for (int i = 0; i < n; ++i) {
        if (equation_idx == 0) {
            rhs[i] = -residuals.conserved(i).mass_mixture;
        } else if (equation_idx == 1) {
            rhs[i] = -residuals.conserved(i).mass_gas;
        } else {
            rhs[i] = -residuals.conserved(i).momentum_mixture;
        }
    }
    
    if (bc.use_pressure_inlet) {
        c[0] = 1.0;
        b[0] = 0.0;
        d[0] = 0.0;
        a[0] = 0.0;
        e[0] = 0.0;
        rhs[0] = bc.inlet_pressure - state.at(0).pressure;
    }
    
    if (bc.use_pressure_outlet) {
        c[n - 1] = 1.0;
        b[n - 1] = 0.0;
        d[n - 1] = 0.0;
        a[n - 1] = 0.0;
        e[n - 1] = 0.0;
        rhs[n - 1] = bc.outlet_pressure - state.at(n - 1).pressure;
    }
}

}
