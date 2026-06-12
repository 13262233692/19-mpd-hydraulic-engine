#include "hydraulic_engine/physics.hpp"
#include <cmath>
#include <algorithm>

namespace hydraulic_engine {

PhysicsModel::PhysicsModel(const FluidProperties& fluid_props,
                           const DriftFluxParameters& drift_params)
    : fluid_props_(fluid_props),
      drift_params_(drift_params) {}

double PhysicsModel::compute_mixture_density(double void_fraction,
                                             double gas_density,
                                             double liquid_density) const {
    return void_fraction * gas_density + (1.0 - void_fraction) * liquid_density;
}

double PhysicsModel::compute_gas_density(double pressure, double temperature) const {
    double Z = 1.0;
    double p = std::max(pressure, fluid_props_.reference_pressure * 0.1);
    return (fluid_props_.molar_mass_gas * p) / 
           (Z * fluid_props_.gas_constant * temperature);
}

double PhysicsModel::compute_liquid_density(double pressure, double temperature) const {
    double dp = pressure - fluid_props_.reference_pressure;
    double dt = temperature - fluid_props_.reference_temperature;
    double compressibility = 4.5e-10;
    double thermal_expansion = 2.1e-4;
    return fluid_props_.liquid_density_ref * (1.0 + compressibility * dp - thermal_expansion * dt);
}

double PhysicsModel::compute_slip_velocity(double void_fraction,
                                           double mixture_velocity,
                                           double mixture_density,
                                           double gas_density,
                                           double liquid_density) const {
    if (void_fraction < 1e-12) {
        return 0.0;
    }
    
    double v_drift = compute_drift_velocity(void_fraction, mixture_density, 
                                            gas_density, liquid_density);
    
    double C0 = drift_params_.distribution_coefficient;
    
    double slip = (C0 - 1.0) * mixture_velocity + v_drift / void_fraction;
    
    return std::max(0.0, slip);
}

double PhysicsModel::compute_drift_velocity(double void_fraction,
                                            double mixture_density,
                                            double gas_density,
                                            double liquid_density) const {
    double g = drift_params_.gravity;
    double sigma = fluid_props_.surface_tension;
    
    double rho_diff = std::max(liquid_density - gas_density, 1.0);
    double rho_m = std::max(mixture_density, 1.0);
    
    double terminal_velocity = 1.53 * std::pow(g * sigma * rho_diff / (rho_m * rho_m), 0.25);
    
    double alpha = std::min(std::max(void_fraction, 0.0), 0.999);
    double hindrance_factor = std::pow(1.0 - alpha, 2.0);
    
    return drift_params_.drift_velocity * hindrance_factor * terminal_velocity;
}

void PhysicsModel::update_phase_velocities(PrimitiveVariables& vars,
                                           double slip_velocity) const {
    double alpha = vars.void_fraction;
    double vm = vars.mixture_velocity;
    
    if (alpha < 1e-12) {
        vars.gas_velocity = vm + slip_velocity;
        vars.liquid_velocity = vm;
    } else if (alpha > 0.999) {
        vars.gas_velocity = vm;
        vars.liquid_velocity = vm - slip_velocity;
    } else {
        vars.gas_velocity = vm + (1.0 - alpha) * slip_velocity;
        vars.liquid_velocity = vm - alpha * slip_velocity;
    }
}

double PhysicsModel::compute_friction_factor(double reynolds_number,
                                             double relative_roughness) const {
    double Re = std::abs(reynolds_number);
    
    if (Re < 1.0) {
        return 16.0;
    }
    
    if (Re < 2300.0) {
        return 16.0 / Re;
    }
    
    double f = 0.02;
    for (int iter = 0; iter < 50; ++iter) {
        double sqrt_f = std::sqrt(f);
        double rhs = -2.0 * std::log10(relative_roughness / 3.7 + 2.51 / (Re * sqrt_f));
        double f_new = 1.0 / (rhs * rhs);
        
        if (std::abs(f_new - f) < 1e-10) {
            f = f_new;
            break;
        }
        f = 0.5 * (f + f_new);
    }
    
    return f;
}

double PhysicsModel::compute_reynolds_number(double mixture_density,
                                             double mixture_velocity,
                                             double hydraulic_diameter,
                                             double mixture_viscosity) const {
    if (mixture_viscosity < eps_) {
        return 0.0;
    }
    return (mixture_density * std::abs(mixture_velocity) * hydraulic_diameter) / 
           mixture_viscosity;
}

double PhysicsModel::compute_mixture_viscosity(double void_fraction,
                                               double gas_viscosity,
                                               double liquid_viscosity) const {
    double alpha = std::min(std::max(void_fraction, 0.0), 0.999);
    double mu_r = gas_viscosity / liquid_viscosity;
    
    double numerator = mu_r * (1.0 + 2.5 * alpha) + (1.0 - alpha);
    double denominator = (1.0 - alpha) + 2.5 * alpha * mu_r;
    
    return liquid_viscosity * numerator / denominator;
}

double PhysicsModel::compute_hydrostatic_gradient(double mixture_density,
                                                  double inclination_angle) const {
    return mixture_density * drift_params_.gravity * std::cos(inclination_angle);
}

double PhysicsModel::compute_friction_gradient(double mixture_density,
                                               double mixture_velocity,
                                               double friction_factor,
                                               double hydraulic_diameter) const {
    if (hydraulic_diameter < eps_) {
        return 0.0;
    }
    return friction_factor * mixture_density * mixture_velocity * 
           std::abs(mixture_velocity) / (2.0 * hydraulic_diameter);
}

double PhysicsModel::compute_acceleration_gradient(double mixture_density,
                                                   double mixture_velocity,
                                                   double dvelocity_dz) const {
    return mixture_density * mixture_velocity * dvelocity_dz;
}

BoundaryConditions::BoundaryConditions()
    : inlet_pressure(101325.0),
      outlet_pressure(101325.0),
      inlet_mass_flow_rate(0.0),
      inlet_void_fraction(0.0),
      inlet_temperature(288.15),
      outlet_temperature(288.15),
      bottom_hole_pressure(1e7),
      kick_flow_rate(0.0),
      kick_void_fraction(0.8),
      kick_depth(3000.0),
      use_pressure_inlet(false),
      use_pressure_outlet(true),
      enable_kick(false),
      kick_start_time(0.0),
      kick_duration(1e6) {}

}
