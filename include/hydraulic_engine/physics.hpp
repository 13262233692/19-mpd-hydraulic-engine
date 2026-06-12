#pragma once

#include <cmath>
#include "variables.hpp"

namespace hydraulic_engine {

struct FluidProperties {
    double liquid_density_ref;
    double gas_density_ref;
    double liquid_viscosity;
    double gas_viscosity;
    double surface_tension;
    double gas_constant;
    double molar_mass_gas;
    double reference_pressure;
    double reference_temperature;
};

struct DriftFluxParameters {
    double distribution_coefficient;
    double drift_velocity;
    double gravity;
};

class PhysicsModel {
public:
    PhysicsModel(const FluidProperties& fluid_props,
                 const DriftFluxParameters& drift_params);
    
    double compute_mixture_density(double void_fraction,
                                   double gas_density,
                                   double liquid_density) const;
    
    double compute_gas_density(double pressure, double temperature) const;
    
    double compute_liquid_density(double pressure, double temperature) const;
    
    double compute_slip_velocity(double void_fraction,
                                 double mixture_velocity,
                                 double mixture_density,
                                 double gas_density,
                                 double liquid_density) const;
    
    double compute_drift_velocity(double void_fraction,
                                  double mixture_density,
                                  double gas_density,
                                  double liquid_density) const;
    
    void update_phase_velocities(PrimitiveVariables& vars,
                                 double slip_velocity) const;
    
    double compute_friction_factor(double reynolds_number,
                                   double relative_roughness) const;
    
    double compute_reynolds_number(double mixture_density,
                                   double mixture_velocity,
                                   double hydraulic_diameter,
                                   double mixture_viscosity) const;
    
    double compute_mixture_viscosity(double void_fraction,
                                     double gas_viscosity,
                                     double liquid_viscosity) const;
    
    double compute_hydrostatic_gradient(double mixture_density,
                                        double inclination_angle = 0.0) const;
    
    double compute_friction_gradient(double mixture_density,
                                     double mixture_velocity,
                                     double friction_factor,
                                     double hydraulic_diameter) const;
    
    double compute_acceleration_gradient(double mixture_density,
                                         double mixture_velocity,
                                         double dvelocity_dz) const;
    
    const FluidProperties& fluid_properties() const { return fluid_props_; }
    const DriftFluxParameters& drift_parameters() const { return drift_params_; }

private:
    FluidProperties fluid_props_;
    DriftFluxParameters drift_params_;
    
    static constexpr double eps_ = 1e-12;
};

class BoundaryConditions {
public:
    BoundaryConditions();
    
    double inlet_pressure;
    double outlet_pressure;
    double inlet_mass_flow_rate;
    double inlet_void_fraction;
    double inlet_temperature;
    double outlet_temperature;
    
    double bottom_hole_pressure;
    double kick_flow_rate;
    double kick_void_fraction;
    double kick_depth;
    
    bool use_pressure_inlet;
    bool use_pressure_outlet;
    bool enable_kick;
    double kick_start_time;
    double kick_duration;
};

}
