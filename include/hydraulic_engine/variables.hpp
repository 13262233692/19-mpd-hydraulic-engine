#pragma once

#include <vector>
#include <cmath>

namespace hydraulic_engine {

struct PrimitiveVariables {
    double pressure;
    double mixture_density;
    double void_fraction;
    double mixture_velocity;
    double gas_density;
    double liquid_density;
    double gas_velocity;
    double liquid_velocity;
    double temperature;
};

struct ConservativeVariables {
    double mass_mixture;
    double mass_gas;
    double momentum_mixture;
};

class FieldState {
public:
    FieldState(int num_cells);
    
    int num_cells() const { return num_cells_; }
    
    PrimitiveVariables& at(int i) { return primitives_[i]; }
    const PrimitiveVariables& at(int i) const { return primitives_[i]; }
    
    ConservativeVariables& conserved(int i) { return conserved_[i]; }
    const ConservativeVariables& conserved(int i) const { return conserved_[i]; }
    
    std::vector<PrimitiveVariables>& primitives() { return primitives_; }
    const std::vector<PrimitiveVariables>& primitives() const { return primitives_; }
    
    std::vector<ConservativeVariables>& conserved() { return conserved_; }
    const std::vector<ConservativeVariables>& conserved() const { return conserved_; }
    
    void initialize(double initial_pressure, double initial_void_fraction,
                    double liquid_density, double gas_density,
                    double initial_temperature);
    
    void primitive_to_conservative(int i, double cell_volume);
    void conservative_to_primitive(int i, double cell_volume);
    
    void update_velocities(int i, double slip_velocity);

private:
    int num_cells_;
    std::vector<PrimitiveVariables> primitives_;
    std::vector<ConservativeVariables> conserved_;
};

class FaceVariables {
public:
    FaceVariables(int num_faces);
    
    int num_faces() const { return num_faces_; }
    
    double& mass_flux_mixture(int i) { return mass_flux_mixture_[i]; }
    double& mass_flux_gas(int i) { return mass_flux_gas_[i]; }
    double& momentum_flux(int i) { return momentum_flux_[i]; }
    
    const double& mass_flux_mixture(int i) const { return mass_flux_mixture_[i]; }
    const double& mass_flux_gas(int i) const { return mass_flux_gas_[i]; }
    const double& momentum_flux(int i) const { return momentum_flux_[i]; }

private:
    int num_faces_;
    std::vector<double> mass_flux_mixture_;
    std::vector<double> mass_flux_gas_;
    std::vector<double> momentum_flux_;
};

}
