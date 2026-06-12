#include "hydraulic_engine/variables.hpp"
#include <cmath>
#include <stdexcept>

namespace hydraulic_engine {

FieldState::FieldState(int num_cells)
    : num_cells_(num_cells) {
    primitives_.resize(num_cells_);
    conserved_.resize(num_cells_);
}

void FieldState::initialize(double initial_pressure, double initial_void_fraction,
                            double liquid_density, double gas_density,
                            double initial_temperature) {
    for (int i = 0; i < num_cells_; ++i) {
        primitives_[i].pressure = initial_pressure;
        primitives_[i].void_fraction = initial_void_fraction;
        primitives_[i].liquid_density = liquid_density;
        primitives_[i].gas_density = gas_density;
        primitives_[i].temperature = initial_temperature;
        primitives_[i].mixture_velocity = 0.0;
        primitives_[i].gas_velocity = 0.0;
        primitives_[i].liquid_velocity = 0.0;
        primitives_[i].mixture_density = initial_void_fraction * gas_density + 
                                        (1.0 - initial_void_fraction) * liquid_density;
    }
}

void FieldState::primitive_to_conservative(int i, double cell_volume) {
    const auto& p = primitives_[i];
    auto& c = conserved_[i];
    
    c.mass_mixture = p.mixture_density * cell_volume;
    c.mass_gas = p.void_fraction * p.gas_density * cell_volume;
    c.momentum_mixture = p.mixture_density * p.mixture_velocity * cell_volume;
}

void FieldState::conservative_to_primitive(int i, double cell_volume) {
    auto& p = primitives_[i];
    const auto& c = conserved_[i];
    
    if (cell_volume <= 0.0) {
        throw std::runtime_error("Cell volume must be positive");
    }
    
    double new_mixture_density = c.mass_mixture / cell_volume;
    double new_mass_gas_per_vol = c.mass_gas / cell_volume;
    
    if (new_mixture_density <= 0.0) {
        new_mixture_density = p.liquid_density * 0.5;
    }
    
    double new_void_fraction = new_mass_gas_per_vol / p.gas_density;
    new_void_fraction = std::max(0.0, std::min(0.999, new_void_fraction));
    
    double new_liquid_density = p.liquid_density;
    if (new_void_fraction < 1.0) {
        new_liquid_density = (new_mixture_density - new_void_fraction * p.gas_density) / 
                             (1.0 - new_void_fraction);
    }
    
    double new_mixture_velocity = 0.0;
    if (std::abs(new_mixture_density) > 1e-12) {
        new_mixture_velocity = c.momentum_mixture / (new_mixture_density * cell_volume);
    }
    
    p.mixture_density = new_mixture_density;
    p.void_fraction = new_void_fraction;
    p.liquid_density = new_liquid_density;
    p.mixture_velocity = new_mixture_velocity;
}

void FieldState::update_velocities(int i, double slip_velocity) {
    auto& p = primitives_[i];
    
    double alpha = p.void_fraction;
    double vm = p.mixture_velocity;
    
    if (alpha < 1e-12) {
        p.gas_velocity = vm + slip_velocity;
        p.liquid_velocity = vm;
    } else if (alpha > 0.999) {
        p.gas_velocity = vm;
        p.liquid_velocity = vm - slip_velocity;
    } else {
        p.gas_velocity = vm + (1.0 - alpha) * slip_velocity;
        p.liquid_velocity = vm - alpha * slip_velocity;
    }
}

FaceVariables::FaceVariables(int num_faces)
    : num_faces_(num_faces) {
    mass_flux_mixture_.resize(num_faces_, 0.0);
    mass_flux_gas_.resize(num_faces_, 0.0);
    momentum_flux_.resize(num_faces_, 0.0);
}

}
