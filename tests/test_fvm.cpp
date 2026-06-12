#include <iostream>
#include <vector>
#include <cmath>
#include "hydraulic_engine/grid.hpp"
#include "hydraulic_engine/variables.hpp"
#include "hydraulic_engine/physics.hpp"
#include "hydraulic_engine/fvm.hpp"

using namespace hydraulic_engine;

int main() {
    std::cout << "Testing FVM Discretization..." << std::endl;
    
    double total_depth = 1000.0;
    int num_cells = 100;
    double wellbore_diameter = 0.216;
    
    Grid grid(total_depth, num_cells, wellbore_diameter);
    
    FluidProperties fluid_props;
    fluid_props.liquid_density_ref = 1000.0;
    fluid_props.gas_density_ref = 100.0;
    fluid_props.liquid_viscosity = 1e-3;
    fluid_props.gas_viscosity = 1.8e-5;
    fluid_props.surface_tension = 0.072;
    fluid_props.gas_constant = 8.314;
    fluid_props.molar_mass_gas = 0.016;
    fluid_props.reference_pressure = 101325.0;
    fluid_props.reference_temperature = 288.15;
    
    DriftFluxParameters drift_params;
    drift_params.distribution_coefficient = 1.2;
    drift_params.drift_velocity = 0.35;
    drift_params.gravity = 9.81;
    
    PhysicsModel physics(fluid_props, drift_params);
    
    FVMDisc fvm(grid, physics, FluxLimiter::FIRST_ORDER_UPWIND);
    
    FieldState state(num_cells);
    state.initialize(5e6, 0.1, 1000.0, 100.0, 300.0);
    
    for (int i = 0; i < num_cells; ++i) {
        double depth = grid.cell(i).z_center;
        state.at(i).pressure = 1e5 + 1000.0 * 9.81 * depth;
        state.at(i).mixture_velocity = 1.0;
        state.at(i).temperature = 300.0 + 0.02 * depth;
        
        double slip = physics.compute_slip_velocity(
            state.at(i).void_fraction,
            state.at(i).mixture_velocity,
            state.at(i).mixture_density,
            state.at(i).gas_density,
            state.at(i).liquid_density
        );
        state.update_velocities(i, slip);
    }
    
    std::cout << "\nTesting first-order upwind flux..." << std::endl;
    
    PrimitiveVariables left, right;
    left.pressure = 5e6;
    left.void_fraction = 0.2;
    left.mixture_density = 800.0;
    left.gas_density = 100.0;
    left.liquid_density = 1000.0;
    left.mixture_velocity = 2.0;
    left.gas_velocity = 2.5;
    left.liquid_velocity = 1.5;
    left.temperature = 300.0;
    
    right = left;
    right.pressure = 4.9e6;
    right.mixture_velocity = 1.8;
    
    double mass_flux_m, mass_flux_g, mom_flux;
    double face_area = grid.face(0).area;
    
    fvm.first_order_upwind(left, right, face_area, mass_flux_m, mass_flux_g, mom_flux);
    
    std::cout << "Mass flux mixture: " << mass_flux_m << " kg/s" << std::endl;
    std::cout << "Mass flux gas: " << mass_flux_g << " kg/s" << std::endl;
    std::cout << "Momentum flux: " << mom_flux << " N" << std::endl;
    
    if (std::abs(mass_flux_m - 800.0 * 1.9 * face_area) > 1e-3) {
        std::cout << "WARNING: Mass flux may be incorrect" << std::endl;
    }
    
    std::cout << "\nTesting face flux computation..." << std::endl;
    
    FaceVariables face_fluxes(grid.num_faces());
    fvm.compute_face_fluxes(state, face_fluxes, 0.01);
    
    double total_mass_in = face_fluxes.mass_flux_mixture(0);
    double total_mass_out = face_fluxes.mass_flux_mixture(grid.num_faces() - 1);
    
    std::cout << "Total mass flux in: " << total_mass_in << " kg/s" << std::endl;
    std::cout << "Total mass flux out: " << total_mass_out << " kg/s" << std::endl;
    
    std::cout << "\nTesting source term computation..." << std::endl;
    
    const auto& cell = grid.cell(50);
    double source = fvm.compute_source_term(state.at(50), cell, physics, 1e-5);
    
    std::cout << "Source term at cell 50: " << source << " N" << std::endl;
    
    double expected_hydro = -1000.0 * 9.81 * cell.area_top * cell.dz;
    std::cout << "Expected hydrostatic contribution: " << expected_hydro << " N" << std::endl;
    
    std::cout << "\nTesting boundary conditions..." << std::endl;
    
    BoundaryConditions bc;
    bc.inlet_pressure = 1e7;
    bc.outlet_pressure = 1e5;
    bc.inlet_void_fraction = 0.0;
    bc.use_pressure_inlet = true;
    bc.use_pressure_outlet = true;
    bc.enable_kick = false;
    
    fvm.apply_boundary_conditions(state, bc, 0.0);
    
    std::cout << "Inlet pressure after BC: " << state.at(0).pressure << " Pa" << std::endl;
    std::cout << "Outlet pressure after BC: " << state.at(num_cells - 1).pressure << " Pa" << std::endl;
    std::cout << "Inlet void fraction after BC: " << state.at(0).void_fraction << std::endl;
    
    if (std::abs(state.at(0).pressure - 1e7) > 1e-3) {
        std::cout << "ERROR: Inlet pressure BC not applied correctly!" << std::endl;
        return 1;
    }
    
    if (std::abs(state.at(num_cells - 1).pressure - 1e5) > 1e-3) {
        std::cout << "ERROR: Outlet pressure BC not applied correctly!" << std::endl;
        return 1;
    }
    
    std::cout << "\nTesting residual computation..." << std::endl;
    
    FieldState residuals(num_cells);
    fvm.compute_residuals(state, face_fluxes, bc, residuals, 0.0);
    
    double max_residual = 0.0;
    for (int i = 0; i < num_cells; ++i) {
        double res_norm = std::sqrt(
            residuals.conserved(i).mass_mixture * residuals.conserved(i).mass_mixture +
            residuals.conserved(i).mass_gas * residuals.conserved(i).mass_gas +
            residuals.conserved(i).momentum_mixture * residuals.conserved(i).momentum_mixture
        );
        max_residual = std::max(max_residual, res_norm);
    }
    
    std::cout << "Maximum residual: " << max_residual << std::endl;
    
    std::cout << "\nTesting kick source term..." << std::endl;
    
    bc.enable_kick = true;
    bc.kick_flow_rate = 10.0;
    bc.kick_void_fraction = 0.8;
    bc.kick_depth = 500.0;
    bc.kick_start_time = 0.0;
    bc.kick_duration = 100.0;
    
    fvm.compute_residuals(state, face_fluxes, bc, residuals, 50.0);
    
    double kick_cell_residual = std::abs(residuals.conserved(50).mass_mixture);
    std::cout << "Kick cell mass residual: " << kick_cell_residual << " kg" << std::endl;
    
    if (kick_cell_residual < 1e-10) {
        std::cout << "WARNING: Kick source term may not be active" << std::endl;
    }
    
    std::cout << "\nFVM tests completed!" << std::endl;
    return 0;
}
