#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>

#include "hydraulic_engine/grid.hpp"
#include "hydraulic_engine/variables.hpp"
#include "hydraulic_engine/physics.hpp"
#include "hydraulic_engine/fvm.hpp"
#include "hydraulic_engine/pentadiagonal_solver.hpp"
#include "hydraulic_engine/time_integrator.hpp"
#include "hydraulic_engine/simulator.hpp"

namespace py = pybind11;
using namespace hydraulic_engine;

PYBIND11_MODULE(hydraulic_engine, m) {
    m.doc() = "Deepwater Drilling Gas-Liquid Two-Phase Flow Simulator";
    
    py::enum_<FluxLimiter::Type>(m, "FluxLimiter")
        .value("FIRST_ORDER_UPWIND", FluxLimiter::FIRST_ORDER_UPWIND)
        .value("MINMOD", FluxLimiter::MINMOD)
        .value("SUPERBEE", FluxLimiter::SUPERBEE)
        .value("VAN_LEER", FluxLimiter::VAN_LEER)
        .value("VAN_ALBADA", FluxLimiter::VAN_ALBADA)
        .value("MUSCL", FluxLimiter::MUSCL)
        .value("QUICK", FluxLimiter::QUICK);
    
    py::enum_<TimeIntegrator::Scheme>(m, "TimeScheme")
        .value("EULER_EXPLICIT", TimeIntegrator::Scheme::EULER_EXPLICIT)
        .value("EULER_IMPLICIT", TimeIntegrator::Scheme::EULER_IMPLICIT)
        .value("CRANK_NICOLSON", TimeIntegrator::Scheme::CRANK_NICOLSON)
        .value("BDF2", TimeIntegrator::Scheme::BDF2);
    
    py::class_<GridCell>(m, "GridCell")
        .def_readonly("z_center", &GridCell::z_center)
        .def_readonly("z_top", &GridCell::z_top)
        .def_readonly("z_bottom", &GridCell::z_bottom)
        .def_readonly("volume", &GridCell::volume)
        .def_readonly("area_top", &GridCell::area_top)
        .def_readonly("area_bottom", &GridCell::area_bottom)
        .def_readonly("dz", &GridCell::dz);
    
    py::class_<GridFace>(m, "GridFace")
        .def_readonly("z", &GridFace::z)
        .def_readonly("area", &GridFace::area);
    
    py::class_<Grid>(m, "Grid")
        .def(py::init<double, int, double>(),
             py::arg("total_depth"),
             py::arg("num_cells"),
             py::arg("wellbore_diameter"))
        .def_property_readonly("num_cells", &Grid::num_cells)
        .def_property_readonly("num_faces", &Grid::num_faces)
        .def_property_readonly("total_depth", &Grid::total_depth)
        .def_property_readonly("wellbore_radius", &Grid::wellbore_radius)
        .def("cell", py::overload_cast<int>(&Grid::cell, py::const_),
             py::return_value_policy::reference_internal)
        .def("face", py::overload_cast<int>(&Grid::face, py::const_),
             py::return_value_policy::reference_internal)
        .def("cells", &Grid::cells,
             py::return_value_policy::reference_internal)
        .def("faces", &Grid::faces,
             py::return_value_policy::reference_internal);
    
    py::class_<PrimitiveVariables>(m, "PrimitiveVariables")
        .def(py::init<>())
        .def_readwrite("pressure", &PrimitiveVariables::pressure)
        .def_readwrite("mixture_density", &PrimitiveVariables::mixture_density)
        .def_readwrite("void_fraction", &PrimitiveVariables::void_fraction)
        .def_readwrite("mixture_velocity", &PrimitiveVariables::mixture_velocity)
        .def_readwrite("gas_density", &PrimitiveVariables::gas_density)
        .def_readwrite("liquid_density", &PrimitiveVariables::liquid_density)
        .def_readwrite("gas_velocity", &PrimitiveVariables::gas_velocity)
        .def_readwrite("liquid_velocity", &PrimitiveVariables::liquid_velocity)
        .def_readwrite("temperature", &PrimitiveVariables::temperature);
    
    py::class_<ConservativeVariables>(m, "ConservativeVariables")
        .def(py::init<>())
        .def_readwrite("mass_mixture", &ConservativeVariables::mass_mixture)
        .def_readwrite("mass_gas", &ConservativeVariables::mass_gas)
        .def_readwrite("momentum_mixture", &ConservativeVariables::momentum_mixture);
    
    py::class_<FieldState>(m, "FieldState")
        .def(py::init<int>(), py::arg("num_cells"))
        .def_property_readonly("num_cells", &FieldState::num_cells)
        .def("at", py::overload_cast<int>(&FieldState::at),
             py::return_value_policy::reference_internal)
        .def("conserved", py::overload_cast<int>(&FieldState::conserved),
             py::return_value_policy::reference_internal)
        .def("primitives", py::overload_cast<>(&FieldState::primitives),
             py::return_value_policy::reference_internal)
        .def("initialize", &FieldState::initialize,
             py::arg("initial_pressure"),
             py::arg("initial_void_fraction"),
             py::arg("liquid_density"),
             py::arg("gas_density"),
             py::arg("initial_temperature"))
        .def("primitive_to_conservative", &FieldState::primitive_to_conservative,
             py::arg("i"), py::arg("cell_volume"))
        .def("conservative_to_primitive", &FieldState::conservative_to_primitive,
             py::arg("i"), py::arg("cell_volume"))
        .def("update_velocities", &FieldState::update_velocities,
             py::arg("i"), py::arg("slip_velocity"));
    
    py::class_<FluidProperties>(m, "FluidProperties")
        .def(py::init<>())
        .def_readwrite("liquid_density_ref", &FluidProperties::liquid_density_ref)
        .def_readwrite("gas_density_ref", &FluidProperties::gas_density_ref)
        .def_readwrite("liquid_viscosity", &FluidProperties::liquid_viscosity)
        .def_readwrite("gas_viscosity", &FluidProperties::gas_viscosity)
        .def_readwrite("surface_tension", &FluidProperties::surface_tension)
        .def_readwrite("gas_constant", &FluidProperties::gas_constant)
        .def_readwrite("molar_mass_gas", &FluidProperties::molar_mass_gas)
        .def_readwrite("reference_pressure", &FluidProperties::reference_pressure)
        .def_readwrite("reference_temperature", &FluidProperties::reference_temperature);
    
    py::class_<DriftFluxParameters>(m, "DriftFluxParameters")
        .def(py::init<>())
        .def_readwrite("distribution_coefficient", &DriftFluxParameters::distribution_coefficient)
        .def_readwrite("drift_velocity", &DriftFluxParameters::drift_velocity)
        .def_readwrite("gravity", &DriftFluxParameters::gravity);
    
    py::class_<BoundaryConditions>(m, "BoundaryConditions")
        .def(py::init<>())
        .def_readwrite("inlet_pressure", &BoundaryConditions::inlet_pressure)
        .def_readwrite("outlet_pressure", &BoundaryConditions::outlet_pressure)
        .def_readwrite("inlet_mass_flow_rate", &BoundaryConditions::inlet_mass_flow_rate)
        .def_readwrite("inlet_void_fraction", &BoundaryConditions::inlet_void_fraction)
        .def_readwrite("inlet_temperature", &BoundaryConditions::inlet_temperature)
        .def_readwrite("outlet_temperature", &BoundaryConditions::outlet_temperature)
        .def_readwrite("bottom_hole_pressure", &BoundaryConditions::bottom_hole_pressure)
        .def_readwrite("kick_flow_rate", &BoundaryConditions::kick_flow_rate)
        .def_readwrite("kick_void_fraction", &BoundaryConditions::kick_void_fraction)
        .def_readwrite("kick_depth", &BoundaryConditions::kick_depth)
        .def_readwrite("use_pressure_inlet", &BoundaryConditions::use_pressure_inlet)
        .def_readwrite("use_pressure_outlet", &BoundaryConditions::use_pressure_outlet)
        .def_readwrite("enable_kick", &BoundaryConditions::enable_kick)
        .def_readwrite("kick_start_time", &BoundaryConditions::kick_start_time)
        .def_readwrite("kick_duration", &BoundaryConditions::kick_duration);
    
    py::class_<PhysicsModel>(m, "PhysicsModel")
        .def(py::init<const FluidProperties&, const DriftFluxParameters&>(),
             py::arg("fluid_props"), py::arg("drift_params"))
        .def("compute_mixture_density", &PhysicsModel::compute_mixture_density,
             py::arg("void_fraction"), py::arg("gas_density"), py::arg("liquid_density"))
        .def("compute_gas_density", &PhysicsModel::compute_gas_density,
             py::arg("pressure"), py::arg("temperature"))
        .def("compute_liquid_density", &PhysicsModel::compute_liquid_density,
             py::arg("pressure"), py::arg("temperature"))
        .def("compute_slip_velocity", &PhysicsModel::compute_slip_velocity,
             py::arg("void_fraction"), py::arg("mixture_velocity"),
             py::arg("mixture_density"), py::arg("gas_density"),
             py::arg("liquid_density"))
        .def("compute_drift_velocity", &PhysicsModel::compute_drift_velocity,
             py::arg("void_fraction"), py::arg("mixture_density"),
             py::arg("gas_density"), py::arg("liquid_density"))
        .def("update_phase_velocities", &PhysicsModel::update_phase_velocities,
             py::arg("vars"), py::arg("slip_velocity"))
        .def("compute_friction_factor", &PhysicsModel::compute_friction_factor,
             py::arg("reynolds_number"), py::arg("relative_roughness"))
        .def("compute_reynolds_number", &PhysicsModel::compute_reynolds_number,
             py::arg("mixture_density"), py::arg("mixture_velocity"),
             py::arg("hydraulic_diameter"), py::arg("mixture_viscosity"))
        .def("compute_mixture_viscosity", &PhysicsModel::compute_mixture_viscosity,
             py::arg("void_fraction"), py::arg("gas_viscosity"),
             py::arg("liquid_viscosity"))
        .def("compute_hydrostatic_gradient", &PhysicsModel::compute_hydrostatic_gradient,
             py::arg("mixture_density"), py::arg("inclination_angle") = 0.0)
        .def("compute_friction_gradient", &PhysicsModel::compute_friction_gradient,
             py::arg("mixture_density"), py::arg("mixture_velocity"),
             py::arg("friction_factor"), py::arg("hydraulic_diameter"))
        .def("compute_acceleration_gradient", &PhysicsModel::compute_acceleration_gradient,
             py::arg("mixture_density"), py::arg("mixture_velocity"),
             py::arg("dvelocity_dz"))
        .def_property_readonly("fluid_properties", &PhysicsModel::fluid_properties)
        .def_property_readonly("drift_parameters", &PhysicsModel::drift_parameters);
    
    py::class_<PentadiagonalSolver>(m, "PentadiagonalSolver")
        .def(py::init<int>(), py::arg("max_size") = 0)
        .def("resize", &PentadiagonalSolver::resize, py::arg("n"))
        .def_property_readonly("size", &PentadiagonalSolver::size)
        .def("solve",
             [](PentadiagonalSolver& self,
                const std::vector<double>& a,
                const std::vector<double>& b,
                const std::vector<double>& c,
                const std::vector<double>& d,
                const std::vector<double>& e,
                const std::vector<double>& rhs,
                std::vector<double>& x) -> py::tuple {
                 bool success = self.solve(a, b, c, d, e, rhs, x);
                 return py::make_tuple(success, x);
             },
             py::arg("a"), py::arg("b"), py::arg("c"),
             py::arg("d"), py::arg("e"), py::arg("rhs"), py::arg("x"))
        .def("solve_modified_thomas",
             [](PentadiagonalSolver& self,
                const std::vector<double>& a,
                const std::vector<double>& b,
                const std::vector<double>& c,
                const std::vector<double>& d,
                const std::vector<double>& e,
                const std::vector<double>& rhs,
                std::vector<double>& x) -> py::tuple {
                 bool success = self.solve_modified_thomas(a, b, c, d, e, rhs, x);
                 return py::make_tuple(success, x);
             },
             py::arg("a"), py::arg("b"), py::arg("c"),
             py::arg("d"), py::arg("e"), py::arg("rhs"), py::arg("x"))
        .def("solve_cyclic_reduction",
             [](PentadiagonalSolver& self,
                const std::vector<double>& a,
                const std::vector<double>& b,
                const std::vector<double>& c,
                const std::vector<double>& d,
                const std::vector<double>& e,
                const std::vector<double>& rhs,
                std::vector<double>& x) -> py::tuple {
                 bool success = self.solve_cyclic_reduction(a, b, c, d, e, rhs, x);
                 return py::make_tuple(success, x);
             },
             py::arg("a"), py::arg("b"), py::arg("c"),
             py::arg("d"), py::arg("e"), py::arg("rhs"), py::arg("x"))
        .def("solve_parallel",
             [](PentadiagonalSolver& self,
                const std::vector<double>& a,
                const std::vector<double>& b,
                const std::vector<double>& c,
                const std::vector<double>& d,
                const std::vector<double>& e,
                const std::vector<double>& rhs,
                std::vector<double>& x) -> py::tuple {
                 bool success = self.solve_parallel(a, b, c, d, e, rhs, x);
                 return py::make_tuple(success, x);
             },
             py::arg("a"), py::arg("b"), py::arg("c"),
             py::arg("d"), py::arg("e"), py::arg("rhs"), py::arg("x"))
        .def("compute_residual_norm", &PentadiagonalSolver::compute_residual_norm,
             py::arg("a"), py::arg("b"), py::arg("c"),
             py::arg("d"), py::arg("e"), py::arg("rhs"), py::arg("x"))
        .def("is_diagonally_dominant", &PentadiagonalSolver::is_diagonally_dominant,
             py::arg("a"), py::arg("b"), py::arg("c"),
             py::arg("d"), py::arg("e"))
        .def_property("tolerance", &PentadiagonalSolver::tolerance,
                      &PentadiagonalSolver::set_tolerance)
        .def_property("max_iterations", &PentadiagonalSolver::max_iterations,
                      &PentadiagonalSolver::set_max_iterations);
    
    py::class_<SimulationConfig>(m, "SimulationConfig")
        .def(py::init<>())
        .def_readwrite("total_depth", &SimulationConfig::total_depth)
        .def_readwrite("num_cells", &SimulationConfig::num_cells)
        .def_readwrite("wellbore_diameter", &SimulationConfig::wellbore_diameter)
        .def_readwrite("initial_pressure", &SimulationConfig::initial_pressure)
        .def_readwrite("initial_void_fraction", &SimulationConfig::initial_void_fraction)
        .def_readwrite("initial_temperature", &SimulationConfig::initial_temperature)
        .def_readwrite("liquid_density_ref", &SimulationConfig::liquid_density_ref)
        .def_readwrite("gas_density_ref", &SimulationConfig::gas_density_ref)
        .def_readwrite("liquid_viscosity", &SimulationConfig::liquid_viscosity)
        .def_readwrite("gas_viscosity", &SimulationConfig::gas_viscosity)
        .def_readwrite("surface_tension", &SimulationConfig::surface_tension)
        .def_readwrite("molar_mass_gas", &SimulationConfig::molar_mass_gas)
        .def_readwrite("distribution_coefficient", &SimulationConfig::distribution_coefficient)
        .def_readwrite("drift_velocity_coeff", &SimulationConfig::drift_velocity_coeff)
        .def_readwrite("gravity", &SimulationConfig::gravity)
        .def_readwrite("dt_init", &SimulationConfig::dt_init)
        .def_readwrite("dt_min", &SimulationConfig::dt_min)
        .def_readwrite("dt_max", &SimulationConfig::dt_max)
        .def_readwrite("cfl", &SimulationConfig::cfl)
        .def_readwrite("max_newton_iter", &SimulationConfig::max_newton_iter)
        .def_readwrite("newton_tol", &SimulationConfig::newton_tol)
        .def_readwrite("adaptive_time_stepping", &SimulationConfig::adaptive_time_stepping)
        .def_readwrite("flux_limiter", &SimulationConfig::flux_limiter);
    
    py::class_<SimulationOutput>(m, "SimulationOutput")
        .def(py::init<>())
        .def_readwrite("time", &SimulationOutput::time)
        .def_readwrite("depth", &SimulationOutput::depth)
        .def_readwrite("pressure", &SimulationOutput::pressure)
        .def_readwrite("void_fraction", &SimulationOutput::void_fraction)
        .def_readwrite("mixture_density", &SimulationOutput::mixture_density)
        .def_readwrite("gas_density", &SimulationOutput::gas_density)
        .def_readwrite("liquid_density", &SimulationOutput::liquid_density)
        .def_readwrite("mixture_velocity", &SimulationOutput::mixture_velocity)
        .def_readwrite("gas_velocity", &SimulationOutput::gas_velocity)
        .def_readwrite("liquid_velocity", &SimulationOutput::liquid_velocity)
        .def_readwrite("temperature", &SimulationOutput::temperature);
    
    py::class_<Simulator>(m, "Simulator")
        .def(py::init<const SimulationConfig&>(), py::arg("config"))
        .def("initialize", &Simulator::initialize)
        .def("run", &Simulator::run,
             py::arg("total_time"), py::arg("output_interval") = 10)
        .def("step", &Simulator::step,
             py::arg("current_time"), py::arg("dt"))
        .def("apply_boundary_conditions", &Simulator::apply_boundary_conditions,
             py::arg("bc"))
        .def("enable_kick", &Simulator::enable_kick,
             py::arg("kick_flow_rate"), py::arg("kick_void_fraction"),
             py::arg("kick_depth"), py::arg("kick_start_time") = 0.0,
             py::arg("kick_duration") = 1e6)
        .def("disable_kick", &Simulator::disable_kick)
        .def("get_output", &Simulator::get_output)
        .def_property_readonly("current_time", &Simulator::get_current_time)
        .def_property_readonly("current_dt", &Simulator::get_current_dt)
        .def_property_readonly("total_steps", &Simulator::get_total_steps)
        .def_property_readonly("last_newton_iterations", &Simulator::get_last_newton_iterations)
        .def_property_readonly("last_newton_residual", &Simulator::get_last_newton_residual)
        .def_property_readonly("grid", &Simulator::grid,
             py::return_value_policy::reference_internal)
        .def_property_readonly("state", &Simulator::state,
             py::return_value_policy::reference_internal)
        .def_property_readonly("boundary_conditions", &Simulator::boundary_conditions)
        .def("set_output_callback", &Simulator::set_output_callback,
             py::arg("callback"));
    
    m.def("create_default_config", []() {
        return SimulationConfig();
    });
}
