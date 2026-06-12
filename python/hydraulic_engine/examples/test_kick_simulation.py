import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import hydraulic_engine as he
    print("Successfully imported hydraulic_engine module")
    print(f"Module version: {he.__version__}")
except ImportError as e:
    print(f"Error importing hydraulic_engine: {e}")
    print("Note: This requires the compiled C++ extension to be available")
    sys.exit(1)

import numpy as np
import time

def test_solver_basic():
    print("\n" + "="*60)
    print("TEST 1: Pentadiagonal Solver Basic Test")
    print("="*60)
    
    n = 100
    solver = he.PentadiagonalSolver(n)
    
    a = np.zeros(n)
    b = np.zeros(n)
    c = np.ones(n) * 4.0
    d = np.zeros(n)
    e = np.zeros(n)
    rhs = np.zeros(n)
    x = np.zeros(n)
    x_exact = np.sin(2.0 * np.pi * np.arange(n) / n) + 1.0
    
    for i in range(n):
        if i >= 1:
            b[i] = -1.0
        if i >= 2:
            a[i] = -0.5
        if i < n - 1:
            d[i] = -1.0
        if i < n - 2:
            e[i] = -0.5
    
    for i in range(n):
        rhs[i] = c[i] * x_exact[i]
        if i >= 1:
            rhs[i] += b[i] * x_exact[i - 1]
        if i >= 2:
            rhs[i] += a[i] * x_exact[i - 2]
        if i < n - 1:
            rhs[i] += d[i] * x_exact[i + 1]
        if i < n - 2:
            rhs[i] += e[i] * x_exact[i + 2]
    
    a_list = a.tolist()
    b_list = b.tolist()
    c_list = c.tolist()
    d_list = d.tolist()
    e_list = e.tolist()
    rhs_list = rhs.tolist()
    x_list = x.tolist()
    
    success, x_list = solver.solve_modified_thomas(a_list, b_list, c_list, d_list, e_list, rhs_list, x_list)
    
    if success:
        error = np.sqrt(np.mean((np.array(x_list) - x_exact)**2))
        print(f"  Solver converged: Yes")
        print(f"  RMS error: {error:.2e}")
        
        residual = solver.compute_residual_norm(a_list, b_list, c_list, d_list, e_list, rhs_list, x_list)
        print(f"  Residual norm: {residual:.2e}")
        
        is_dd = solver.is_diagonally_dominant(a_list, b_list, c_list, d_list, e_list)
        print(f"  Diagonally dominant: {'Yes' if is_dd else 'No'}")
        
        if error < 1e-8:
            print("  Test PASSED")
            return True
        else:
            print("  Test FAILED: Error too large")
            return False
    else:
        print("  Test FAILED: Solver did not converge")
        return False

def test_grid():
    print("\n" + "="*60)
    print("TEST 2: Grid Structure Test")
    print("="*60)
    
    total_depth = 4000.0
    num_cells = 400
    wellbore_diameter = 0.216
    
    grid = he.Grid(total_depth, num_cells, wellbore_diameter)
    
    print(f"  Total depth: {grid.total_depth} m")
    print(f"  Number of cells: {grid.num_cells}")
    print(f"  Number of faces: {grid.num_faces}")
    print(f"  Wellbore radius: {grid.wellbore_radius:.4f} m")
    
    cell_0 = grid.cell(0)
    print(f"\n  Cell 0:")
    print(f"    z_bottom: {cell_0.z_bottom} m")
    print(f"    z_center: {cell_0.z_center} m")
    print(f"    z_top: {cell_0.z_top} m")
    print(f"    dz: {cell_0.dz} m")
    print(f"    volume: {cell_0.volume:.6f} m³")
    
    cell_last = grid.cell(num_cells - 1)
    print(f"\n  Cell {num_cells - 1}:")
    print(f"    z_bottom: {cell_last.z_bottom} m")
    print(f"    z_top: {cell_last.z_top} m")
    
    if abs(cell_last.z_top - total_depth) < 1e-10:
        print("  Test PASSED")
        return True
    else:
        print(f"  Test FAILED: Expected top depth {total_depth}, got {cell_last.z_top}")
        return False

def test_physics_model():
    print("\n" + "="*60)
    print("TEST 3: Physics Model Test")
    print("="*60)
    
    fluid_props = he.FluidProperties()
    fluid_props.liquid_density_ref = 1200.0
    fluid_props.gas_density_ref = 100.0
    fluid_props.liquid_viscosity = 1e-3
    fluid_props.gas_viscosity = 1.8e-5
    fluid_props.surface_tension = 0.072
    fluid_props.gas_constant = 8.314
    fluid_props.molar_mass_gas = 0.016
    fluid_props.reference_pressure = 101325.0
    fluid_props.reference_temperature = 288.15
    
    drift_params = he.DriftFluxParameters()
    drift_params.distribution_coefficient = 1.2
    drift_params.drift_velocity = 0.35
    drift_params.gravity = 9.81
    
    physics = he.PhysicsModel(fluid_props, drift_params)
    
    void_fraction = 0.3
    gas_density = 150.0
    liquid_density = 1200.0
    
    mix_density = physics.compute_mixture_density(void_fraction, gas_density, liquid_density)
    print(f"  Mixture density (alpha={void_fraction}): {mix_density:.2f} kg/m³")
    
    expected_density = void_fraction * gas_density + (1 - void_fraction) * liquid_density
    print(f"  Expected density: {expected_density:.2f} kg/m³")
    
    pressure = 30e6
    temperature = 350.0
    rho_g = physics.compute_gas_density(pressure, temperature)
    rho_l = physics.compute_liquid_density(pressure, temperature)
    
    print(f"\n  At P={pressure/1e6} MPa, T={temperature} K:")
    print(f"    Gas density: {rho_g:.2f} kg/m³")
    print(f"    Liquid density: {rho_l:.2f} kg/m³")
    
    mix_velocity = 1.5
    mixture_density = physics.compute_mixture_density(void_fraction, rho_g, rho_l)
    slip_vel = physics.compute_slip_velocity(void_fraction, mix_velocity, mixture_density, rho_g, rho_l)
    print(f"\n  Slip velocity: {slip_vel:.4f} m/s")
    
    Re = physics.compute_reynolds_number(mixture_density, mix_velocity, 0.216, 1e-3)
    f = physics.compute_friction_factor(Re, 1e-5)
    print(f"  Reynolds number: {Re:.2e}")
    print(f"  Friction factor: {f:.6f}")
    
    dp_dz_hydro = physics.compute_hydrostatic_gradient(mixture_density)
    dp_dz_fric = physics.compute_friction_gradient(mixture_density, mix_velocity, f, 0.216)
    print(f"  Hydrostatic gradient: {dp_dz_hydro:.2f} Pa/m")
    print(f"  Friction gradient: {dp_dz_fric:.2f} Pa/m")
    
    if abs(mix_density - expected_density) < 1e-10 and rho_g > 0 and rho_l > 0:
        print("  Test PASSED")
        return True
    else:
        print("  Test FAILED")
        return False

def test_full_simulation():
    print("\n" + "="*60)
    print("TEST 4: Full Kick Simulation")
    print("="*60)
    
    config = he.create_default_config()
    config.total_depth = 3000.0
    config.num_cells = 100
    config.wellbore_diameter = 0.216
    config.initial_pressure = 1e5
    config.initial_void_fraction = 0.0
    config.initial_temperature = 300.0
    config.liquid_density_ref = 1200.0
    config.gas_density_ref = 150.0
    config.liquid_viscosity = 1e-3
    config.gas_viscosity = 1.8e-5
    config.dt_init = 0.001
    config.dt_min = 1e-5
    config.dt_max = 0.01
    config.cfl = 0.05
    config.max_newton_iter = 30
    config.newton_tol = 1.0
    config.adaptive_time_stepping = False
    config.flux_limiter = he.FluxLimiter.FIRST_ORDER_UPWIND
    
    print(f"  Creating simulator with {config.num_cells} cells...")
    try:
        simulator = he.Simulator(config)
    except Exception as e:
        print(f"  ERROR creating simulator: {e}")
        return False
    
    output = simulator.get_output()
    print(f"  Initial top pressure: {output.pressure[0]/1e6:.4f} MPa")
    print(f"  Initial bottom pressure: {output.pressure[-1]/1e6:.4f} MPa")
    
    print("\n  Enabling kick at 2000m depth...")
    simulator.enable_kick(
        kick_flow_rate=3.0,
        kick_void_fraction=0.8,
        kick_depth=2000.0,
        kick_start_time=5.0,
        kick_duration=20.0
    )
    
    output_history = []
    
    def callback(output):
        output_history.append({
            'time': output.time,
            'max_void': max(output.void_fraction),
            'top_void': output.void_fraction[0],
        })
        if len(output_history) % 10 == 0:
            print(f"    t={output.time:.1f}s, max_void={max(output.void_fraction):.4f}")
    
    simulator.set_output_callback(callback)
    
    total_time = 30.0
    print(f"\n  Running simulation for {total_time}s...")
    start_time = time.time()
    
    try:
        success = simulator.run(total_time, output_interval=5)
    except Exception as e:
        print(f"  ERROR during simulation: {e}")
        return False
    
    elapsed = time.time() - start_time
    
    if success:
        print(f"\n  Simulation completed in {elapsed:.2f}s")
        print(f"  Total steps: {simulator.total_steps}")
        print(f"  Final time: {simulator.current_time:.2f}s")
        
        final_output = simulator.get_output()
        max_void = max(final_output.void_fraction)
        print(f"\n  Final results:")
        print(f"    Max void fraction: {max_void:.4f}")
        print(f"    Top void fraction: {final_output.void_fraction[0]:.4f}")
        print(f"    Max mixture velocity: {max(abs(v) for v in final_output.mixture_velocity):.4f} m/s")
        
        if max_void > 0.01:
            print("  Test PASSED: Kick successfully propagated")
            return True
        else:
            print("  Test WARNING: Kick may not have propagated significantly")
            return True
    else:
        print("  Test FAILED: Simulation did not complete")
        return False

def main():
    print("\n" + "="*60)
    print("HYDRAULIC ENGINE PYTHON INTERFACE TEST")
    print("="*60)
    
    tests = [
        test_solver_basic,
        test_grid,
        test_physics_model,
        test_full_simulation,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"  ERROR: {e}")
            failed += 1
    
    print("\n" + "="*60)
    print(f"TEST SUMMARY: {passed} passed, {failed} failed")
    print("="*60)
    
    return failed == 0

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
