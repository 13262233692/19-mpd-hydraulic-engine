#pragma once

#include <vector>
#include "grid.hpp"
#include "variables.hpp"
#include "physics.hpp"
#include "fvm.hpp"
#include "pentadiagonal_solver.hpp"

namespace hydraulic_engine {

struct TimeIntegratorParams {
    double dt_init;
    double dt_min;
    double dt_max;
    double cfl;
    int max_newton_iter;
    double newton_tol;
    int max_time_steps;
    double total_time;
    bool adaptive_time_stepping;
};

class TimeIntegrator {
public:
    enum class Scheme {
        EULER_EXPLICIT = 0,
        EULER_IMPLICIT,
        CRANK_NICOLSON,
        BDF2
    };
    
    TimeIntegrator(const Grid& grid,
                 const PhysicsModel& physics,
                 const FVMDisc& fvm,
                 const TimeIntegratorParams& params,
                 Scheme scheme = Scheme::EULER_IMPLICIT);
    
    bool step(FieldState& state,
              const BoundaryConditions& bc,
              double& current_time,
              double& dt);
    
    bool solve_nonlinear(FieldState& state,
                      const BoundaryConditions& bc,
                      double current_time,
                      double dt);
    
    bool solve_linear_system(FieldState& state,
                             const BoundaryConditions& bc,
                             double current_time,
                             double dt);
    
    double compute_suggested_dt(const FieldState& state) const;
    
    void update_physical_properties(FieldState& state);
    
    const TimeIntegratorParams& params() const { return params_; }
    Scheme scheme() const { return scheme_; }
    
    int last_newton_iterations() const { return last_newton_iter_; }
    double last_newton_residual() const { return last_newton_residual_; }
    int total_time_steps() const { return total_time_steps_; }

private:
    const Grid& grid_;
    const PhysicsModel& physics_;
    const FVMDisc& fvm_;
    TimeIntegratorParams params_;
    Scheme scheme_;
    
    PentadiagonalSolver solver_;
    
    std::vector<double> a_, b_, c_, d_, e_, rhs_;
    std::vector<double> delta_rho_, delta_alpha_, delta_p_;
    
    FieldState state_old_;
    FieldState state_prev_;
    FieldState residuals_;
    FaceVariables face_fluxes_;
    
    int last_newton_iter_;
    double last_newton_residual_;
    int total_time_steps_;
    
    static constexpr double eps_ = 1e-12;
    
    void compute_newton_update(const FieldState& current_state,
                       FieldState& updated_state,
                       const std::vector<double>& delta_rho,
                       const std::vector<double>& delta_alpha,
                       const std::vector<double>& delta_p);
    
    double compute_residual_norm(const FieldState& residuals) const;
    
    void apply_under_relaxation(FieldState& new_state,
                              const FieldState& old_state,
                              double relaxation_factor = 0.8);
};

}
