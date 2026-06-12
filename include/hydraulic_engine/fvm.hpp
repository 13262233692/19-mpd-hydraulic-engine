#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <vector>
#include <cmath>
#include "grid.hpp"
#include "variables.hpp"
#include "physics.hpp"

namespace hydraulic_engine {

struct FluxLimiter {
    enum Type {
        FIRST_ORDER_UPWIND = 0,
        MINMOD,
        SUPERBEE,
        VAN_LEER,
        VAN_ALBADA,
        MUSCL,
        QUICK
    };
};

class FVMDisc {
public:
    FVMDisc(const Grid& grid,
            const PhysicsModel& physics,
            FluxLimiter::Type limiter_type = FluxLimiter::FIRST_ORDER_UPWIND);
    
    void compute_face_fluxes(const FieldState& state,
                             FaceVariables& face_fluxes,
                             double dt) const;
    
    void compute_residuals(const FieldState& state,
                           const FaceVariables& face_fluxes,
                           const BoundaryConditions& bc,
                           FieldState& residuals,
                           double current_time) const;
    
    void apply_boundary_conditions(FieldState& state,
                                   const BoundaryConditions& bc,
                                   double current_time) const;
    
    void first_order_upwind(const PrimitiveVariables& left,
                            const PrimitiveVariables& right,
                            const double face_area,
                            double& mass_flux_mixture,
                            double& mass_flux_gas,
                            double& momentum_flux) const;
    
    void muscl_reconstruction(const FieldState& state,
                              int face_idx,
                              PrimitiveVariables& left,
                              PrimitiveVariables& right) const;
    
    double compute_source_term(const PrimitiveVariables& vars,
                               const GridCell& cell,
                               const PhysicsModel& physics,
                               double relative_roughness) const;
    
    void compute_jacobian_vectors(const FieldState& state,
                                  const BoundaryConditions& bc,
                                  std::vector<double>& a,
                                  std::vector<double>& b,
                                  std::vector<double>& c,
                                  std::vector<double>& d,
                                  std::vector<double>& e,
                                  std::vector<double>& rhs,
                                  double dt,
                                  double current_time,
                                  int equation_idx) const;

    void set_relative_roughness(double roughness) { relative_roughness_ = roughness; }
    double relative_roughness() const { return relative_roughness_; }

private:
    const Grid& grid_;
    const PhysicsModel& physics_;
    FluxLimiter::Type limiter_type_;
    double relative_roughness_;
    
    static constexpr double eps_ = 1e-12;
    
    double minmod(double a, double b) const;
    double minmod3(double a, double b, double c) const;
    double flux_limiter(double r, FluxLimiter::Type type) const;
    
    void apply_kick_source(FieldState& residuals,
                           const FieldState& state,
                           const BoundaryConditions& bc,
                           double current_time) const;
};

}
