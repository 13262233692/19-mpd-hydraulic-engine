#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include <vector>
#include <cmath>
#include "hydraulic_engine/pentadiagonal_solver.hpp"

using namespace hydraulic_engine;

int main() {
    std::cout << "Testing Pentadiagonal Solver..." << std::endl;
    
    int n = 100;
    PentadiagonalSolver solver(n);
    
    std::vector<double> a(n, 0.0);
    std::vector<double> b(n, 0.0);
    std::vector<double> c(n, 0.0);
    std::vector<double> d(n, 0.0);
    std::vector<double> e(n, 0.0);
    std::vector<double> rhs(n, 0.0);
    std::vector<double> x(n, 0.0);
    std::vector<double> x_exact(n, 0.0);
    
    for (int i = 0; i < n; ++i) {
        x_exact[i] = std::sin(2.0 * M_PI * i / n) + 1.0;
    }
    
    for (int i = 0; i < n; ++i) {
        c[i] = 4.0;
        if (i >= 1) b[i] = -1.0;
        if (i >= 2) a[i] = -0.5;
        if (i < n - 1) d[i] = -1.0;
        if (i < n - 2) e[i] = -0.5;
    }
    
    for (int i = 0; i < n; ++i) {
        rhs[i] = c[i] * x_exact[i];
        if (i >= 1) rhs[i] += b[i] * x_exact[i - 1];
        if (i >= 2) rhs[i] += a[i] * x_exact[i - 2];
        if (i < n - 1) rhs[i] += d[i] * x_exact[i + 1];
        if (i < n - 2) rhs[i] += e[i] * x_exact[i + 2];
    }
    
    std::cout << "Testing modified Thomas algorithm..." << std::endl;
    bool success = solver.solve_modified_thomas(a, b, c, d, e, rhs, x);
    
    if (!success) {
        std::cout << "ERROR: Solver failed!" << std::endl;
        return 1;
    }
    
    double error = 0.0;
    for (int i = 0; i < n; ++i) {
        error += (x[i] - x_exact[i]) * (x[i] - x_exact[i]);
    }
    error = std::sqrt(error / n);
    
    std::cout << "RMS error: " << error << std::endl;
    
    if (error > 1e-8) {
        std::cout << "ERROR: Solution error too large!" << std::endl;
        return 1;
    }
    
    double residual = solver.compute_residual_norm(a, b, c, d, e, rhs, x);
    std::cout << "Residual norm: " << residual << std::endl;
    
    if (residual > 1e-10) {
        std::cout << "ERROR: Residual too large!" << std::endl;
        return 1;
    }
    
    bool dd = solver.is_diagonally_dominant(a, b, c, d, e);
    std::cout << "Diagonally dominant: " << (dd ? "Yes" : "No") << std::endl;
    
    std::cout << "\nTesting cyclic reduction..." << std::endl;
    std::vector<double> x2(n, 0.0);
    success = solver.solve_cyclic_reduction(a, b, c, d, e, rhs, x2);
    
    if (!success) {
        std::cout << "ERROR: Cyclic reduction failed!" << std::endl;
        return 1;
    }
    
    error = 0.0;
    for (int i = 0; i < n; ++i) {
        error += (x2[i] - x_exact[i]) * (x2[i] - x_exact[i]);
    }
    error = std::sqrt(error / n);
    std::cout << "Cyclic reduction RMS error: " << error << std::endl;
    
    if (error > 1e-8) {
        std::cout << "ERROR: Cyclic reduction error too large!" << std::endl;
        return 1;
    }
    
    std::cout << "\nTesting parallel solver..." << std::endl;
    std::vector<double> x3(n, 0.0);
    success = solver.solve_parallel(a, b, c, d, e, rhs, x3);
    
    if (!success) {
        std::cout << "ERROR: Parallel solver failed!" << std::endl;
        return 1;
    }
    
    error = 0.0;
    for (int i = 0; i < n; ++i) {
        error += (x3[i] - x_exact[i]) * (x3[i] - x_exact[i]);
    }
    error = std::sqrt(error / n);
    std::cout << "Parallel solver RMS error: " << error << std::endl;
    
    if (error > 1e-6) {
        std::cout << "ERROR: Parallel solver error too large!" << std::endl;
        return 1;
    }
    
    std::cout << "\nTesting larger problem (n=1000)..." << std::endl;
    int n_large = 1000;
    solver.resize(n_large);
    
    std::vector<double> a_l(n_large, 0.0);
    std::vector<double> b_l(n_large, 0.0);
    std::vector<double> c_l(n_large, 0.0);
    std::vector<double> d_l(n_large, 0.0);
    std::vector<double> e_l(n_large, 0.0);
    std::vector<double> rhs_l(n_large, 0.0);
    std::vector<double> x_l(n_large, 0.0);
    std::vector<double> x_exact_l(n_large, 0.0);
    
    for (int i = 0; i < n_large; ++i) {
        x_exact_l[i] = std::cos(4.0 * M_PI * i / n_large) * std::exp(-0.001 * i);
    }
    
    for (int i = 0; i < n_large; ++i) {
        c_l[i] = 6.0;
        if (i >= 1) b_l[i] = -2.0;
        if (i >= 2) a_l[i] = -1.0;
        if (i < n_large - 1) d_l[i] = -2.0;
        if (i < n_large - 2) e_l[i] = -1.0;
    }
    
    for (int i = 0; i < n_large; ++i) {
        rhs_l[i] = c_l[i] * x_exact_l[i];
        if (i >= 1) rhs_l[i] += b_l[i] * x_exact_l[i - 1];
        if (i >= 2) rhs_l[i] += a_l[i] * x_exact_l[i - 2];
        if (i < n_large - 1) rhs_l[i] += d_l[i] * x_exact_l[i + 1];
        if (i < n_large - 2) rhs_l[i] += e_l[i] * x_exact_l[i + 2];
    }
    
    success = solver.solve_modified_thomas(a_l, b_l, c_l, d_l, e_l, rhs_l, x_l);
    
    if (!success) {
        std::cout << "ERROR: Large problem solver failed!" << std::endl;
        return 1;
    }
    
    error = 0.0;
    for (int i = 0; i < n_large; ++i) {
        error += (x_l[i] - x_exact_l[i]) * (x_l[i] - x_exact_l[i]);
    }
    error = std::sqrt(error / n_large);
    std::cout << "Large problem RMS error: " << error << std::endl;
    
    if (error > 1e-8) {
        std::cout << "ERROR: Large problem error too large!" << std::endl;
        return 1;
    }
    
    std::cout << "\nAll tests passed!" << std::endl;
    return 0;
}
