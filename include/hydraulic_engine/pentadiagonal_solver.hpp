#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <omp.h>

namespace hydraulic_engine {

class PentadiagonalSolver {
public:
    PentadiagonalSolver(int max_size = 0);
    
    void resize(int n);
    
    int size() const { return n_; }
    
    bool solve(const std::vector<double>& a,
               const std::vector<double>& b,
               const std::vector<double>& c,
               const std::vector<double>& d,
               const std::vector<double>& e,
               const std::vector<double>& rhs,
               std::vector<double>& x);
    
    bool solve_modified_thomas(const std::vector<double>& a,
                               const std::vector<double>& b,
                               const std::vector<double>& c,
                               const std::vector<double>& d,
                               const std::vector<double>& e,
                               const std::vector<double>& rhs,
                               std::vector<double>& x);
    
    bool solve_cyclic_reduction(const std::vector<double>& a,
                                const std::vector<double>& b,
                                const std::vector<double>& c,
                                const std::vector<double>& d,
                                const std::vector<double>& e,
                                const std::vector<double>& rhs,
                                std::vector<double>& x);
    
    bool solve_parallel(const std::vector<double>& a,
                        const std::vector<double>& b,
                        const std::vector<double>& c,
                        const std::vector<double>& d,
                        const std::vector<double>& e,
                        const std::vector<double>& rhs,
                        std::vector<double>& x);
    
    double compute_residual_norm(const std::vector<double>& a,
                                 const std::vector<double>& b,
                                 const std::vector<double>& c,
                                 const std::vector<double>& d,
                                 const std::vector<double>& e,
                                 const std::vector<double>& rhs,
                                 const std::vector<double>& x) const;
    
    bool is_diagonally_dominant(const std::vector<double>& a,
                                const std::vector<double>& b,
                                const std::vector<double>& c,
                                const std::vector<double>& d,
                                const std::vector<double>& e) const;
    
    void set_tolerance(double tol) { tolerance_ = tol; }
    double tolerance() const { return tolerance_; }
    
    void set_max_iterations(int iter) { max_iter_ = iter; }
    int max_iterations() const { return max_iter_; }

private:
    int n_;
    int max_size_;
    double tolerance_;
    int max_iter_;
    
    std::vector<double> alpha_;
    std::vector<double> beta_;
    std::vector<double> gamma_;
    std::vector<double> delta_;
    std::vector<double> mu_;
    std::vector<double> z_;
    std::vector<double> work1_;
    std::vector<double> work2_;
    std::vector<double> work3_;
    
    static constexpr double eps_ = 1e-30;
    
    void forward_elimination(const std::vector<double>& a,
                             const std::vector<double>& b,
                             const std::vector<double>& c,
                             const std::vector<double>& d,
                             const std::vector<double>& e,
                             std::vector<double>& rhs);
    
    void back_substitution(const std::vector<double>& a,
                           const std::vector<double>& b,
                           const std::vector<double>& c,
                           const std::vector<double>& d,
                           const std::vector<double>& e,
                           std::vector<double>& x);
};

}
