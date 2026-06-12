#include "hydraulic_engine/pentadiagonal_solver.hpp"
#include <iostream>
#include <algorithm>

namespace hydraulic_engine {

PentadiagonalSolver::PentadiagonalSolver(int max_size)
    : n_(0),
      max_size_(max_size),
      tolerance_(1e-12),
      max_iter_(100) {
    if (max_size_ > 0) {
        resize(max_size_);
    }
}

void PentadiagonalSolver::resize(int n) {
    if (n <= 0) {
        throw std::invalid_argument("Size must be positive");
    }
    n_ = n;
    max_size_ = std::max(max_size_, n);
    
    alpha_.resize(n);
    beta_.resize(n);
    gamma_.resize(n);
    delta_.resize(n);
    mu_.resize(n);
    z_.resize(n);
    work1_.resize(n);
    work2_.resize(n);
    work3_.resize(n);
}

bool PentadiagonalSolver::is_diagonally_dominant(const std::vector<double>& a,
                                                 const std::vector<double>& b,
                                                 const std::vector<double>& c,
                                                 const std::vector<double>& d,
                                                 const std::vector<double>& e) const {
    for (int i = 0; i < n_; ++i) {
        double sum = 0.0;
        if (i >= 2) sum += std::abs(a[i]);
        if (i >= 1) sum += std::abs(b[i]);
        if (i < n_ - 1) sum += std::abs(d[i]);
        if (i < n_ - 2) sum += std::abs(e[i]);
        
        if (std::abs(c[i]) < sum - eps_) {
            return false;
        }
    }
    return true;
}

double PentadiagonalSolver::compute_residual_norm(const std::vector<double>& a,
                                                  const std::vector<double>& b,
                                                  const std::vector<double>& c,
                                                  const std::vector<double>& d,
                                                  const std::vector<double>& e,
                                                  const std::vector<double>& rhs,
                                                  const std::vector<double>& x) const {
    double norm = 0.0;
    for (int i = 0; i < n_; ++i) {
        double residual = rhs[i];
        if (i >= 2) residual -= a[i] * x[i - 2];
        if (i >= 1) residual -= b[i] * x[i - 1];
        residual -= c[i] * x[i];
        if (i < n_ - 1) residual -= d[i] * x[i + 1];
        if (i < n_ - 2) residual -= e[i] * x[i + 2];
        norm += residual * residual;
    }
    return std::sqrt(norm / n_);
}

bool PentadiagonalSolver::solve(const std::vector<double>& a,
                                const std::vector<double>& b,
                                const std::vector<double>& c,
                                const std::vector<double>& d,
                                const std::vector<double>& e,
                                const std::vector<double>& rhs,
                                std::vector<double>& x) {
    return solve_modified_thomas(a, b, c, d, e, rhs, x);
}

bool PentadiagonalSolver::solve_modified_thomas(const std::vector<double>& a,
                                                const std::vector<double>& b,
                                                const std::vector<double>& c,
                                                const std::vector<double>& d,
                                                const std::vector<double>& e,
                                                const std::vector<double>& rhs,
                                                std::vector<double>& x) {
    if (static_cast<int>(a.size()) != n_ || static_cast<int>(b.size()) != n_ ||
        static_cast<int>(c.size()) != n_ || static_cast<int>(d.size()) != n_ ||
        static_cast<int>(e.size()) != n_ || static_cast<int>(rhs.size()) != n_) {
        throw std::invalid_argument("Vector size mismatch");
    }
    
    x.assign(n_, 0.0);
    
    std::vector<double> cp = c;
    std::vector<double> dp = d;
    std::vector<double> ep = e;
    std::vector<double> rhsp = rhs;
    
    if (std::abs(cp[0]) < eps_) {
        return false;
    }
    
    dp[0] /= cp[0];
    ep[0] /= cp[0];
    rhsp[0] /= cp[0];
    
    if (n_ > 1) {
        double m = b[1];
        if (std::abs(cp[1] - m * dp[0]) < eps_) {
            return false;
        }
        cp[1] = cp[1] - m * dp[0];
        dp[1] = (dp[1] - m * ep[0]) / cp[1];
        ep[1] /= cp[1];
        rhsp[1] = (rhsp[1] - m * rhsp[0]) / cp[1];
    }
    
    for (int i = 2; i < n_; ++i) {
        double m1 = b[i] - a[i] * dp[i - 2];
        double m2 = a[i];
        
        double denom = cp[i] - m1 * dp[i - 1] - m2 * ep[i - 2];
        if (std::abs(denom) < eps_) {
            return false;
        }
        
        ep[i] /= denom;
        dp[i] = (dp[i] - m1 * ep[i - 1]) / denom;
        rhsp[i] = (rhsp[i] - m1 * rhsp[i - 1] - m2 * rhsp[i - 2]) / denom;
        cp[i] = 1.0;
    }
    
    x[n_ - 1] = rhsp[n_ - 1];
    if (n_ > 1) {
        x[n_ - 2] = rhsp[n_ - 2] - dp[n_ - 2] * x[n_ - 1];
    }
    
    for (int i = n_ - 3; i >= 0; --i) {
        x[i] = rhsp[i] - dp[i] * x[i + 1] - ep[i] * x[i + 2];
    }
    
    double residual = compute_residual_norm(a, b, c, d, e, rhs, x);
    return residual < tolerance_;
}

bool PentadiagonalSolver::solve_cyclic_reduction(const std::vector<double>& a,
                                                 const std::vector<double>& b,
                                                 const std::vector<double>& c,
                                                 const std::vector<double>& d,
                                                 const std::vector<double>& e,
                                                 const std::vector<double>& rhs,
                                                 std::vector<double>& x) {
    if (static_cast<int>(a.size()) != n_ || static_cast<int>(b.size()) != n_ ||
        static_cast<int>(c.size()) != n_ || static_cast<int>(d.size()) != n_ ||
        static_cast<int>(e.size()) != n_ || static_cast<int>(rhs.size()) != n_) {
        throw std::invalid_argument("Vector size mismatch");
    }
    
    x.assign(n_, 0.0);
    
    std::vector<double> a_curr = a;
    std::vector<double> b_curr = b;
    std::vector<double> c_curr = c;
    std::vector<double> d_curr = d;
    std::vector<double> e_curr = e;
    std::vector<double> rhs_curr = rhs;
    
    std::vector<double> a_next(n_), b_next(n_), c_next(n_);
    std::vector<double> d_next(n_), e_next(n_), rhs_next(n_);
    
    int n_curr = n_;
    int level = 0;
    const int max_levels = 20;
    
    std::vector<std::vector<double>> x_levels(max_levels);
    std::vector<int> n_levels(max_levels);
    
    while (n_curr > 2 && level < max_levels) {
        x_levels[level].resize(n_curr);
        n_levels[level] = n_curr;
        
        int n_next = (n_curr + 1) / 2;
        
        a_next.assign(n_next, 0.0);
        b_next.assign(n_next, 0.0);
        c_next.assign(n_next, 0.0);
        d_next.assign(n_next, 0.0);
        e_next.assign(n_next, 0.0);
        rhs_next.assign(n_next, 0.0);
        
        for (int i = 0; i < n_curr; i += 2) {
            int ii = i / 2;
            
            double diag = c_curr[i];
            if (std::abs(diag) < eps_) return false;
            
            double factor_b = (i >= 1) ? -b_curr[i] / diag : 0.0;
            double factor_a = (i >= 2) ? -a_curr[i] / diag : 0.0;
            double factor_d = (i < n_curr - 1) ? -d_curr[i] / diag : 0.0;
            double factor_e = (i < n_curr - 2) ? -e_curr[i] / diag : 0.0;
            
            if (i >= 2 && i - 2 >= 0) {
                a_next[ii] += factor_a * a_curr[i - 2];
            }
            if (i >= 1 && i - 1 >= 0) {
                if (ii - 1 >= 0) {
                    b_next[ii] += factor_b * b_curr[i - 1];
                }
                a_next[ii] += factor_b * a_curr[i - 1];
            }
            
            c_next[ii] = 1.0;
            if (i >= 2) c_next[ii] += factor_a * e_curr[i - 2];
            if (i >= 1) c_next[ii] += factor_b * d_curr[i - 1];
            if (i < n_curr - 1) c_next[ii] += factor_d * b_curr[i + 1];
            if (i < n_curr - 2) c_next[ii] += factor_e * a_curr[i + 2];
            
            if (i < n_curr - 1 && i + 1 < n_curr) {
                if (ii + 1 < n_next) {
                    d_next[ii] += factor_d * d_curr[i + 1];
                }
                e_next[ii] += factor_d * e_curr[i + 1];
            }
            if (i < n_curr - 2 && i + 2 < n_curr) {
                if (ii + 1 < n_next) {
                    d_next[ii] += factor_e * b_curr[i + 2];
                }
                e_next[ii] += factor_e * c_curr[i + 2];
            }
            
            rhs_next[ii] = rhs_curr[i] / diag;
            if (i >= 2) rhs_next[ii] += factor_a * rhs_curr[i - 2];
            if (i >= 1) rhs_next[ii] += factor_b * rhs_curr[i - 1];
            if (i < n_curr - 1) rhs_next[ii] += factor_d * rhs_curr[i + 1];
            if (i < n_curr - 2) rhs_next[ii] += factor_e * rhs_curr[i + 2];
        }
        
        a_curr.swap(a_next);
        b_curr.swap(b_next);
        c_curr.swap(c_next);
        d_curr.swap(d_next);
        e_curr.swap(e_next);
        rhs_curr.swap(rhs_next);
        n_curr = n_next;
        level++;
    }
    
    if (n_curr <= 2) {
        if (n_curr == 1) {
            x_levels[level - 1][0] = rhs_curr[0] / c_curr[0];
        } else {
            double det = c_curr[0] * c_curr[1] - d_curr[0] * b_curr[1];
            if (std::abs(det) < eps_) return false;
            x_levels[level - 1][0] = (rhs_curr[0] * c_curr[1] - rhs_curr[1] * d_curr[0]) / det;
            x_levels[level - 1][1] = (c_curr[0] * rhs_curr[1] - b_curr[1] * rhs_curr[0]) / det;
        }
    }
    
    for (int l = level - 1; l > 0; --l) {
        int n_prev = n_levels[l - 1];
        int n_curr_l = n_levels[l];
        
        for (int i = 0; i < n_prev; i += 2) {
            int ii = i / 2;
            if (ii < n_curr_l) {
                x_levels[l - 1][i] = x_levels[l][ii];
            }
        }
        
        for (int i = 1; i < n_prev; i += 2) {
            double val = rhs_curr[i];
            if (i >= 2) val -= a_curr[i] * x_levels[l - 1][i - 2];
            if (i >= 1) val -= b_curr[i] * x_levels[l - 1][i - 1];
            if (i < n_prev - 1) val -= d_curr[i] * x_levels[l - 1][i + 1];
            if (i < n_prev - 2) val -= e_curr[i] * x_levels[l - 1][i + 2];
            if (std::abs(c_curr[i]) < eps_) return false;
            x_levels[l - 1][i] = val / c_curr[i];
        }
    }
    
    for (int i = 0; i < n_; ++i) {
        x[i] = x_levels[0][i];
    }
    
    double residual = compute_residual_norm(a, b, c, d, e, rhs, x);
    return residual < tolerance_;
}

bool PentadiagonalSolver::solve_parallel(const std::vector<double>& a,
                                         const std::vector<double>& b,
                                         const std::vector<double>& c,
                                         const std::vector<double>& d,
                                         const std::vector<double>& e,
                                         const std::vector<double>& rhs,
                                         std::vector<double>& x) {
    if (static_cast<int>(a.size()) != n_ || static_cast<int>(b.size()) != n_ ||
        static_cast<int>(c.size()) != n_ || static_cast<int>(d.size()) != n_ ||
        static_cast<int>(e.size()) != n_ || static_cast<int>(rhs.size()) != n_) {
        throw std::invalid_argument("Vector size mismatch");
    }
    
    x.assign(n_, 0.0);
    
    int num_threads = omp_get_max_threads();
    int block_size = (n_ + num_threads - 1) / num_threads;
    block_size = std::max(block_size, 4);
    num_threads = (n_ + block_size - 1) / block_size;
    
    std::vector<int> block_start(num_threads);
    std::vector<int> block_end(num_threads);
    
    for (int t = 0; t < num_threads; ++t) {
        block_start[t] = t * block_size;
        block_end[t] = std::min((t + 1) * block_size, n_);
    }
    
    std::vector<std::vector<double>> x_blocks(num_threads);
    std::vector<bool> success(num_threads, true);
    
    #pragma omp parallel for
    for (int t = 0; t < num_threads; ++t) {
        int start = block_start[t];
        int end = block_end[t];
        int size = end - start;
        
        if (size <= 0) {
            success[t] = true;
            continue;
        }
        
        std::vector<double> a_local(size, 0.0);
        std::vector<double> b_local(size, 0.0);
        std::vector<double> c_local(size);
        std::vector<double> d_local(size, 0.0);
        std::vector<double> e_local(size, 0.0);
        std::vector<double> rhs_local(size);
        std::vector<double> x_local(size, 0.0);
        
        for (int i = 0; i < size; ++i) {
            int global_i = start + i;
            c_local[i] = c[global_i];
            rhs_local[i] = rhs[global_i];
            
            if (i >= 2) a_local[i] = a[global_i];
            if (i >= 1) b_local[i] = b[global_i];
            if (i < size - 1) d_local[i] = d[global_i];
            if (i < size - 2) e_local[i] = e[global_i];
        }
        
        PentadiagonalSolver local_solver(size);
        local_solver.set_tolerance(tolerance_);
        success[t] = local_solver.solve_modified_thomas(
            a_local, b_local, c_local, d_local, e_local, rhs_local, x_local
        );
        
        x_blocks[t] = x_local;
    }
    
    for (int t = 0; t < num_threads; ++t) {
        if (!success[t]) return false;
        int start = block_start[t];
        for (int i = 0; i < static_cast<int>(x_blocks[t].size()); ++i) {
            x[start + i] = x_blocks[t][i];
        }
    }
    
    for (int iter = 0; iter < 3; ++iter) {
        std::vector<double> x_new = x;
        
        #pragma omp parallel for
        for (int i = 0; i < n_; ++i) {
            double sum = rhs[i];
            if (i >= 2) sum -= a[i] * x[i - 2];
            if (i >= 1) sum -= b[i] * x[i - 1];
            if (i < n_ - 1) sum -= d[i] * x[i + 1];
            if (i < n_ - 2) sum -= e[i] * x[i + 2];
            if (std::abs(c[i]) < eps_) {
                x_new[i] = 0.0;
            } else {
                x_new[i] = sum / c[i];
            }
        }
        
        x.swap(x_new);
    }
    
    double residual = compute_residual_norm(a, b, c, d, e, rhs, x);
    return residual < tolerance_;
}

void PentadiagonalSolver::forward_elimination(const std::vector<double>& a,
                                              const std::vector<double>& b,
                                              const std::vector<double>& c,
                                              const std::vector<double>& d,
                                              const std::vector<double>& e,
                                              std::vector<double>& rhs) {
    work1_[0] = c[0];
    work2_[0] = d[0] / work1_[0];
    work3_[0] = e[0] / work1_[0];
    rhs[0] = rhs[0] / work1_[0];
    
    if (n_ > 1) {
        double m = b[1];
        work1_[1] = c[1] - m * work2_[0];
        work2_[1] = (d[1] - m * work3_[0]) / work1_[1];
        work3_[1] = e[1] / work1_[1];
        rhs[1] = (rhs[1] - m * rhs[0]) / work1_[1];
    }
    
    for (int i = 2; i < n_; ++i) {
        double m1 = b[i] - a[i] * work2_[i - 2];
        double m2 = a[i];
        work1_[i] = c[i] - m1 * work2_[i - 1] - m2 * work3_[i - 2];
        work2_[i] = (d[i] - m1 * work3_[i - 1]) / work1_[i];
        work3_[i] = e[i] / work1_[i];
        rhs[i] = (rhs[i] - m1 * rhs[i - 1] - m2 * rhs[i - 2]) / work1_[i];
    }
}

void PentadiagonalSolver::back_substitution(const std::vector<double>& a,
                                            const std::vector<double>& b,
                                            const std::vector<double>& c,
                                            const std::vector<double>& d,
                                            const std::vector<double>& e,
                                            std::vector<double>& x) {
    x[n_ - 1] = work1_[n_ - 1];
    if (n_ > 1) {
        x[n_ - 2] = work1_[n_ - 2] - work2_[n_ - 2] * x[n_ - 1];
    }
    
    for (int i = n_ - 3; i >= 0; --i) {
        x[i] = work1_[i] - work2_[i] * x[i + 1] - work3_[i] * x[i + 2];
    }
}

}
