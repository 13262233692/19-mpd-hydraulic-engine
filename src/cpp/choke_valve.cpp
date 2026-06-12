#include "hydraulic_engine/choke_valve.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

namespace hydraulic_engine {

ChokeValveConfig::ChokeValveConfig()
    : discharge_coefficient(0.85),
      max_area(0.01),
      min_area(1e-6),
      initial_area(0.005),
      max_area_change_rate(0.1),
      target_bhp(40e6),
      bhp_tolerance(1e5),
      kp(1e-7),
      ki(5e-6),
      kd(1e-8),
      integral_min(-1.0),
      integral_max(1.0),
      derivative_filter_tau(0.1) {}

ChokeValveController::ChokeValveController(const ChokeValveConfig& config)
    : config_(config),
      current_area_(config.initial_area),
      integral_error_(0.0),
      prev_error_(0.0),
      prev_time_(0.0),
      filtered_derivative_(0.0),
      initialized_(false) {
    if (config_.max_area <= 0.0) {
        throw std::invalid_argument("Max choke area must be positive");
    }
    if (config_.min_area <= 0.0) {
        throw std::invalid_argument("Min choke area must be positive");
    }
    if (config_.discharge_coefficient <= 0.0) {
        throw std::invalid_argument("Discharge coefficient must be positive");
    }
}

void ChokeValveController::reset() {
    current_area_ = config_.initial_area;
    integral_error_ = 0.0;
    prev_error_ = 0.0;
    prev_time_ = 0.0;
    filtered_derivative_ = 0.0;
    initialized_ = false;
    error_history_.clear();
    control_log_.clear();
}

double ChokeValveController::get_choke_opening_percent() const {
    return (current_area_ - config_.min_area) / 
           (config_.max_area - config_.min_area) * 100.0;
}

double ChokeValveController::compute_backpressure(double flow_rate, 
                                                   double mixture_density, 
                                                   double void_fraction) const {
    if (current_area_ < eps_ || mixture_density < eps_) {
        return 1e9;
    }
    
    double effective_density = mixture_density;
    if (void_fraction > eps_ && void_fraction < 1.0 - eps_) {
        effective_density = mixture_density;
    }
    
    double velocity = flow_rate / current_area_;
    double dynamic_head = 0.5 * effective_density * velocity * velocity;
    
    double delta_p = dynamic_head / (config_.discharge_coefficient * config_.discharge_coefficient);
    
    return delta_p;
}

double ChokeValveController::compute_flow_rate(double delta_p, 
                                                double mixture_density, 
                                                double void_fraction) const {
    if (delta_p < 0.0) delta_p = 0.0;
    if (mixture_density < eps_) return 0.0;
    
    double effective_density = mixture_density;
    
    double velocity = std::sqrt(2.0 * delta_p * config_.discharge_coefficient * config_.discharge_coefficient 
                                / effective_density);
    
    return current_area_ * velocity;
}

double ChokeValveController::compute_pid(double error, double dt) {
    if (dt <= eps_) dt = eps_;
    
    double p_term = config_.kp * error;
    
    integral_error_ += error * dt;
    integral_error_ = std::max(config_.integral_min, 
                               std::min(config_.integral_max, integral_error_));
    double i_term = config_.ki * integral_error_;
    
    double raw_derivative = (error - prev_error_) / dt;
    
    double alpha = dt / (config_.derivative_filter_tau + dt);
    filtered_derivative_ = alpha * raw_derivative + (1.0 - alpha) * filtered_derivative_;
    
    double d_term = config_.kd * filtered_derivative_;
    
    prev_error_ = error;
    
    return p_term + i_term + d_term;
}

double ChokeValveController::saturate_area(double area) const {
    return std::max(config_.min_area, std::min(config_.max_area, area));
}

double ChokeValveController::limit_area_rate(double new_area, double old_area, double dt) const {
    if (dt <= eps_) return old_area;
    
    double max_change = config_.max_area_change_rate * config_.max_area * dt;
    double change = new_area - old_area;
    
    if (change > max_change) {
        return old_area + max_change;
    } else if (change < -max_change) {
        return old_area - max_change;
    }
    
    return new_area;
}

void ChokeValveController::update(double current_time,
                                   double bhp,
                                   double wellhead_flow_rate,
                                   double wellhead_void_fraction,
                                   double wellhead_pressure,
                                   double& out_choke_area,
                                   double& out_backpressure) {
    double dt = 0.0;
    if (initialized_) {
        dt = current_time - prev_time_;
    } else {
        initialized_ = true;
        dt = config_.derivative_filter_tau * 0.1;
    }
    
    prev_time_ = current_time;
    
    double error = config_.target_bhp - bhp;
    
    double p_term = config_.kp * error;
    double i_term = 0.0;
    double d_term = 0.0;
    double control_output = 0.0;
    
    double desired_area = current_area_;
    
    if (std::abs(error) > config_.bhp_tolerance) {
        control_output = compute_pid(error, dt);
        
        p_term = config_.kp * error;
        i_term = config_.ki * integral_error_;
        d_term = config_.kd * filtered_derivative_;
        
        double area_change = -control_output * config_.max_area;
        desired_area = current_area_ + area_change;
        
        desired_area = limit_area_rate(desired_area, current_area_, dt);
        desired_area = saturate_area(desired_area);
        
        current_area_ = desired_area;
    }
    
    double backpressure = compute_backpressure(wellhead_flow_rate, 
                                               1000.0 * (1.0 - wellhead_void_fraction) + 
                                               150.0 * wellhead_void_fraction,
                                               wellhead_void_fraction);
    
    ControlLogEntry entry;
    entry.time = current_time;
    entry.bhp = bhp;
    entry.bhp_target = config_.target_bhp;
    entry.bhp_error = error;
    entry.choke_area = current_area_;
    entry.choke_opening_pct = get_choke_opening_percent();
    entry.wellhead_flow_rate = wellhead_flow_rate;
    entry.wellhead_void_fraction = wellhead_void_fraction;
    entry.wellhead_pressure = wellhead_pressure;
    entry.p_term = p_term;
    entry.i_term = i_term;
    entry.d_term = d_term;
    entry.control_output = control_output;
    
    control_log_.push_back(entry);
    
    out_choke_area = current_area_;
    out_backpressure = backpressure;
}

void ChokeValveController::print_report() const {
    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "  CHOKE VALVE CONTROL REPORT" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;
    
    if (control_log_.empty()) {
        std::cout << "  No control data available." << std::endl;
        std::cout << std::endl;
        return;
    }
    
    const auto& first = control_log_.front();
    const auto& last = control_log_.back();
    
    std::cout << "  Simulation time: " << std::fixed << std::setprecision(3) 
              << last.time << " s" << std::endl;
    std::cout << "  Target BHP: " << config_.target_bhp / 1e6 
              << " MPa" << std::endl;
    std::cout << "  BHP tolerance: " << config_.bhp_tolerance / 1e6 
              << " MPa" << std::endl;
    std::cout << std::endl;
    
    std::cout << "  Initial conditions:" << std::endl;
    std::cout << "    BHP: " << first.bhp / 1e6 << " MPa" << std::endl;
    std::cout << "    Choke opening: " << first.choke_opening_pct << " %" << std::endl;
    std::cout << "    Choke area: " << first.choke_area << " m^2" << std::endl;
    std::cout << std::endl;
    
    std::cout << "  Final conditions:" << std::endl;
    std::cout << "    BHP: " << last.bhp / 1e6 << " MPa" << std::endl;
    std::cout << "    BHP error: " << last.bhp_error / 1e6 << " MPa" << std::endl;
    std::cout << "    Choke opening: " << last.choke_opening_pct << " %" << std::endl;
    std::cout << "    Choke area: " << last.choke_area << " m^2" << std::endl;
    std::cout << "    Wellhead flow rate: " << last.wellhead_flow_rate << " m^3/s" << std::endl;
    std::cout << "    Wellhead void fraction: " << last.wellhead_void_fraction << std::endl;
    std::cout << "    Wellhead pressure: " << last.wellhead_pressure / 1e6 << " MPa" << std::endl;
    std::cout << std::endl;
    
    double max_bhp = -1e99;
    double min_bhp = 1e99;
    double max_error = 0.0;
    double avg_error = 0.0;
    
    for (const auto& entry : control_log_) {
        max_bhp = std::max(max_bhp, entry.bhp);
        min_bhp = std::min(min_bhp, entry.bhp);
        max_error = std::max(max_error, std::abs(entry.bhp_error));
        avg_error += std::abs(entry.bhp_error);
    }
    avg_error /= control_log_.size();
    
    std::cout << "  Control performance:" << std::endl;
    std::cout << "    Max BHP: " << max_bhp / 1e6 << " MPa" << std::endl;
    std::cout << "    Min BHP: " << min_bhp / 1e6 << " MPa" << std::endl;
    std::cout << "    Max absolute error: " << max_error / 1e6 << " MPa" << std::endl;
    std::cout << "    Average absolute error: " << avg_error / 1e6 << " MPa" << std::endl;
    std::cout << "    Number of control steps: " << control_log_.size() << std::endl;
    std::cout << std::endl;
    
    int within_tolerance = 0;
    for (const auto& entry : control_log_) {
        if (std::abs(entry.bhp_error) <= config_.bhp_tolerance) {
            within_tolerance++;
        }
    }
    double pct_within = 100.0 * within_tolerance / control_log_.size();
    
    std::cout << "  Time within tolerance: " << pct_within << " %" << std::endl;
    std::cout << std::endl;
    
    std::cout << "  Controller gains:" << std::endl;
    std::cout << "    Kp: " << config_.kp << std::endl;
    std::cout << "    Ki: " << config_.ki << std::endl;
    std::cout << "    Kd: " << config_.kd << std::endl;
    std::cout << std::endl;
    
    std::cout << "  Time history (every 5th entry):" << std::endl;
    std::cout << "  " << std::setw(10) << "Time (s)"
              << std::setw(15) << "BHP (MPa)"
              << std::setw(15) << "Error (MPa)"
              << std::setw(15) << "Choke (%)"
              << std::setw(15) << "Flow (m3/s)" << std::endl;
    std::cout << "  " << std::string(70, '-') << std::endl;
    
    size_t step = std::max(size_t(1), control_log_.size() / 20);
    for (size_t i = 0; i < control_log_.size(); i += step) {
        const auto& e = control_log_[i];
        std::cout << "  " << std::setw(10) << std::fixed << std::setprecision(3) << e.time
                  << std::setw(15) << std::setprecision(4) << e.bhp / 1e6
                  << std::setw(15) << std::setprecision(4) << e.bhp_error / 1e6
                  << std::setw(15) << std::setprecision(2) << e.choke_opening_pct
                  << std::setw(15) << std::setprecision(4) << e.wellhead_flow_rate 
                  << std::endl;
    }
    
    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;
}

}
