#pragma once

#include <vector>
#include <string>
#include <deque>

namespace hydraulic_engine {

struct ChokeValveConfig {
    double discharge_coefficient;
    double max_area;
    double min_area;
    double initial_area;
    double max_area_change_rate;
    
    double target_bhp;
    double bhp_tolerance;
    
    double kp;
    double ki;
    double kd;
    
    double integral_min;
    double integral_max;
    
    double derivative_filter_tau;
    
    ChokeValveConfig();
};

struct ControlLogEntry {
    double time;
    double bhp;
    double bhp_target;
    double bhp_error;
    double choke_area;
    double choke_opening_pct;
    double wellhead_flow_rate;
    double wellhead_void_fraction;
    double wellhead_pressure;
    double p_term;
    double i_term;
    double d_term;
    double control_output;
};

class ChokeValveController {
public:
    explicit ChokeValveController(const ChokeValveConfig& config);
    
    void reset();
    
    double compute_backpressure(double flow_rate, double mixture_density, double void_fraction) const;
    
    double compute_flow_rate(double delta_p, double mixture_density, double void_fraction) const;
    
    void update(double current_time,
                double bhp,
                double wellhead_flow_rate,
                double wellhead_void_fraction,
                double wellhead_pressure,
                double& out_choke_area,
                double& out_backpressure);
    
    void set_target_bhp(double target) { config_.target_bhp = target; }
    double get_target_bhp() const { return config_.target_bhp; }
    
    double get_choke_area() const { return current_area_; }
    double get_choke_opening_percent() const;
    
    const std::vector<ControlLogEntry>& get_log() const { return control_log_; }
    
    void clear_log() { control_log_.clear(); }
    
    void print_report() const;
    
    const ChokeValveConfig& config() const { return config_; }

private:
    ChokeValveConfig config_;
    
    double current_area_;
    double integral_error_;
    double prev_error_;
    double prev_time_;
    double filtered_derivative_;
    
    std::deque<double> error_history_;
    
    std::vector<ControlLogEntry> control_log_;
    
    bool initialized_;
    
    static constexpr double eps_ = 1e-12;
    
    double compute_pid(double error, double dt);
    
    double saturate_area(double area) const;
    
    double limit_area_rate(double new_area, double old_area, double dt) const;
};

}
