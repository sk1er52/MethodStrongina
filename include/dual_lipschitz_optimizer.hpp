#ifndef DUAL_LIPSCHITZ_OPTIMIZER_HPP
#define DUAL_LIPSCHITZ_OPTIMIZER_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <stdexcept>
#include <iomanip>

#include "common_types.hpp"
#include "peano_mapper.hpp"

template <typename TProblemType> 
class DualLipschitzOptimizer {
public:
    DualLipschitzOptimizer(const std::vector<double>& domain_left,
                           const std::vector<double>& domain_right, 
                           double epsilon, double r_loc, double r_glob, 
                           const TProblemType& problem, std::ostream& logger, int max_iterations)
        : search_left_bound_(domain_left), search_right_bound_(domain_right),
          epsilon_(epsilon), r_loc_(r_loc), r_glob_(r_glob),
          r_selected_for_next_trial_(r_glob), problem_(problem), logger_(logger),
          max_iterations_(max_iterations), iteration_counter_(0),
          exit_main_criteria_count_(0), exit_test_criteria_count_(0),
          use_peano_mapping_(false), problem_dimension_(1), current_max_mu_estimate_(1.0) {
        if (search_left_bound_.empty() || search_right_bound_.empty()) {
            throw std::invalid_argument("Границы задачи (1D) не могут быть пустыми.");
        }
    }

    DualLipschitzOptimizer(double peano_domain_left, double peano_domain_right,
                           double epsilon, double r_loc, double r_glob, 
                           const TProblemType& problem, int peano_m_order, 
                           int original_dimension, int peano_key_type, 
                           std::ostream& logger, int max_iterations)
        : search_left_bound_({peano_domain_left}), search_right_bound_({peano_domain_right}),
          epsilon_(epsilon), r_loc_(r_loc), r_glob_(r_glob),
          r_selected_for_next_trial_(r_glob), problem_(problem), logger_(logger),
          max_iterations_(max_iterations), iteration_counter_(0),
          exit_main_criteria_count_(0), exit_test_criteria_count_(0),
          use_peano_mapping_(true), peano_m_order_(peano_m_order),
          problem_dimension_(original_dimension), peano_key_type_(peano_key_type),
          current_max_mu_estimate_(1.0) {}

    std::vector<double> Optimize() {
        trial_points_.clear();
        iteration_counter_ = 0;
        exit_main_criteria_count_ = 0;
        exit_test_criteria_count_ = 0;

        PerformFirstIteration();

        if (trial_points_.size() < 2) {
            return trial_points_.empty() ? std::vector<double>{} : trial_points_[0].y_coords;
        }

        for (int iter = 0; iter < max_iterations_; ++iter) {
            iteration_counter_ = iter + 1;

            UpdateLipschitzEstimate();

            size_t interval_idx = FindIntervalIndexWithMaxR();

            if (interval_idx + 1 >= trial_points_.size()) break;

            if (CheckStoppingConditions(trial_points_[interval_idx], trial_points_[interval_idx + 1])) {
                break;
            }

            TrialPoint new_point = GenerateNewTrialPoint(interval_idx);
            InsertPointUnique(new_point);
        }

        return GetBestPointCoords();
    }

    int GetIterationCount() const { return iteration_counter_; }
    int GetExitMainCount() const { return exit_main_criteria_count_; }
    int GetExitTestCount() const { return exit_test_criteria_count_; }
    const std::vector<TrialPoint>& GetTrialPointsHistory() const { return trial_points_; }

private:
    PeanoMapper peano_mapper_;
    std::vector<double> search_left_bound_;
    std::vector<double> search_right_bound_;
    double epsilon_;
    double r_loc_;
    double r_glob_;
    double r_selected_for_next_trial_;
    double current_max_mu_estimate_;

    const TProblemType& problem_;
    std::ostream& logger_;

    int max_iterations_;
    int iteration_counter_;
    int exit_main_criteria_count_;
    int exit_test_criteria_count_;

    std::vector<TrialPoint> trial_points_;

    bool use_peano_mapping_;
    int peano_m_order_;
    int problem_dimension_;
    int peano_key_type_;

    double GetHolderDelta(double x_left, double x_right) const {
        double dx = std::abs(x_right - x_left);
        if (!use_peano_mapping_ || problem_dimension_ <= 1) return dx;
        return dx < 1e-15 ? 1e-15 : std::pow(dx, 1.0 / problem_dimension_);
    }

    void PerformFirstIteration() {
        TrialPoint p_left = CreatePoint(search_left_bound_[0]);
        TrialPoint p_right = CreatePoint(search_right_bound_[0]);
        trial_points_.push_back(p_left);
        trial_points_.push_back(p_right);

        if (use_peano_mapping_ && search_left_bound_[0] <= 0.5 && search_right_bound_[0] >= 0.5) {
            TrialPoint center = CreatePoint(0.5);
            trial_points_.push_back(center);
            std::sort(trial_points_.begin(), trial_points_.end());
        }
    }

    TrialPoint CreatePoint(double x_param) {
        TrialPoint p;
        p.x_param = x_param;
        if (!use_peano_mapping_) {
            p.y_coords = {x_param};
        } else {
            p.y_coords = peano_mapper_.Map(x_param, peano_m_order_, problem_dimension_, peano_key_type_);
        }
        p.z_value = problem_.ComputeFunction(p.y_coords);
        return p;
    }

    void UpdateLipschitzEstimate() {
        current_max_mu_estimate_ = 0.0;
        for (size_t i = 0; i < trial_points_.size() - 1; ++i) {
            double delta = GetHolderDelta(trial_points_[i].x_param, trial_points_[i + 1].x_param);
            if (delta > 1e-12) {
                double mu = std::abs(trial_points_[i + 1].z_value - trial_points_[i].z_value) / delta;
                if (mu > current_max_mu_estimate_) {
                    current_max_mu_estimate_ = mu;
                }
            }
        }
        if (current_max_mu_estimate_ < 1e-9) {
            current_max_mu_estimate_ = 1.0;
        }
    }

    // Вспомогательный метод для получения z* (текущего минимума)
    double GetCurrentZStar() const {
        double z_star = std::numeric_limits<double>::infinity();
        for (const auto& p : trial_points_) {
            if (!std::isnan(p.z_value) && p.z_value < z_star) {
                z_star = p.z_value;
            }
        }
        return z_star;
    }

    // ИСПРАВЛЕННАЯ ФОРМУЛА ИЗ СТАТЬИ
    double ComputeIntervalCharacteristicR(const TrialPoint& p1, const TrialPoint& p2, double z_star) const {
        double delta_i = GetHolderDelta(p1.x_param, p2.x_param);
        if (delta_i < 1e-12) return -std::numeric_limits<double>::infinity();
        
        double dz = p2.z_value - p1.z_value;
        double mu_v = std::max(current_max_mu_estimate_, 1e-9);

        // Строго по формулам со страницы 5
        double r_mu_glob = r_glob_ * mu_v;
        double r_glob_val = delta_i + (dz * dz) / (r_mu_glob * r_mu_glob * delta_i) 
                            - 2.0 * (p1.z_value + p2.z_value - 2.0 * z_star) / r_mu_glob;

        double r_mu_loc = r_loc_ * mu_v;
        double r_loc_val = delta_i + (dz * dz) / (r_mu_loc * r_mu_loc * delta_i) 
                           - 2.0 * (p1.z_value + p2.z_value - 2.0 * z_star) / r_mu_loc;

        double num = 1.0 - 1.0 / r_glob_;
        double den = 1.0 - 1.0 / r_loc_;
        double rho = (num * num) / (den * den);

        return std::max(rho * r_loc_val, r_glob_val);
    }

    size_t FindIntervalIndexWithMaxR() {
        double z_star = GetCurrentZStar(); // Находим z* для текущей итерации
        
        double max_r_val = -std::numeric_limits<double>::infinity();
        size_t max_r_idx = 0;

        for (size_t i = 0; i < trial_points_.size() - 1; ++i) {
            double current_r_val = ComputeIntervalCharacteristicR(trial_points_[i], trial_points_[i + 1], z_star);
            if (current_r_val > max_r_val) {
                max_r_val = current_r_val;
                max_r_idx = i;
            }
        }

        // Определяем, какой r победил для генерации новой точки
        const TrialPoint& p1 = trial_points_[max_r_idx];
        const TrialPoint& p2 = trial_points_[max_r_idx + 1];
        double delta_i = GetHolderDelta(p1.x_param, p2.x_param);
        double dz = p2.z_value - p1.z_value;
        double mu_v = std::max(current_max_mu_estimate_, 1e-9);

        double r_mu_glob = r_glob_ * mu_v;
        double r_glob_val = delta_i + (dz * dz) / (r_mu_glob * r_mu_glob * delta_i) 
                            - 2.0 * (p1.z_value + p2.z_value - 2.0 * z_star) / r_mu_glob;

        double r_mu_loc = r_loc_ * mu_v;
        double r_loc_val = delta_i + (dz * dz) / (r_mu_loc * r_mu_loc * delta_i) 
                           - 2.0 * (p1.z_value + p2.z_value - 2.0 * z_star) / r_mu_loc;
        
        double rho = std::pow((1.0 - 1.0 / r_glob_) / (1.0 - 1.0 / r_loc_), 2);

        r_selected_for_next_trial_ = (rho * r_loc_val > r_glob_val) ? r_loc_ : r_glob_;

        return max_r_idx;
    }

    TrialPoint GenerateNewTrialPoint(size_t interval_idx) {
        const TrialPoint& p_left = trial_points_[interval_idx];
        const TrialPoint& p_right = trial_points_[interval_idx + 1];
        
        double x_l = p_left.x_param;
        double x_r = p_right.x_param;
        double z_l = p_left.z_value;
        double z_r = p_right.z_value;
        double new_x = 0.0;

        // Формула (8) из статьи
        if (use_peano_mapping_ && problem_dimension_ > 0) {
            double mu_v = std::max(current_max_mu_estimate_, 1e-9);
            double term_dz_mu = std::abs(z_r - z_l) / mu_v;
            double power_term = std::pow(term_dz_mu, static_cast<double>(problem_dimension_));
            int sign_dz = (z_r - z_l > 0) ? 1 : ((z_r - z_l < 0) ? -1 : 0);

            new_x = 0.5 * (x_l + x_r) - sign_dz * (1.0 / (2.0 * r_selected_for_next_trial_)) * power_term;
        } else {
            double mu_v = std::max(current_max_mu_estimate_, 1e-9);
            double reliable_r_mu = r_selected_for_next_trial_ * mu_v;
            new_x = 0.5 * (x_l + x_r) - (z_r - z_l) / (2.0 * reliable_r_mu);
        }

        double interval_width = x_r - x_l;
        double step_from_boundary = std::max(1e-10, interval_width * 0.0001);

        if (new_x <= x_l + 1e-10) new_x = x_l + step_from_boundary;
        if (new_x >= x_r - 1e-10) new_x = x_r - step_from_boundary;
        if (new_x <= x_l || new_x >= x_r) new_x = x_l + interval_width / 2.0;

        new_x = std::max(search_left_bound_[0], std::min(search_right_bound_[0], new_x));

        return CreatePoint(new_x);
    }

    void InsertPointUnique(const TrialPoint& new_point) {
        for (const auto& p : trial_points_) {
            if (std::abs(p.x_param - new_point.x_param) < 1e-10) return; 
        }
        auto it = std::lower_bound(trial_points_.begin(), trial_points_.end(), new_point);
        trial_points_.insert(it, new_point);
    }

    bool CheckStoppingConditions(const TrialPoint& p_left, const TrialPoint& p_right) {
        double holder_length = GetHolderDelta(p_left.x_param, p_right.x_param);
        if (holder_length <= epsilon_) {
            exit_main_criteria_count_++;
            return true;
        }
        return false;
    }

    std::vector<double> GetBestPointCoords() const {
        if (trial_points_.empty()) return {};
        auto best_it = std::min_element(trial_points_.begin(), trial_points_.end(),
            [](const TrialPoint& a, const TrialPoint& b) {
                if (std::isnan(a.z_value)) return false;
                if (std::isnan(b.z_value)) return true;
                return a.z_value < b.z_value;
            });
        return best_it->y_coords;
    }
};

#endif // DUAL_LIPSCHITZ_OPTIMIZER_HPP