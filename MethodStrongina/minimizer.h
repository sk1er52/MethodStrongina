// minimizer.h (ВЕРСИЯ С РАСШИРЕННЫМ ЛОГГИРОВАНИЕМ В ФАЙЛ)
#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

// Структуры и объявления функций остаются без изменений
struct Point {
  double x_param;
  double z_value;
  std::vector<double> y_coords;
};

class FunctionInterface {
public:
  virtual double ComputeFunction(const std::vector<double> &x) const = 0;
  virtual std::vector<double> GetOptimumPoint() const = 0;
  virtual double GetOptimumValue() const = 0;
  virtual ~FunctionInterface() {}
};

extern int n1, nexp, l, iq, iu[10], iv[10];
void mapd(double x, int m, double *y, int n, int key);
void node(int is);
void resetMappingGlobals(int n);
std::vector<double> peanoMapping(double x, int m, int n, int key);

template <typename T_ProblemType> class Minimizer {
private:
  // ... все поля класса ...
  std::vector<double> searchLeftBound_param;
  std::vector<double> searchRightBound_param;
  double epsilon_val;
  double r_strongin_parameter;
  const T_ProblemType &function_instance;
  std::ostream &logger; // Основной лог для UI

  // [ИЗМЕНЕНИЕ 1]: Добавляем отдельный файловый поток для трассировки
  std::ofstream trace_logger;

  int max_iterations;
  int iteration_counter;
  int exit_main_criteria_count;
  int exit_test_criteria_count;
  std::vector<Point> trial_points_list;
  bool use_peano_mapping_flag;
  int peano_mapping_m_order;
  int problem_actual_dimension;
  int peano_mapping_key_type;
  double current_m_phi_estimate;

  // --- Приватные методы ---
  // Большинство методов остаются без изменений, изменения в основном в
  // findMinimum и findIntervalIndexWithMaxR

  void performFirstIteration() {
    // ... этот метод без изменений ...
    logger << "Первая итерация:\n";
    trace_logger << "--- Initial State ---\n";
    Point p_left;
    p_left.x_param = searchLeftBound_param[0];
    if (!use_peano_mapping_flag)
      p_left.y_coords = {p_left.x_param};
    else
      p_left.y_coords =
          peanoMapping(p_left.x_param, peano_mapping_m_order,
                       problem_actual_dimension, peano_mapping_key_type);
    p_left.z_value = function_instance.ComputeFunction(p_left.y_coords);
    trial_points_list.push_back(p_left);
    logTrialPointToFile(p_left, "Начальная точка (левая): ");
    logPointToTraceFile(p_left, "Initial Left");

    Point p_right;
    p_right.x_param = searchRightBound_param[0];
    if (!use_peano_mapping_flag)
      p_right.y_coords = {p_right.x_param};
    else
      p_right.y_coords =
          peanoMapping(p_right.x_param, peano_mapping_m_order,
                       problem_actual_dimension, peano_mapping_key_type);
    p_right.z_value = function_instance.ComputeFunction(p_right.y_coords);
    if (p_right.x_param < p_left.x_param)
      trial_points_list.insert(trial_points_list.begin(), p_right);
    else
      trial_points_list.push_back(p_right);
    logTrialPointToFile(p_right, "Начальная точка (правая): ");
    logPointToTraceFile(p_right, "Initial Right");
  }

  // [ИЗМЕНЕНИЕ 2]: Добавляем логгирование в computeIntervalCharacteristicR
  double computeIntervalCharacteristicR(const Point &p1, const Point &p2,
                                        size_t interval_idx) {
    trace_logger << "  Interval " << interval_idx << ": x_left=" << p1.x_param
                 << ", z_left=" << p1.z_value << " | x_right=" << p2.x_param
                 << ", z_right=" << p2.z_value << "\n";

    double dx_param_raw = p2.x_param - p1.x_param;
    if (std::abs(dx_param_raw) < 1e-12) {
      trace_logger << "    -> Interval too small. R = -inf\n";
      return -std::numeric_limits<double>::infinity();
    }

    double delta_i =
        use_peano_mapping_flag
            ? pow(dx_param_raw,
                  1.0 / static_cast<double>(problem_actual_dimension))
            : dx_param_raw;

    trace_logger << "    dx_raw=" << dx_param_raw
                 << ", delta_i(1/N)=" << delta_i << "\n";

    double dz = p2.z_value - p1.z_value;
    double reliable_m_phi = std::max(current_m_phi_estimate, 1e-9);

    double term1 = reliable_m_phi * delta_i;
    double term2 = (dz * dz) / (reliable_m_phi * delta_i);
    double term3 = -2.0 * (p1.z_value + p2.z_value);
    double r_val = term1 + term2 + term3;

    trace_logger << "    m_phi=" << reliable_m_phi << ", R = " << term1 << " + "
                 << term2 << " + (" << term3 << ") = " << r_val << "\n";

    return r_val;
  }

  size_t findIntervalIndexWithMaxR() {
    if (trial_points_list.size() < 2)
      return 0;
    double max_R_val = -std::numeric_limits<double>::infinity();
    size_t max_R_idx = 0;
    trace_logger << "  Calculating R for all intervals:\n";
    for (size_t i = 0; i < trial_points_list.size() - 1; ++i) {
      double current_R_val = computeIntervalCharacteristicR(
          trial_points_list[i], trial_points_list[i + 1], i + 1);
      if (current_R_val > max_R_val) {
        max_R_val = current_R_val;
        max_R_idx = i;
      }
    }
    logger << "Макс. характеристика R = " << max_R_val
           << " на интервале x_param=[" << trial_points_list[max_R_idx].x_param
           << ", " << trial_points_list[max_R_idx + 1].x_param << "]\n";
    trace_logger << "  -> CHOSEN interval " << max_R_idx + 1
                 << " with max R = " << max_R_val << "\n";
    return max_R_idx;
  }

  Point generateNewTrialPoint(size_t interval_index_left) {
    // ... этот метод без изменений ...
    const Point &p_left = trial_points_list[interval_index_left];
    const Point &p_right = trial_points_list[interval_index_left + 1];
    Point new_p;
    double reliable_m_phi_for_new_point =
        std::max(current_m_phi_estimate, 1e-9);
    new_p.x_param = 0.5 * (p_left.x_param + p_right.x_param) -
                    (p_right.z_value - p_left.z_value) /
                        (2.0 * reliable_m_phi_for_new_point);

    double interval_width = p_right.x_param - p_left.x_param;
    double min_rel_step = 0.001, abs_min_step = 1e-9;
    double step_from_boundary =
        std::max(abs_min_step, interval_width * min_rel_step);
    if (interval_width <= 0)
      new_p.x_param = p_left.x_param;
    else {
      if (new_p.x_param <= p_left.x_param)
        new_p.x_param = p_left.x_param + step_from_boundary;
      if (new_p.x_param >= p_right.x_param)
        new_p.x_param = p_right.x_param - step_from_boundary;
      if (new_p.x_param <= p_left.x_param || new_p.x_param >= p_right.x_param)
        new_p.x_param = p_left.x_param + interval_width / 2.0;
    }

    if (!use_peano_mapping_flag)
      new_p.y_coords = {new_p.x_param};
    else
      new_p.y_coords =
          peanoMapping(new_p.x_param, peano_mapping_m_order,
                       problem_actual_dimension, peano_mapping_key_type);
    new_p.z_value = function_instance.ComputeFunction(new_p.y_coords);
    trace_logger << "  Generated new point:\n";
    logPointToTraceFile(new_p, "    New Point");
    return new_p;
  }

  // Вспомогательная функция для логгирования в trace_logger
  void logPointToTraceFile(const Point &p, const std::string &prefix = "") {
    trace_logger << prefix << ": x=" << p.x_param << ", z=" << p.z_value
                 << ", y=[";
    for (size_t i = 0; i < p.y_coords.size(); ++i)
      trace_logger << p.y_coords[i] << (i == p.y_coords.size() - 1 ? "" : ",");
    trace_logger << "]\n";
  }

  // Остальные методы (checkStoppingConditions, logTrialPointToFile) без
  // изменений из финальной версии
  bool checkStoppingConditions(const Point &, const Point &) { /* ... */
    return false;
  } // Пока заглушка, основная логика в findMinimum
  void logTrialPointToFile(const Point &p,
                           const std::string &prefix = "") { /* ... */
    logger << prefix;
    logger << "x_param: " << p.x_param << ", z_value: " << p.z_value
           << ", y_coords: [";
    for (size_t i = 0; i < p.y_coords.size(); ++i) {
      logger << p.y_coords[i] << (i == p.y_coords.size() - 1 ? "" : ", ");
    }
    logger << "]" << std::endl;
  }

public:
  // [ИЗМЕНЕНИЕ 3]: Конструкторы теперь открывают trace_logger
  Minimizer(const std::vector<double> &problem_domain_a,
            const std::vector<double> &problem_domain_b, double eps,
            double r_val, const T_ProblemType &func, std::ostream &log_stream,
            int max_iter)
      : epsilon_val(eps), r_strongin_parameter(r_val), function_instance(func),
        logger(log_stream), max_iterations(max_iter),
        use_peano_mapping_flag(false), problem_actual_dimension(1) {
    searchLeftBound_param = problem_domain_a;
    searchRightBound_param = problem_domain_b;
    trace_logger.open("minimizer_trace_log.txt",
                      std::ios::out | std::ios::trunc);
    trace_logger << std::fixed << std::setprecision(15);
  }

  Minimizer(double peano_search_param_a, double peano_search_param_b,
            double eps, double r_val, const T_ProblemType &func,
            int mapping_m_order, int original_problem_dimension,
            int mapping_key, std::ostream &log_stream, int max_iter)
      : epsilon_val(eps), r_strongin_parameter(r_val), function_instance(func),
        logger(log_stream), max_iterations(max_iter),
        use_peano_mapping_flag(true), peano_mapping_m_order(mapping_m_order),
        problem_actual_dimension(original_problem_dimension),
        peano_mapping_key_type(mapping_key) {
    searchLeftBound_param = {peano_search_param_a};
    searchRightBound_param = {peano_search_param_b};
    trace_logger.open("minimizer_trace_log.txt",
                      std::ios::out | std::ios::trunc);
    trace_logger << std::fixed << std::setprecision(15);
  }

  ~Minimizer() {
    if (trace_logger.is_open()) {
      trace_logger.close();
    }
  }

  // Основной цикл с логгированием
  std::vector<double> findMinimum() {
    trial_points_list.clear();
    iteration_counter = 0;
    exit_main_criteria_count = 0;
    exit_test_criteria_count = 0;
    performFirstIteration();

    for (int current_iter = 0; current_iter < this->max_iterations;
         ++current_iter) {
      iteration_counter = current_iter + 1;
      trace_logger << "\n--- Iteration " << iteration_counter << " ---\n";

      double max_abs_slope = 0.0;
      for (size_t i = 0; i < trial_points_list.size() - 1; ++i) {
        double dx_param_raw =
            trial_points_list[i + 1].x_param - trial_points_list[i].x_param;
        if (dx_param_raw > 1e-12) {
          double delta_i =
              use_peano_mapping_flag
                  ? pow(dx_param_raw,
                        1.0 / static_cast<double>(problem_actual_dimension))
                  : dx_param_raw;
          if (delta_i > 1e-12) {
            double slope = std::abs((trial_points_list[i + 1].z_value -
                                     trial_points_list[i].z_value) /
                                    delta_i);
            if (slope > max_abs_slope)
              max_abs_slope = slope;
          }
        }
      }
      current_m_phi_estimate =
          (max_abs_slope > 1e-9) ? r_strongin_parameter * max_abs_slope : 1.0;
      trace_logger << "Global state: max_abs_slope=" << max_abs_slope
                   << ", current_m_phi_estimate=" << current_m_phi_estimate
                   << "\n";

      size_t interval_idx = findIntervalIndexWithMaxR();

      // ... (логика останова, как в предыдущем ответе) ...
      Point left_p = trial_points_list[interval_idx];
      Point right_p = trial_points_list[interval_idx + 1];
      double metric_y_check =
          use_peano_mapping_flag
              ? pow(right_p.x_param - left_p.x_param,
                    1.0 / static_cast<double>(problem_actual_dimension))
              : (right_p.x_param - left_p.x_param);
      if (metric_y_check <= epsilon_val) {
        exit_main_criteria_count++;
        trace_logger << "STOP condition (by interval length) MET.\n";
        logger << "УСЛОВИЕ ОСТАНОВА (по длине интервала) СРАБОТАЛО.\n";
        break;
      }

      Point new_point = generateNewTrialPoint(interval_idx);

      bool already_exists = false;
      for (const auto &p : trial_points_list)
        if (std::abs(p.x_param - new_point.x_param) < 1e-10)
          already_exists = true;
      if (!already_exists) {
        auto it = std::lower_bound(
            trial_points_list.begin(), trial_points_list.end(), new_point,
            [](auto &p1, auto &p2) { return p1.x_param < p2.x_param; });
        trial_points_list.insert(it, new_point);
      } else {
        trace_logger << "  New point is a duplicate, skipping insertion.\n";
      }

      // ... (логика останова по Y*, как в предыдущем ответе) ...
      auto best_it = std::min_element(
          trial_points_list.begin(), trial_points_list.end(),
          [](auto &a, auto &b) { return a.z_value < b.z_value; });
      auto opt_y = function_instance.GetOptimumPoint();
      bool all_close = true;
      for (size_t i = 0; i < opt_y.size(); ++i)
        if (std::abs(best_it->y_coords[i] - opt_y[i]) > epsilon_val)
          all_close = false;
      if (all_close) {
        exit_test_criteria_count++;
        trace_logger << "STOP condition (by Y* proximity) MET.\n";
        logger << "УСЛОВИЕ ОСТАНОВА (ЦЕЛЕВОЕ, по Y*) СРАБОТАЛО.\n";
        break;
      }
    }

    auto min_it = std::min_element(
        trial_points_list.begin(), trial_points_list.end(),
        [](const Point &a, const Point &b) { return a.z_value < b.z_value; });
    return min_it->y_coords;
  }

  // Getters без изменений
  int GetIterationCount() const { return iteration_counter; }
  int GetExitMainCount() const { return exit_main_criteria_count; }
  int GetExitTestCount() const { return exit_test_criteria_count; }
  const std::vector<Point> &GetTrialPointsHistory() const {
    return trial_points_list;
  }
};

#endif // MINIMIZER_H