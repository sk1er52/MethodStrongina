// minimizer.h
#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <algorithm> // Для std::max_element, std::max, std::min, std::sort, std::lower_bound
#include <cmath>
#include <fstream> // Для std::ofstream в main.cpp, но не в Minimizer
#include <iomanip> // Для std::fixed, std::setprecision в логировании, если нужно
#include <iostream> // Для std::cerr и std::ostream
#include <limits>
#include <numeric> // Для std::distance (не используется сейчас, но может пригодиться)
#include <stdexcept> // Для std::invalid_argument
#include <string>
#include <vector>

// Используем std:: для краткости в этом файле
using namespace std;

// Структура для представления точки
struct Point {
  double x_param; // Одномерная координата для поиска (параметр t для кривой
                  // Пеано в nD)
  double z_value; // Значение функции в точке: f(y_coords)
  std::vector<double>
      y_coords; // Координаты точки в исходном пространстве (n-мерном)
};

// Интерфейс для целевой функции
class FunctionInterface {
public:
  virtual double ComputeFunction(const std::vector<double> &x) const = 0;
  virtual std::vector<double> GetOptimumPoint() const = 0;
  virtual double GetOptimumValue() const = 0; // Убедитесь, что этот метод есть
  virtual ~FunctionInterface() {}
};

//
// Объявления функций кривой Пеано
//
extern int n1, nexp, l, iq, iu[10], iv[10];

void mapd(double x, int m, double *y, int n, int key);
void node(int is);
void resetMappingGlobals(int n); // n - размерность
std::vector<double> peanoMapping(double x, int m, int n, int key);

//
// Шаблонный класс Minimizer
//
template <typename T_ProblemType> class Minimizer {
private:
  std::vector<double> searchLeftBound_param;
  std::vector<double> searchRightBound_param;
  double epsilon_val; // Используется для критерия останова по координатам Y
  double r_strongin_parameter;
  const T_ProblemType &function_instance;
  std::ostream &logger;

  int max_iterations_limit_member; // Максимальное количество итераций
  int iteration_counter;
  int exit_main_criteria_count; // Счетчик для останова по x_param интервалу
  int exit_test_criteria_count; // Счетчик для останова по Y координатам

  std::vector<Point> trial_points_list;

  bool use_peano_mapping_flag;
  int peano_mapping_m_order;
  int problem_actual_dimension;
  int peano_mapping_key_type;

  double current_m_phi_estimate;

  // --- Приватные методы ---

  void performFirstIteration() {
    logger << "First iteration (Minimizer):\n";

    Point p_left;
    p_left.x_param = searchLeftBound_param[0];
    if (!use_peano_mapping_flag) {
      p_left.y_coords = {p_left.x_param};
    } else {
      p_left.y_coords =
          peanoMapping(p_left.x_param, peano_mapping_m_order,
                       problem_actual_dimension, peano_mapping_key_type);
    }
    p_left.z_value = function_instance.ComputeFunction(p_left.y_coords);
    trial_points_list.push_back(p_left);
    logTrialPointToFile(p_left, "  First trial(left): ");

    Point p_right;
    p_right.x_param = searchRightBound_param[0];
    if (!use_peano_mapping_flag) {
      p_right.y_coords = {p_right.x_param};
    } else {
      p_right.y_coords =
          peanoMapping(p_right.x_param, peano_mapping_m_order,
                       problem_actual_dimension, peano_mapping_key_type);
    }
    p_right.z_value = function_instance.ComputeFunction(p_right.y_coords);
    trial_points_list.push_back(p_right);
    logTrialPointToFile(p_right, "  First trial (right): ");

    if (use_peano_mapping_flag &&
        (searchLeftBound_param[0] <= 0.5 &&
         searchRightBound_param[0] >=
             0.5) &&                    // Убедимся, что 0.5 внутри диапазона
        problem_actual_dimension > 0) { // И что есть размерность

      double x_param_for_center = 0.5;
      bool center_is_new = true;
      for (const auto &p : trial_points_list) {
        if (std::abs(p.x_param - x_param_for_center) < 1e-9) {
          center_is_new = false;
          break;
        }
      }

      if (center_is_new) {
        logger << "  Center trial" << x_param_for_center << std::endl;
        Point center_point;
        center_point.x_param = x_param_for_center;
        center_point.y_coords =
            peanoMapping(center_point.x_param, peano_mapping_m_order,
                         problem_actual_dimension, peano_mapping_key_type);
        center_point.z_value =
            function_instance.ComputeFunction(center_point.y_coords);

        trial_points_list.push_back(center_point);
        std::sort(trial_points_list.begin(), trial_points_list.end(),
                  [](const Point &a, const Point &b) {
                    return a.x_param < b.x_param;
                  });
        logTrialPointToFile(center_point, "  ПCenter trial (x_param=0.5): ");
      }
    }
    logger << "First iteration (Minimizer) ENDED. Trials in list: "
           << trial_points_list.size() << std::endl;
  }

  double computeIntervalCharacteristicR(const Point &p1,
                                        const Point &p2) const {
    double dx_param = p2.x_param - p1.x_param;
    if (std::abs(dx_param) < 1e-12) {
      return -std::numeric_limits<double>::infinity();
    }
    double dz = p2.z_value - p1.z_value;
    double reliable_m_phi = std::max(current_m_phi_estimate, 1e-9);
    return reliable_m_phi * dx_param + (dz * dz) / (reliable_m_phi * dx_param) -
           2.0 * (p1.z_value + p2.z_value);
  }

  size_t findIntervalIndexWithMaxR() {
    if (trial_points_list.size() < 2)
      return 0;
    double max_R_val = -std::numeric_limits<double>::infinity();
    size_t max_R_idx = 0;
    for (size_t i = 0; i < trial_points_list.size() - 1; ++i) {
      double current_R_val = computeIntervalCharacteristicR(
          trial_points_list[i], trial_points_list[i + 1]);
      if (current_R_val > max_R_val) {
        max_R_val = current_R_val;
        max_R_idx = i;
      }
    }
    logger << " R = " << max_R_val << " in interval x_param=["
           << trial_points_list[max_R_idx].x_param << ", "
           << trial_points_list[max_R_idx + 1].x_param << "]\n";
    return max_R_idx;
  }

  Point generateNewTrialPoint(size_t interval_index_left) {
    const Point &p_left = trial_points_list[interval_index_left];
    const Point &p_right = trial_points_list[interval_index_left + 1];
    Point new_p;
    double reliable_m_phi_for_new_point =
        std::max(current_m_phi_estimate, 1e-9);
    new_p.x_param = 0.5 * (p_left.x_param + p_right.x_param) -
                    (p_right.z_value - p_left.z_value) /
                        (2.0 * reliable_m_phi_for_new_point);

    double interval_width = p_right.x_param - p_left.x_param;
    double min_rel_step = 0.001;
    double abs_min_step = 1e-9;
    double step_from_boundary =
        std::max(abs_min_step, interval_width * min_rel_step);

    if (interval_width <= 1e-12) { // Если интервал очень мал или нулевой
      new_p.x_param =
          p_left.x_param +
          interval_width / 2.0; // Середина (может быть = p_left.x_param)
    } else {
      if (new_p.x_param <=
          p_left.x_param + abs_min_step) { // Сдвигаем от левой границы
        new_p.x_param = p_left.x_param + step_from_boundary;
      }
      if (new_p.x_param >=
          p_right.x_param - abs_min_step) { // Сдвигаем от правой границы
        new_p.x_param = p_right.x_param - step_from_boundary;
      }
      // Если после коррекции точка "перепрыгнула" или осталась на границе,
      // ставим в середину
      if (new_p.x_param <= p_left.x_param || new_p.x_param >= p_right.x_param) {
        new_p.x_param = p_left.x_param + interval_width / 2.0;
      }
    }
    // Гарантируем, что точка не выходит за глобальные границы поиска по x_param
    new_p.x_param =
        std::max(searchLeftBound_param[0],
                 std::min(searchRightBound_param[0], new_p.x_param));

    if (!use_peano_mapping_flag) {
      new_p.y_coords = {new_p.x_param};
    } else {
      new_p.y_coords =
          peanoMapping(new_p.x_param, peano_mapping_m_order,
                       problem_actual_dimension, peano_mapping_key_type);
    }
    new_p.z_value = function_instance.ComputeFunction(new_p.y_coords);
    return new_p;
  }

  bool checkStoppingConditions(const Point &p_left_of_interval,
                               const Point & /* new_p_ref_not_used */,
                               const Point &p_right_of_interval) {
    logger << "\n--- CheckStoppingConditions (interation: " << iteration_counter
           << ") ---" << std::endl;

    Point y_current_best_point;
    bool best_point_found_for_check = false;

    if (!trial_points_list.empty()) {
      size_t best_idx = 0;
      double min_z_val = std::numeric_limits<double>::infinity();
      for (size_t i = 0; i < trial_points_list.size(); ++i) {
        if (std::isnan(trial_points_list[i].z_value))
          continue;
        if (!best_point_found_for_check ||
            trial_points_list[i].z_value < min_z_val) {
          min_z_val = trial_points_list[i].z_value;
          best_idx = i;
          best_point_found_for_check = true;
        }
      }
      if (best_point_found_for_check) {
        y_current_best_point = trial_points_list[best_idx];
      }
    }

    // --- Критерий 1: Проверка близости КООРДИНАТ Y лучшей точки к известному
    // оптимуму Y* (для nD) ---
    if (use_peano_mapping_flag && best_point_found_for_check) {
      std::vector<double> known_optimum_y_coords =
          function_instance.GetOptimumPoint();

      if (!known_optimum_y_coords.empty() &&
          known_optimum_y_coords.size() == problem_actual_dimension &&
          y_current_best_point.y_coords.size() == problem_actual_dimension) {

        bool all_coords_close_enough = true;
        logger << "  Check coordinates Y:" << std::endl;
        logger << "    Best find Y: [";
        for (size_t j = 0; j < y_current_best_point.y_coords.size(); ++j)
          logger << y_current_best_point.y_coords[j]
                 << (j == y_current_best_point.y_coords.size() - 1 ? "" : ", ");
        logger << "], Z_best=" << y_current_best_point.z_value << std::endl;
        logger << "    Known optimum Y*: [";
        for (size_t j = 0; j < known_optimum_y_coords.size(); ++j)
          logger << known_optimum_y_coords[j]
                 << (j == known_optimum_y_coords.size() - 1 ? "" : ", ");
        logger << "], Z*=" << function_instance.GetOptimumValue() << std::endl;
        logger << "    Epsilon Y: " << epsilon_val << std::endl;

        for (int j = 0; j < problem_actual_dimension; ++j) {
          double coord_diff = std::abs(y_current_best_point.y_coords[j] -
                                       known_optimum_y_coords[j]);
          logger << "      Coordinate " << j << ": |"
                 << y_current_best_point.y_coords[j] << " - "
                 << known_optimum_y_coords[j] << "| = " << coord_diff
                 << std::endl;
          if (coord_diff >
              epsilon_val) { // epsilon_val используется для координат Y
            all_coords_close_enough = false;
            // break; // Можно выйти раньше, если важна только общая оценка
          }
        }

        if (all_coords_close_enough) {
          logger << "  First condition worked\n";
          exit_test_criteria_count++;
          return true;
        } else {
          logger << "  First condition dont worked" << std::endl;
        }
      } else {
      }
    }

    // // --- Критерий 2: Основной по длине интервала x_param ---
    // double chosen_interval_length =
    //     p_right_of_interval.x_param - p_left_of_interval.x_param;
    // double epsilon_for_x_param = 1e-6; // Фиксированное очень малое значение.

    // logger << "  Основной критерий: длина интервала x_param = "
    //        << chosen_interval_length
    //        << ", epsilon_for_x_param = " << epsilon_for_x_param << std::endl;

    // if (chosen_interval_length <= epsilon_for_x_param) {
    //   logger << "  Условие останова (ОСНОВНОЕ, по X_PARAM) СРАБОТАЛО.\n";
    //   exit_main_criteria_count++;
    //   return true;
    // }

    logger << "--- End CheckStoppingConditions ---" << std::endl;
    return false;
  }

  void logTrialPointToFile(const Point &p, const std::string &prefix = "") {
    logger << prefix;
    logger << "x_param: " << p.x_param << ", z_value: " << p.z_value
           << ", y_coords: [";
    for (size_t i = 0; i < p.y_coords.size(); ++i) {
      logger << p.y_coords[i] << (i == p.y_coords.size() - 1 ? "" : ", ");
    }
    logger << "]" << std::endl;
  }

public:
  Minimizer(const std::vector<double> &problem_domain_a,
            const std::vector<double> &problem_domain_b, double eps,
            double r_val, const T_ProblemType &func, std::ostream &log_stream,
            int max_iter)
      : epsilon_val(eps), r_strongin_parameter(r_val), function_instance(func),
        logger(log_stream), max_iterations_limit_member(max_iter),
        iteration_counter(0), exit_main_criteria_count(0),
        exit_test_criteria_count(0), use_peano_mapping_flag(false),
        problem_actual_dimension(1), current_m_phi_estimate(1.0) {
    if (problem_domain_a.empty() || problem_domain_b.empty()) {
      // logger
      //     << "ОШИБКА в конструкторе 1D: Границы задачи не могут быть
      //     пустыми."
      //     << std::endl;
      // throw std::invalid_argument("Границы задачи (1D) не могут быть
      // пустыми.");
    }
    searchLeftBound_param = problem_domain_a;
    searchRightBound_param = problem_domain_b;
    // logger << ">>>> 1D КОНСТРУКТОР MINIMIZER ВЫЗВАН (лог через ссылку
    // ostream) "
    //           "<<<<"
    //        << std::endl;
  }

  Minimizer(double peano_search_param_a, double peano_search_param_b,
            double eps, double r_val, const T_ProblemType &func,
            int mapping_m_order, int original_problem_dimension,
            int mapping_key, std::ostream &log_stream, int max_iter)
      : epsilon_val(eps), r_strongin_parameter(r_val), function_instance(func),
        logger(log_stream), max_iterations_limit_member(max_iter),
        iteration_counter(0), exit_main_criteria_count(0),
        exit_test_criteria_count(0), use_peano_mapping_flag(true),
        peano_mapping_m_order(mapping_m_order),
        problem_actual_dimension(original_problem_dimension),
        peano_mapping_key_type(mapping_key), current_m_phi_estimate(1.0) {
    searchLeftBound_param = {peano_search_param_a};
    searchRightBound_param = {peano_search_param_b};
    // logger << ">>>> nD КОНСТРУКТОР MINIMIZER ВЫЗВАН (с Пеано, лог через
    // ссылку "
    //           "ostream) <<<<"
    //        << std::endl;
  }

  ~Minimizer() {
    // logger << ">>>> MINIMIZER ДЕСТРУКТОР ВЫЗВАН <<<<" << std::endl;
  }

  std::vector<double> findMinimum() {
    logger << ">>>> findMinimum() START <<<<" << std::endl;
    trial_points_list.clear();
    iteration_counter = 0;
    exit_main_criteria_count = 0;
    exit_test_criteria_count = 0;

    logger << "  call performFirstIteration()..." << std::endl;
    performFirstIteration();
    logger << "  performFirstIteration() ended. Trials in list: "
           << trial_points_list.size() << std::endl;

    if (trial_points_list.size() < 2) {
      // logger << "Ошибка: Первая итерация не создала достаточно точек для "
      //           "начала основного цикла."
      //        << std::endl;
      if (!trial_points_list.empty())
        return trial_points_list[0]
            .y_coords; // Возвращаем хоть что-то если есть
      return {};       // Пустой вектор если совсем ничего
    }

    for (int current_iter_loop = 0;
         current_iter_loop < this->max_iterations_limit_member;
         ++current_iter_loop) {
      iteration_counter = current_iter_loop + 1;

      if (trial_points_list.size() < 2) {
        // logger << "  В trial_points_list меньше 2 точек, невозможно "
        //           "продолжить. Итерация: "
        //        << iteration_counter << std::endl;
        break;
      }

      double max_abs_slope = 0.0;
      for (size_t i = 0; i < trial_points_list.size() - 1; ++i) {
        double dx_param =
            trial_points_list[i + 1].x_param - trial_points_list[i].x_param;
        if (dx_param > 1e-12) {
          double slope = std::abs((trial_points_list[i + 1].z_value -
                                   trial_points_list[i].z_value) /
                                  dx_param);
          if (slope > max_abs_slope)
            max_abs_slope = slope;
        }
      }
      current_m_phi_estimate =
          (max_abs_slope > 1e-9) ? r_strongin_parameter * max_abs_slope : 1.0;

      size_t interval_to_split_idx = findIntervalIndexWithMaxR();

      if (interval_to_split_idx + 1 >= trial_points_list.size()) {
        // logger << "ОШИБКА: interval_to_split_idx (" << interval_to_split_idx
        //        << ") указывает за пределы массива точек (размер "
        //        << trial_points_list.size()
        //        << "). Выход. Итерация: " << iteration_counter << std::endl;
        break;
      }

      if (checkStoppingConditions(
              trial_points_list[interval_to_split_idx],
              trial_points_list[interval_to_split_idx],
              trial_points_list[interval_to_split_idx + 1])) {
        break;
      }

      Point new_point = generateNewTrialPoint(interval_to_split_idx);

      bool already_exists_nearby = false;
      for (const auto &existing_p : trial_points_list) {
        if (std::abs(existing_p.x_param - new_point.x_param) <
            1e-10) { // Очень жесткий порог для "той же точки"
          already_exists_nearby = true;
          break;
        }
      }

      if (already_exists_nearby) {
        logger << "    New trial x_param=" << new_point.x_param
               << " (z=" << new_point.z_value
               << ") more close. Iteration: " << iteration_counter << std::endl;
        // Если мы постоянно генерируем одну и ту же точку в очень маленьком
        // интервале, основной критерий останова по длине интервала должен это
        // поймать.
      } else {
        auto it_insert =
            std::lower_bound(trial_points_list.begin(), trial_points_list.end(),
                             new_point, [](const Point &p1, const Point &p2) {
                               return p1.x_param < p2.x_param;
                             });
        trial_points_list.insert(it_insert, new_point);
        logTrialPointToFile(new_point, "  New trial added: ");
      }
    }

    size_t min_z_idx = 0;
    bool valid_points_exist = false;
    if (!trial_points_list.empty()) {
      double current_min_z = std::numeric_limits<double>::infinity();
      for (size_t i = 0; i < trial_points_list.size(); ++i) {
        if (std::isnan(trial_points_list[i].z_value))
          continue;
        if (!valid_points_exist ||
            trial_points_list[i].z_value < current_min_z) {
          current_min_z = trial_points_list[i].z_value;
          min_z_idx = i;
          valid_points_exist = true;
        }
      }
    }

    logger << ">>>> findMinimum() END <<<<" << std::endl;
    if (valid_points_exist) {
      logger << "Min finded in trial (x_param="
             << trial_points_list[min_z_idx].x_param
             << ", z=" << trial_points_list[min_z_idx].z_value << ")"
             << std::endl;
      logTrialPointToFile(trial_points_list[min_z_idx], "Best trial: ");
      return trial_points_list[min_z_idx].y_coords;
    }

    // logger
    //     << "Список пробных точек пуст или все точки NaN, оптимум не
    //     найден.\n";
    return {};
  }

  int GetIterationCount() const { return iteration_counter; }
  int GetExitMainCount() const { return exit_main_criteria_count; }
  int GetExitTestCount() const { return exit_test_criteria_count; }
  const std::vector<Point> &GetTrialPointsHistory() const {
    return trial_points_list;
  }
};

#endif // MINIMIZER_H