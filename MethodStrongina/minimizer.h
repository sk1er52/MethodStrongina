// minimizer.h
#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <algorithm> // Для std::max_element, std::max, std::min, std::sort, std::lower_bound
#include <cmath>
#include <fstream>
#include <iomanip> // Для std::fixed, std::setprecision в логировании, если нужно
#include <iostream>
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
  virtual double GetOptimumValue() const = 0;
  virtual ~FunctionInterface() {}
};

// Объявления функций кривой Пеано
extern int n1, nexp, l, iq, iu[10], iv[10];
void mapd(double x, int m, double *y, int n, int key);
void node(int is);
void resetMappingGlobals(int n);
std::vector<double> peanoMapping(
    double x, int m, int n,
    int key); // Изменено: убран логгер отсюда, т.к. mapd его не принимает

// Шаблонный класс Minimizer
template <typename T_ProblemType> class Minimizer {
private:
  std::vector<double> searchLeftBound_param;
  std::vector<double> searchRightBound_param;
  double epsilon_val; // Используется для критерия останова по координатам Y и
                      // для Гёльдеровской длины интервала
  double r_strongin_parameter;
  const T_ProblemType &function_instance;
  std::ostream &logger;

  int max_iterations_limit_member;
  int iteration_counter;
  int exit_main_criteria_count;
  int exit_test_criteria_count;

  std::vector<Point> trial_points_list;

  bool use_peano_mapping_flag;
  int peano_mapping_m_order;
  int problem_actual_dimension; // N из теории
  int peano_mapping_key_type;

  double current_r_mu_estimate;   // r_v * mu_v из книги (аналог вашего
                                  // current_m_phi_estimate)
  double current_max_mu_estimate; // M или mu_v из книги (аналог вашего
                                  // max_abs_slope)

  // --- Приватные методы ---

  // Метод для вычисления гёльдеровской длины интервала (delta_i из книги)
  double getHolderDelta(double x_left, double x_right) const {
    double dx_param_abs = std::abs(x_right - x_left);
    if (!use_peano_mapping_flag || problem_actual_dimension <= 0) {
      return dx_param_abs;
    }
    if (problem_actual_dimension == 1 &&
        use_peano_mapping_flag) { // Если Пеано для 1D (редко, но возможно)
      return dx_param_abs;
    }
    if (dx_param_abs < 1e-15)
      return 1e-15; // Избегаем pow(0,...) и очень малых чисел, которые могут
                    // дать 0 Возвращаем малое положительное число, чтобы
                    // избежать деления на ноль в R
    return std::pow(dx_param_abs,
                    1.0 / static_cast<double>(problem_actual_dimension));
  }

  void performFirstIteration() {
    logger << "Первая итерация (Minimizer):\n";

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
    logTrialPointToFile(p_left, "  Начальная точка (левая): ");

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
    logTrialPointToFile(p_right, "  Начальная точка (правая): ");

    if (use_peano_mapping_flag &&
        (searchLeftBound_param[0] <= 0.5 && searchRightBound_param[0] >= 0.5) &&
        problem_actual_dimension > 0) {

      double x_param_for_center = 0.5;
      bool center_is_new = true;
      for (const auto &p : trial_points_list) {
        if (std::abs(p.x_param - x_param_for_center) < 1e-9) {
          center_is_new = false;
          break;
        }
      }

      if (center_is_new) {
        logger << "  Принудительное добавление 'центральной' точки для x_param="
               << x_param_for_center << std::endl;
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
        logTrialPointToFile(
            center_point,
            "  Принудительная 'центральная' точка (x_param=0.5): ");
      }
    }
    logger << "Первая итерация (Minimizer) ЗАВЕРШЕНА. Точек в списке: "
           << trial_points_list.size() << std::endl;
  }

  double computeIntervalCharacteristicR(const Point &p1,
                                        const Point &p2) const {
    double delta_i = getHolderDelta(p1.x_param, p2.x_param);

    if (delta_i < 1e-12) { // Если гёльдеровская длина очень мала
      return -std::numeric_limits<double>::infinity();
    }
    double dz = p2.z_value - p1.z_value;
    // current_r_mu_estimate это r_v * mu_v из книги
    double reliable_r_mu =
        std::max(current_r_mu_estimate,
                 r_strongin_parameter); // Если m_phi=0, используем r (т.к. M=1)
                                        // или просто r_strongin_parameter, если
                                        // mu_v=0, то r_v*mu_v -> r_v*1

    // Стандартная формула Стронгина, но с delta_i
    return reliable_r_mu * delta_i + (dz * dz) / (reliable_r_mu * delta_i) -
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
    logger << "  Макс. характеристика R = " << max_R_val
           << " на интервале x_param=[" << trial_points_list[max_R_idx].x_param
           << ", " << trial_points_list[max_R_idx + 1].x_param << "]\n";
    return max_R_idx;
  }

  Point generateNewTrialPoint(size_t interval_index_left) {
    const Point &p_left = trial_points_list[interval_index_left];
    const Point &p_right = trial_points_list[interval_index_left + 1];
    Point new_p;

    double xL = p_left.x_param;
    double xR = p_right.x_param;
    double zL = p_left.z_value;
    double zR = p_right.z_value;

    if (use_peano_mapping_flag &&
        problem_actual_dimension >
            0) { // Используем формулу из книги для N-мерного случая
      // x^{k+1} = (x_t + x_{t-1})/2 - sign(z_t - z_{t-1}) * (1/(2*r_v)) * [|z_t
      // - z_{t-1}| / mu_v]^N current_max_mu_estimate это mu_v (M)
      // r_strongin_parameter это r_v
      double mu_v_calc =
          std::max(current_max_mu_estimate, 1e-9); // Защита от деления на ноль
      if (current_max_mu_estimate < 1e-9)
        mu_v_calc = 1.0; // Если наклон 0, mu_v=1 по книге

      double term_dz_mu = std::abs(zR - zL) / mu_v_calc;
      double power_term =
          std::pow(term_dz_mu, static_cast<double>(problem_actual_dimension));

      int sign_dz = (zR - zL > 0) ? 1 : ((zR - zL < 0) ? -1 : 0);
      if (std::abs(zR - zL) < 1e-9)
        sign_dz = 0; // Если dz=0, то и поправка 0

      new_p.x_param = 0.5 * (xL + xR) -
                      static_cast<double>(sign_dz) *
                          (1.0 / (2.0 * r_strongin_parameter)) * power_term;
      logger << "    Новая точка (nD по книге): x_L=" << xL << ", x_R=" << xR
             << ", z_L=" << zL << ", z_R=" << zR << ", sign_dz=" << sign_dz
             << ", r_v=" << r_strongin_parameter << ", mu_v=" << mu_v_calc
             << ", |zR-zL|/mu_v=" << term_dz_mu
             << ", N=" << problem_actual_dimension
             << ", power_term=" << power_term
             << ", x_new_calc=" << new_p.x_param << std::endl;

    } else { // 1D случай (или если problem_actual_dimension некорректна)
      double reliable_r_mu_for_new_point =
          std::max(current_r_mu_estimate, r_strongin_parameter);
      new_p.x_param =
          0.5 * (xL + xR) - (zR - zL) / (2.0 * reliable_r_mu_for_new_point);
    }

    // Ограничение новой точки
    double interval_width = xR - xL;
    double min_rel_step = 0.0001;
    double abs_min_step = 1e-10;
    double step_from_boundary =
        std::max(abs_min_step, interval_width * min_rel_step);

    if (interval_width <= 1e-12) {
      new_p.x_param = xL + interval_width / 2.0;
    } else {
      if (new_p.x_param <= xL + abs_min_step) {
        new_p.x_param = xL + step_from_boundary;
      }
      if (new_p.x_param >= xR - abs_min_step) {
        new_p.x_param = xR - step_from_boundary;
      }
      if (new_p.x_param <= xL || new_p.x_param >= xR) {
        new_p.x_param = xL + interval_width / 2.0;
      }
    }
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
    logger << "\n--- CheckStoppingConditions (Итерация: " << iteration_counter
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
    // оптимуму Y* ---
    if (use_peano_mapping_flag &&
        best_point_found_for_check) { // Добавил use_peano_mapping_flag
      std::vector<double> known_optimum_y_coords =
          function_instance.GetOptimumPoint();

      if (!known_optimum_y_coords.empty() &&
          known_optimum_y_coords.size() == problem_actual_dimension &&
          y_current_best_point.y_coords.size() == problem_actual_dimension) {

        bool all_coords_close_enough = true;
        logger << "  Проверка критерия по КООРДИНАТАМ Y:" << std::endl;
        // ... (логирование y_best, y*, epsilon_val) ...
        for (int j = 0; j < problem_actual_dimension; ++j) {
          double coord_diff = std::abs(y_current_best_point.y_coords[j] -
                                       known_optimum_y_coords[j]);
          if (coord_diff > epsilon_val) {
            all_coords_close_enough = false;
            break;
          }
        }
        if (all_coords_close_enough) {
          logger << "  Условие останова (доп. nD по ВСЕМ КООРД. РАЗНИЦАМ Y) "
                    "СРАБОТАЛО.\n";
          exit_test_criteria_count++;
          return true;
        } else {
          logger << "  Условие останова (доп. nD по ВСЕМ КООРД. РАЗНИЦАМ Y) НЕ "
                    "СРАБОТАЛО."
                 << std::endl;
        }
      } else { /* ... */
      }
    }

    // --- Критерий 2: Основной по ГЁЛЬДЕРОВСКОЙ длине интервала x_param ---
    // (Xt - Xt-1)^(1/N) <= epsilon (из книги стр. 219)
    double holder_interval_length =
        getHolderDelta(p_left_of_interval.x_param, p_right_of_interval.x_param);

    // epsilon_val из UI используется как "покоординатная точность решения
    // задачи" (epsilon из книги)
    double epsilon_for_holder_length_stop = epsilon_val;

    logger << "  Основной критерий: гёльдеровская длина интервала x_param "
              "(dx^(1/N))="
           << holder_interval_length
           << ", epsilon_for_stop=" << epsilon_for_holder_length_stop
           << std::endl;

    if (holder_interval_length <= epsilon_for_holder_length_stop) {
      logger << "  Условие останова (ОСНОВНОЕ, по Гёльдеровской длине X_PARAM) "
                "СРАБОТАЛО.\n";
      exit_main_criteria_count++;
      return true;
    }

    logger << "--- Конец CheckStoppingConditions (никакое условие не "
              "сработало) ---"
           << std::endl;
    return false;
  }

  void logTrialPointToFile(const Point &p, const std::string &prefix = "") {
    logger << prefix;
    logger << std::fixed
           << std::setprecision(15); // Устанавливаем точность для логгера
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
        problem_actual_dimension(1), current_r_mu_estimate(r_val),
        current_max_mu_estimate(1.0) // Начальные значения
  {
    if (problem_domain_a.empty() || problem_domain_b.empty()) {
      logger
          << "ОШИБКА в конструкторе 1D: Границы задачи не могут быть пустыми."
          << std::endl;
      throw std::invalid_argument("Границы задачи (1D) не могут быть пустыми.");
    }
    searchLeftBound_param = problem_domain_a;
    searchRightBound_param = problem_domain_b;
    logger << ">>>> 1D КОНСТРУКТОР MINIMIZER ВЫЗВАН (лог через ссылку ostream) "
              "<<<<"
           << std::endl;
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
        peano_mapping_key_type(mapping_key), current_r_mu_estimate(r_val),
        current_max_mu_estimate(1.0) // Начальные значения
  {
    searchLeftBound_param = {peano_search_param_a};
    searchRightBound_param = {peano_search_param_b};
    logger << ">>>> nD КОНСТРУКТОР MINIMIZER ВЫЗВАН (с Пеано, лог через ссылку "
              "ostream) <<<<"
           << std::endl;
  }

  ~Minimizer() {
    logger << ">>>> MINIMIZER ДЕСТРУКТОР ВЫЗВАН <<<<" << std::endl;
  }

  std::vector<double> findMinimum() {
    logger << ">>>> findMinimum() НАЧАЛО <<<<" << std::endl;
    trial_points_list.clear();
    iteration_counter = 0;
    exit_main_criteria_count = 0;
    exit_test_criteria_count = 0;

    logger << "  Вызов performFirstIteration()..." << std::endl;
    performFirstIteration();
    logger << "  performFirstIteration() ЗАВЕРШЕН. Точек в списке: "
           << trial_points_list.size() << std::endl;

    if (trial_points_list.size() < 2) {
      logger << "Ошибка: Первая итерация не создала достаточно точек для "
                "начала основного цикла."
             << std::endl;
      if (!trial_points_list.empty())
        return trial_points_list[0].y_coords;
      return {};
    }

    for (int current_iter_loop = 0;
         current_iter_loop < this->max_iterations_limit_member;
         ++current_iter_loop) {
      iteration_counter = current_iter_loop + 1;

      if (trial_points_list.size() < 2) {
        logger << "  В trial_points_list меньше 2 точек, невозможно "
                  "продолжить. Итерация: "
               << iteration_counter << std::endl;
        break;
      }

      // 1. Обновить оценку current_max_mu_estimate (M или mu_v) и
      // current_r_mu_estimate (r*M)
      current_max_mu_estimate = 0.0;
      size_t idx_max_slope_p1 = 0, idx_max_slope_p2 = 0;

      for (size_t i = 0; i < trial_points_list.size() - 1; ++i) {
        // dx_param_raw используется для getHolderDelta, которое само возьмет
        // модуль
        double dx_param_raw =
            trial_points_list[i + 1].x_param - trial_points_list[i].x_param;
        double holder_delta_for_slope = getHolderDelta(
            trial_points_list[i].x_param, trial_points_list[i + 1].x_param);

        if (holder_delta_for_slope >
            1e-12) { // Знаменатель не должен быть слишком мал
          double current_z_diff = std::abs(trial_points_list[i + 1].z_value -
                                           trial_points_list[i].z_value);
          double mu_val = current_z_diff / holder_delta_for_slope;
          if (mu_val > current_max_mu_estimate) {
            current_max_mu_estimate = mu_val;
            idx_max_slope_p1 = i;
            idx_max_slope_p2 = i + 1;
          }
        }
      }
      // Если наклон 0 (или очень мал), mu_v=1 по книге. r_v*mu_v = r_v.
      if (current_max_mu_estimate < 1e-9) {
        current_max_mu_estimate = 1.0; // Это mu_v
      }
      current_r_mu_estimate =
          r_strongin_parameter * current_max_mu_estimate; // Это r_v * mu_v

      logger << "  Итерация " << iteration_counter
             << ": current_max_mu_estimate (M) = " << current_max_mu_estimate
             << " (на интервале x=["
             << trial_points_list[idx_max_slope_p1].x_param << ", "
             << trial_points_list[idx_max_slope_p2].x_param
             << "]), current_r_mu_estimate (r*M) = " << current_r_mu_estimate
             << std::endl;

      size_t interval_to_split_idx = findIntervalIndexWithMaxR();

      if (interval_to_split_idx + 1 >= trial_points_list.size()) {
        logger << "ОШИБКА: interval_to_split_idx (" << interval_to_split_idx
               << ") указывает за пределы массива точек (размер "
               << trial_points_list.size()
               << "). Выход. Итерация: " << iteration_counter << std::endl;
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
        if (std::abs(existing_p.x_param - new_point.x_param) < 1e-10) {
          already_exists_nearby = true;
          break;
        }
      }
      if (already_exists_nearby) {
        logger
            << "    Новая точка x_param=" << new_point.x_param
            << " (z=" << new_point.z_value
            << ") очень близка к существующей. Пропуск добавления. Итерация: "
            << iteration_counter << std::endl;
      } else {
        auto it_insert =
            std::lower_bound(trial_points_list.begin(), trial_points_list.end(),
                             new_point, [](const Point &p1, const Point &p2) {
                               return p1.x_param < p2.x_param;
                             });
        trial_points_list.insert(it_insert, new_point);
        logTrialPointToFile(new_point, "  Новая точка добавлена: ");
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

    logger << ">>>> findMinimum() КОНЕЦ <<<<" << std::endl;
    if (valid_points_exist) {
      logger << "Минимум найден в точке (x_param="
             << trial_points_list[min_z_idx].x_param
             << ", z=" << trial_points_list[min_z_idx].z_value << ")"
             << std::endl;
      logTrialPointToFile(trial_points_list[min_z_idx], "Лучшая точка: ");
      return trial_points_list[min_z_idx].y_coords;
    }

    logger
        << "Список пробных точек пуст или все точки NaN, оптимум не найден.\n";
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