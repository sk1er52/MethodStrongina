// Source.cpp
#include <chrono>
#include <clocale>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept> // Для std::invalid_argument, std::stod, std::stoi
#include <string>
#include <vector>

#include "GrishaginProblemFamily.hpp"
#include "HillProblem.hpp"
#include "HillProblemFamily.hpp"
#include "ShekelProblem.hpp"
#include "ShekelProblemFamily.hpp"
#include "grishagin_function.hpp"
#include "minimizer.h"

// Вспомогательная функция для вывода JSON данных минимизации
void print_minimization_json_to_stdout(
    const std::string &problem_family, int problem_idx_cpp, // 0-based
    double found_z, const std::vector<double> &found_y, int iterations,
    double time_ms, int exit_main_criteria_count, int exit_test_criteria_count,
    double known_opt_z, const std::vector<double> &known_opt_y,
    const std::vector<Point> &trial_history) {

  char *original_locale_numeric = std::setlocale(LC_NUMERIC, nullptr);
  std::string original_locale_numeric_str =
      (original_locale_numeric) ? original_locale_numeric : "";
  std::setlocale(LC_NUMERIC, "C");

  std::cout << std::fixed << std::setprecision(15);
  std::cout << "{" << std::endl;
  std::cout << "  \"type\": \"minimization_result\","
            << std::endl; // Добавляем тип ответа
  std::cout << "  \"problem_family\": \"" << problem_family << "\","
            << std::endl;
  std::cout << "  \"problem_index_cpp\": " << problem_idx_cpp << ","
            << std::endl;
  std::cout << "  \"found_minimum_z\": "
            << (std::isnan(found_z) ? "null" : std::to_string(found_z)) << ","
            << std::endl;
  std::cout << "  \"found_minimum_y\": [";
  for (size_t i = 0; i < found_y.size(); ++i) {
    std::cout << (std::isnan(found_y[i]) ? "null" : std::to_string(found_y[i]))
              << (i == found_y.size() - 1 ? "" : ", ");
  }
  std::cout << "]," << std::endl;
  std::cout << "  \"iterations\": " << iterations << "," << std::endl;
  std::cout << "  \"time_ms\": " << time_ms << "," << std::endl;
  std::cout << "  \"exit_main_count\": " << exit_main_criteria_count << ","
            << std::endl;
  std::cout << "  \"exit_test_count\": " << exit_test_criteria_count << ","
            << std::endl;
  std::cout << "  \"known_optimum_z\": "
            << (std::isnan(known_opt_z) ? "null" : std::to_string(known_opt_z))
            << "," << std::endl;
  std::cout << "  \"known_optimum_y\": [";
  for (size_t i = 0; i < known_opt_y.size(); ++i) {
    std::cout << (std::isnan(known_opt_y[i]) ? "null"
                                             : std::to_string(known_opt_y[i]))
              << (i == known_opt_y.size() - 1 ? "" : ", ");
  }
  std::cout << "]," << std::endl;

  std::cout << "  \"trial_history_x_param\": [";
  for (size_t i = 0; i < trial_history.size(); ++i) {
    std::cout << trial_history[i].x_param
              << (i == trial_history.size() - 1 ? "" : ", ");
  }
  std::cout << "]," << std::endl;

  std::cout << "  \"trial_history_z_value\": [";
  for (size_t i = 0; i < trial_history.size(); ++i) {
    std::cout << (std::isnan(trial_history[i].z_value)
                      ? "null"
                      : std::to_string(trial_history[i].z_value))
              << (i == trial_history.size() - 1 ? "" : ", ");
  }
  std::cout << "]," << std::endl;

  std::cout << "  \"trial_history_y_coords\": [";
  for (size_t i = 0; i < trial_history.size(); ++i) {
    std::cout << "[";
    for (size_t j = 0; j < trial_history[i].y_coords.size(); ++j) {
      std::cout << (std::isnan(trial_history[i].y_coords[j])
                        ? "null"
                        : std::to_string(trial_history[i].y_coords[j]))
                << (j == trial_history[i].y_coords.size() - 1 ? "" : ", ");
    }
    std::cout << "]" << (i == trial_history.size() - 1 ? "" : ", ");
  }
  std::cout << "]" << std::endl;
  std::cout << "}" << std::endl;

  if (!original_locale_numeric_str.empty()) {
    std::setlocale(LC_NUMERIC, original_locale_numeric_str.c_str());
  } else {
    std::setlocale(LC_NUMERIC, "C");
  }
}

// Новая функция для вывода JSON данных линий уровня
void print_level_lines_json_to_stdout(
    const std::vector<double> &x_grid, const std::vector<double> &y_grid,
    const std::vector<std::vector<double>> &z_grid_values) {

  char *original_locale_numeric = std::setlocale(LC_NUMERIC, nullptr);
  std::string original_locale_numeric_str =
      (original_locale_numeric) ? original_locale_numeric : "";
  std::setlocale(LC_NUMERIC, "C");

  std::cout << std::fixed << std::setprecision(15);
  std::cout << "{" << std::endl;
  std::cout << "  \"type\": \"level_lines_data\"," << std::endl;
  std::cout << "  \"x_coords_grid\": [";
  for (size_t i = 0; i < x_grid.size(); ++i) {
    std::cout << x_grid[i] << (i == x_grid.size() - 1 ? "" : ", ");
  }
  std::cout << "]," << std::endl;

  std::cout << "  \"y_coords_grid\": [";
  for (size_t i = 0; i < y_grid.size(); ++i) {
    std::cout << y_grid[i] << (i == y_grid.size() - 1 ? "" : ", ");
  }
  std::cout << "]," << std::endl;

  std::cout << "  \"z_values_grid\": [";
  for (size_t i = 0; i < z_grid_values.size(); ++i) {
    std::cout << "[";
    for (size_t j = 0; j < z_grid_values[i].size(); ++j) {
      std::cout << (std::isnan(z_grid_values[i][j])
                        ? "null"
                        : std::to_string(z_grid_values[i][j]))
                << (j == z_grid_values[i].size() - 1 ? "" : ", ");
    }
    std::cout << "]" << (i == z_grid_values.size() - 1 ? "" : ", ");
  }
  std::cout << "]" << std::endl;
  std::cout << "}" << std::endl;

  if (!original_locale_numeric_str.empty()) {
    std::setlocale(LC_NUMERIC, original_locale_numeric_str.c_str());
  } else {
    std::setlocale(LC_NUMERIC, "C");
  }
}

int main(int argc, char *argv[]) {
  // setlocale(LC_ALL, "Russian");

  // --- Проверка на запрос данных для линий уровня ---
  if (argc > 1 && std::string(argv[1]) == "--get-level-lines") {
    if (argc < 4) { // Ожидаем: ./exe --get-level-lines <problem_idx_0_based>
                    // <grid_resolution>
      std::cerr << "{\"error\": \"Not enough arguments for --get-level-lines. "
                   "Expected: <problem_idx_0_based> <grid_resolution>\"}"
                << std::endl;
      return 1;
    }
    try {
      int problem_idx_ll = std::stoi(argv[2]);
      int grid_resolution_ll = std::stoi(argv[3]);

      if (grid_resolution_ll <= 1 || grid_resolution_ll > 500) {
        std::cerr << "{\"error\": \"Invalid grid_resolution for level lines. "
                     "Expected 2-500.\"}"
                  << std::endl;
        return 1;
      }

      TGrishaginProblemFamily family_ll;
      if (problem_idx_ll < 0 || problem_idx_ll >= family_ll.GetFamilySize()) {
        std::cerr << "{\"error\": \"Invalid problem_idx for level lines. Index "
                     "out of bounds.\"}"
                  << std::endl;
        return 1;
      }
      TGrishaginProblem *instance_ll =
          dynamic_cast<TGrishaginProblem *>(family_ll[problem_idx_ll]);
      if (!instance_ll) {
        std::cerr << "{\"error\": \"Failed to create Grishagin problem "
                     "instance for level lines.\"}"
                  << std::endl;
        return 1;
      }

      std::vector<double> x_grid(grid_resolution_ll);
      std::vector<double> y_grid(grid_resolution_ll);
      std::vector<std::vector<double>> z_values_grid(
          grid_resolution_ll, std::vector<double>(grid_resolution_ll));

      double x_min = 0.0, x_max = 1.0; // Для Гришагина
      double y_min = 0.0, y_max = 1.0;
      double x_step = (grid_resolution_ll > 1)
                          ? (x_max - x_min) / (grid_resolution_ll - 1)
                          : 0;
      double y_step = (grid_resolution_ll > 1)
                          ? (y_max - y_min) / (grid_resolution_ll - 1)
                          : 0;

      for (int i = 0; i < grid_resolution_ll; ++i) {
        y_grid[i] = y_min + i * y_step;
        for (int j = 0; j < grid_resolution_ll; ++j) {
          if (i == 0)
            x_grid[j] = x_min + j * x_step;
          std::vector<double> current_point_coords = {x_grid[j], y_grid[i]};
          z_values_grid[i][j] =
              instance_ll->ComputeFunction(current_point_coords);
        }
      }
      print_level_lines_json_to_stdout(x_grid, y_grid, z_values_grid);
      return 0;

    } catch (const std::exception &e) {
      std::cerr
          << "{\"error\": \"Error processing --get-level-lines arguments: "
          << e.what() << "\"}" << std::endl;
      return 1;
    }
  }

  // --- Основной режим работы (минимизация) ---
  std::ofstream cpp_internal_log("minimization_log_cpp.txt", std::ios::out);
  if (!cpp_internal_log.is_open()) {
    std::cerr << "{\"error\": \"CRITICAL CPP ERROR: Failed to open C++ "
                 "internal log file 'minimization_log_cpp.txt'!'\"}"
              << std::endl;
  } else {
    cpp_internal_log
        << "--- C++ Internal Log Session Start (Minimization Mode) ---"
        << std::endl;
  }

  if (argc < 7) {
    if (cpp_internal_log.is_open())
      cpp_internal_log << "Error: Not enough arguments for minimization. "
                          "Expected 6 parameters."
                       << std::endl;
    std::cerr
        << "{\"error\": \"Not enough arguments for minimization. Expected: "
           "<taskType> <subChoice/idx> <eps> <r> <idx_1D/peano_m> <maxIter>\"}"
        << std::endl;
    if (cpp_internal_log.is_open())
      cpp_internal_log.close();
    return 1;
  }

  int taskType_choice = 0;
  int subChoice_val = 0;
  int problem_idx_cpp_val = 0;
  double epsilon_param = 0.001;
  double r_param_strongin = 3.0;
  int peano_m_order_val = 10;
  int max_iterations_val = 10000;

  try {
    taskType_choice = std::stoi(argv[1]);
    epsilon_param = std::stod(argv[3]);
    r_param_strongin = std::stod(argv[4]);
    max_iterations_val = std::stoi(argv[6]);

    if (taskType_choice == 1) {
      subChoice_val = std::stoi(argv[2]);
      problem_idx_cpp_val = std::stoi(argv[5]);
    } else if (taskType_choice == 2) {
      problem_idx_cpp_val = std::stoi(argv[2]);
      peano_m_order_val = std::stoi(argv[5]);
    } else {
      throw std::invalid_argument("Invalid taskType_choice for minimization.");
    }
  } catch (const std::exception &e) {
    if (cpp_internal_log.is_open())
      cpp_internal_log << "Error parsing arguments for minimization: "
                       << e.what() << std::endl;
    std::cerr << "{\"error\": \"Error parsing arguments for minimization: "
              << e.what() << "\"}" << std::endl;
    if (cpp_internal_log.is_open())
      cpp_internal_log.close();
    return 1;
  }

  if (cpp_internal_log.is_open()) {
    cpp_internal_log << "Parsed Args: TaskType=" << taskType_choice
                     << ", Epsilon=" << epsilon_param
                     << ", R=" << r_param_strongin
                     << ", MaxIter=" << max_iterations_val;
    if (taskType_choice == 1) {
      cpp_internal_log << ", SubChoice=" << subChoice_val
                       << ", ProblemIdx(0-based)=" << problem_idx_cpp_val;
    } else {
      cpp_internal_log << ", ProblemIdx(0-based)=" << problem_idx_cpp_val
                       << ", PeanoOrder(m)=" << peano_m_order_val;
    }
    cpp_internal_log << std::endl;
  }

  std::ofstream dataFile("plot_data.txt", std::ios::out);
  if (!dataFile.is_open() && cpp_internal_log.is_open()) {
    cpp_internal_log << "Warning: Could not open plot_data.txt for writing."
                     << std::endl;
  }

  std::string problem_family_name_str;
  double found_z = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> found_y;
  int iterations = 0;
  double time_ms = 0;
  int exit_main = 0;
  int exit_test = 0;
  double known_opt_z_val = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> known_opt_y_coords;
  std::vector<Point> history_of_trials;

  if (taskType_choice == 1) {
    if (subChoice_val == 1) {
      problem_family_name_str = "Hill";
      THillProblemFamily family;
      if (problem_idx_cpp_val < 0 ||
          problem_idx_cpp_val >= family.GetFamilySize())
        problem_idx_cpp_val = 0;
      THillProblem *instance =
          dynamic_cast<THillProblem *>(family[problem_idx_cpp_val]);
      if (!instance) {
        if (cpp_internal_log.is_open())
          cpp_internal_log.close();
        return 1;
      }

      known_opt_z_val = instance->GetOptimumValue();
      known_opt_y_coords = instance->GetOptimumPoint();
      std::vector<double> domain_a = {0.0}, domain_b = {1.0};
      Minimizer<THillProblem> minimizer(domain_a, domain_b, epsilon_param,
                                        r_param_strongin, *instance,
                                        cpp_internal_log, max_iterations_val);

      auto start_t = std::chrono::high_resolution_clock::now();
      found_y = minimizer.findMinimum();
      auto end_t = std::chrono::high_resolution_clock::now();
      time_ms =
          std::chrono::duration<double, std::milli>(end_t - start_t).count();
      if (!found_y.empty())
        found_z = instance->ComputeFunction(found_y);
      iterations = minimizer.GetIterationCount();
      exit_main = minimizer.GetExitMainCount();
      exit_test = minimizer.GetExitTestCount();
      history_of_trials = minimizer.GetTrialPointsHistory();

    } else if (subChoice_val == 2) {
      problem_family_name_str = "Shekel";
      TShekelProblemFamily family;
      if (problem_idx_cpp_val < 0 ||
          problem_idx_cpp_val >= family.GetFamilySize())
        problem_idx_cpp_val = 0;
      TShekelProblem *instance =
          dynamic_cast<TShekelProblem *>(family[problem_idx_cpp_val]);
      if (!instance) {
        if (cpp_internal_log.is_open())
          cpp_internal_log.close();
        return 1;
      }

      known_opt_z_val = instance->GetOptimumValue();
      known_opt_y_coords = instance->GetOptimumPoint();
      std::vector<double> domain_a = {0.0}, domain_b = {10.0};
      Minimizer<TShekelProblem> minimizer(domain_a, domain_b, epsilon_param,
                                          r_param_strongin, *instance,
                                          cpp_internal_log, max_iterations_val);

      auto start_t = std::chrono::high_resolution_clock::now();
      found_y = minimizer.findMinimum();
      auto end_t = std::chrono::high_resolution_clock::now();
      time_ms =
          std::chrono::duration<double, std::milli>(end_t - start_t).count();
      if (!found_y.empty())
        found_z = instance->ComputeFunction(found_y);
      iterations = minimizer.GetIterationCount();
      exit_main = minimizer.GetExitMainCount();
      exit_test = minimizer.GetExitTestCount();
      history_of_trials = minimizer.GetTrialPointsHistory();
    } else {
      if (cpp_internal_log.is_open())
        cpp_internal_log << "Error: Invalid subChoice for 1D task: "
                         << subChoice_val << std::endl;
      std::cerr << "{\"error\": \"Invalid subChoice for 1D task.\"}"
                << std::endl;
      if (dataFile.is_open())
        dataFile.close();
      if (cpp_internal_log.is_open())
        cpp_internal_log.close();
      return 1;
    }
  } else if (taskType_choice == 2) {
    problem_family_name_str = "Grishagin";
    TGrishaginProblemFamily family;
    if (problem_idx_cpp_val < 0 ||
        problem_idx_cpp_val >= family.GetFamilySize())
      problem_idx_cpp_val = 0;
    TGrishaginProblem *instance =
        dynamic_cast<TGrishaginProblem *>(family[problem_idx_cpp_val]);
    if (!instance) {
      if (cpp_internal_log.is_open())
        cpp_internal_log.close();
      return 1;
    }

    known_opt_z_val = instance->GetOptimumValue();
    known_opt_y_coords = instance->GetOptimumPoint();
    double p_min = 0.0, p_max = 1.0;
    int p_dim = 2, p_key = 1;
    Minimizer<TGrishaginProblem> minimizer(
        p_min, p_max, epsilon_param, r_param_strongin, *instance,
        peano_m_order_val, p_dim, p_key, cpp_internal_log, max_iterations_val);

    auto start_t = std::chrono::high_resolution_clock::now();
    found_y = minimizer.findMinimum();
    auto end_t = std::chrono::high_resolution_clock::now();
    time_ms =
        std::chrono::duration<double, std::milli>(end_t - start_t).count();
    if (!found_y.empty())
      found_z = instance->ComputeFunction(found_y);
    iterations = minimizer.GetIterationCount();
    exit_main = minimizer.GetExitMainCount();
    exit_test = minimizer.GetExitTestCount();
    history_of_trials = minimizer.GetTrialPointsHistory();
  } else {
    if (cpp_internal_log.is_open())
      cpp_internal_log << "Error: Invalid taskType_choice (final check): "
                       << taskType_choice << std::endl;
    std::cerr << "{\"error\": \"Invalid taskType_choice (final check).\"}"
              << std::endl;
    if (dataFile.is_open())
      dataFile.close();
    if (cpp_internal_log.is_open())
      cpp_internal_log.close();
    return 1;
  }

  if (dataFile.is_open()) {
    if (iterations > 0)
      dataFile << 1 << " " << iterations << std::endl;
    dataFile.close();
  }

  print_minimization_json_to_stdout(
      problem_family_name_str, problem_idx_cpp_val, found_z, found_y,
      iterations, time_ms, exit_main, exit_test, known_opt_z_val,
      known_opt_y_coords, history_of_trials);

  if (cpp_internal_log.is_open()) {
    cpp_internal_log
        << "--- C++ Internal Log Session End (Minimization Mode) ---"
        << std::endl;
    cpp_internal_log.close();
  }
  return 0;
}