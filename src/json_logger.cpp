#include "json_logger.hpp"

#include <iostream>
#include <iomanip>
#include <clocale>
#include <cmath>

namespace JsonLogger {

void PrintMinimizationResult(
    const std::string& problem_family, int problem_idx_cpp, 
    double found_z, const std::vector<double>& found_y, int iterations,
    double time_ms, int exit_main_criteria_count, int exit_test_criteria_count,
    double known_opt_z, const std::vector<double>& known_opt_y,
    const std::vector<TrialPoint>& trial_history) {

    char* original_locale_numeric = std::setlocale(LC_NUMERIC, nullptr);
    std::string original_locale_numeric_str = (original_locale_numeric) ? original_locale_numeric : "";
    std::setlocale(LC_NUMERIC, "C");

    std::cout << std::fixed << std::setprecision(15);
    std::cout << "{" << std::endl;
    std::cout << "  \"type\": \"minimization_result\"," << std::endl;
    std::cout << "  \"problem_family\": \"" << problem_family << "\"," << std::endl;
    std::cout << "  \"problem_index_cpp\": " << problem_idx_cpp << "," << std::endl;
    std::cout << "  \"found_minimum_z\": " << (std::isnan(found_z) ? "null" : std::to_string(found_z)) << "," << std::endl;
    
    std::cout << "  \"found_minimum_y\": [";
    for (size_t i = 0; i < found_y.size(); ++i) {
        std::cout << (std::isnan(found_y[i]) ? "null" : std::to_string(found_y[i])) << (i == found_y.size() - 1 ? "" : ", ");
    }
    std::cout << "]," << std::endl;
    
    std::cout << "  \"iterations\": " << iterations << "," << std::endl;
    std::cout << "  \"time_ms\": " << time_ms << "," << std::endl;
    std::cout << "  \"exit_main_count\": " << exit_main_criteria_count << "," << std::endl;
    std::cout << "  \"exit_test_count\": " << exit_test_criteria_count << "," << std::endl;
    
    std::cout << "  \"known_optimum_z\": " << (std::isnan(known_opt_z) ? "null" : std::to_string(known_opt_z)) << "," << std::endl;
    std::cout << "  \"known_optimum_y\": [";
    for (size_t i = 0; i < known_opt_y.size(); ++i) {
        std::cout << (std::isnan(known_opt_y[i]) ? "null" : std::to_string(known_opt_y[i])) << (i == known_opt_y.size() - 1 ? "" : ", ");
    }
    std::cout << "]," << std::endl;

    std::cout << "  \"trial_history_x_param\": [";
    for (size_t i = 0; i < trial_history.size(); ++i) {
        std::cout << trial_history[i].x_param << (i == trial_history.size() - 1 ? "" : ", ");
    }
    std::cout << "]," << std::endl;

    std::cout << "  \"trial_history_z_value\": [";
    for (size_t i = 0; i < trial_history.size(); ++i) {
        std::cout << (std::isnan(trial_history[i].z_value) ? "null" : std::to_string(trial_history[i].z_value)) << (i == trial_history.size() - 1 ? "" : ", ");
    }
    std::cout << "]," << std::endl;

    std::cout << "  \"trial_history_y_coords\": [";
    for (size_t i = 0; i < trial_history.size(); ++i) {
        std::cout << "[";
        for (size_t j = 0; j < trial_history[i].y_coords.size(); ++j) {
            std::cout << (std::isnan(trial_history[i].y_coords[j]) ? "null" : std::to_string(trial_history[i].y_coords[j])) << (j == trial_history[i].y_coords.size() - 1 ? "" : ", ");
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

void PrintLevelLines(
    const std::vector<double>& x_grid, const std::vector<double>& y_grid,
    const std::vector<std::vector<double>>& z_grid_values) {

    char* original_locale_numeric = std::setlocale(LC_NUMERIC, nullptr);
    std::string original_locale_numeric_str = (original_locale_numeric) ? original_locale_numeric : "";
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
            std::cout << (std::isnan(z_grid_values[i][j]) ? "null" : std::to_string(z_grid_values[i][j])) << (j == z_grid_values[i].size() - 1 ? "" : ", ");
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

} // namespace JsonLogger