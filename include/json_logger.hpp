#ifndef JSON_LOGGER_HPP
#define JSON_LOGGER_HPP

#include <string>
#include <vector>
#include "common_types.hpp"

/**
 * @brief Пространство имен для вывода данных в формате JSON.
 * Используется для передачи результатов из C++ в Python.
 */
namespace JsonLogger {

    /**
     * @brief Выводит результаты минимизации в формате JSON в стандартный поток вывода (stdout).
     */
    void PrintMinimizationResult(
        const std::string& problem_family, 
        int problem_idx_cpp, 
        double found_z, 
        const std::vector<double>& found_y, 
        int iterations,
        double time_ms, 
        int exit_main_criteria_count, 
        int exit_test_criteria_count,
        double known_opt_z, 
        const std::vector<double>& known_opt_y,
        const std::vector<TrialPoint>& trial_history
    );

    /**
     * @brief Выводит данные для построения линий уровня (тепловой карты) в формате JSON.
     */
    void PrintLevelLines(
        const std::vector<double>& x_grid, 
        const std::vector<double>& y_grid,
        const std::vector<std::vector<double>>& z_grid_values
    );

} // namespace JsonLogger

#endif // JSON_LOGGER_HPP