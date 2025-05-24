// minimizer.h
#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <cmath>
#include <vector>
#include <limits>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm> // Для std::max_element, std::max, std::min, std::lower_bound
#include <numeric>   // Для std::distance (не используется сейчас, но может пригодиться)
#include <stdexcept> // Для std::invalid_argument

// Используем std:: для краткости в этом файле
using namespace std;

// Модифицированная структура для представления точки
struct Point {
    double x_param;               // Одномерная координата для поиска (параметр t для кривой Пеано в nD)
    double z_value;               // Значение функции в точке: f(y_coords)
    std::vector<double> y_coords; // Координаты точки в исходном пространстве (n-мерном)
};

// Интерфейс для целевой функции (без изменений)
class FunctionInterface {
public:
    virtual double ComputeFunction(const std::vector<double>& x) const = 0;
    virtual std::vector<double> GetOptimumPoint() const = 0;
    // Добавляем GetOptimumValue, если его нет в вашем FunctionInterface,
    // но он есть в IOptProblem, который вы, вероятно, используете.
    // Если он уже есть, эта строка не нужна.
    // virtual double GetOptimumValue() const = 0; // Раскомментируйте, если нужно объявить здесь
    virtual ~FunctionInterface() {}
};

//
// Объявления функций кривой Пеано
//
extern int n1, nexp, l, iq, iu[10], iv[10];

void mapd( double x, int m, double* y, int n, int key );
void node ( int is );
void resetMappingGlobals(int n); // n - размерность
std::vector<double> peanoMapping(double x, int m, int n, int key);


//
// Шаблонный класс Minimizer
//
template <typename T_ProblemType>
class Minimizer {
private:
    std::vector<double> searchLeftBound_param;
    std::vector<double> searchRightBound_param;
    double epsilon_val;
    double r_strongin_parameter;
    const T_ProblemType& function_instance;
    std::ofstream logFile_stream;

    int iteration_counter;
    int exit_main_criteria_count;
    int exit_test_criteria_count;
    std::vector<Point> trial_points_list;

    bool use_peano_mapping_flag;
    int peano_mapping_m_order;
    int problem_actual_dimension;
    int peano_mapping_key_type;

    double current_m_phi_estimate;

    // --- Приватные методы с определениями прямо здесь (inline) ---
    
    void performFirstIteration() {
        logFile_stream << "Первая итерация:\n";
        Point p_left;
        p_left.x_param = searchLeftBound_param[0];
        if (!use_peano_mapping_flag) {
            p_left.y_coords = {p_left.x_param};
        } else {
            p_left.y_coords = peanoMapping(p_left.x_param, peano_mapping_m_order, problem_actual_dimension, peano_mapping_key_type);
        }
        p_left.z_value = function_instance.ComputeFunction(p_left.y_coords);
        trial_points_list.push_back(p_left);
        logTrialPointToFile(p_left, "Начальная точка (левая): ");

        Point p_right;
        p_right.x_param = searchRightBound_param[0];
         if (!use_peano_mapping_flag) {
            p_right.y_coords = {p_right.x_param};
        } else {
            p_right.y_coords = peanoMapping(p_right.x_param, peano_mapping_m_order, problem_actual_dimension, peano_mapping_key_type);
        }
        p_right.z_value = function_instance.ComputeFunction(p_right.y_coords);
        
        if (p_right.x_param < p_left.x_param) {
            trial_points_list.insert(trial_points_list.begin(), p_right);
            logTrialPointToFile(p_right, "Начальная точка (правая, вставлена первой): ");
        } else {
            trial_points_list.push_back(p_right);
            logTrialPointToFile(p_right, "Начальная точка (правая): ");
        }
    }

    double computeIntervalCharacteristicR(const Point& p1, const Point& p2) const {
        double dx_param = p2.x_param - p1.x_param;
        if (std::abs(dx_param) < 1e-12) { 
            return -std::numeric_limits<double>::infinity();
        }
        double dz = p2.z_value - p1.z_value;
        double reliable_m_phi = std::max(current_m_phi_estimate, 1e-9); // Защита от нуля или очень малого m_phi
        return reliable_m_phi * dx_param + (dz * dz) / (reliable_m_phi * dx_param) - 2.0 * (p1.z_value + p2.z_value);
    }

    size_t findIntervalIndexWithMaxR() {
        if (trial_points_list.size() < 2) return 0; // Не должно случиться после performFirstIteration
        double max_R_val = -std::numeric_limits<double>::infinity();
        size_t max_R_idx = 0; 
        for (size_t i = 0; i < trial_points_list.size() - 1; ++i) {
            double current_R_val = computeIntervalCharacteristicR(trial_points_list[i], trial_points_list[i+1]);
            if (current_R_val > max_R_val) {
                max_R_val = current_R_val;
                max_R_idx = i;
            }
        }
        logFile_stream << "Макс. характеристика R = " << max_R_val << " на интервале x_param=["
                       << trial_points_list[max_R_idx].x_param << ", "
                       << trial_points_list[max_R_idx+1].x_param << "]\n";
        return max_R_idx;
    }

    Point generateNewTrialPoint(size_t interval_index_left) {
        const Point& p_left = trial_points_list[interval_index_left];
        const Point& p_right = trial_points_list[interval_index_left + 1];
        Point new_p;
        double reliable_m_phi_for_new_point = std::max(current_m_phi_estimate, 1e-9);
        new_p.x_param = 0.5 * (p_left.x_param + p_right.x_param) -
                        (p_right.z_value - p_left.z_value) / (2.0 * reliable_m_phi_for_new_point);

        double interval_width = p_right.x_param - p_left.x_param;
        // Минимальный шаг от границы, чтобы избежать точного совпадения
        // и обеспечить некоторое продвижение.
        double min_rel_step = 0.001; // 0.1% от ширины интервала
        double abs_min_step = 1e-9;  // Абсолютный минимальный шаг
        
        double step_from_boundary = std::max(abs_min_step, interval_width * min_rel_step);

        if (interval_width <= 0) { // Нулевой или инвертированный интервал
            new_p.x_param = p_left.x_param; // Просто возвращаем левую точку (или середину, если они разные)
        } else {
            if (new_p.x_param <= p_left.x_param) {
                new_p.x_param = p_left.x_param + step_from_boundary;
            }
            if (new_p.x_param >= p_right.x_param) {
                new_p.x_param = p_right.x_param - step_from_boundary;
            }
            // Если после коррекции точка "перепрыгнула" на другую сторону (очень маленький интервал)
            // или осталась на границе, ставим ее в середину.
            if (new_p.x_param <= p_left.x_param || new_p.x_param >= p_right.x_param) {
                 new_p.x_param = p_left.x_param + interval_width / 2.0;
            }
        }
        
        if (!use_peano_mapping_flag) {
            new_p.y_coords = {new_p.x_param};
        } else {
            new_p.y_coords = peanoMapping(new_p.x_param, peano_mapping_m_order, problem_actual_dimension, peano_mapping_key_type);
        }
        new_p.z_value = function_instance.ComputeFunction(new_p.y_coords);
        return new_p;
    }
    
    double getHolderAlphaValue() const { // Не используется в текущей версии checkStoppingConditions
        return use_peano_mapping_flag ? (1.0 / static_cast<double>(problem_actual_dimension)) : 1.0;
    }
    
    bool checkStoppingConditions(const Point& p_left_of_interval, const Point& /* new_p_ref - не используется */, const Point& p_right_of_interval) {
        // double chosen_interval_length = p_right_of_interval.x_param - p_left_of_interval.x_param;
        
        // // 1. Основное условие останова по Стронгину (длина выбранного интервала)
        // if (chosen_interval_length <= epsilon_val) {
        //      logFile_stream << "Условие останова (основное): длина выбранного интервала (" << chosen_interval_length 
        //                << ") <= epsilon_abs (" << epsilon_val << ")\n";
        //      exit_main_criteria_count++;
        //      return true;
        // }

        // 2. Дополнительное условие для МНОГОМЕРНОГО случая (если используется Пеано)
        //    Остановка, если максимум модуля разности |z_i - z*| <= epsilon
        if (use_peano_mapping_flag) {
            // Предполагаем, что T_ProblemType имеет метод GetOptimumValue()
            // Это должно быть частью вашего интерфейса IOptProblem или аналогичного.
            double known_optimum_z_value = function_instance.GetOptimumValue();

            if (!std::isnan(known_optimum_z_value) && !trial_points_list.empty()) {
                double max_abs_diff_z = 0.0;

                for (const auto& p : trial_points_list) {
                    double current_diff = std::abs(p.z_value - known_optimum_z_value);
                    if (current_diff > max_abs_diff_z) {
                        max_abs_diff_z = current_diff;
                    }
                }

                logFile_stream << "Проверка доп. критерия (nD): max|z_i - z*| = " << max_abs_diff_z 
                               << ", epsilon = " << epsilon_val << "\n";

                if (max_abs_diff_z <= epsilon_val) {
                    logFile_stream << "Условие останова (дополнительное, nD): max|z_i - z*| (" << max_abs_diff_z 
                                   << ") <= epsilon (" << epsilon_val << ").\n";
                    exit_test_criteria_count++; // Используем другой счетчик
                    return true;
                }
            } else if (std::isnan(known_optimum_z_value)) {
                 logFile_stream << "Проверка доп. критерия (nD): известное значение оптимума z* недоступно (NAN).\n";
            }
        }
        // Дополнительное условие для 1D, если оно нужно, можно добавить здесь.
        // Например, сравнение x_param лучшей найденной точки с x_param известного оптимума.
        // if (!use_peano_mapping_flag) { ... }

        return false;
    }

    // Объявление с аргументом по умолчанию
    void logTrialPointToFile(const Point& p, const std::string& prefix = "") {
        if (!logFile_stream.is_open()) return;
        logFile_stream << prefix;
        logFile_stream << "x_param: " << p.x_param << ", z_value: " << p.z_value << ", y_coords: [";
        for (size_t i = 0; i < p.y_coords.size(); ++i) {
            logFile_stream << p.y_coords[i] << (i == p.y_coords.size() - 1 ? "" : ", ");
        }
        logFile_stream << "]" << std::endl;
    }

public:
    Minimizer(const std::vector<double>& problem_domain_a, const std::vector<double>& problem_domain_b,
              double eps, double r_val, const T_ProblemType& func)
        : epsilon_val(eps), r_strongin_parameter(r_val), function_instance(func),
          iteration_counter(0), exit_main_criteria_count(0), exit_test_criteria_count(0),
          use_peano_mapping_flag(false), problem_actual_dimension(1), current_m_phi_estimate(1.0)
    {
        if (problem_domain_a.empty() || problem_domain_b.empty()) {
            throw std::invalid_argument("Границы задачи (1D) не могут быть пустыми.");
        }
        searchLeftBound_param = problem_domain_a;
        searchRightBound_param = problem_domain_b;
        logFile_stream.open("minimization_log.txt", std::ios::out);
        if (!logFile_stream.is_open()) {
            std::cerr << "Ошибка открытия файла журнала!" << std::endl;
        }
    }

    Minimizer(double peano_search_param_a, double peano_search_param_b,
              double eps, double r_val, const T_ProblemType& func,
              int mapping_m_order, int original_problem_dimension, int mapping_key)
        : epsilon_val(eps), r_strongin_parameter(r_val), function_instance(func),
          iteration_counter(0), exit_main_criteria_count(0), exit_test_criteria_count(0),
          use_peano_mapping_flag(true), peano_mapping_m_order(mapping_m_order),
          problem_actual_dimension(original_problem_dimension), peano_mapping_key_type(mapping_key), current_m_phi_estimate(1.0)
    {
        searchLeftBound_param = {peano_search_param_a};
        searchRightBound_param = {peano_search_param_b};
        logFile_stream.open("minimization_log.txt", std::ios::out);
        if (!logFile_stream.is_open()) {
            std::cerr << "Ошибка открытия файла журнала!" << std::endl;
        }
    }

    ~Minimizer() {
        if (logFile_stream.is_open()) {
            logFile_stream.close();
        }
    }

    std::vector<double> findMinimum() {
        trial_points_list.clear();
        iteration_counter = 0;
        exit_main_criteria_count = 0;
        exit_test_criteria_count = 0;
        const int MAX_ITERATIONS_LIMIT = 10000; // Можно сделать настраиваемым

        performFirstIteration();

        if (trial_points_list.size() < 2) {
            logFile_stream << "Ошибка: Первая итерация не создала достаточно точек для начала основного цикла." << std::endl;
            if (!trial_points_list.empty()) return trial_points_list[0].y_coords;
            return {}; 
        }
        
        for (int current_iter = 0; current_iter < MAX_ITERATIONS_LIMIT; ++current_iter) {
            iteration_counter = current_iter + 1;

            // 1. Обновить оценку m_phi_estimate (M*r)
            double max_abs_slope = 0.0;
            for (size_t i = 0; i < trial_points_list.size() - 1; ++i) {
                double dx_param = trial_points_list[i+1].x_param - trial_points_list[i].x_param; // p[i+1] всегда правее p[i]
                if (dx_param > 1e-12) { // Избегаем деления на очень малое число
                    double slope = std::abs((trial_points_list[i+1].z_value - trial_points_list[i].z_value) / dx_param);
                    if (slope > max_abs_slope) {
                        max_abs_slope = slope;
                    }
                }
            }
            current_m_phi_estimate = (max_abs_slope > 1e-9) ? r_strongin_parameter * max_abs_slope : 1.0;

            // 2. Найти интервал с максимальной характеристикой R
            size_t interval_to_split_idx = findIntervalIndexWithMaxR();
            
            // 3. Проверить условия останова для ВЫБРАННОГО интервала *перед* генерацией новой точки
            if (checkStoppingConditions(trial_points_list[interval_to_split_idx], 
                                        trial_points_list[interval_to_split_idx], // Фиктивная new_p, не используется для основного условия
                                        trial_points_list[interval_to_split_idx+1])) {
                break; 
            }

            // 4. Сгенерировать новую пробную точку
            Point new_point = generateNewTrialPoint(interval_to_split_idx);
            
            // 5. Проверка, не совпадает ли новая точка с существующими (слишком близко)
            bool already_exists_nearby = false;
            for(const auto& existing_p : trial_points_list) {
                if (std::abs(existing_p.x_param - new_point.x_param) < 1e-10) { // Порог близости
                    already_exists_nearby = true;
                    logFile_stream << "Итерация " << iteration_counter << ": Новая точка x_param=" << new_point.x_param 
                                   << " слишком близка к существующей. Попытка сдвига или пропуск.\n";
                    // Если точка очень близка, пытаемся немного сдвинуть или обработать иначе,
                    // чтобы избежать застревания. В generateNewTrialPoint уже есть логика отступа.
                    // Если и после этого точка близка, возможно, достигнут предел точности.
                    // Здесь можно было бы пропустить итерацию, но это может привести к зацикливанию,
                    // если всегда выбирается один и тот же "застрявший" интервал.
                    // Вместо continue, лучше положиться на то, что generateNewTrialPoint пытается найти место.
                    // Если интервал настолько мал, что некуда ставить, условие останова должно сработать.
                    break; 
                }
            }
             if (already_exists_nearby && (trial_points_list[interval_to_split_idx+1].x_param - trial_points_list[interval_to_split_idx].x_param <= epsilon_val * 0.1 )) {
                logFile_stream << "Интервал очень мал и новая точка близка к существующей. Вероятна остановка.\n";
                // Не добавляем точку, пусть сработает условие останова по длине интервала.
                // Это более безопасный выход, чем continue, который может зациклить.
             } else {
                // 6. Вставить новую точку, сохраняя сортировку по x_param
                auto it_insert = std::lower_bound(trial_points_list.begin(), trial_points_list.end(), new_point,
                                        [](const Point& p1, const Point& p2){ return p1.x_param < p2.x_param; });
                trial_points_list.insert(it_insert, new_point);
                logTrialPointToFile(new_point, "Итерация " + std::to_string(iteration_counter) + ": Новая точка ");
             }
        }

        // 7. Поиск лучшей точки среди всех исследованных
        size_t min_z_idx = 0;
        if (!trial_points_list.empty()) {
            for (size_t i = 1; i < trial_points_list.size(); ++i) {
                if (trial_points_list[i].z_value < trial_points_list[min_z_idx].z_value) {
                    min_z_idx = i;
                }
            }
            logFile_stream << "Минимум найден в точке (x_param=" << trial_points_list[min_z_idx].x_param
                           << ", z=" << trial_points_list[min_z_idx].z_value << ")\n";
            logTrialPointToFile(trial_points_list[min_z_idx], "Лучшая точка: ");
            return trial_points_list[min_z_idx].y_coords;
        }
        
        logFile_stream << "Список пробных точек пуст или произошла ошибка, оптимум не найден.\n";
        return {};
    }

    int GetIterationCount() const { return iteration_counter; }
    int GetExitMainCount() const { return exit_main_criteria_count; }
    int GetExitTestCount() const { return exit_test_criteria_count; }
};

#endif // MINIMIZER_H