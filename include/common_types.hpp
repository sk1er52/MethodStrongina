// common_types.hpp
#ifndef COMMON_TYPES_HPP
#define COMMON_TYPES_HPP

#include <vector>

/**
 * @brief Структура для представления пробной точки в процессе оптимизации.
 */
struct TrialPoint {
    /// Одномерная координата для поиска (параметр t для кривой Пеано в nD или x в 1D).
    double x_param; 
    
    /// Значение целевой функции в точке: f(y_coords).
    double z_value; 
    
    /// Координаты точки в исходном n-мерном пространстве.
    std::vector<double> y_coords; 

    /**
     * @brief Оператор сравнения для сортировки точек по координате x.
     */
    bool operator<(const TrialPoint& other) const {
        return x_param < other.x_param;
    }
};

#endif // COMMON_TYPES_HPP