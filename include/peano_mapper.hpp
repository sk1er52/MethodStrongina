#ifndef PEANO_MAPPER_HPP
#define PEANO_MAPPER_HPP

#include <vector>

/**
 * @brief Класс для отображения одномерного отрезка [0, 1] в n-мерный гиперкуб 
 * с помощью развертки (кривой) Пеано.
 */
class PeanoMapper {
public:
    PeanoMapper();

    /**
     * @brief Выполняет развертку Пеано.
     * 
     * @param x_param Одномерный параметр на отрезке [0, 1].
     * @param m_order Плотность развертки (порядок).
     * @param n_dim Размерность исходного пространства.
     * @param key_type Тип развертки (1, 2 или 3).
     * @return std::vector<double> Координаты в n-мерном пространстве.
     */
    std::vector<double> Map(double x_param, int m_order, int n_dim, int key_type);

private:
    void Node(int is_param);
    void Mapd(double x_input, int m_order, double* y_output, int n_dim, int key_type);
    void ResetState(int n_dim);

    // Внутреннее состояние (ранее это были глобальные переменные)
    int n1_;
    int nexp_;
    int l_;
    int iq_;
    int iu_[10];
    int iv_[10];
};

#endif // PEANO_MAPPER_HPP