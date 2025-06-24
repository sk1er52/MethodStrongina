// minimizer_functions.cpp (или minimizer.cpp, если вы его так назвали)

#include "minimizer.h" // Для объявлений функций Пеано, структуры Point и extern глобальных переменных
#include <algorithm> // Для std::max, std::min (если понадобятся где-то еще)
#include <cmath>     // Для fmod, floor
#include <numeric> // Для std::abs в С++17 (хотя cmath::abs тоже работает для double)
#include <vector> // Для std::vector (используется в peanoMapping)


// Определения глобальных переменных, используемых в mapd и node
// (они объявлены как extern в minimizer.h)
int n1, nexp, l, iq, iu[10], iv[10];

// Вспомогательная функция для mapd
// (код node остается таким же, как вы предоставили)
void node(int is_param) {
  /* calculate iu, iv, l by is_param */
  int current_n_dim, i, j_idx, k1, k2, iff_val;

  current_n_dim = n1 + 1;
  if (is_param == 0) {
    l = n1;
    for (i = 0; i < current_n_dim; i++) {
      iu[i] = -1;
      iv[i] = -1;
    }
  } else if (is_param == (nexp - 1)) {
    l = n1;
    iu[0] = 1;
    iv[0] = 1;
    for (i = 1; i < current_n_dim; i++) {
      iu[i] = -1;
      iv[i] = -1;
    }
    iv[n1] = 1;
  } else {
    iff_val = nexp;
    k1 = -1;
    for (i = 0; i < current_n_dim; i++) {
      iff_val = iff_val / 2;
      if (is_param >= iff_val) {
        if ((is_param == iff_val) && (is_param != 1)) {
          l = i;
          iq = -1;
        }
        is_param = is_param - iff_val;
        k2 = 1;
      } else {
        k2 = -1;
        if ((is_param == (iff_val - 1)) && (is_param != 0)) {
          l = i;
          iq = 1;
        }
      }
      j_idx = -k1 * k2;
      iv[i] = j_idx;
      iu[i] = j_idx;
      k1 = k2;
    }
    iv[l] = iv[l] * iq;
    iv[n1] = -iv[n1];
  }
}

// Реализация кривой Пеано с изменениями для x_input = 1.0
void mapd(double x_input, int m_order, double *y_output, int n_dim,
          int key_type) {
  double mne_val, dd_val, dr_val;
  double p_calc, r_calc; // Используем double
  int iw_arr[11];        // Размер массива iw, убедитесь, что 11 достаточно
  int it_val, is_val, i_loop, j_loop, k_idx = 0;

  p_calc = 0.0;
  n1 = n_dim - 1;
  if (n_dim <= 0 ||
      n_dim >= 10) { // Защита, т.к. iu/iv/iw_arr имеют размер 10/11
    for (i_loop = 0; i_loop < n_dim; i_loop++)
      y_output[i_loop] = 0.5; // Возвращаем центр по умолчанию
    return;
  }
  for (nexp = 1, i_loop = 0; i_loop < n_dim; nexp *= 2, i_loop++)
    ;
  if (nexp <= 0 && n_dim > 0)
    nexp = 1; // Избегаем nexp=0 или отрицательного

  double d_for_key_logic =
      x_input; // Копия x_input для логики ключей key_type 2 и >2

  r_calc = 0.5;
  it_val = 0;
  dr_val = static_cast<double>(nexp);
  for (mne_val = 1.0, i_loop = 0; i_loop < m_order; mne_val *= dr_val, i_loop++)
    ;

  for (i_loop = 0; i_loop < n_dim; i_loop++) {
    iw_arr[i_loop] = 1;
    y_output[i_loop] = 0.0;
  }

  // Предварительная обработка d_for_key_logic, если key_type != 1
  // Эта логика может изменять значение, которое будет использоваться для
  // извлечения "цифр" Если она должна применяться к x_input ДО извлечения цифр,
  // то d_for_digit_extraction ниже нужно будет инициализировать этим измененным
  // значением. Оригинальный код, похоже, менял 'd' (которое было и для цифр, и
  // для key_type==2). Это сложный момент, так как 'd' имело двойное назначение.
  // Пока оставим d_for_key_logic для использования в финальной коррекции, а для
  // цифр - чистый x_input.

  if (key_type == 2) {
    if (std::abs(mne_val) > 1e-12) {
      d_for_key_logic = d_for_key_logic * (1.0 - 1.0 / mne_val);
    } else {
      d_for_key_logic = (x_input == 1.0) ? 1.0 : 0.0;
    }
    k_idx = 0;
  } else if (key_type > 2) {
    // Логика из вашего оригинального кода для key > 2, применяемая к
    // d_for_key_logic Это d_for_key_logic будет использоваться для извлечения
    // is_val в цикле. Оригинальное 'd' модифицировалось здесь.
    double temp_d_for_key_gt2 = x_input; // Работаем с копией x_input

    dr_val = mne_val / nexp;
    dr_val = dr_val - fmod(dr_val, 1.0);
    dd_val = mne_val - dr_val;
    dr_val = temp_d_for_key_gt2 * dd_val;
    dd_val = dr_val - fmod(dr_val, 1.0);
    if (nexp - 1 != 0)
      dr_val = dd_val + (dd_val - 1.0) / (nexp - 1.0);
    else
      dr_val = dd_val;
    dd_val = dr_val - fmod(dr_val, 1.0);
    if (std::abs(mne_val) > 1e-12)
      d_for_key_logic = dd_val * (1.0 / mne_val);
    else
      d_for_key_logic = 0.0;
  }
  // Теперь d_for_key_logic содержит значение, которое будет использоваться для
  // итеративного извлечения цифр. А x_input нам нужен для корректного
  // определения случая x_input == 1.0

  double d_iter_for_digits = (key_type == 1) ? x_input : d_for_key_logic;
  double final_d_remainder =
      0.0; // Остаток после всех m_order итераций, для key_type==2

  for (j_loop = 0; j_loop < m_order; j_loop++) {
    iq = 0;
    // Логика извлечения is_val
    if (std::abs(x_input - 1.0) < 1e-9 && key_type != 2) {
      // Для x_input=1.0 (и key_type не 2, где d_iter уже специально обработан),
      // is_val всегда nexp-1
      is_val = nexp - 1;
      // d_iter_for_digits не меняем, чтобы is_val оставался nexp-1
      // final_d_remainder будет 0.0, если это последняя итерация
      if (j_loop == m_order - 1)
        final_d_remainder = 0.0;
    } else {
      d_iter_for_digits = d_iter_for_digits * static_cast<double>(nexp);
      is_val = static_cast<int>(floor(d_iter_for_digits));
      if (is_val >= nexp)
        is_val = nexp - 1;
      if (is_val < 0)
        is_val = 0;
      d_iter_for_digits = d_iter_for_digits - static_cast<double>(is_val);
      if (j_loop == m_order - 1)
        final_d_remainder = d_iter_for_digits;
    }

    node(is_val);

    int temp_i;
    temp_i = iu[0];
    iu[0] = iu[it_val];
    iu[it_val] = temp_i;
    temp_i = iv[0];
    iv[0] = iv[it_val];
    iv[it_val] = temp_i;

    if (l == 0 && it_val != 0)
      l = it_val;
    else if (l == it_val)
      l = 0;

    if ((iq > 0) || ((iq == 0) && (is_val == 0)))
      k_idx = l;
    else if (iq < 0)
      k_idx = (it_val == n1) ? 0 : n1;

    r_calc = r_calc * 0.5;
    it_val = l;
    for (i_loop = 0; i_loop < n_dim; i_loop++) {
      iu[i_loop] = iu[i_loop] * iw_arr[i_loop];
      iw_arr[i_loop] = -iv[i_loop] * iw_arr[i_loop];
      p_calc = r_calc * static_cast<double>(iu[i_loop]);
      p_calc = p_calc + y_output[i_loop];
      y_output[i_loop] = p_calc;
    }
  }

  // Финальная коррекция, используя final_d_remainder (который является
  // последним d_iter_for_digits) is_val здесь - это is_val от *последней*
  // итерации j_loop
  if (key_type == 2) {
    int sign_val = (is_val == (nexp - 1)) ? -1 : 1;
    if (k_idx < 0 || k_idx >= 10)
      k_idx = 0;
    p_calc = 2.0 * static_cast<double>(sign_val) *
             static_cast<double>(iu[k_idx]) * r_calc * final_d_remainder;
    y_output[k_idx] = y_output[k_idx] - p_calc;
  } else if (key_type == 3) {
    for (i_loop = 0; i_loop < n_dim; i_loop++) {
      if (i_loop >= 10)
        continue;
      p_calc = r_calc * static_cast<double>(iu[i_loop]);
      p_calc = p_calc + y_output[i_loop];
      y_output[i_loop] = p_calc;
    }
  }
}

// Функция для сброса глобальных переменных
void resetMappingGlobals(int /*dim_unused*/) {
  n1 = 0;
  nexp = 1;
  l = 0;
  iq = 0;
  for (int i = 0; i < 10; ++i) {
    iu[i] = 0;
    iv[i] = 0;
  }
}

// Обёртка вокруг кривой Пеано
std::vector<double> peanoMapping(double x_param, int m_order, int n_dim,
                                 int key_type) {
  if (n_dim <= 0 ||
      n_dim >= ((sizeof(iu) / sizeof(iu[0])) - 1)) { // Проверка на размер iu/iv
    std::cerr << "Ошибка Peano: некорректная размерность n_dim = " << n_dim
              << " (должна быть >0 и <" << (sizeof(iu) / sizeof(iu[0])) - 1
              << ")" << std::endl;
    std::vector<double> error_vec(n_dim > 0 ? n_dim : 1,
                                  0.5); // Возвращаем центр
    // Можно бросить исключение: throw std::out_of_range("Invalid dimension for
    // Peano curve");
    return error_vec;
  }
  resetMappingGlobals(n_dim);
  std::vector<double> y_result(n_dim); // По умолчанию инициализируется нулями
  mapd(x_param, m_order, y_result.data(), n_dim, key_type);
  for (int i = 0; i < n_dim; ++i) {
    y_result[i] += 0.5;
  }
  return y_result;
}