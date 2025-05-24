#include "minimizer.h" 

// Определения глобальных переменных, используемых в mapd и node
int n1, nexp, l, iq, iu[10], iv[10];

// Вспомогательная функция для mapd
void node ( int is_param )
{
    /* calculate iu, iv, l by is_param */
    int current_n_dim, i, j_idx, k1, k2, iff_val;

    current_n_dim = n1 + 1; // n1 - глобальная
    if ( is_param == 0 ) {
        l = n1; // l - глобальная
        for ( i = 0; i < current_n_dim; i++ ) {
            iu[i] = -1; iv[i] = -1; // iu, iv - глобальные
        }
    } else if ( is_param == (nexp - 1) ) { // nexp - глобальная
        l = n1;
        iu[0] = 1;
        iv[0] = 1;
        for ( i = 1; i < current_n_dim; i++ ) {
            iu[i] = -1; iv[i] = -1;
        }
        iv[n1] = 1;
    } else {
        iff_val = nexp;
        k1 = -1;
        for ( i = 0; i < current_n_dim; i++ ) {
            iff_val = iff_val / 2;
            if ( is_param >= iff_val ) {
                if ( (is_param == iff_val) && (is_param != 1) ) { l = i; iq = -1; } // iq - глобальная
                is_param = is_param - iff_val;
                k2 = 1;
            } else {
                k2 = -1;
                if ( (is_param == (iff_val - 1)) && (is_param != 0) ) { l = i; iq = 1; }
            }
            j_idx = -k1 * k2;
            iv[i] = j_idx;
            iu[i] = j_idx;
            k1 = k2;
        }
        iv[l] = iv[l] * iq; // l - глобальная, обновлена в цикле выше
        iv[n1] = -iv[n1];   // n1 - глобальная
    }
}

// Реализация кривой Пеано
void mapd( double x_input, int m_order, double* y_output, int n_dim, int key_type )
{
    /* mapping y_output(x_input) : 1 - center, 2 - line, 3 - node */
    double d_current, mne_val, dd_val, dr_val;
    double p_calc, r_calc;
    int iw_arr[11];
    int it_val, is_val, i_loop, j_loop, k_idx = 0;

    p_calc = 0.0;
    n1 = n_dim - 1;
    for ( nexp = 1, i_loop = 0; i_loop < n_dim; nexp *= 2, i_loop++ );
    d_current = x_input;
    r_calc = 0.5;
    it_val = 0;
    dr_val = static_cast<double>(nexp);
    for ( mne_val = 1.0, i_loop = 0; i_loop < m_order; mne_val *= dr_val, i_loop++ );

    for ( i_loop = 0; i_loop < n_dim; i_loop++ ) {
        iw_arr[i_loop] = 1; y_output[i_loop] = 0.0;
    }

    if ( key_type == 2 ) {
        if (mne_val > 1e-9) { // Защита от деления на ноль или очень маленькое mne_val
            d_current = d_current * (1.0 - 1.0 / mne_val);
        } else {
            d_current = 0.0;
        }
        k_idx = 0;
    } else if ( key_type > 2 ) {
        if (nexp > 0) dr_val = mne_val / nexp; // nexp должно быть > 0
        else dr_val = mne_val;
        dr_val = dr_val - fmod(dr_val, 1.0);
        dd_val = mne_val - dr_val;
        if (abs(dd_val) > 1e-9) { // Защита от деления на ноль
            d_current = d_current * dd_val;
        } else {
            d_current = 0.0;
        }

        dd_val = d_current - fmod(d_current, 1.0);
        d_current = floor(d_current);
        if (nexp -1 > 0) { // Защита от деления на ноль
            d_current = d_current + (d_current - 1.0) / (nexp - 1.0);
        }
        d_current = floor(d_current);
        if (mne_val > 1e-9) d_current = d_current * (1.0 / mne_val);
        else d_current = 0.0;
    }

    for ( j_loop = 0; j_loop < m_order; j_loop++ ) {
        iq = 0; // iq - глобальная
        if ( x_input == 1.0 && j_loop == 0 ) {
            is_val = nexp - 1;
            d_current = 0.0;
        } else {
            d_current = d_current * nexp;
            is_val = static_cast<int>(floor(d_current));
            if (is_val >= nexp) is_val = nexp -1; // Ограничение сверху, если d_current было близко к 1.0 до умножения
            if (is_val < 0) is_val = 0;          // Ограничение снизу
            d_current = d_current - static_cast<double>(is_val);
        }

        node(is_val); // iu, iv, l, iq обновляются в node

        int temp_iu = iu[0];
        iu[0] = iu[it_val];
        iu[it_val] = temp_iu;

        int temp_iv = iv[0];
        iv[0] = iv[it_val];
        iv[it_val] = temp_iv;

        // Обновление l (глобальная)
        if ( l == 0 && it_val != 0) { // Если l было 0, и it_val не 0, l становится it_val
            l = it_val;
        } else if ( l == it_val ) { // Если l совпадает с it_val, l обнуляется
            l = 0;
        }
        // Если l != 0 и l != it_val, l остается без изменений - это неявное условие

        // Обновление k_idx в зависимости от iq и is_val
        if ( (iq > 0) || ((iq == 0) && (is_val == 0)) ) {
            k_idx = l;
        } else if ( iq < 0 ) {
            k_idx = ( it_val == n1 ) ? 0 : n1;
        }

        r_calc = r_calc * 0.5;
        it_val = l;

        for ( i_loop = 0; i_loop < n_dim; i_loop++ ) {
            iu[i_loop] = iu[i_loop] * iw_arr[i_loop];
            iw_arr[i_loop] = -iv[i_loop] * iw_arr[i_loop];
            p_calc = r_calc * static_cast<double>(iu[i_loop]);
            p_calc = p_calc + y_output[i_loop];
            y_output[i_loop] = p_calc;
        }
    }

    if ( key_type == 2 ) {
        int sign_val = (is_val == (nexp - 1)) ? -1 : 1;
        if (k_idx < 0 || k_idx >= 10) k_idx = 0;

        p_calc = 2.0 * static_cast<double>(sign_val) * static_cast<double>(iu[k_idx]) * r_calc * d_current;
        p_calc = y_output[k_idx] - p_calc;
        y_output[k_idx] = p_calc;
    } else if ( key_type == 3 ) {
        for ( i_loop = 0; i_loop < n_dim; i_loop++ ) {
            p_calc = r_calc * static_cast<double>(iu[i_loop]);
            p_calc = p_calc + y_output[i_loop];
            y_output[i_loop] = p_calc;
        }
    }
}

// Функция для сброса глобальных переменных, используемых в mapd
void resetMappingGlobals(int dim) {
    n1 = 0;
    nexp = 1;
    l = 0;
    iq = 0;
    for (int i = 0; i < 10; ++i) { // 10 - размер массивов iu/iv
        iu[i] = 0;
        iv[i] = 0;
    }
}

// Обёртка вокруг кривой Пеано: принимает одномерное значение t и возвращает вектор размерности n.
vector<double> peanoMapping(double x_param, int m_order, int n_dim, int key_type) {
    resetMappingGlobals(n_dim); // Сбрасываем глобальные переменные
    vector<double> y_result(n_dim); // Инициализируем нулями по умолчанию
    mapd(x_param, m_order, y_result.data(), n_dim, key_type);
    for(int i = 0; i < n_dim; ++i ) {
        y_result[i] += 0.5; // Сдвиг для получения значений в диапазоне [0,1], если mapd выдает [-0.5, 0.5]
    }
    return y_result;
}