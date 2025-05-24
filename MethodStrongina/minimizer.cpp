// minimizer_functions.cpp (или minimizer.cpp, если вы его так назвали)

#include "minimizer.h" // Для объявлений функций Пеано, структуры Point и extern глобальных переменных
#include <cmath>       // Для fmod

// Определения глобальных переменных, используемых в mapd и node
// (они объявлены как extern в minimizer.h)
int n1, nexp, l, iq, iu[10], iv[10];

// Вспомогательная функция для mapd
void node ( int is_param )
{
    /* calculate iu, iv, l by is_param */
    int current_n_dim, i, j_idx, k1, k2, iff_val;

    current_n_dim = n1 + 1;
    if ( is_param == 0 ) {
        l = n1;
        for ( i = 0; i < current_n_dim; i++ ) {
            iu[i] = -1; iv[i] = -1;
        }
    } else if ( is_param == (nexp - 1) ) {
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
                if ( (is_param == iff_val) && (is_param != 1) ) { l = i; iq = -1; }
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
        iv[l] = iv[l] * iq;
        iv[n1] = -iv[n1];
    }
}

// Реализация кривой Пеано
void mapd( double x_input, int m_order, double* y_output, int n_dim, int key_type )
{
    double d_current, mne_val, dd_val, dr_val;
    double p_calc, r_calc;
    int iw_arr[11];
    int it_val, is_val, i_loop, j_loop, k_idx = 0;

    p_calc = 0.0;
    n1 = n_dim - 1;
    if (n_dim <= 0) return; // Защита от некорректной размерности
    for ( nexp = 1, i_loop = 0; i_loop < n_dim; nexp *= 2, i_loop++ );
    if (nexp == 0 && n_dim > 0) nexp = 1; // Избегаем nexp=0 если n_dim > 0 но мал

    d_current = x_input;
    r_calc = 0.5;
    it_val = 0;
    dr_val = static_cast<double>(nexp);
    for ( mne_val = 1.0, i_loop = 0; i_loop < m_order; mne_val *= dr_val, i_loop++ );

    for ( i_loop = 0; i_loop < n_dim; i_loop++ ) {
        iw_arr[i_loop] = 1; y_output[i_loop] = 0.0;
    }

    if ( key_type == 2 ) {
        if (std::abs(mne_val) > 1e-12) {
             d_current = d_current * (1.0 - 1.0 / mne_val);
        } else {
             d_current = (x_input == 1.0) ? 1.0 : 0.0; // Специальная обработка для x=1
        }
        k_idx = 0;
    } else if ( key_type > 2 ) {
        if (nexp > 0) dr_val = mne_val / nexp;
        else dr_val = mne_val; 
        dr_val = floor(dr_val); // dr = dr - fmod(dr, 1.0)
        
        dd_val = mne_val - dr_val;
        d_current = x_input * dd_val; // В оригинале было dr = d * dd
        
        d_current = floor(d_current); // dd_val = dr - fmod(dr,1.0) -> это floor(dr)
        if (nexp -1 != 0) { // Защита от деления на ноль
            d_current = d_current + (d_current - 1.0) / (nexp - 1.0);
        }
        d_current = floor(d_current);
        if (std::abs(mne_val) > 1e-12) d_current = d_current * (1.0 / mne_val);
        else d_current = 0.0;
    }

    for ( j_loop = 0; j_loop < m_order; j_loop++ ) {
        iq = 0;
        if ( x_input == 1.0 && j_loop == 0 && key_type != 2 ) { // Условие для x_input=1.0 (кроме key_type=2, где d_current уже обработан)
            is_val = nexp - 1;
            d_current = 0.0; 
        } else {
            d_current = d_current * nexp;
            is_val = static_cast<int>(floor(d_current));
            if (is_val >= nexp) is_val = nexp -1; 
            if (is_val < 0) is_val = 0;          
            d_current = d_current - static_cast<double>(is_val);
        }

        node(is_val); 

        int temp_iu = iu[0]; iu[0] = iu[it_val]; iu[it_val] = temp_iu;
        int temp_iv = iv[0]; iv[0] = iv[it_val]; iv[it_val] = temp_iv;

        if ( l == 0 && it_val != 0) l = it_val;
        else if ( l == it_val ) l = 0;
        
        if ( (iq > 0) || ((iq == 0) && (is_val == 0)) ) k_idx = l;
        else if ( iq < 0 ) k_idx = ( it_val == n1 ) ? 0 : n1;
        
        r_calc *= 0.5; 
        it_val = l; 

        for ( i_loop = 0; i_loop < n_dim; i_loop++ ) {
            // Проверка границ для iu, iv, iw_arr, если n_dim может быть больше 10
            if (i_loop >= 10 || k_idx >= 10) { /* обработка ошибки или пропуск */ continue; }

            iu[i_loop] = iu[i_loop] * iw_arr[i_loop];
            iw_arr[i_loop] = -iv[i_loop] * iw_arr[i_loop];
            p_calc = r_calc * static_cast<double>(iu[i_loop]);
            p_calc = p_calc + y_output[i_loop];
            y_output[i_loop] = p_calc;
        }
    }

    if ( key_type == 2 ) {
        // is_val здесь будет последним вычисленным значением is_val из цикла по j_loop
        int sign_val = (is_val == (nexp - 1)) ? -1 : 1;
        if (k_idx < 0 || k_idx >= 10) k_idx = 0; 

        p_calc = 2.0 * static_cast<double>(sign_val) * static_cast<double>(iu[k_idx]) * r_calc * d_current;
        y_output[k_idx] = y_output[k_idx] - p_calc; // Изменение было y[k] = y[k] - p
    } else if ( key_type == 3 ) {
        for ( i_loop = 0; i_loop < n_dim; i_loop++ ) {
            if (i_loop >= 10) continue;
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
std::vector<double> peanoMapping(double x_param, int m_order, int n_dim, int key_type) {
    if (n_dim <= 0 || n_dim > 10) { // Ограничение на размерность из-за iu/iv
        // Можно бросить исключение или вернуть пустой вектор
        std::cerr << "Ошибка Peano: некорректная размерность n_dim = " << n_dim << std::endl;
        return std::vector<double>(n_dim > 0 ? n_dim : 1, 0.0); // Возвращаем нули, если размерность некорректна, но >0
    }
    resetMappingGlobals(n_dim);
    std::vector<double> y_result(n_dim);
    mapd(x_param, m_order, y_result.data(), n_dim, key_type);
    for(int i = 0; i < n_dim; ++i ) {
        y_result[i] += 0.5;
    }
    return y_result;
}