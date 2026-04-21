#include "peano_mapper.hpp"
#include <cmath>
#include <iostream>

PeanoMapper::PeanoMapper() {
    ResetState(1);
}

void PeanoMapper::ResetState(int n_dim) {
    n1_ = 0;
    nexp_ = 1;
    l_ = 0;
    iq_ = 0;
    for (int i = 0; i < 10; ++i) {
        iu_[i] = 0;
        iv_[i] = 0;
    }
}

void PeanoMapper::Node(int is_param) {
    int current_n_dim = n1_ + 1;
    int i, j_idx, k1, k2, iff_val;

    if (is_param == 0) {
        l_ = n1_;
        for (i = 0; i < current_n_dim; i++) {
            iu_[i] = -1;
            iv_[i] = -1;
        }
    } else if (is_param == (nexp_ - 1)) {
        l_ = n1_;
        iu_[0] = 1;
        iv_[0] = 1;
        for (i = 1; i < current_n_dim; i++) {
            iu_[i] = -1;
            iv_[i] = -1;
        }
        iv_[n1_] = 1;
    } else {
        iff_val = nexp_;
        k1 = -1;
        for (i = 0; i < current_n_dim; i++) {
            iff_val = iff_val / 2;
            if (is_param >= iff_val) {
                if ((is_param == iff_val) && (is_param != 1)) {
                    l_ = i;
                    iq_ = -1;
                }
                is_param = is_param - iff_val;
                k2 = 1;
            } else {
                k2 = -1;
                if ((is_param == (iff_val - 1)) && (is_param != 0)) {
                    l_ = i;
                    iq_ = 1;
                }
            }
            j_idx = -k1 * k2;
            iv_[i] = j_idx;
            iu_[i] = j_idx;
            k1 = k2;
        }
        iv_[l_] = iv_[l_] * iq_;
        iv_[n1_] = -iv_[n1_];
    }
}

void PeanoMapper::Mapd(double x_input, int m_order, double* y_output, int n_dim, int key_type) {
    double mne_val, dd_val, dr_val;
    double p_calc, r_calc;
    int iw_arr[11];
    int it_val, is_val, i_loop, j_loop, k_idx = 0;

    p_calc = 0.0;
    n1_ = n_dim - 1;
    if (n_dim <= 0 || n_dim >= 10) {
        for (i_loop = 0; i_loop < n_dim; i_loop++) y_output[i_loop] = 0.5;
        return;
    }
    for (nexp_ = 1, i_loop = 0; i_loop < n_dim; nexp_ *= 2, i_loop++);
    if (nexp_ <= 0 && n_dim > 0) nexp_ = 1;

    double d_for_key_logic = x_input;
    r_calc = 0.5;
    it_val = 0;
    dr_val = static_cast<double>(nexp_);
    for (mne_val = 1.0, i_loop = 0; i_loop < m_order; mne_val *= dr_val, i_loop++);

    for (i_loop = 0; i_loop < n_dim; i_loop++) {
        iw_arr[i_loop] = 1;
        y_output[i_loop] = 0.0;
    }

    if (key_type == 2) {
        if (std::abs(mne_val) > 1e-12) {
            d_for_key_logic = d_for_key_logic * (1.0 - 1.0 / mne_val);
        } else {
            d_for_key_logic = (x_input == 1.0) ? 1.0 : 0.0;
        }
        k_idx = 0;
    } else if (key_type > 2) {
        double temp_d_for_key_gt2 = x_input;
        dr_val = mne_val / nexp_;
        dr_val = dr_val - std::fmod(dr_val, 1.0);
        dd_val = mne_val - dr_val;
        dr_val = temp_d_for_key_gt2 * dd_val;
        dd_val = dr_val - std::fmod(dr_val, 1.0);
        if (nexp_ - 1 != 0)
            dr_val = dd_val + (dd_val - 1.0) / (nexp_ - 1.0);
        else
            dr_val = dd_val;
        dd_val = dr_val - std::fmod(dr_val, 1.0);
        if (std::abs(mne_val) > 1e-12)
            d_for_key_logic = dd_val * (1.0 / mne_val);
        else
            d_for_key_logic = 0.0;
    }

    double d_iter_for_digits = (key_type == 1) ? x_input : d_for_key_logic;
    double final_d_remainder = 0.0;

    for (j_loop = 0; j_loop < m_order; j_loop++) {
        iq_ = 0;
        if (std::abs(x_input - 1.0) < 1e-9 && key_type != 2) {
            is_val = nexp_ - 1;
            if (j_loop == m_order - 1) final_d_remainder = 0.0;
        } else {
            d_iter_for_digits = d_iter_for_digits * static_cast<double>(nexp_);
            is_val = static_cast<int>(std::floor(d_iter_for_digits));
            if (is_val >= nexp_) is_val = nexp_ - 1;
            if (is_val < 0) is_val = 0;
            d_iter_for_digits = d_iter_for_digits - static_cast<double>(is_val);
            if (j_loop == m_order - 1) final_d_remainder = d_iter_for_digits;
        }

        Node(is_val);

        int temp_i;
        temp_i = iu_[0];
        iu_[0] = iu_[it_val];
        iu_[it_val] = temp_i;
        temp_i = iv_[0];
        iv_[0] = iv_[it_val];
        iv_[it_val] = temp_i;

        if (l_ == 0 && it_val != 0)
            l_ = it_val;
        else if (l_ == it_val)
            l_ = 0;

        if ((iq_ > 0) || ((iq_ == 0) && (is_val == 0)))
            k_idx = l_;
        else if (iq_ < 0)
            k_idx = (it_val == n1_) ? 0 : n1_;

        r_calc = r_calc * 0.5;
        it_val = l_;
        for (i_loop = 0; i_loop < n_dim; i_loop++) {
            iu_[i_loop] = iu_[i_loop] * iw_arr[i_loop];
            iw_arr[i_loop] = -iv_[i_loop] * iw_arr[i_loop];
            p_calc = r_calc * static_cast<double>(iu_[i_loop]);
            p_calc = p_calc + y_output[i_loop];
            y_output[i_loop] = p_calc;
        }
    }

    if (key_type == 2) {
        int sign_val = (is_val == (nexp_ - 1)) ? -1 : 1;
        if (k_idx < 0 || k_idx >= 10) k_idx = 0;
        p_calc = 2.0 * static_cast<double>(sign_val) * static_cast<double>(iu_[k_idx]) * r_calc * final_d_remainder;
        y_output[k_idx] = y_output[k_idx] - p_calc;
    } else if (key_type == 3) {
        for (i_loop = 0; i_loop < n_dim; i_loop++) {
            if (i_loop >= 10) continue;
            p_calc = r_calc * static_cast<double>(iu_[i_loop]);
            p_calc = p_calc + y_output[i_loop];
            y_output[i_loop] = p_calc;
        }
    }
}

std::vector<double> PeanoMapper::Map(double x_param, int m_order, int n_dim, int key_type) {
    if (n_dim <= 0 || n_dim >= 10) {
        std::cerr << "Ошибка Peano: некорректная размерность n_dim = " << n_dim << std::endl;
        return std::vector<double>(n_dim > 0 ? n_dim : 1, 0.5);
    }
    ResetState(n_dim);
    std::vector<double> y_result(n_dim, 0.0);
    Mapd(x_param, m_order, y_result.data(), n_dim, key_type);
    for (int i = 0; i < n_dim; ++i) {
        y_result[i] += 0.5;
    }
    return y_result;
}