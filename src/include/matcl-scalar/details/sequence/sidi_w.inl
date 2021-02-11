/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#pragma once

#include "matcl-scalar/lib_functions/sequence/sidi_w.h"
#include "matcl-scalar/details/sequence/accelerators_helpers.h"

#include <limits>
#include <algorithm>

namespace matcl
{

template<class Float>
sidi_w<Float>::sidi_w(int precision)
    :seq_transform_base(precision)
{
    m_points.reserve(50);
    m_table_den.reserve(50);
    m_table_res.reserve(50);
};

template<class Float>
void sidi_w<Float>::clear(int precision)
{
    seq_transform_base::clear(precision);
}

template<class Float>
void sidi_w<Float>::eval(const Float& x, const Float& omega, const Float& y, 
                            Float& res, Float& abs_err)
{
    ++m_size;

    Float x_m           = seq_helpers::prepare_value<Float>::eval(x, m_precision);
    Float y_m           = seq_helpers::prepare_value<Float>::eval(y, m_precision);
    Float om_m          = seq_helpers::prepare_value<Float>::eval(omega, m_precision);
        
    m_table.resize(m_size);
    m_table_den.resize(m_size);
    m_table_res.resize(m_size);
    m_points.resize(m_size);

    if (matcl::abs(om_m) <= m_tiny_value)
    {
        details::seq_error_omega_is_zero();
    };

    // initialize table with y_m/omega and table_den with 1/omega
    Float om_inv        = Float(1) / om_m;
    Float om_inv_abs    = matcl::abs(om_inv);
    Float om_inv_e      = Float(0.5) * m_epsilon * om_inv_abs;

    y_m                 = y_m * om_inv;

    // 1 ulp in y + 0.5 ulp in om_inv + 0.5 ulp in y_m * om_inv
    Float y_e           = Float(2.0) * matcl::abs(y_m) * m_epsilon;    

    m_table[m_size -1]      = val_err(y_m, y_e);
    m_table_den[m_size -1]  = val_err(om_inv, om_inv_e);
    m_points[m_size-1]      = x_m;

    make();

    res                 = m_result;
    abs_err             = m_error_limit + m_error_numeric;
}

template<class Float>
void sidi_w<Float>::eval(const Float& x, const Float& omega, const Float& omega_err0, 
                        const Float& y, const Float& y_err0, Float& res, Float& abs_err)
{
    Float y_err         = y_err0;
    Float omega_err     = omega_err0;

    if (y_err <= Float(0))
        y_err           = m_epsilon * matcl::abs(y);

    if (omega_err <= Float(0))
        omega_err       = m_epsilon * matcl::abs(omega);

    ++m_size;

    Float x_m           = seq_helpers::prepare_value<Float>::eval(x, m_precision);
    Float y_m           = seq_helpers::prepare_value<Float>::eval(y, m_precision);
    Float y_e           = seq_helpers::prepare_value<Float>::eval(y_err, m_precision);
    Float om_m          = seq_helpers::prepare_value<Float>::eval(omega, m_precision);
    Float om_e          = seq_helpers::prepare_value<Float>::eval(omega_err, m_precision);

    m_table.resize(m_size);
    m_table_den.resize(m_size);
    m_points.resize(m_size);
    m_table_res.resize(m_size);

    if (matcl::abs(om_m) <= m_tiny_value)
    {
        details::seq_error_omega_is_zero();
    };

    // initialize table with y_m/omega and table_den with 1/omega
    Float om_inv        = Float(1) / om_m;
    Float om_inv_abs    = matcl::abs(om_inv);
    Float om_inv_e      = om_inv_abs * om_inv_abs * om_e;

    y_m                 = y_m * om_inv;
    y_e                 = om_inv_abs * y_e + matcl::abs(y_m) * m_epsilon;

    m_table[m_size - 1]     = val_err(y_m, y_e);
    m_table_den[m_size -1]  = val_err(om_inv, om_inv_e);
    m_points[m_size- 1]     = x_m;

    make();

    res                 = m_result;
    abs_err             = m_error_limit + m_error_numeric;
}

template<class Float>
void sidi_w<Float>::make()
{
    using seq_helpers::sqr;

    using const_ref     = const Float&;

    Float zero          = Float(0);
    Float one           = Float(1);
    Float three_half    = Float(1.5);
    const int min_size  = 5;

    m_num_calls         = m_num_calls+1;

    val_err* table_nom  = m_table.data();
    val_err* table_den  = m_table_den.data();
    Float* point_arr    = m_points.data();
    val_err* table_res  = m_table_res.data();

    {
        Float res_n     = table_nom[m_size - 1].first;
        Float res_d     = table_den[m_size - 1].first;
        Float res_ne    = table_nom[m_size - 1].second;
        Float res_de    = table_den[m_size - 1].second;

        div_nom_den(res_n, res_ne, res_d, res_de, m_result, m_error_numeric);

        table_res[m_size - 1] = val_err(m_result, m_error_numeric);
    };

    if(m_size < 2)
    {
        m_error_limit   = m_huge_value;
        m_best_position = m_size;   

        update_last_results(m_result);
        return;
    }

    m_error_limit       = m_huge_value;
    m_error_numeric     = m_huge_value;
    m_best_position     = 0;    

    int old_n           = m_size;

    // arbitrary parameters
    const int num_keep  = 4;    
    const int min_it    = 2;

    for(int j = 1; j < m_size; ++j)
    {   
        Float M1_v          = table_nom[m_size - j].first;
        Float M2_v          = table_nom[m_size - j - 1].first;
        Float M1_e          = table_nom[m_size - j].second;
        Float M2_e          = table_nom[m_size - j - 1].second;

        Float N1_v          = table_den[m_size - j].first;
        Float N2_v          = table_den[m_size - j - 1].first;
        Float N1_e          = table_den[m_size - j].second;
        Float N2_e          = table_den[m_size - j - 1].second;

        Float MN1_v         = table_res[m_size - j].first;
        Float MN2_v         = table_res[m_size - j - 1].first;

        Float x1_v          = point_arr[m_size - j - 1];
        Float x2_v          = point_arr[m_size - 1];

        Float M_num         = M1_v - M2_v;
        Float N_num         = N1_v - N2_v;

        Float den           = x1_v - x2_v;
        Float den_abs       = matcl::abs(den);

        // if den is close to zero, omit a part of the table by adjusting
        // the value of m_size

        if (den_abs <= m_tiny_value)
        {
            details::seq_error_two_equal_points
                        (seq_helpers::cast_double<Float>::eval(x1_v));

            m_size          = j;
            break;
        };

        Float den_inv       = one / den;
        Float den_inv_abs   = matcl::abs(den_inv);

        Float M_res_v       = M_num * den_inv;
        Float M_res_abs     = matcl::abs(M_res_v);

        Float N_res_v       = N_num * den_inv;
        Float N_res_abs     = matcl::abs(N_res_v);

        // compute error        
        Float M_res_e       = zero;
        Float N_res_e       = zero;
        {
            // contribution of errors in M1, M2, N1, N2
            Float err_M     = den_inv_abs * (M1_e + M2_e);
            Float err_N     = den_inv_abs * (N1_e + N2_e);

            // total numerical error
            M_res_e         = err_M + M_res_abs * m_epsilon * three_half;
            N_res_e         = err_N + N_res_abs * m_epsilon * three_half;
        }

        table_nom[m_size-j-1]= val_err(M_res_v, M_res_e);
        table_den[m_size-j-1]= val_err(N_res_v, N_res_e);

        Float res_v, res_e;
        div_nom_den(M_res_v, M_res_e, N_res_v, N_res_e, res_v, res_e);

        table_res[m_size-j-1]= val_err(res_v, res_e);

        Float abs_d0        = matcl::abs(MN1_v - MN2_v);
        Float abs_cor       = matcl::abs(res_v - MN2_v);
        Float lim_error     = abs_d0 + abs_cor;

        update_result_and_error(lim_error, res_e, res_v, j);
    };

    // compute error estimate

    m_error_limit       = m_error_limit + get_lim_error(m_result);

    m_size              = std::max(m_size, min_size);
    m_size              = std::min(m_size, old_n);

    // if approximations at low positions are not good, omit a part of the
    // table by adjusting the value of m_size 
    
    if (m_best_position + 1 + num_keep < m_size)
    {
        // keep num_keep additional approximations
        int n2          = m_best_position + 1 + num_keep;
        m_size          = n2;
    };

    shift_table(old_n, m_table);
    shift_table(old_n, m_table_den);
    shift_table(old_n, m_table_res);
    shift_table(old_n, m_points);

    update_last_results(m_result);
    return;
};

template<class Float>
void sidi_w<Float>::div_nom_den(const Float& res_n, const Float& res_ne, 
                    const Float& res_d, const Float& res_de, Float& result, 
                    Float& error)
{
    Float den_inv       = Float(1) / res_d;
    Float den_inv_abs   = matcl::abs(den_inv);

    result              = res_n / res_d;
    Float res_abs       = matcl::abs(result);

    error               = den_inv_abs * res_ne + res_abs * den_inv_abs * res_de
                        + res_abs * m_epsilon * Float(0.5);
}

template<class Float>
const Float& sidi_w<Float>::last_result() const
{
    return m_result;
}

template<class Float>
const Float& sidi_w<Float>::numerical_error() const
{
    return m_error_numeric;
}

template<class Float>
const Float& sidi_w<Float>::limit_error() const
{
    return m_error_limit;
}

template<class Float>
Float sidi_w<Float>::total_error() const
{
    return numerical_error() + limit_error();
}

}
