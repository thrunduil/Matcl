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

#include "matcl-scalar/lib_functions/sequence/richardson.h"
#include "matcl-scalar/details/sequence/accelerators_helpers.h"

#include <limits>
#include <algorithm>

namespace matcl
{

template<class Float>
richardson<Float>::richardson(int precision)
    :seq_transform_base(precision)
{
    m_points.reserve(50);
};

template<class Float>
void richardson<Float>::clear(int precision)
{
    seq_transform_base::clear(precision);
}

template<class Float>
void richardson<Float>::eval(const Float& x, const Float& y, 
                            Float& res, Float& abs_err)
{
    ++m_size;

    Float x_m           = seq_helpers::prepare_value<Float>::eval(x, m_precision);
    Float y_m           = seq_helpers::prepare_value<Float>::eval(y, m_precision);
    
    m_table.resize(m_size);
    m_points.resize(m_size);

    m_table[m_size -1]  = val_err(y_m, m_epsilon * matcl::abs(y_m));
    m_points[m_size-1]  = x_m;

    make();

    res                 = m_result;
    abs_err             = m_error_limit + m_error_numeric;
}

template<class Float>
void richardson<Float>::eval(const Float& x, const Float& y, const Float& y_err0, 
                            Float& res, Float& abs_err)
{
    Float y_err         = y_err0;

    if (y_err <= Float(0))
        y_err           = m_epsilon * matcl::abs(y);

    ++m_size;

    Float x_m           = seq_helpers::prepare_value<Float>::eval(x, m_precision);
    Float y_m           = seq_helpers::prepare_value<Float>::eval(y, m_precision);
    Float y_e           = seq_helpers::prepare_value<Float>::eval(y_err, m_precision);

    m_table.resize(m_size);
    m_points.resize(m_size);

    m_table[m_size - 1] = val_err(y_m, y_e);
    m_points[m_size- 1] = x_m;

    make();

    res                 = m_result;
    abs_err             = m_error_limit + m_error_numeric;
}

template<class Float>
void richardson<Float>::make()
{
    using seq_helpers::sqr;

    using const_ref     = const Float&;

    Float zero          = Float(0);
    Float one           = Float(1);
    Float half          = Float(0.5);
    Float three         = Float(3.0);
    const int min_size  = 5;

    m_num_calls         = m_num_calls+1;

    val_err* table_arr  = m_table.data();
    Float* point_arr    = m_points.data();

    if(m_size < 2)
    {
        m_error_limit   = m_huge_value;
        m_error_numeric = table_arr[m_size - 1].second;
        m_result        = table_arr[m_size - 1].first;
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
        Float N1_v          = table_arr[m_size - j].first;
        Float N2_v          = table_arr[m_size - j - 1].first;
        Float N1_e          = table_arr[m_size - j].second;
        Float N2_e          = table_arr[m_size - j - 1].second;

        Float x1_v          = point_arr[m_size - j - 1];
        Float x2_v          = point_arr[m_size - 1];

        Float num1          = x1_v * N1_v;
        Float num2          = x2_v * N2_v;

        Float num1_abs      = matcl::abs(num1);
        Float num2_abs      = matcl::abs(num2);

        Float num           = num1 - num2;
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

        Float res_v         = num * den_inv;
        Float res_abs       = matcl::abs(res_v);        

        // compute error        
        Float res_e         = zero;
        {
            // contribution of errors in N1, N2
            Float x1_den    = matcl::abs(x1_v * den_inv);
            Float x2_den    = matcl::abs(x2_v * den_inv);
            Float err_N     = x1_den * N1_e + x2_den * N2_e;

            // error in computing x1 * N1, x2 * N2
            Float err_xN    = (num1_abs + num2_abs) * den_inv_abs;

            // error in computing num, den and num/den
            Float err_numden= res_abs;

            // total numerical error
            res_e           = err_N
                            + (err_xN + err_numden * three) * m_epsilon * half;
        }

        table_arr[m_size-j-1] = val_err(res_v, res_e);

        Float abs_d0        = matcl::abs(N1_v - N2_v);
        Float abs_cor       = matcl::abs(res_v - N2_v);
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
    shift_table(old_n, m_points);

    update_last_results(m_result);
    return;
};

template<class Float>
const Float& richardson<Float>::last_result() const
{
    return m_result;
}

template<class Float>
const Float& richardson<Float>::numerical_error() const
{
    return m_error_numeric;
}

template<class Float>
const Float& richardson<Float>::limit_error() const
{
    return m_error_limit;
}

template<class Float>
Float richardson<Float>::total_error() const
{
    return numerical_error() + limit_error();
}

}
