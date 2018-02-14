/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017-2018
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

#include "matcl-scalar/lib_functions/sequence/aitken_delta.h"
#include "matcl-scalar/details/sequence/accelerators_helpers.h"

#include <limits>
#include <algorithm>

namespace matcl
{

template<class Float>
aitken_delta<Float>::aitken_delta(int precision)
    :seq_transform_base(precision)
{};

template<class Float>
void aitken_delta<Float>::clear(int precision)
{
    seq_transform_base::clear(precision);
}

template<class Float>
void aitken_delta<Float>::eval(const Float& elem, Float& res, Float& abs_err)
{
    ++m_size;

    m_table.resize(m_size);

    Float elem_m        = seq_helpers::prepare_value<Float>::eval(elem, m_precision);
    m_table[m_size - 1] = val_err(elem_m, m_epsilon * matcl::abs(elem_m));

    make();

    res                 = m_result;
    abs_err             = m_error_limit + m_error_numeric;
}

template<class Float>
void aitken_delta<Float>::eval(const Float& elem, const Float& elem_err, Float& res, Float& abs_err)
{
    if (elem_err <= Float(0))
        return eval(elem, res, abs_err);

    ++m_size;

    m_table.resize(m_size);

    Float elem_m        = seq_helpers::prepare_value<Float>::eval(elem, m_precision);
    Float elem_e        = seq_helpers::prepare_value<Float>::eval(elem_err, m_precision);
    m_table[m_size - 1] = val_err(elem_m, elem_e);

    make();

    res                 = m_result;
    abs_err             = m_error_limit + m_error_numeric;
}

template<class Float>
void aitken_delta<Float>::make()
{
    using seq_helpers::sqr;

    using const_ref = const Float&;

    const Float half    = Float(0.5);
    const Float two     = Float(2);
    const Float three   = Float(3);
    const Float one     = Float(1);
    const int min_size  = 5;

    m_num_calls         = m_num_calls+1;
    val_err* table_arr  = m_table.data();

    if(m_size < 3)
    {
        m_error_limit   = m_huge_value;
        m_error_numeric = table_arr[m_size - 1].second;
        m_result        = table_arr[m_size - 1].first;
        m_best_position = 1;   

        update_last_results(m_result);
        return;
    }

    m_error_limit       = m_huge_value;
    m_error_numeric     = m_huge_value;
    m_best_position     = 0;   
    int old_n           = m_size;

    // arbitrary parameters
    const int num_keep      = 4;    
    const int min_it        = 2;

    int iter_limit      = (m_size - 1) / 2;

    for(int i = 1; i <= iter_limit; ++i)
    {     
        int M           = m_size - 2 * i;

        const_ref A2_v  = table_arr[M + 2 - 1].first;
        const_ref A1_v  = table_arr[M + 1 - 1].first;
        const_ref A0_v  = table_arr[M - 1].first;

        Float delta_21  = A2_v - A1_v;
        Float delta_01  = A0_v - A1_v;
        Float delta_02  = A0_v - A2_v;

        Float abs_d21   = matcl::abs(delta_21);
        Float abs_d01   = matcl::abs(delta_01);
        Float abs_d02   = matcl::abs(delta_02);
        
        Float den       = delta_21 + delta_01;
        Float abs_den   = matcl::abs(den);

        // if two elements are very close to each other, omit
        // a part of the table by adjusting the value of m_size

        if (abs_den <= m_tiny_value)
        {
            // in this case A0, A1, A2 are good limit estimations,
            // try use them as the final estimation

            Float res_e     = table_arr[M - 1].second;
            Float res_v     = A0_v;
            Float lim_error = abs_d21 + abs_d01 + abs_d02;

            if (i > 1)
            {
                const_ref A3_v  = table_arr[M + 3 - 1].first;
                Float abs_d32   = matcl::abs(A3_v - A2_v);
                lim_error       = lim_error + abs_d32;
            };

            update_result_and_error(lim_error, res_e, res_v, i);

            m_size           = 2 * i;
            break;
        };

        Float num       = delta_01;
        Float den_inv   = one / den;
        Float num_den1  = num * den_inv;
        Float cor       = num * num_den1;

        Float abs_den_inv   = matcl::abs(den_inv);
        Float abs_cor       = matcl::abs(cor);

        const_ref A2_e  = table_arr[M + 2 - 1].second;
        const_ref A1_e  = table_arr[M + 1 - 1].second;
        const_ref A0_e  = table_arr[M - 1].second;            
        
        // compute a new element and eventually adjust the value of 
        // result

        Float res_v         = A0_v - cor;
        Float res_abs       = matcl::abs(res_v);

        // compute error        
    
        Float res_e;
        {
            Float num_den2  = num_den1 * num_den1;

            // contribution of errors in A0, A1, A2
            Float erprop_A0, erprop_A1, erprop_A2;
    
            // (A1 - A2)^2 / (A0 - 2A1 + A2)^2
            erprop_A0       = sqr(delta_21 * den_inv);

            // 1/2 - 1/2 * (A0 - A2)^2/(A0 - 2A1 + A2)^2
            erprop_A1       = matcl::abs(half - half * sqr(delta_02 * den_inv));

            // (A0 - A1)^2 / (A0 - 2A1 + A2)^2;
            erprop_A2       = num_den2;

            Float err_A     = erprop_A0 * A0_e + erprop_A1 * A1_e + erprop_A2 * A2_e;

            Float err_A012  = A0_e + two * A1_e + A2_e;
            Float err_rel   = err_A012 * abs_den_inv;

            if (err_rel > Float(0.01))
            {                
                if (err_rel <= Float(0.99999999))
                {
                    Float err_den   = one / (abs_den - err_A012) - abs_den_inv;
                    err_A           = err_A + (err_den * abs_d01) * abs_d01;
                }
                else
                {
                    // in this case algorithm cannot proceed further, A0 - A2 are
                    // numerically equal

                    res_e           = A0_e;
                    res_v           = A0_v;

                    Float lim_error = abs_d21 + abs_d01 + abs_d02;

                    if (i > 1)
                    {
                        const_ref A3_v  = table_arr[M + 3 - 1].first;
                        Float abs_d32   = matcl::abs(A3_v - A2_v);
                        lim_error       = lim_error + abs_d32;
                    };

                    update_result_and_error(lim_error, res_e, res_v, i);

                    m_size           = 2 * i;
                    break;
                }
            }

            //error in computing den
            Float err_D2    = seq_helpers::error_minus(A2_v, A1_v, abs_d21);
            Float err_D0    = seq_helpers::error_minus(A0_v, A1_v, abs_d01);
            Float err_den   = seq_helpers::error_plus(delta_21, delta_01, abs_den) + err_D0 + err_D2;

            //error in computing num
            Float err_num   = err_D0;

            Float erprop_num, erprop_den;

            erprop_num      = matcl::abs(num_den1) * two;
            erprop_den      = num_den2;

            // contribution of errors in num, den
            Float err_nd    = erprop_num * err_num + erprop_den * err_den;            

            // total error
            res_e           = err_A                                 // inaccuracy of A0, A1, A2
                            + (err_nd                               // error in  num, den
                                + abs_cor * three                   // error in num * (num * 1/den)
                                + res_abs                           // error in (A0) + (cor)
                              ) * m_epsilon * half;                 // floating point errors in
                                                                    // u ( = eps/2)
        };

        table_arr[M - 1]    = val_err(res_v, res_e);
        
        // select best approximation in case of absence of round-off errors
        // m_table[0 or 1], however due to round-off errors the best approximation
        // can be somewhere in the middle of the sequence

        Float lim_err       = abs_d01 + abs_d21 + abs_d02 + abs_cor;
        update_result_and_error(lim_err, res_e, res_v, i);
    };

    // compute error estimate

    m_error_limit       = m_error_limit + get_lim_error(m_result);

    m_size              = std::max(m_size, min_size);
    m_size              = std::min(m_size, old_n);

    // if approximations at low positions are not good, omit a part of the
    // table by adjusting the value of m_size 

    if (m_best_position + num_keep < iter_limit)
    {
        // keep num_keep additional approximations
        int i           = m_best_position + num_keep;
        m_size          = std::min(2 * i - 1, m_size);
    };

    shift_table(old_n, m_table);

    update_last_results(m_result);
    return;
};

template<class Float>
const Float& aitken_delta<Float>::last_result() const
{
    return m_result;
}

template<class Float>
const Float& aitken_delta<Float>::numerical_error() const
{
    return m_error_numeric;
}

template<class Float>
const Float& aitken_delta<Float>::limit_error() const
{
    return m_error_limit;
}

template<class Float>
Float aitken_delta<Float>::total_error() const
{
    return numerical_error() + limit_error();
}

}
