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

#include "matcl-scalar/lib_functions/sequence/wynn_epsilon.h"
#include "matcl-scalar/details/sequence/accelerators_helpers.h"

#include <limits>
#include <algorithm>

namespace matcl
{

template<class Float>
wynn_epsilon<Float>::wynn_epsilon(int precision)
    :seq_transform_base(precision)
{};

template<class Float>
void wynn_epsilon<Float>::clear(int precision)
{
    seq_transform_base::clear(precision);
}

template<class Float>
void wynn_epsilon<Float>::eval(const Float& elem, Float& res, 
                                    Float& abs_err)
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
void wynn_epsilon<Float>::eval(const Float& elem, const Float& elem_err, Float& res, Float& abs_err)
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
void wynn_epsilon<Float>::make()
{
    using seq_helpers::sqr;

    using const_ref     = const Float&;

    const Float zero    = Float(0);
    const Float one     = Float(1);
    const Float half    = Float(0.5);
    const int min_size  = 5;

    m_num_calls         = m_num_calls+1;
    val_err* table_arr  = m_table.data();

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
    const int min_it    = 3;

    val_err aux2        = val_err(zero, zero);
    bool approx_pos     = false;

    for(int j = m_size, num_it = 1; j >= 2; --j, ++num_it)
    {   
        Float E0_v      = aux2.first;
        Float E0_e      = aux2.second;

        Float E1_v      = table_arr[j - 1].first;
        Float E2_v      = table_arr[j - 2].first;

        Float E1_e      = table_arr[j - 1].second;
        Float E2_e      = table_arr[j - 2].second;        

        aux2            = table_arr[j - 2];

        Float delta     = E1_v - E2_v;
        Float absd      = matcl::abs(delta);

        // if two elements are very close to each other, omit
        // a part of the table by adjusting the value of m_size

        if (absd <= m_tiny_value)
        {
            if (approx_pos == false)
            {
                // in this case E1 and E2 are good limit estimations,
                // try use them as the final estimation

                Float res_e     = E1_e;
                Float res_v     = E1_v;
                Float lim_err   = absd;                      

                if (num_it > 1)
                {
                    const_ref E3_v  = table_arr[j + 1].first;
                    Float abs_d31   = matcl::abs(E1_v - E3_v);
                    lim_err         = lim_err + abs_d31;
                };

                update_result_and_error(lim_err, res_e, res_v, j);
            }
            else
            {
                // try use E0 as the final estimation

                Float res_e     = E0_e;
                Float res_v     = E0_v;

                const_ref E3_v  = table_arr[j].first;
                Float abs_d30   = matcl::abs(E0_v - E3_v);
                Float lim_err   = abs_d30;

                if (num_it > 2)
                {
                    const_ref E4_v  = table_arr[j + 2].first;
                    Float abs_d40   = matcl::abs(E0_v - E4_v);
                    lim_err         = lim_err + abs_d40;
                };

                update_result_and_error(lim_err, res_e, res_v, j);
            }

            m_size      = m_size - j + 1;
            break;
        };

        // compute a new element and eventually adjust
        // the value of result.        

        Float del_inv   = one / delta;
        Float dinv_abs  = matcl::abs(del_inv);

        Float res_v     = E0_v + del_inv;
        Float res_abs   = matcl::abs(res_v);

        // compute error        

        Float res_e;
        {   
            Float E12_err   = E1_e + E2_e;
            Float err_E;

            Float rel_err   = E12_err *  dinv_abs;

            if (rel_err <= Float(0.00001))
            {
                // use first order approximation
                err_E       = rel_err *  dinv_abs;
            }
            else if (rel_err <= Float(0.99999))
            {
                // exact error
                err_E       = one / (absd - E12_err) - dinv_abs;
            }
            else
            {
                // in this case algorithm cannot proceed further, E1 and E2 are
                // numerically equal

                if (approx_pos == false)
                {
                    // in this case E1 and E2 are good limit estimations,
                    // try use them as the final estimation

                    res_e           = E1_e;
                    res_v           = E1_v;

                    Float lim_err   = absd;

                    if (num_it > 1)
                    {
                        const_ref E3_v  = table_arr[j + 1].first;
                        Float abs_d31   = matcl::abs(E1_v - E3_v);
                        lim_err         = lim_err + abs_d31;
                    };

                    update_result_and_error(lim_err, res_e, res_v, j);
                }
                else
                {
                    // try use E0 as the final estimation

                    res_e           = E0_e;
                    res_v           = E0_v;

                    const_ref E3_v  = table_arr[j].first;
                    Float abs_d30   = matcl::abs(E0_v - E3_v);
                    Float lim_err   = abs_d30;

                    if (num_it > 2)
                    {
                        const_ref E4_v  = table_arr[j + 2].first;
                        Float abs_d40   = matcl::abs(E0_v - E4_v);
                        lim_err         = lim_err + abs_d40;
                    };

                    update_result_and_error(lim_err, res_e, res_v, j);
                }

                m_size      = m_size - j + 1;
                break;
            }

            err_E           = err_E + E0_e;

            Float del_inv2  = del_inv * del_inv;

            //error in computing delta
            Float err_d     = seq_helpers::error_minus(E1_v, E2_v, absd) * del_inv2; 

            //error in computing del_inv
            Float err_dinv  = dinv_abs;

            //error is computing res_v
            Float err_res   = res_abs;

            // total error
            res_e           = err_E
                            + (err_d + err_dinv + err_res) * m_epsilon * half;
        };

        table_arr[j - 2]    = val_err(res_v, res_e);

        if (approx_pos == true)
        {
            Float E3_v      = table_arr[j].first;
            Float abs_d0    = matcl::abs(E0_v - E3_v);
            Float abs_cor   = dinv_abs;
            Float lim_err   = abs_d0 + abs_cor;

            update_result_and_error(lim_err, res_e, res_v, j);
        };

        approx_pos          = !approx_pos;
    };

    // compute error estimate

    m_error_limit       = m_error_limit + get_lim_error(m_result);

    m_size              = std::max(m_size, min_size);
    m_size              = std::min(m_size, old_n);

    // if approximations at low positions are not good, omit a part of the
    // table by adjusting the value of m_size 

    if (m_best_position > num_keep + 1)
    {
        // keep num_keep additional approximations
        int n2          = m_size - m_best_position + 1  + num_keep;
        m_size          = std::min(n2, m_size);
    };

    shift_table(old_n, m_table);

    update_last_results(m_result);
    return;
};

template<class Float>
const Float& wynn_epsilon<Float>::last_result() const
{
    return m_result;
}

template<class Float>
const Float& wynn_epsilon<Float>::numerical_error() const
{
    return m_error_numeric;
}

template<class Float>
const Float& wynn_epsilon<Float>::limit_error() const
{
    return m_error_limit;
}

template<class Float>
Float wynn_epsilon<Float>::total_error() const
{
    return numerical_error() + limit_error();
}

}
