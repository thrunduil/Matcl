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

#include "matcl-scalar/details/sequence/seq_transform_base.h"
#include "matcl-scalar/details/sequence/accelerators_helpers.h"

namespace matcl { namespace details
{

template<class Float>
seq_transform_base<Float>::seq_transform_base(int precision)    
{    
    m_table.reserve(50);
    clear(precision);
};

template<class Float>
void seq_transform_base<Float>::clear(int precision)
{
    m_num_calls     = 0;
    m_size          = 0;
    m_precision     = precision;

    m_huge_value    = seq_helpers::max_value<Float>::value(precision) * Float(0.1);
    m_tiny_value    = seq_helpers::min_value<Float>::value(precision) * Float(10.0);
    m_epsilon       = seq_helpers::epsilon<Float>::value(precision);

    m_result        = Float(0);
    m_error_numeric = m_huge_value;
    m_error_limit   = m_huge_value;
    m_best_position = 0;

    init_last_results();
}

template<class Float>
void seq_transform_base<Float>::init_last_results()
{
    m_last_results[0]   = m_huge_value;
    m_last_results[1]   = m_huge_value;
    m_last_results[2]   = m_huge_value;
}

template<class Float>
void seq_transform_base<Float>::update_last_results(const Float& res)
{
    m_last_results[0]   = m_last_results[1];
    m_last_results[1]   = m_last_results[2];
    m_last_results[2]   = res;
};

template<class Float>
template<class Arg_type>
void seq_transform_base<Float>::shift_table(int old_size, std::vector<Arg_type>& table)
{
    m_size              = std::max(m_size, 2);

    if (old_size > m_size)
    {
        Arg_type* table_arr = table.data();

        int index0          = old_size - m_size + 1;

        for (int i = 1, index = index0; i <= m_size; ++i, ++index)
            table_arr[i-1]  = table_arr[index - 1];
    };

    table.resize(m_size);
}

template<class Float>
void seq_transform_base<Float>::update_result_and_error(const Float& lim_err_it, 
                const Float& num_err_it, const Float& res_v, int pos)
{
    // select best approximation; in case of absence of round-off errors this
    // is m_table[0 or 1], however due to round-off errors the best approximation
    // can be somewhere in the middle of the sequence

    if (lim_err_it + num_err_it <= m_error_numeric + m_error_limit) 
    {
        m_error_limit   = lim_err_it;
        m_error_numeric = num_err_it;
        m_result        = res_v;
        m_best_position = pos;
    };            
}

template<class Float>
Float seq_transform_base<Float>::get_lim_error(const Float& res)
{
    Float err0          = matcl::abs(res - m_last_results[0]);
    Float err1          = matcl::abs(res - m_last_results[1]);
    Float err2          = matcl::abs(res - m_last_results[2]);

    Float dif_err       = err0 + err1 + err2;

    return dif_err;
}

}}
