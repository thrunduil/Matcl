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

#include "matcl-core/config.h"
#include "matcl-scalar/config.h"

#include <utility>
#include <limits>
#include <vector>

namespace matcl { namespace details
{
    
template<class Float>
class seq_transform_base
{
    protected:
        using val_err           = std::pair<Float, Float>;
        using vector_estim      = std::vector<val_err>;
        using vector_points     = std::vector<Float>;

    protected:
        //  number of calls to the routine
        int             m_num_calls;

        // vector containing the elements of the two lower diagonals of the
        // condensated triangular epsilon table and interpolation points
        vector_estim    m_table;

        // vector containing the last 3 results
        Float           m_last_results[3];

        // floating point type properties
        Float           m_huge_value;
        Float           m_tiny_value;
        Float           m_epsilon;

        // m_table[m_size - 1] contains the new element in the first column 
        // of the epsilon table
        int             m_size;

        // working precision; has effect only when Float is a multiprecision
        // type
        int             m_precision;

        // computed limit
        Float           m_result;

        // estimated errors
        Float           m_error_numeric;
        Float           m_error_limit;
        int             m_best_position;

    protected:
        seq_transform_base(int precision);

    protected:
        void            init_last_results();
        void            update_last_results(const Float& res);
        void            clear(int precision);

        template<class Arg_type>
        void            shift_table(int old_size, std::vector<Arg_type>& table);

        Float           get_lim_error(const Float& res);
        void            update_result_and_error(const Float& lim_err_it, 
                            const Float& num_err_it, const Float& res_v, int pos);
};

void MATCL_SCALAR_EXPORT seq_error_omega_is_zero();
void MATCL_SCALAR_EXPORT seq_error_two_equal_points(double x);

}}

#include "matcl-scalar/details/sequence/seq_transform_base.inl"
