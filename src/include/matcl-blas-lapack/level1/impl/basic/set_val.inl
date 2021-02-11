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

#include "matcl-blas-lapack/level1/level1_basic.h"
#include "matcl-blas-lapack/level1/details/eval_row_functors.h"

namespace matcl { namespace level1
{

template<class T1, Integer Rows>
struct set_val<T1, Rows, false, details::true_t>
{
    force_inline
    static void eval(T1* Y1, Integer rows, const T1& a)
    {
        Integer size =  get_int<Rows>::eval(rows);

        T1* MATCL_RESTRICTED Y  = Y1;

        for (Integer pos = 0; pos < size; ++pos)
            Y[pos] = a;
    };

    force_inline
    static void eval(T1* Y1, Integer Y_step, Integer rows, const T1& a)
    {
        if (Y_step == 1)
            return eval(Y1, rows, a);

        Integer size =  get_int<Rows>::eval(rows);

        T1* MATCL_RESTRICTED Y = Y1;

        for (Integer pos = 0; pos < size; ++pos)
            Y[pos*Y_step]   = a;
    };
};

template<class T, Integer Rows>
struct set_val<T, Rows, true, details::true_t>
{
    force_inline
    static void eval(T* Y1, Integer rows, const T& a)
    {
        using simd_type = typename simd::default_simd_type<T>::type;
        static const int vec_size   = simd_type::vector_size;

        Integer size   =  get_int<Rows>::eval(rows);
        Integer size_1 =  get_int<Rows>::eval(rows) / vec_size;

        T* MATCL_RESTRICTED Y= Y1;

        simd_type sa(a);

		for (Integer i = 0; i < size_1; ++i)
            sa.store(Y+i*vec_size, std::false_type());

        for (Integer i = size_1*vec_size; i < size; ++i)
            Y[i]        = a;
    };

    force_inline
    static void eval(T* Y1, Integer Y_step, Integer rows, const T& a)
    {
        if (Y_step == 1)
            return eval(Y1, rows, a);

        return set_val<T,Rows,false,details::true_t>::eval(Y1, Y_step, rows, a); 
    };
};

template<bool Select, class T1, Integer Rows, Integer Cols, Integer Continuous>
struct set_val_mat<Select, T1, Rows, Cols, Continuous, details::true_t>
{
    static void eval(T1* Y, Integer Y_ld, Integer rows, Integer cols, const T1& a)
    {
        return eval_mat_func_Y<Select, T1, Rows, Cols, Continuous, details::func_set_val<T1>>
                    ::eval(Y, Y_ld, rows, cols, details::func_set_val<T1>(a));
    };
};

}}