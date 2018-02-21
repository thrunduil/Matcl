/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-blas/level1/level1_other.h"

namespace matcl { namespace level1
{

template<class T1, Integer Rows>
struct dot<T1, Rows, false, details::true_t>
{
    force_inline
    static T1 eval(const T1* Y1, const T1* X1, Integer rows)
    {
        Integer size =  get_int<Rows>::eval(rows);

        const T1* MATCL_RESTRICTED Y  = Y1;
        const T1* MATCL_RESTRICTED X  = X1;

        T1 sum  = T1();

		for (Integer i = 0; i < size; ++i)
			sum = sum + Y[i] * X[i];

        return sum;
    };

    force_inline
    static T1 eval_step1(const T1* Y1, Integer Y_step, const T1* X1, Integer rows)
    {
        Integer size =  get_int<Rows>::eval(rows);

        const T1* MATCL_RESTRICTED Y  = Y1;
        const T1* MATCL_RESTRICTED X  = X1;

        T1 sum  = T1();

		for (Integer i = 0; i < size; ++i)
			sum = sum + Y[i*Y_step] * X[i];

        return sum;
    };

    force_inline
    static T1 eval_step2(const T1* Y1, const T1* X1, Integer X_step, Integer rows)
    {
        Integer size =  get_int<Rows>::eval(rows);

        const T1* MATCL_RESTRICTED Y  = Y1;
        const T1* MATCL_RESTRICTED X  = X1;

        T1 sum  = T1();

		for (Integer i = 0; i < size; ++i)
			sum = sum + Y[i] * X[i*X_step];

        return sum;
    };

    force_inline
    static T1 eval(const T1* Y1, Integer Y_step, const T1* X1, Integer X_step, Integer rows)
    {
        if (Y_step == 1)
        {
            if (X_step == 1)
                return eval(Y1, X1, rows);
            else
                return eval_step2(Y1, X1, X_step, rows);
        }
        else if (X_step == 1)                
        {
            return eval_step1(Y1, Y_step, X1, rows);
        };

        Integer size =  get_int<Rows>::eval(rows);

        const T1* MATCL_RESTRICTED Y  = Y1;
        const T1* MATCL_RESTRICTED X  = X1;

        T1 sum  = T1();

		for (Integer i = 0; i < size; ++i)
			sum = sum + Y[i*Y_step] * X[i*X_step];

        return sum;
    };
};

template<class T, Integer Rows>
struct dot<T, Rows, true, details::true_t>
{
    force_inline
    static T eval(const T* Y1, const T* X1, Integer rows)
    {
        using simd_type = typename simd::default_simd_type<T>::type;
        static const int vec_size   = simd_type::vector_size;

        Integer size   =  get_int<Rows>::eval(rows);
        Integer size_1 =  get_int<Rows>::eval(rows) / vec_size;

        const T* MATCL_RESTRICTED Y = Y1;
        const T* MATCL_RESTRICTED X = X1;

        simd_type res_1 = simd_type::zero();

		for (Integer i = 0; i < size_1; ++i)
        {
            simd_type sX_1  = simd_type::load(X+i*vec_size, std::false_type());
            simd_type sY_1  = simd_type::load(Y+i*vec_size, std::false_type());
			res_1           = fma_f(sX_1,sY_1, res_1);
        };

        T res               = horizontal_sum(res_1);
		
        for (Integer i = size_1*vec_size; i < size; ++i)
            res             += X[i] * Y[i];

        return res;
    };

    force_inline
    static T eval(const T* Y1, Integer Y_step, const T* X1, Integer X_step, Integer rows)
    {
        if (Y_step == 1 && X_step == 1)
            return eval(Y1, X1, rows);

        return dot<T,Rows,false,details::true_t>::eval(Y1, Y_step, X1, X_step, rows);
    };
};

}}