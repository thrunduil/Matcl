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

namespace matcl { namespace level1
{

template<class T1, Integer Rows>
struct swap<T1, Rows, false, details::true_t>
{
    force_inline
    static void eval(T1* Y1, T1* X1, Integer rows)
    {
        Integer size =  get_int<Rows>::eval(rows);

        T1* MATCL_RESTRICTED Y    = Y1;
        T1* MATCL_RESTRICTED X    = X1;

        for (Integer pos = 0; pos < size; ++pos)
            std::swap(Y[pos], X[pos]);
    };

    force_inline
    static void eval_step1(T1* Y1, Integer Y_step, T1* X1, Integer rows)
    {
        Integer size =  get_int<Rows>::eval(rows);

        T1* MATCL_RESTRICTED Y    = Y1;
        T1* MATCL_RESTRICTED X    = X1;

        for (Integer pos = 0; pos < size; ++pos)
            std::swap(Y[pos*Y_step], X[pos]);
    };

    force_inline
    static void eval_step2(T1* Y1, T1* X1, Integer X_step, Integer rows)
    {
        Integer size =  get_int<Rows>::eval(rows);

        T1* MATCL_RESTRICTED Y    = Y1;
        T1* MATCL_RESTRICTED X    = X1;

        for (Integer pos = 0; pos < size; ++pos)
            std::swap(Y[pos], X[pos*X_step]);
    };

    force_inline
    static void eval(T1* Y1, Integer Y_step, T1* X1, Integer X_step, Integer rows)
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

        T1* MATCL_RESTRICTED Y    = Y1;
        T1* MATCL_RESTRICTED X    = X1;

        for (Integer pos = 0; pos < size; ++pos)
            std::swap(Y[pos*Y_step], X[pos*X_step]);
    };
};

template<class T, Integer Rows>
struct swap<T, Rows, true, details::true_t>
{
    force_inline
    static void eval(T* Y1, T* X1, Integer rows)
    {
        using simd_type = typename simd::default_simd_type<T>::type;
        static const int vec_size   = simd_type::vector_size;

        Integer size   =  get_int<Rows>::eval(rows);
        Integer size_1 =  get_int<Rows>::eval(rows) / vec_size;

        T* MATCL_RESTRICTED Y    = Y1;
        T* MATCL_RESTRICTED X    = X1;

		for (Integer i = 0; i < size_1; ++i)
        {
            simd_type sX_1  = simd_type::load(X+i*vec_size, std::false_type());
            simd_type sY_1  = simd_type::load(Y+i*vec_size, std::false_type());

            sX_1.store(Y+i*vec_size, std::false_type());
            sY_1.store(X+i*vec_size, std::false_type());
        };

        for (Integer i = size_1*vec_size; i < size; ++i)
            std::swap(X[i], Y[i]);
    };

    force_inline
    static void eval(T* Y1, Integer Y_step, T* X1, Integer X_step, Integer rows)
    {
        if (Y_step == 1 && X_step == 1)
            return eval(Y1, X1, rows);

        return swap<T,Rows,false,details::true_t>::eval(Y1, Y_step, X1, X_step, rows);
    };
};

}}