/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2018
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

#include "matcl-blas/level1/level1_axpby.h"

#pragma warning(push)
#pragma warning(disable: 4127)  //conditional expression is constant

namespace matcl { namespace level1
{

template<class TY, class TX, Integer Rows>
struct ypxm<TY, TX, Rows, false, details::true_t>
{
    force_inline
    static void eval(TY* Y1, const TX* X1, Integer rows)
    {
        Integer size =  get_int<Rows>::eval(rows);

        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

		for (Integer i = 0; i < size; ++i)
			Y[i] = -(Y[i] + X[i]);
    };

    force_inline
    static void eval_step1(TY* Y1, Integer Y_step, const TX* X1, Integer rows)
    {
        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        Integer size =  get_int<Rows>::eval(rows);

		for (Integer i = 0; i < size; ++i)
			Y[i*Y_step] = -(Y[i*Y_step] + X[i]);
    };

    force_inline
    static void eval_step2(TY* Y1, const TX* X1, Integer X_step, Integer rows)
    {
        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        Integer size =  get_int<Rows>::eval(rows);

		for (Integer i = 0; i < size; ++i)
			Y[i] = -(Y[i] + X[i*X_step]);
    };

    force_inline
    static void eval(TY* Y1, Integer Y_step, const TX* X1, Integer X_step, Integer rows)
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

        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        Integer size =  get_int<Rows>::eval(rows);

		for (Integer i = 0; i < size; ++i)
			Y[i*Y_step] = -(Y[i*Y_step] + X[i*X_step]);
    };
};

template<class TY, Integer Rows>
struct ypxm<TY, TY, Rows, true, details::true_t>
{
    force_inline
    static void eval(TY* Y1, const TY* X1, Integer rows)
    {
        using simd_type = typename simd::default_simd_type<TY>::type;
        static const int vec_size   = simd_type::vector_size;

        //do not use fms for complex types; this is expanded to two
        //operations including complex multiplication, which is not needed here
        static const int use_fma    = details::is_complex<TY>::value == false;

        Integer size   =  get_int<Rows>::eval(rows);
        Integer size_1 =  get_int<Rows>::eval(rows) / vec_size;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TY* MATCL_RESTRICTED X  = X1;

        if (use_fma)
        {
            simd_type mone      = simd_type(TY(-1.0));

		    for (Integer i = 0; i < size_1; ++i)
            {
                simd_type sX_1  = simd_type::load(X+i*vec_size, std::false_type());
                simd_type sY_1  = simd_type::load(Y+i*vec_size, std::false_type());
			    simd_type res_1 = fms_f(mone, sY_1, sX_1);

                res_1.store(Y+i*vec_size, std::false_type());
            };
        }
        else
        {
		    for (Integer i = 0; i < size_1; ++i)
            {
                simd_type sX_1  = simd_type::load(X+i*vec_size, std::false_type());
                simd_type sY_1  = simd_type::load(Y+i*vec_size, std::false_type());
			    simd_type res_1 = -(sY_1 + sX_1);

                res_1.store(Y+i*vec_size, std::false_type());
            };
        };

        for (Integer i = size_1*vec_size; i < size; ++i)
        {
            Y[i]        = -(Y[i] + X[i]);
        };
    };

    force_inline
    static void eval(TY* Y1, Integer Y_step, const TY* X1, Integer X_step, Integer rows)
    {
        if (Y_step == 1 && X_step == 1)
            return eval(Y1,X1,rows);

        return ypxm<TY,TY,Rows,false,details::true_t>::eval(Y1, Y_step, X1, X_step, rows);
    };
};

}}

#pragma warning(pop)