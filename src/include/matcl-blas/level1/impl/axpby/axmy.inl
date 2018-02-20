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

#include "matcl-blas/level1/level1_axpby.h"

namespace matcl { namespace level1
{

template<class TY, class TX, class TA, Integer Rows>
struct axmy<TY, TX, TA, Rows, false, details::true_t>
{
    force_inline
    static void eval(TY* Y1, const TX* X1, Integer rows, const TA& a)
    {
        Integer size    =  get_int<Rows>::eval(rows);
        using TAP       = typename details::promote_scalar<TA, TY>::type;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        TAP ap                  = TAP(a);

		for (Integer i = 0; i < size; ++i)
			Y[i] = ap * X[i] - Y[i];
    };

    force_inline
    static void eval_step1(TY* Y1, Integer Y_step, const TX* X1, Integer rows, const TA& a)
    {
        Integer size    =  get_int<Rows>::eval(rows);
        using TAP       = typename details::promote_scalar<TA, TY>::type;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        TAP ap                  = TAP(a);

		for (Integer i = 0; i < size; ++i)
			Y[i*Y_step] = ap * X[i] - Y[i*Y_step];
    };

    force_inline
    static void eval_step2(TY* Y1, const TX* X1, Integer X_step, Integer rows, const TA& a)
    {
        Integer size    =  get_int<Rows>::eval(rows);
        using TAP       = typename details::promote_scalar<TA, TY>::type;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        TAP ap                  = TAP(a);

		for (Integer i = 0; i < size; ++i)
			Y[i] = ap * X[i*X_step] - Y[i];
    };

    force_inline
    static void eval(TY* Y1, Integer Y_step, const TX* X1, Integer X_step, Integer rows, const TA& a)
    {
        if (Y_step == 1)
        {
            if (X_step == 1)
                return eval(Y1, X1, rows, a);
            else
                return eval_step2(Y1, X1, X_step, rows, a);
        }
        else if (X_step == 1)                
        {
            return eval_step1(Y1, Y_step, X1, rows, a);
        };

        Integer size    =  get_int<Rows>::eval(rows);
        using TAP       = typename details::promote_scalar<TA, TY>::type;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        TAP ap                  = TAP(a);

		for (Integer i = 0; i < size; ++i)
			Y[i*Y_step] = ap * X[i*X_step] - Y[i*Y_step];
    };
};

template<class TY, class TA, Integer Rows>
struct axmy<TY, TY, TA, Rows, true, details::true_t>
{
    force_inline
    static void eval(TY* Y1, const TY* X1, Integer rows, const TA& a)
    {
        using simd_type = typename simd::default_simd_type<TY>::type;
        using TAP       = typename details::promote_scalar<TA, TY>::type;
        using simd_a    = typename simd::default_simd_type<TAP>::type;

        static const int vec_size   = simd_type::vector_size;

        Integer size   =  get_int<Rows>::eval(rows);
        Integer size_1 =  get_int<Rows>::eval(rows) / vec_size;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TY* MATCL_RESTRICTED X  = X1;

        TAP ap          = TAP(a);
        simd_a sa(ap);

		for (Integer i = 0; i < size_1; ++i)
        {
            simd_type sX_1  = simd_type::load(X+i*vec_size, std::false_type());
            simd_type sY_1  = simd_type::load(Y+i*vec_size, std::false_type());
			simd_type res_1 = fms_f(sa,sX_1,sY_1);

            res_1.store(Y+i*vec_size, std::false_type());
        };

        for (Integer i = size_1*vec_size; i < size; ++i)
            Y[i]        = ap * X[i]  - Y[i];
    };

    force_inline
    static void eval(TY* Y1, Integer Y_step, const TY* X1, Integer X_step, Integer rows, const TA& a)
    {
        if (Y_step == 1 && X_step == 1)
            return eval(Y1, X1, rows, a);

        return axmy<TY,TY,TA,Rows,false,details::true_t>::eval(Y1, Y_step, X1, X_step, rows, a);
    };
};

}}