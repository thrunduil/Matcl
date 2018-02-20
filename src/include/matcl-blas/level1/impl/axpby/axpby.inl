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

namespace matcl { namespace level1
{

template<class TX, class TY, class TA, class TB, Integer Rows>
struct axpby<TY, TX, TA, TB, Rows, false, details::true_t>
{
    force_inline
    static void eval(TY* Y1, const TX* X1, Integer rows, const TA& a, const TB& b)
    {
        Integer size    =  get_int<Rows>::eval(rows);
        using TAP       = typename details::promote_scalar<TA, TY>::type;
        using TBP       = typename details::promote_scalar<TB, TY>::type;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        TAP ap          = TAP(a);
        TBP bp          = TBP(b);

		for (Integer i = 0; i < size; ++i)
			Y[i] = ap * X[i] + bp * Y[i];
    };

    force_inline
    static void eval_step1(TY* Y1, Integer Y_step, const TX* X1, Integer rows, const TA& a, const TB& b)
    {
        Integer size    =  get_int<Rows>::eval(rows);
        using TAP       = typename details::promote_scalar<TA, TY>::type;
        using TBP       = typename details::promote_scalar<TB, TY>::type;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        TAP ap          = TAP(a);
        TBP bp          = TBP(b);

		for (Integer i = 0; i < size; ++i)
			Y[i*Y_step] = ap * X[i] + bp * Y[i*Y_step];
    };

    force_inline
    static void eval_step2(TY* Y1, const TX* X1, Integer X_step, Integer rows, const TA& a, const TB& b)
    {
        Integer size    =  get_int<Rows>::eval(rows);
        using TAP       = typename details::promote_scalar<TA, TY>::type;
        using TBP       = typename details::promote_scalar<TB, TY>::type;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        TAP ap          = TAP(a);
        TBP bp          = TBP(b);

		for (Integer i = 0; i < size; ++i)
			Y[i] = ap * X[i*X_step] + bp * Y[i];
    };

    force_inline
    static void eval(TY* Y1, Integer Y_step, const TX* X1, Integer X_step, Integer rows, const TA& a, const TB& b)
    {
        if (Y_step == 1)
        {
            if (X_step == 1)
                return eval(Y1, X1, rows, a, b);
            else
                return eval_step2(Y1, X1, X_step, rows, a, b);
        }
        else if (X_step == 1)                
        {
            return eval_step1(Y1, Y_step, X1, rows, a, b);
        };

        Integer size    =  get_int<Rows>::eval(rows);
        using TAP       = typename details::promote_scalar<TA, TY>::type;
        using TBP       = typename details::promote_scalar<TB, TY>::type;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        TAP ap          = TAP(a);
        TBP bp          = TBP(b);

		for (Integer i = 0; i < size; ++i)
			Y[i*Y_step] = ap * X[i*X_step] + bp * Y[i*Y_step];
    };
};

template<class TY, class TA, class TB, Integer Rows>
struct axpby<TY, TY, TA, TB, Rows, true, details::true_t>
{
    force_inline
    static void eval(TY* Y1, const TY* X1, Integer rows, const TA& a, const TB& b)
    {
        using simd_type = typename simd::default_simd_type<TY>::type;
        using TAP       = typename details::promote_scalar<TA, TY>::type;
        using TBP       = typename details::promote_scalar<TB, TY>::type;
        using simd_a    = typename simd::default_simd_type<TAP>::type;
        using simd_b    = typename simd::default_simd_type<TBP>::type;
        static const int vec_size   = simd_type::vector_size;

        Integer size   =  get_int<Rows>::eval(rows);
        Integer size_1 =  get_int<Rows>::eval(rows) / vec_size;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TY* MATCL_RESTRICTED X  = X1;

        TAP ap          = TAP(a);
        TBP bp          = TBP(b);

        simd_a sa(ap);
        simd_b sb(bp);

		for (Integer i = 0; i < size_1; ++i)
        {
            simd_type sX_1  = simd_type::load(X+i*vec_size, std::false_type());
            simd_type sY_1  = simd_type::load(Y+i*vec_size, std::false_type());
			simd_type res_1 = sa * sX_1 + sb * sY_1;

            res_1.store(Y+i*vec_size, std::false_type());
        };

        for (Integer i = size_1*vec_size; i < size; ++i)
            Y[i]        = ap * X[i]  + bp * Y[i];
    };

    force_inline
    static void eval(TY* Y1, Integer Y_step, const TY* X1, Integer X_step, Integer rows, const TA& a, const TB& b)
    {
        if (Y_step == 1 && X_step == 1)
            return eval(Y1, X1, rows, a, b);

        return axpby<TY,TY,TA,TB,Rows,false,details::true_t>::eval(Y1, Y_step, X1, X_step, rows, a, b);
    };
};

}}