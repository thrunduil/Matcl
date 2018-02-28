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

#include "matcl-blas-lapack/level1/level1_axpby.h"

namespace matcl { namespace level1
{

template<class TY, class TX, class TA, Integer Rows>
struct ax<TY, TX, TA, Rows, false, details::true_t>
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
			Y[i] = ap * X[i];
    };

    force_inline
    static void eval_stepx(TY* Y1, const TX* X1, Integer X_step, Integer rows, const TA& a)
    {
        Integer size    =  get_int<Rows>::eval(rows);
        using TAP       = typename details::promote_scalar<TA, TY>::type;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        TAP ap                  = TAP(a);

		for (Integer i = 0; i < size; ++i)
			Y[i] = ap * X[i*X_step];
    };

    force_inline
    static void eval_stepy(TY* Y1, Integer Y_step, const TX* X1, Integer rows, const TA& a)
    {
        Integer size    =  get_int<Rows>::eval(rows);
        using TAP       = typename details::promote_scalar<TA, TY>::type;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        TAP ap                  = TAP(a);

		for (Integer i = 0; i < size; ++i)
			Y[i*Y_step] = ap * X[i];
    };

    force_inline
    static void eval(TY* Y1, Integer Y_step, const TX* X1, Integer X_step, Integer rows, const TA& a)
    {
        Integer size    =  get_int<Rows>::eval(rows);
        using TAP       = typename details::promote_scalar<TA, TY>::type;

        if (Y_step == 1)
        {
            if (X_step == 1)
                return eval(Y1,X1,rows,a);
            else
                return eval_stepx(Y1,X1,X_step,rows,a);
        }
        else if (X_step == 1)
        {
            return eval_stepy(Y1,Y_step,X1,rows,a);
        };

        TY* MATCL_RESTRICTED Y        = Y1;
        const TX* MATCL_RESTRICTED X  = X1;

        TAP ap                  = TAP(a);

		for (Integer i = 0; i < size; ++i)
			Y[i*Y_step] = ap * X[i*X_step];
    };
};

template<class TY, class TA, Integer Rows>
struct ax<TY, TY, TA, Rows, true, details::true_t>
{
    force_inline
    static void eval(TY* Y1, const TY* X1, Integer rows, const TA& a)
    {
        using simd_y    = typename simd::default_simd_type<TY>::type;
        using TAP       = typename details::promote_scalar<TA, TY>::type;
        using simd_a    = typename simd::default_simd_type<TAP>::type;

        static const int vec_size   = simd_y::vector_size;

        Integer size   =  get_int<Rows>::eval(rows);
        Integer size_1 =  get_int<Rows>::eval(rows) / vec_size;

        TY* MATCL_RESTRICTED Y        = Y1;
        const TY* MATCL_RESTRICTED X  = X1;

        TAP ap          = TAP(a);
        simd_a sa(ap);

		for (Integer i = 0; i < size_1; ++i)
        {
            simd_y sX_1  = simd_y::load(X+i*vec_size, std::false_type());
			simd_y res_1 = sa * sX_1;

            res_1.store(Y+i*vec_size, std::false_type());
        };

        for (Integer i = size_1*vec_size; i < size; ++i)
            Y[i]        = ap * X[i];
    };

    force_inline
    static void eval(TY* Y1, Integer Y_step, const TY* X1, Integer X_step, Integer rows, 
                     const TA& a)
    {
        if (Y_step == 1 && X_step == 1)
            return eval(Y1, X1, rows, a);

        return ax<TY,TY,TA,Rows,false,details::true_t>::eval(Y1, Y_step, X1, X_step, rows, a);
    };
};

template<bool Select, class TY, class TX, class TA, Integer Rows, Integer Cols, Integer Continuous>
struct ax_test_mat<Select, TY, TX, TA, Rows, Cols, Continuous, details::true_t>
{
    static void eval(TY* Y, Integer Y_ld, const TX* X, Integer X_ld, 
                     Integer rows, Integer cols, const TA& a)
    {
        if (a == TA(0.0))
        {
            return eval_mat_func_Y<Select, TY, Rows, Cols, Continuous, details::func_set_val_zero<TY>>
                        ::eval(Y, Y_ld, rows, cols, details::func_set_val_zero<TY>());
        }
        else if ( a == TA(1.0) )
        {
            //Y = X
            return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_copy<TY,TX>>
                        ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_copy<TY,TX>());
        }
        else if ( a == TA(-1.0) )
        {
            //Y = -X
            return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_mx<TY,TX>>
                        ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_mx<TY,TX>());
        }
        else if (details::is_real<TA>::eval(a))
        {
            using TAR   = typename details::real_type<TA>::type;
            //Y = real(a)*X
            return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_ax<TY,TX,TAR>>
                        ::eval(Y, Y_ld, X, X_ld, rows, cols, 
                                details::func_ax<TY,TX,TAR>(details::eval_real<TA>::eval(a)));
        }
        else
        {
            //Y = a*X
            return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_ax<TY,TX,TA>>
                        ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_ax<TY,TX,TA>(a));
        };
    };
};

}}