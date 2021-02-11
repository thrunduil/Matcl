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

#include "matcl-blas-lapack/level1/level1_other.h"

namespace matcl { namespace level1
{

template<class TZ, class TX, class TY, Integer Rows>
struct xy<TZ, TX, TY, Rows, false, details::true_t>
{
    force_inline
    static void eval(TZ* Z1, const TX* X1, const TY* Y1, Integer rows)
    {
        Integer size    =  get_int<Rows>::eval(rows);
        using TXP       = typename details::promote_scalar<TX, TZ>::type;

        TZ* MATCL_RESTRICTED Z        = Z1;
        const TX* MATCL_RESTRICTED X  = X1;
        const TY* MATCL_RESTRICTED Y  = Y1;

		for (Integer i = 0; i < size; ++i)
			Z[i] = TXP(X[i]) * Y[i];
    };
    
    template <Integer ZS, Integer YS, Integer XS>
    force_inline
    static void eval_step(TZ* Z1, Integer Z_step, const TX* X1, Integer X_step, const TY* Y1, Integer Y_step, 
                     Integer rows)
    {
        Integer size    = get_int<Rows>::eval(rows);
        Integer zs      = get_int<ZS>::eval(Z_step);
        Integer ys      = get_int<YS>::eval(Y_step);
        Integer xs      = get_int<XS>::eval(X_step);

        using TXP       = typename details::promote_scalar<TX, TZ>::type;

        TZ* MATCL_RESTRICTED Z        = Z1;
        const TX* MATCL_RESTRICTED X  = X1;
        const TY* MATCL_RESTRICTED Y  = Y1;

		for (Integer i = 0; i < size; ++i)
			Z[i*zs] = TXP(X[i*xs]) * Y[i*ys];
    };

    static void eval(TZ* Z1, Integer Z_step, const TX* X1, Integer X_step, const TY* Y1, Integer Y_step, 
                     Integer rows)
    {
        if (Z_step == 1)
        {
            if (Y_step == 1)
            {
                 if (X_step == 1)
                    return eval(Z1,X1,Y1,rows);
                 else
                     return eval_step<1,1,0>(Z1, Z_step, X1, X_step, Y1, Y_step, rows);
            }
            else
            {
                 if (X_step == 1)
                    return eval_step<1,0,1>(Z1, Z_step, X1, X_step, Y1, Y_step, rows);
                 else
                     return eval_step<1,0,0>(Z1, Z_step, X1, X_step, Y1, Y_step, rows);
            };
        }
        else
        {
            if (Y_step == 1)
            {
                 if (X_step == 1)
                    return eval_step<0,1,1>(Z1, Z_step, X1, X_step, Y1, Y_step, rows);
                 else
                     return eval_step<0,1,0>(Z1, Z_step, X1, X_step, Y1, Y_step, rows);
            }
            else
            {
                 if (X_step == 1)
                    return eval_step<0,0,1>(Z1, Z_step, X1, X_step, Y1, Y_step, rows);
                 else
                     return eval_step<0,0,0>(Z1, Z_step, X1, X_step, Y1, Y_step, rows);
            };
        };
    };
};

template<class TZ, Integer Rows>
struct xy<TZ, TZ, TZ, Rows, true, details::true_t>
{
    force_inline
    static void eval(TZ* Z, const TZ* X, const TZ* Y, Integer rows)
    {
        using simd_type = typename simd::default_simd_type<TZ>::type;
        static const int vec_size   = simd_type::vector_size;

        Integer size   =  get_int<Rows>::eval(rows);
        Integer size_1 =  get_int<Rows>::eval(rows) / vec_size;

		for (Integer i = 0; i < size_1; ++i)
        {
            simd_type sX_1  = simd_type::load(X+i*vec_size, std::false_type());
            simd_type sY_1  = simd_type::load(Y+i*vec_size, std::false_type());
			simd_type res_1 = sX_1 * sY_1;

            res_1.store(Z+i*vec_size, std::false_type());
        };

        for (Integer i = size_1*vec_size; i < size; ++i)
            Z[i] = X[i] * Y[i];
    };

    static void eval(TZ* Z1, Integer Z_step, const TZ* X1, Integer X_step, const TZ* Y1, 
                     Integer Y_step, Integer rows)
    {
        if (Z_step == 1 && Y_step == 1 && X_step == 1)
            return eval(Z1, X1, Y1, rows);

        return xy<TZ,TZ,TZ,Rows,false,details::true_t>::eval(Z1, Z_step, X1, X_step, Y1, Y_step, rows);
    };
};

}}