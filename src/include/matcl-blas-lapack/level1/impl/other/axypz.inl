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
#include "matcl-simd/simd.h"
#include "matcl-simd/simd_complex.h"

#pragma warning(push)
#pragma warning(disable: 4127) // conditional expression is constant

namespace matcl { namespace level1
{

namespace details
{

//how to evaluate expression a*x*y
enum class axy
{
    ax_y,       //form (a * x) * y
    a_xy,       //form a * (x * y)
    x_ay,       //form x * (a * y)
};

template<bool Is_compl_a, bool Is_compl_x, bool Is_compl_y, bool Req_pr_X, bool Req_pr_Y>
struct get_axy_impl
{
    static const axy value = axy::ax_y;
};

template<bool Req_pr_X, bool Req_pr_Y>
struct get_axy_impl<true, false, false, Req_pr_X, Req_pr_Y>
{
    static const axy value = (Req_pr_X == false || Req_pr_Y == false)? axy::a_xy : axy::ax_y;
};

template<bool Req_pr_X, bool Req_pr_Y>
struct get_axy_impl<false, true, false, Req_pr_X, Req_pr_Y>
{
    static const axy value = axy::x_ay;
};

//how to evaluate expression a*x*y
//ax_y and x_ay are always possible
//a_xy is possible if x or y has the same precision as a
//a_xy is preferred over ax_y if x and y are real and a is complex
//x_ay is preferred over ax_y if a and y are real and x is complex
template<class TA, class TX, class TY>
struct get_axy
{
    using TXP       = typename details::promote_scalar<TX, TA>::type;
    using TYP       = typename details::promote_scalar<TY, TA>::type;

    static const bool X_req_pr   = std::is_same<TX,TXP>::value == false;
    static const bool Y_req_pr   = std::is_same<TY,TYP>::value == false;
    static const bool is_compl_a = details::is_complex<TA>::value;
    static const bool is_compl_x = details::is_complex<TX>::value;
    static const bool is_compl_y = details::is_complex<TY>::value;

    static const axy value = get_axy_impl<is_compl_a,is_compl_x,is_compl_y,X_req_pr,Y_req_pr>
                                ::value;
};

};

template<class TZ, class TX, class TY, class TA, Integer Rows>
struct axypz<TZ, TX, TY, TA, Rows, false, details::true_t>
{
    force_inline
    static void eval(TZ* Z1, const TX* X1, const TY* Y1, Integer rows, const TA& a)
    {
        Integer size =  get_int<Rows>::eval(rows);

        using TAP       = typename details::promote_scalar<TA, TZ>::type;
        using axy       = details::axy;
        static const axy a_XY  = details::get_axy<TAP,TX,TY>::value;

        TZ* MATCL_RESTRICTED Z        = Z1;
        const TX* MATCL_RESTRICTED X  = X1;
        const TY* MATCL_RESTRICTED Y  = Y1;

        TAP ap          = TAP(a);        

        if (a_XY == axy::ax_y)
        {
		    for (Integer i = 0; i < size; ++i)
			    Z[i] = (ap * X[i]) * Y[i] + Z[i];
        }
        else if (a_XY == axy::a_xy)
        {
		    for (Integer i = 0; i < size; ++i)
			    Z[i] = ap * (X[i] * Y[i]) + Z[i];
        }
        else
        {
		    for (Integer i = 0; i < size; ++i)
			    Z[i] = X[i] * (ap *Y[i]) + Z[i];
        }
    };

    template<Integer ZS, Integer XS, Integer YS>
    force_inline
    static void eval_step(TZ* Z1, Integer Z_step, const TX* X1, Integer X_step, const TY* Y1, 
                     Integer Y_step, Integer rows, const TA& a)
    {
        Integer size = get_int<Rows>::eval(rows);
        Integer zs   = get_int<ZS>::eval(Z_step);
        Integer ys   = get_int<YS>::eval(Y_step);
        Integer xs   = get_int<XS>::eval(X_step);

        using TAP       = typename details::promote_scalar<TA, TZ>::type;
        using axy       = details::axy;
        static const axy a_XY  = details::get_axy<TAP,TX,TY>::value;

        TZ* MATCL_RESTRICTED Z        = Z1;
        const TX* MATCL_RESTRICTED X  = X1;
        const TY* MATCL_RESTRICTED Y  = Y1;

        TAP ap          = TAP(a);

        if (a_XY == axy::ax_y)
        {
		    for (Integer i = 0; i < size; ++i)
			    Z[i*zs] = (ap * X[i*xs]) * Y[i*ys] + Z[i*zs];
        }
        else if (a_XY == axy::a_xy)
        {
		    for (Integer i = 0; i < size; ++i)
			    Z[i*zs] = ap * (X[i*xs] * Y[i*ys]) + Z[i*zs];
        }
        else
        {
		    for (Integer i = 0; i < size; ++i)
			    Z[i*zs] = X[i*xs] * (ap * Y[i*ys]) + Z[i*zs];
        }
    };

    static void eval(TZ* Z1, Integer Z_step, const TX* X1, Integer X_step, const TY* Y1, Integer Y_step, 
                     Integer rows, const TA& a)
    {
        if (Z_step == 1)
        {
            if (Y_step == 1)
            {
                 if (X_step == 1)
                    return eval(Z1,X1,Y1,rows,a);
                 else
                     return eval_step<1,1,0>(Z1, Z_step, X1, X_step, Y1, Y_step, rows,a);
            }
            else
            {
                 if (X_step == 1)
                    return eval_step<1,0,1>(Z1, Z_step, X1, X_step, Y1, Y_step, rows,a);
                 else
                     return eval_step<1,0,0>(Z1, Z_step, X1, X_step, Y1, Y_step, rows,a);
            };
        }
        else
        {
            if (Y_step == 1)
            {
                 if (X_step == 1)
                    return eval_step<0,1,1>(Z1, Z_step, X1, X_step, Y1, Y_step, rows,a);
                 else
                     return eval_step<0,1,0>(Z1, Z_step, X1, X_step, Y1, Y_step, rows,a);
            }
            else
            {
                 if (X_step == 1)
                    return eval_step<0,0,1>(Z1, Z_step, X1, X_step, Y1, Y_step, rows,a);
                 else
                     return eval_step<0,0,0>(Z1, Z_step, X1, X_step, Y1, Y_step, rows,a);
            };
        };
    };
};

template<class TZ, class TA, Integer Rows>
struct axypz<TZ, TZ, TZ, TA, Rows, true, details::true_t>
{
    force_inline
    static void eval(TZ* Z, const TZ* X, const TZ* Y, Integer rows, const TA& a)
    {
        using simd_type = typename simd::default_simd_type<TZ>::type;
        using TAP       = typename details::promote_scalar<TA, TZ>::type;
        using simd_a    = typename simd::default_simd_type<TAP>::type;

        static const int vec_size   = simd_type::vector_size;

        Integer size   =  get_int<Rows>::eval(rows);
        Integer size_1 =  get_int<Rows>::eval(rows) / vec_size;

        TAP ap          = TAP(a);
        simd_a sa(ap);

		for (Integer i = 0; i < size_1; ++i)
        {
            simd_type sX_1  = simd_type::load(X+i*vec_size, std::false_type());
            simd_type sY_1  = simd_type::load(Y+i*vec_size, std::false_type());
            simd_type sZ_1  = simd_type::load(Z+i*vec_size, std::false_type());
			simd_type res_1 = fma_f(sa,sX_1*sY_1,sZ_1);

            res_1.store(Z+i*vec_size, std::false_type());
        };

        for (Integer i = size_1*vec_size; i < size; ++i)
            Z[i]        = ap * X[i]*Y[i]  + Z[i];
    };

    static void eval(TZ* Z1, Integer Z_step, const TZ* X1, Integer X_step, 
                     const TZ* Y1, Integer Y_step, Integer rows, const TA& a)
    {
        if (Z_step == 1 && X_step == 1 && Y_step == 1)
            return eval(Z1, X1, Y1, rows, a);

        return axypz<TZ,TZ,TZ,TA,Rows,false,details::true_t>
            ::eval(Z1, Z_step, X1, X_step, Y1, Y_step, rows, a);
    };
};

}}

#pragma warning(pop)