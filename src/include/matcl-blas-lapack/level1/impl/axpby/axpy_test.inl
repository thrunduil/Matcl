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

#include "matcl-blas-lapack/level1/level1_axpby.h"

namespace matcl { namespace level1
{

template<bool Select, class TY, class TX, class TA, Integer Rows, 
        Integer Cols, Integer Continuous>
struct axpy_test_mat<Select, TY, TX, TA, Rows, Cols, Continuous, details::true_t>
{
    static void eval(TY* Y, Integer Y_ld, const TX* X, Integer X_ld, Integer rows, Integer cols, 
                     const TA& a)
    {
        if (a == TA(0.0))
        {
            //Y = Y
            return;
        }
        else if ( a == TA(1.0) )
        {
            //Y = X + Y
            return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_ypx<TY,TX>>
                        ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_ypx<TY,TX>());
        }
        else if (a == TA(-1.0))
        {
            //Y = Y-X
            return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_ymx<TY,TX>>
                        ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_ymx<TY,TX>());                
        }
        else if (details::is_real<TA>::eval(a))
        {
            using TAR   = typename details::real_type<TA>::type;
            //Y = real(a)*X + Y
            return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_axpy<TY,TX,TAR>>
                        ::eval(Y, Y_ld, X, X_ld, rows, cols, 
                                details::func_axpy<TY,TX,TAR>(details::eval_real<TA>::eval(a)));
        }
        else
        {
            //Y = a*X + Y
            return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_axpy<TY,TX,TA>>
                        ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_axpy<TY,TX,TA>(a));
        };
    };
};

}}