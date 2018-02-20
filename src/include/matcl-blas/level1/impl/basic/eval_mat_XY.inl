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

#include "matcl-blas/level1/level1_basic.h"

#pragma warning(push)
#pragma warning(disable: 4127)  //conditional expression is constant

namespace matcl { namespace level1
{

//-------------------------------------------------------------
//              remove Continuous == 0
//-------------------------------------------------------------
template<bool Select, class T1, class T2, Integer Rows, Integer Cols, class Func>
struct eval_mat_func_XY<Select, T1, T2, Rows, Cols, 0, Func, details::true_t>
{
    static void eval(T1* Y, Integer Y_ld, const T2* X, Integer X_ld, 
                     Integer rows, Integer cols, const Func& fun)
    {
        if (Cols == 1)
        {
            return eval_mat_func_XY<Select, T1, T2, Rows, 1, 1, Func, details::true_t>
                    ::eval(Y, Y_ld, X, X_ld, rows, cols, fun);
        };

        if (Select == false)
        {
            return eval_mat_func_XY<false, T1, T2, Rows, Cols, 2, Func, details::true_t>
                ::eval(Y, Y_ld, X, X_ld, rows, cols, fun);
        };

        if (Y_ld == get_int<Rows>::eval(rows) && X_ld == get_int<Rows>::eval(rows))
        {
            return eval_mat_func_XY<true, T1, T2, Rows, Cols, 1, Func, details::true_t>
                    ::eval(Y, Y_ld, X, X_ld, rows, cols, fun);
        }
        else
        {
            return eval_mat_func_XY<true, T1, T2, Rows, Cols, 2, Func, details::true_t>
                    ::eval(Y, Y_ld, X, X_ld, rows, cols, fun);
        } 
    };
};

//-------------------------------------------------------------
//              Continuous == 1
//-------------------------------------------------------------
template<bool Select, class T1, class T2, Integer Rows, Integer Cols, class Func>
struct eval_mat_func_XY<Select, T1, T2, Rows, Cols, 1, Func, details::true_t>
{
    static void eval(T1* Y, Integer Y_ld, const T2* X, Integer X_ld, 
                     Integer rows, Integer cols, const Func& func)
    {
        (void)Y_ld;
        (void)X_ld;
        Integer size =  get_int<Rows>::eval(rows) * get_int<Cols>::eval(cols);
        func.eval<Rows*Cols>(Y, X, size);
    };
};

//-------------------------------------------------------------
//              Continuous == 2
//-------------------------------------------------------------
template<class T1, class T2, Integer Cols, class Func>
struct eval_mat_func_XY<true, T1, T2, 0, Cols, 2, Func, details::true_t>
{
    static void eval(T1* Y, Integer Y_ld, const T2* X, Integer X_ld, 
                     Integer rows, Integer cols, const Func& func)
    {
        switch(rows)
        {
            case 1:
                return eval_mat_func_XY<true, T1, T2, 1, Cols, 2, Func, details::true_t>
                    ::eval(Y, Y_ld, X, X_ld, rows, cols, func);
            case 2:
                return eval_mat_func_XY<true, T1, T2, 2, Cols, 2, Func, details::true_t>
                    ::eval(Y, Y_ld, X, X_ld, rows, cols, func);
            case 3:
                return eval_mat_func_XY<true, T1, T2, 3, Cols, 2, Func, details::true_t>
                    ::eval(Y, Y_ld, X, X_ld, rows, cols, func);
            case 4:
                return eval_mat_func_XY<true, T1, T2, 4, Cols, 2, Func, details::true_t>
                    ::eval(Y, Y_ld, X, X_ld, rows, cols, func);
            default:
                return eval_mat_func_XY<false, T1, T2, 0, Cols, 2, Func, details::true_t>
                    ::eval(Y, Y_ld, X, X_ld, rows, cols, func);
        };
    };
};

template<bool Select, class T1, class T2, Integer Rows, Integer Cols, class Func>
struct eval_mat_func_XY<Select, T1, T2, Rows, Cols, 2, Func, details::true_t>
{
    static void eval(T1* Y, Integer Y_ld, const T2* X, Integer X_ld, 
                     Integer rows, Integer cols, const Func& func)
    {
        for (Integer i = 0 ; i < get_int<Cols>::eval(cols); ++i)
        {
            func.eval<Rows>(Y, X, rows);

            Y           += Y_ld;
            X           += X_ld;
        }
    };
};

}}

#pragma warning(pop)