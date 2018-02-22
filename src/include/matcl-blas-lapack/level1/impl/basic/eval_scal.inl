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

#pragma warning(push)
#pragma warning(disable:4100)   //unreferenced formal parameter

namespace matcl { namespace level1
{

namespace details
{

template<class TY, Integer Rows, class Func, bool Use_simd, class ... Args>
struct eval_scal_func_Y_impl;

template<class TY, Integer Rows, class Func, class ... Args>
struct eval_scal_func_Y_impl<TY, Rows, Func, false, Args...>
{
    force_inline
    static void eval(TY* Y1, Integer rows, const Func& func, Args&& ... args)
    {
        Integer size    =  get_int<Rows>::eval(rows);

        TY* MATCL_RESTRICTED Y  = Y1;

		for (Integer i = 0; i < size; ++i)
			Y[i]        = func.eval(Y[i],std::forward<Args>(args)...);
    };

    force_inline
    static void eval(TY* Y1, Integer Y_step, Integer rows, const Func& func, Args&& ... args)
    {
        if (Y_step == 1)
            return eval(Y1,rows,func,std::forward<Args>(args)...);

        Integer size    =  get_int<Rows>::eval(rows);

        TY* MATCL_RESTRICTED Y  = Y1;

		for (Integer i = 0; i < size; ++i)
			Y[i*Y_step] = func.eval(Y[i*Y_step],std::forward<Args>(args)...);
    };
};

template<class TY, Integer Rows, class Func, class ... Args>
struct eval_scal_func_Y_impl<TY, Rows, Func, true, Args...>
{
    force_inline
    static void eval(TY* Y1, Integer rows, const Func& func, Args&& ... args)
    {
        using simd_y    = typename simd::default_simd_type<TY>::type;

        static const int vec_size   = simd_y::vector_size;

        Integer size    =  get_int<Rows>::eval(rows);
        Integer size_1  =  get_int<Rows>::eval(rows) / vec_size;

        TY* MATCL_RESTRICTED Y = Y1;

		for (Integer i = 0; i < size_1; ++i)
        {
            simd_y sY_1     = simd_y::load(Y+i*vec_size, std::false_type());
			simd_y res_1    = func.eval_simd(sY_1,std::forward<Args>(args)...);

            res_1.store(Y+i*vec_size, std::false_type());
        };        

        for (Integer i = size_1*vec_size; i < size; ++i)
            Y[i]        = func.eval(Y[i],std::forward<Args>(args)...);
    };

    force_inline
    static void eval(TY* Y1, Integer Y_step, Integer rows, const Func& fun, Args&& ... args)
    {
        if (Y_step == 1)
            return eval(Y1,rows,fun,std::forward<Args>(args)...);

        return eval_scal_func_Y_impl<TY, Rows, Func, false, Args...>
                ::eval(Y1, Y_step, rows, fun,std::forward<Args>(args)...);
    };
};

};

//-----------------------------------------------------------------------
//                  eval_vec_func_Y
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = func(Y, args) for every elements in vector Y, where args
*  denote zero or more scalar arguments
*
*  Template arguments
*  =======
*  TY           - type of elements in Y matrix
*  Rows         - number of rows of the matrix, 0 if unknown statically
*  Func         - class with member function implementing:
*
*                       template<class TY, Args>
*                       TY eval(T& y, Args args) const
*
*                       template<class TS, Args>
*                       TS eval_simd(T& y, Args args) const
*
*                 that evaluates the function func for every element y of
*                 type TY in case of eval function or one of the simd types TS
*                 is case of eval_simd function, Args denotes types
*                 of additional arguments and args are additional scalar arguments
*                 (possibly zero).
*
*                  Additionally one must define a template
*
*                       template<class T, Args>
*                       using enable_simd   = [some type ES];
*
*                  such that if ES::value = true, then simd version can be used 
*                  for given type of array T and additional scalar arguments Args
*
*               Example:
*                   struct Fun
*                   {
*                       template<class TY, class TA>
*                       using enable_simd   = details::check_simd<TY>;
*
*                       template<class T, class TA>
*                       static T eval(const T& in, const TA& a)
*                       {
*                           return a * in;
*                       };
*
*                       template<class T, class TA>
*                       static T eval_simd(const T& in, const TA& a)
*                       {
*                           return T(a) * in;
*                       };
*                   };
*
*  Args         - types of additional arguments
*   
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  rows         - number of elements
*  fun          - instance of the class Func
*  args         - additional arguments passed to the function fun
*/
template<class TY, Integer Rows, class Func, class ... Args>
struct eval_scal_func_Y
{
    static void eval(TY* Y, Integer rows, const Func& fun, Args&& ... args);
    static void eval(TY* Y, Integer Y_step, Integer rows, const Func& fun, Args&& ... args);
};

template<class TY, Integer Rows, class Func, class ... Args>
force_inline
void eval_scal_func_Y<TY, Rows, Func, Args...>
        ::eval(TY* Y1, Integer Y_step, Integer rows, const Func& fun, Args&& ... args)
{
    static const bool use_simd  = Func::enable_simd<TY, Args...>::value;

    return details::eval_scal_func_Y_impl<TY, Rows, Func, use_simd, Args...>
                ::eval(Y1, Y_step, rows, fun, std::forward<Args>(args)...);
};

template<class TY, Integer Rows, class Func, class ... Args>
force_inline
void eval_scal_func_Y<TY, Rows, Func, Args...>
        ::eval(TY* Y1, Integer rows, const Func& func, Args&& ... args)
{
    static const bool use_simd  = Func::enable_simd<TY, Args...>::value;

    return details::eval_scal_func_Y_impl<TY, Rows, Func, use_simd, Args...>
                ::eval(Y1, rows, func, std::forward<Args>(args)...);
};

}}

#pragma warning(pop)