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

#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/details/extract_type_switch.h"

//------------------------------------------------------------------
//                  eval_scalar_func
//------------------------------------------------------------------
namespace matcl { namespace details
{

template<class Mat>
struct get_value_type_raw_mat{};

template<class Val, class Struct>
struct get_value_type_raw_mat<mr::Matrix<Val,Struct>>
{
    using type = Val;
};

template<class Ret, class In>
class virtual_functor_base
{
    private:
        using ti_type   = typename ti::get_ti_type<Ret>::type;

    public:
        virtual Ret eval(const In&) const = 0;
        virtual Ret test_zero(const In& zero) const = 0;

    public:
        bool        is_special_case() const 
                        { return false; };

        template<class ret_type,class in_type>
        void        eval_special_case(matcl::Matrix& ret, ti_type,const in_type& ) const
                        { ret = 0; };
};

template<class Ret, class In, class Func, class Test>
class virtual_functor : public virtual_functor_base<Ret, In>
{
    private:
        const Func& fun;
        const Test& test;

    public:
        virtual_functor(const Func& f, const Test& t)       :fun(f),test(t) {};

        virtual Ret eval(const In& v) const override        { return fun.eval(v); };
        virtual Ret test_zero(const In& v) const override   { return test.eval(v); };

        virtual_functor(const virtual_functor&) = delete;        
};

template<class Ret, class Mat>
struct MATCL_MATREP_EXPORT function_eval_func
{
    using In = typename get_value_type_raw_mat<Mat>::type;
    static void eval(matcl::Matrix& ret, const Mat& in, const virtual_functor_base<Ret, In>& vf);
};

struct eval_func_scal : public extract_type_switch<void,eval_func_scal,true> 
{
    template<class T, class scalar_f, class test_f>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret, const scalar_f& f,
                     const test_f& t)
    {
        using value_type    = typename get_value_type_raw_mat<T>::type;
        using ret_type0     = decltype(f.eval(*(value_type*)nullptr));
        using ret_type      = typename md::promote_scalar<ret_type0>::type;

        virtual_functor<ret_type,value_type,scalar_f, test_f> vf(f,t); 

        function_eval_func<ret_type,T>::eval(ret, mat, vf);
    };

    template<class T, class scalar_f, class test_f>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret, 
                            const scalar_f& f, const test_f&)
    {
        using ret_type0     = decltype(f.eval(*(T*)nullptr));
        using ret_type      = typename md::promote_scalar<ret_type0>::type;

        ret = ret_type(f.eval(mat));
    };
};

}}

namespace matcl
{

template<class scalar_fun>
Matrix matcl::eval_scalar_func(const Matrix& A0, const scalar_fun& func)
{
    Matrix A(A0);

    Matrix ret;
    details::eval_func_scal::make<const Matrix&>(A,ret, func, func);
    return ret;
};

template<class scalar_fun>
Matrix matcl::eval_scalar_func(Matrix&& A0, const scalar_fun& func)
{
    Matrix A(std::move(A0));

    Matrix ret;
    details::eval_func_scal::make<const Matrix&>(A,ret, func, func);
    return ret;
};

template<class scalar_fun, class test_fun>
Matrix matcl::eval_scalar_func(const Matrix& A0, const scalar_fun& func,
                               const test_fun& t)
{
    Matrix A(A0);

    Matrix ret;
    details::eval_func_scal::make<const Matrix&>(A,ret, func, t);
    return ret;
};

template<class scalar_fun, class test_fun>
Matrix matcl::eval_scalar_func(Matrix&& A0, const scalar_fun& func,
                               const test_fun& t)
{
    Matrix A(std::move(A0));

    Matrix ret;
    details::eval_func_scal::make<const Matrix&>(A,ret, func, t);
    return ret;
};

};