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

#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-matrep/details/extract_type2_switch.h"
#include "matcl-matrep/details/func_unary.inl"

namespace matcl { namespace details
{

namespace mr = matcl::raw; 

template<class Ret, class V1, class V2>
struct op_virtual_base;

template<class Ret, class M1, class M2>
struct MATCL_MATREP_EXPORT function_eval_bin_func
{
    using In1   = typename get_value_type_raw_mat<M1>::type;
    using In2   = typename get_value_type_raw_mat<M2>::type;
    using func  = op_virtual_base<Ret, In1, In2>;

    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const func& vf);
};

template<class Ret, class V1, class V2>
struct op_virtual_base
{
    using ret_type  = matcl::raw::Matrix<Ret,struct_dense>;

    template<bool zero_on_right>
    using is_eval_zero_id   = std::false_type;

    virtual Ret eval(const V1& v1, const V2& v2) const = 0;
    virtual Ret test(const V1& v1, const V2& v2) const = 0;

    template<class ret, bool is_inv, class TI1, class TI2>
    ti::ti_type<ret> return_type(TI1 t1, TI2 t2) const
    {
        V1 Z1   = default_value<V1>(t1);
        V2 Z2   = default_value<V2>(t2);

        Ret ret = this->test(Z1,Z2);
        return ti::get_ti(ret);
    };

    template<class Ret, class In1, class In2>
    Ret eval(const V1& a, const V2& b) const
    {
        return this->eval(a,b);
    };
};


template<class Ret, class V1, class V2, class Bin_f, class Test_f>
class op_virtual : public op_virtual_base<Ret,V1,V2>
{
    public:
        using base_type = op_virtual_base<Ret,V1,V2>;

    private:
        const Bin_f&    m_fun;
        const Test_f&   m_test;

    public:
        op_virtual(const Bin_f& f, const Test_f& t) : m_fun(f), m_test(t) {};

        virtual Ret eval(const V1& v1, const V2& v2) const override
        {
            return m_fun.eval(v1,v2);
        };

        virtual Ret test(const V1& v1, const V2& v2) const override
        {
            return m_test.eval(v1,v2);
        };

        op_virtual(const op_virtual&) = delete;
};

struct eval_FUNC_BIN : public extract_type2_switch<void,eval_FUNC_BIN, mr::val_type_corrector_dense>
{
    template<class T1, class T2, class scalar_f, class test_f>
    static void eval_mat_mat(const T1& A, const T2& B, matcl::Matrix& ret, const scalar_f& f, 
                            const test_f& t)
    {
        using V1            = typename get_value_type_raw_mat<T1>::type;
        using V2            = typename get_value_type_raw_mat<T2>::type;
        using ret_type0     = decltype(f.eval(*(V1*)nullptr,*(V2*)nullptr));
        using ret_type      = typename md::promote_scalar<ret_type0>::type;

        using op_type       = op_virtual<ret_type,V1,V2,scalar_f,test_f>;
        using op_base       = typename op_type::base_type;
        op_type op(f,t);

        return function_eval_bin_func<ret_type, T1, T2>::eval(ret,A,B,op);        
    };

    template<class T1, class T2, class scalar_f, class test_f>
    static void eval_scal_scal(const T1& A, const T2& B, matcl::Matrix& ret, const scalar_f& f,
                               const test_f&)
    {
        (void)f;
        using ret_type0     = decltype(f.eval(*(T1*)nullptr,*(T2*)nullptr));
        using ret_type      = typename md::promote_scalar<ret_type0>::type;

        ret = ret_type(f.eval(A,B));
    };

    template<class T1, class T2, class scalar_f, class test_f>
    static void eval_mat_scal(const T1& mat, const T2& scal, matcl::Matrix& ret, const scalar_f& f,
                              const test_f& t)
    {
        using FullMatrix = raw::Matrix<T2,struct_dense>;
        FullMatrix m_scal(ti::get_ti(scal),scal,1,1);
        return eval_mat_mat<T1,FullMatrix>(mat,m_scal, ret, f, t);
    };

    template<class T1, class T2, class scalar_f, class test_f>
    static void eval_scal_mat(const T1& scal, const T2& mat, matcl::Matrix& ret, const scalar_f& f,
                              const test_f& t)
    {
        using FullMatrix = raw::Matrix<T1,struct_dense>;
        FullMatrix m_scal(ti::get_ti(scal),scal,1,1);
        return eval_mat_mat<FullMatrix,T2>(m_scal,mat, ret, f, t);
    };
};

}};

namespace matcl
{

template<class binary_function>
Matrix matcl::eval_binary_func(const Matrix& A10, const Matrix& A20, const binary_function& f)
{
    //make copy
    Matrix A1(A10);
    Matrix A2(A20);

    Matrix ret;
    details::eval_FUNC_BIN::make(A1, A2, ret, f, f);
    return ret;
};

template<class binary_function>
Matrix matcl::eval_binary_func(Matrix&& A10, const Matrix& A20, const binary_function& f)
{
    //make copy
    Matrix A1(std::move(A10));
    Matrix A2(A20);

    Matrix ret;
    details::eval_FUNC_BIN::make(A1, A2, ret, f, f);
    return ret;
};

template<class binary_function>
Matrix matcl::eval_binary_func(const Matrix& A10, Matrix&& A20, const binary_function& f)
{
    //make copy
    Matrix A1(A10);
    Matrix A2(std::move(A20));

    Matrix ret;
    details::eval_FUNC_BIN::make(A1, A2, ret, f, f);
    return ret;
};

template<class binary_function>
Matrix matcl::eval_binary_func(Matrix&& A10, Matrix&& A20, const binary_function& f)
{
    //make copy
    Matrix A1(std::move(A10));
    Matrix A2(std::move(A20));

    Matrix ret;
    details::eval_FUNC_BIN::make(A1, A2, ret, f, f);
    return ret;
};


template<class binary_function, class test_function>
Matrix matcl::eval_binary_func(const Matrix& A10, const Matrix& A20, const binary_function& f,
                               const test_function& t)
{
    //make copy
    Matrix A1(A10);
    Matrix A2(A20);

    Matrix ret;
    details::eval_FUNC_BIN::make(A1, A2, ret, f, t);
    return ret;
};

template<class binary_function, class test_function>
Matrix matcl::eval_binary_func(Matrix&& A10, const Matrix& A20, const binary_function& f,
                               const test_function& t)
{
    //make copy
    Matrix A1(std::move(A10));
    Matrix A2(A20);

    Matrix ret;
    details::eval_FUNC_BIN::make(A1, A2, ret, f, t);
    return ret;
};

template<class binary_function, class test_function>
Matrix matcl::eval_binary_func(const Matrix& A10, Matrix&& A20, const binary_function& f,
                               const test_function& t)
{
    //make copy
    Matrix A1(A10);
    Matrix A2(std::move(A20));

    Matrix ret;
    details::eval_FUNC_BIN::make(A1, A2, ret, f, t);
    return ret;
};

template<class binary_function, class test_function>
Matrix matcl::eval_binary_func(Matrix&& A10, Matrix&& A20, const binary_function& f,
                               const test_function& t)
{
    //make copy
    Matrix A1(std::move(A10));
    Matrix A2(std::move(A20));

    Matrix ret;
    details::eval_FUNC_BIN::make(A1, A2, ret, f, t);
    return ret;
};

};
