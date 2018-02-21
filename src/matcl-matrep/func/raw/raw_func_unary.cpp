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

#include "matcl-matrep/func/raw/raw_func_unary.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-matrep/func/raw/eval_functor.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl { namespace raw { namespace details
{

template<class ret_type,class M,bool is_compl>
struct eval_real_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        ret = matcl::Matrix(m,false);
    };
};

template<class ret_type,class M>
struct eval_real_helper_impl<ret_type,M,true>
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        eval_functor_impl<ret_type,M>::eval(ret, m,details::real_helper<typename M::value_type>());

        struct_flag so  = md::predefined_struct::get_real(m.get_struct());
        so.add(ret.get_struct());
        ret.set_struct(so);
    };
};

template<class ret_type,class M,bool is_compl>
struct eval_imag_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        using val_type      = typename M::value_type;
        using sparse_mat    = Matrix<val_type,struct_sparse>;        
        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;

        sparse_mat out(ti::get_return_ti<ret_ti_type>(mdyf::imag::eval(), ti::get_ti(m)),
                         m.rows(),m.cols());
        ret = matcl::Matrix(out,false);
    };
};

template<class ret_type,class M>
struct eval_imag_helper_impl<ret_type,M,true>
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        return eval_functor_impl<ret_type,M>::eval(ret, m,details::imag_helper<typename M::value_type>());
    };
};

template<class ret_type,class M,bool is_compl>
struct eval_conj_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        ret = matcl::Matrix(m, false);
    };
};

template<class ret_type,class M>
struct eval_conj_helper_impl<ret_type,M,true>
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        eval_functor_impl<ret_type,M>::eval(ret, m,details::conj_helper<typename M::value_type>());

        struct_flag so  = md::predefined_struct::get_conj(m.get_struct());
        so.add(ret.get_struct());
        ret.set_struct(so);
    };
};

template<class ret_type,class M,class value_type>
struct eval_abs_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        eval_functor_impl<ret_type,M>::eval(ret, m,details::abs_helper<typename M::value_type>());

        bool is_square  = m.rows() == m.cols();
        struct_flag so  = md::predefined_struct::get_abs(m.get_struct(), is_square);
        so.add(ret.get_struct());
        ret.set_struct(so);
    };
};

template<class ret_type,class M,class value_type>
struct eval_abs2_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        eval_functor_impl<ret_type,M>::eval(ret, m,details::abs2_helper<typename M::value_type>());

        bool is_square  = m.rows() == m.cols();
        struct_flag so  = md::predefined_struct::get_abs(m.get_struct(), is_square);
        so.add(ret.get_struct());
        ret.set_struct(so);
    };
};

template<class ret_type,class M>
struct eval_arg_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        return eval_functor_impl<ret_type,M>::eval(ret, m,details::arg_helper<typename M::value_type>());
    };
};

template<class ret_type,class M>
struct eval_eps_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        return eval_functor_impl<ret_type,M>::eval(ret, m,details::eps_helper<typename M::value_type>());
    };
};

template<class ret_type,class M,bool is_int>
struct eval_is_finite_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        return eval_functor_impl<ret_type,M>::eval(ret, m,details::isfinite_helper<typename M::value_type>());
    };
};

template<class ret_type,class M>
struct eval_is_finite_helper_impl<ret_type,M,true>
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        ret = matcl::Matrix(integer_dense(ti::ti_empty(),1,m.rows(),m.cols()), false);
    };
};

template<class ret_type,class M,bool is_int>
struct eval_is_inf_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        return eval_functor_impl<ret_type,M>::eval(ret, m,details::isinf_helper<typename M::value_type>());
    };
};

template<class ret_type,class M>
struct eval_is_inf_helper_impl<ret_type,M,true>
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        integer_sparse out(ti::ti_empty(),m.rows(),m.cols());
        ret = matcl::Matrix(out, false);
    };
};

template<class ret_type, class M, bool is_int>
struct eval_is_nan_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        return eval_functor_impl<ret_type,M>::eval(ret, m,details::isnan_helper<typename M::value_type>());
    };
};

template<class ret_type,class M>
struct eval_is_nan_helper_impl<ret_type,M,true>
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        integer_sparse out(ti::ti_empty(),m.rows(),m.cols());
        ret = matcl::Matrix(out, false);
    };
};

template<class ret_type,class M>
struct eval_is_regular_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        return eval_functor_impl<ret_type,M>::eval(ret, m,details::isregular_helper<typename M::value_type>());
    };
};

template<class ret_type,class M>
struct eval_is_normal_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        return eval_functor_impl<ret_type,M>::eval(ret, m,details::isnormal_helper<typename M::value_type>());
    };
};

template<class ret_type,class M>
struct eval_is_zero_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        return eval_functor_impl<ret_type,M>::eval(ret, m,details::iszero_helper<typename M::value_type>());
    };
};

template<class ret_type,class M>
struct eval_is_one_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        return eval_functor_impl<ret_type,M>::eval(ret, m,details::isone_helper<typename M::value_type>());
    };
};

template<class ret_type,class M,bool is_int>
struct eval_is_int_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        return eval_functor_impl<ret_type,M>::eval(ret, m,details::isint_helper<typename M::value_type>());
    };
};

template<class ret_type,class M>
struct eval_is_int_helper_impl<ret_type,M,true>
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        ret = matcl::Matrix(integer_dense(ti::ti_empty(),1,m.rows(),m.cols()), false);
    };
};

template<class ret_type,class M,bool is_int>
struct eval_is_real_helper_impl
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        return eval_functor_impl<ret_type,M>::eval(ret, m,details::isreal_helper<typename M::value_type>());
    };
};

template<class ret_type,class M>
struct eval_is_real_helper_impl<ret_type,M,true>
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        ret = matcl::Matrix(integer_dense(ti::ti_empty(),1,m.rows(),m.cols()), false);
    };
};

template<class ret_type,class in_type,class value_type>
struct eval_floor_helper_impl
{
    static void eval(matcl::Matrix& ret, const in_type& m)
    {
        return eval_functor_impl<ret_type,in_type>::eval(ret, m,details::floor_helper<value_type>());
    };
};

template<class ret_type,class in_type>
struct eval_floor_helper_impl<ret_type,in_type,Integer>
{
    static void eval(matcl::Matrix& ret, const in_type& mat)
    {
        ret = matcl::Matrix(mat,false);
    };
};

template<class ret_type,class in_type,class value_type>
struct eval_ceil_helper_impl
{
    static void eval(matcl::Matrix& ret, const in_type& m)
    {
        return eval_functor_impl<ret_type,in_type>::eval(ret, m,details::ceil_helper<value_type>());
    };
};

template<class ret_type,class in_type>
struct eval_ceil_helper_impl<ret_type,in_type,Integer>
{
    static void eval(matcl::Matrix& ret, const in_type& mat)
    {
        ret = matcl::Matrix(mat,false);
    };
};

template<class ret_type,class in_type,class value_type>
struct eval_round_helper_impl
{
    static void eval(matcl::Matrix& ret, const in_type& m)
    {
        eval_functor_impl<ret_type,in_type>::eval(ret, m,details::round_helper<value_type>());
    };
};

template<class ret_type,class in_type>
struct eval_round_helper_impl<ret_type,in_type,Integer>
{
    static void eval(matcl::Matrix& ret, const in_type& mat)
    {
        ret = matcl::Matrix(mat,false);
    };
};

template<class ret_type,class in_type,class value_type>
struct eval_trunc_helper_impl
{
    static void eval(matcl::Matrix& ret, const in_type& m)
    {
        return eval_functor_impl<ret_type,in_type>::eval(ret, m,details::trunc_helper<value_type>());
    };
};

template<class ret_type,class in_type>
struct eval_trunc_helper_impl<ret_type,in_type,Integer>
{
    static void eval(matcl::Matrix& ret, const in_type& mat)
    {
        ret = matcl::Matrix(mat,false);
    };
};

template<class ret_type,class in_type,class value_type>
struct eval_ifloor_helper_impl
{
    static void eval(matcl::Matrix& ret, const in_type& m)
    {
        return eval_functor_impl<ret_type,in_type>::eval(ret, m,details::ifloor_helper<value_type>());
    };
};

template<class ret_type,class in_type>
struct eval_ifloor_helper_impl<ret_type,in_type,Integer>
{
    static void eval(matcl::Matrix& ret, const in_type& mat)
    {
        ret = matcl::Matrix(mat,false);
    };
};

template<class ret_type,class in_type,class value_type>
struct eval_iceil_helper_impl
{
    static void eval(matcl::Matrix& ret, const in_type& m)
    {
        return eval_functor_impl<ret_type,in_type>::eval(ret, m,details::iceil_helper<value_type>());
    };
};

template<class ret_type,class in_type>
struct eval_iceil_helper_impl<ret_type,in_type,Integer>
{
    static void eval(matcl::Matrix& ret, const in_type& mat)
    {
        ret = matcl::Matrix(mat,false);
    };
};

template<class ret_type,class in_type,class value_type>
struct eval_iround_helper_impl
{
    static void eval(matcl::Matrix& ret, const in_type& m)
    {
        return eval_functor_impl<ret_type,in_type>::eval(ret, m,details::iround_helper<value_type>());
    };
};

template<class ret_type,class in_type>
struct eval_iround_helper_impl<ret_type,in_type,Integer>
{
    static void eval(matcl::Matrix& ret, const in_type& mat)
    {
        ret = matcl::Matrix(mat,false);
    };
};

template<class ret_type,class in_type,class value_type>
struct eval_itrunc_helper_impl
{
    static void eval(matcl::Matrix& ret, const in_type& m)
    {
        return eval_functor_impl<ret_type,in_type>::eval(ret, m,details::itrunc_helper<value_type>());
    };
};

template<class ret_type,class in_type>
struct eval_itrunc_helper_impl<ret_type,in_type,Integer>
{
    static void eval(matcl::Matrix& ret, const in_type& mat)
    {
        ret = matcl::Matrix(mat,false);
    };
};

template<class M>
void details::scalfunc_isa_helper<M>::eval_is_nan(matcl::Matrix& ret, const M& m)
{
    static const bool is_integer = std::is_same<typename M::value_type,Integer>::value;
    return eval_is_nan_helper_impl<ret_type_nan,M,is_integer>::eval(ret, m);
};

template<class M>
void details::scalfunc_isa_helper<M>::eval_is_inf(matcl::Matrix& ret, const M& m)
{
    static const bool is_integer = std::is_same<typename M::value_type,Integer>::value;
    return eval_is_inf_helper_impl<ret_type_inf,M,is_integer>::eval(ret, m);
};

template<class M>
void details::scalfunc_isa_helper<M>::eval_is_finite(matcl::Matrix& ret, const M& m)
{
    static const bool is_integer = std::is_same<typename M::value_type,Integer>::value;
    return eval_is_finite_helper_impl<ret_type_finite,M,is_integer>::eval(ret, m);
};

template<class M>
void details::scalfunc_isa_helper<M>::eval_is_regular(matcl::Matrix& ret, const M& m)
{
    return eval_is_regular_helper_impl<ret_type_nan,M>::eval(ret, m);
};

template<class M>
void details::scalfunc_isa_helper<M>::eval_is_normal(matcl::Matrix& ret, const M& m)
{
    return eval_is_normal_helper_impl<ret_type_nan,M>::eval(ret, m);
};

template<class M>
void details::scalfunc_isa_helper<M>::eval_is_int(matcl::Matrix& ret, const M& m)
{
    static const bool is_integer = std::is_same<typename M::value_type,Integer>::value;
    return eval_is_int_helper_impl<ret_type_nan,M,is_integer>::eval(ret, m);
};

template<class M>
void details::scalfunc_isa_helper<M>::eval_is_real(matcl::Matrix& ret, const M& m)
{
    static const bool is_integer = std::is_same<typename M::value_type,Integer>::value;
    return eval_is_real_helper_impl<ret_type_nan,M,is_integer>::eval(ret, m);
};

//-----------------------------------------------------------------------------
//                      scalar_func_helper
//-----------------------------------------------------------------------------
template<class MP, class Val = typename MP::value_type>
struct eval_c_func_helper
{
    static void eval_sqrt_c(matcl::Matrix&, const MP&){};
    static void eval_sqrt1pm1_c(matcl::Matrix&, const MP&){};
    static void eval_log_c(matcl::Matrix&, const MP&){};
    static void eval_log1p_c(matcl::Matrix&, const MP&){};
    static void eval_log2_c(matcl::Matrix&, const MP&){};
    static void eval_log10_c(matcl::Matrix&, const MP&){};
    static void eval_asin_c(matcl::Matrix&, const MP&){};
    static void eval_acos_c(matcl::Matrix&, const MP&){};
    static void eval_asec_c(matcl::Matrix&, const MP&){};
    static void eval_acsc_c(matcl::Matrix&, const MP&){};
    static void eval_acosh_c(matcl::Matrix&, const MP&){};
    static void eval_atanh_c(matcl::Matrix&, const MP&){};
    static void eval_acoth_c(matcl::Matrix&, const MP&){};
    static void eval_asech_c(matcl::Matrix&, const MP&){};
};

template<class MP>
struct eval_c_func_helper<MP, Object>
{
    static void eval_sqrt_c(matcl::Matrix& ret, const MP& m);
    static void eval_sqrt1pm1_c(matcl::Matrix& ret, const MP& m);
    static void eval_log_c(matcl::Matrix& ret, const MP& m);
    static void eval_log1p_c(matcl::Matrix& ret, const MP& m);
    static void eval_log2_c(matcl::Matrix& ret, const MP& m);
    static void eval_log10_c(matcl::Matrix& ret, const MP& m);
    static void eval_asin_c(matcl::Matrix& ret, const MP& m);
    static void eval_acos_c(matcl::Matrix& ret, const MP& m);
    static void eval_asec_c(matcl::Matrix& ret, const MP& m);
    static void eval_acsc_c(matcl::Matrix& ret, const MP& m);
    static void eval_acosh_c(matcl::Matrix& ret, const MP& m);
    static void eval_atanh_c(matcl::Matrix& ret, const MP& m);
    static void eval_acoth_c(matcl::Matrix& ret, const MP& m);
    static void eval_asech_c(matcl::Matrix& ret, const MP& m);
};

template<class MP>
void scalar_func_helper<MP>::eval_sqrt_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_sqrt_c(ret, m);
};

template<class MP>
void scalar_func_helper<MP>::eval_sqrt1pm1_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_sqrt1pm1_c(ret, m);
};

template<class MP>
void scalar_func_helper<MP>::eval_log_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_log_c(ret, m);
};

template<class MP>
void scalar_func_helper<MP>::eval_log1p_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_log1p_c(ret, m);
};

template<class MP>
void scalar_func_helper<MP>::eval_log2_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_log2_c(ret, m);
};

template<class MP>
void scalar_func_helper<MP>::eval_log10_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_log10_c(ret, m);
};

template<class MP>
void scalar_func_helper<MP>::eval_asin_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_asin_c(ret, m);
};

template<class MP>
void scalar_func_helper<MP>::eval_acos_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_acos_c(ret, m);
};

template<class MP>
void scalar_func_helper<MP>::eval_asec_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_asec_c(ret, m);
};

template<class MP>
void scalar_func_helper<MP>::eval_acsc_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_acsc_c(ret, m);
};

template<class MP>
void scalar_func_helper<MP>::eval_acosh_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_acosh_c(ret, m);
};

template<class MP>
void scalar_func_helper<MP>::eval_atanh_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_atanh_c(ret, m);
};

template<class MP>
void scalar_func_helper<MP>::eval_acoth_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_acoth_c(ret, m);
};

template<class MP>
void scalar_func_helper<MP>::eval_asech_c(matcl::Matrix& ret, const MP& m)
{
    eval_c_func_helper<MP>::eval_asech_c(ret, m);
};

template<class M>
void details::scalar_func_helper<M>::eval_sqrt(matcl::Matrix& ret, const M& m)
{
    return eval_functor_impl<ret_type_sqrt,M>::eval(ret, m,sqrt_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_sqrt_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_sqrt = typename scalar_func_helper<M>::ret_type_sqrt;
    return eval_functor_impl<ret_type_sqrt,M>::eval(ret, m,sqrt_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_sqrt1pm1(matcl::Matrix& ret, const M& m)
{
    return eval_functor_impl<ret_type_sqrt,M>::eval(ret, m,sqrt1pm1_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_sqrt1pm1_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_sqrt = typename scalar_func_helper<M>::ret_type_sqrt;
    return eval_functor_impl<ret_type_sqrt,M>::eval(ret, m,sqrt1pm1_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_log(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_log,M>::eval(ret, m, log_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_log_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_log = typename scalar_func_helper<M>::ret_type_log;
    eval_functor_impl<ret_type_log,M>::eval(ret, m, log_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_log1p(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_log,M>::eval(ret, m, log1p_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_log1p_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_log = typename scalar_func_helper<M>::ret_type_log;
    eval_functor_impl<ret_type_log,M>::eval(ret, m, log1p_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_log2(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_log2,M>::eval(ret, m,log2_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_log2_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_log2 = typename scalar_func_helper<M>::ret_type_log2;
    eval_functor_impl<ret_type_log2,M>::eval(ret, m,log2_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_log10(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_log10,M>::eval(ret, m, log10_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_log10_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_log10 = typename scalar_func_helper<M>::ret_type_log10;
    eval_functor_impl<ret_type_log10,M>::eval(ret, m, log10_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_exp(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_exp,M>::eval(ret, m,exp_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_sin(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_sin,M>::eval(ret, m, sin_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_cos(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_cos,M>::eval(ret, m, cos_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_tan(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_tan,M>::eval(ret, m, tan_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_cot(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_cot,M>::eval(ret, m, cot_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_sec(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_sec,M>::eval(ret, m,sec_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_csc(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_csc,M>::eval(ret, m, csc_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_asin(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_asin,M>::eval(ret, m, asin_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_asin_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_asin = typename scalar_func_helper<M>::ret_type_asin;
    eval_functor_impl<ret_type_asin,M>::eval(ret, m, asin_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_acos(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_acos,M>::eval(ret, m, acos_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_acos_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_acos = typename scalar_func_helper<M>::ret_type_acos;
    eval_functor_impl<ret_type_acos,M>::eval(ret, m, acos_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_atan(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_atan,M>::eval(ret, m, atan_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_acot(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_acot,M>::eval(ret, m, acot_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_asec(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_asec,M>::eval(ret, m, asec_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_asec_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_asec = typename scalar_func_helper<M>::ret_type_asec;
    eval_functor_impl<ret_type_asec,M>::eval(ret, m, asec_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_acsc(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_acsc,M>::eval(ret, m, acsc_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_acsc_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_acsc = typename scalar_func_helper<M>::ret_type_acsc;
    eval_functor_impl<ret_type_acsc,M>::eval(ret, m, acsc_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_sinh(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_sinh,M>::eval(ret, m, sinh_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_cosh(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_cosh,M>::eval(ret, m, cosh_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_tanh(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_tanh,M>::eval(ret, m,tanh_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_coth(matcl::Matrix& ret, const M& m)
{
    return eval_functor_impl<ret_type_coth,M>::eval(ret, m,coth_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_sech(matcl::Matrix& ret, const M& m)
{
    return eval_functor_impl<ret_type_sech,M>::eval(ret, m, sech_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_csch(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_csch,M>::eval(ret, m, csch_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_asinh(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_asinh,M>::eval(ret, m, asinh_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_acosh(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_acosh,M>::eval(ret, m, acosh_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_acosh_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_acosh = typename scalar_func_helper<M>::ret_type_acosh;
    eval_functor_impl<ret_type_acosh,M>::eval(ret, m, acosh_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_atanh(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_atanh,M>::eval(ret, m, atanh_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_atanh_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_atanh = typename scalar_func_helper<M>::ret_type_atanh;
    eval_functor_impl<ret_type_atanh,M>::eval(ret, m, atanh_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_acoth(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_acoth,M>::eval(ret, m, acoth_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_acoth_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_acoth = typename scalar_func_helper<M>::ret_type_acoth;
    eval_functor_impl<ret_type_acoth,M>::eval(ret, m, acoth_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_asech(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_asech,M>::eval(ret, m, asech_helper<typename M::value_type>());
};

template<class M>
void eval_c_func_helper<M, Object>::eval_asech_c(matcl::Matrix& ret, const M& m)
{
    using ret_type_asech = typename scalar_func_helper<M>::ret_type_asech;
    eval_functor_impl<ret_type_asech,M>::eval(ret, m, asech_c_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_acsch(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_acsch,M>::eval(ret, m, acsch_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_floor(matcl::Matrix& ret, const M& m)
{
    eval_floor_helper_impl<ret_type_floor,M, typename M::value_type>::eval(ret, m);
};

template<class M>
void details::scalar_func_helper<M>::eval_ceil(matcl::Matrix& ret, const M& m)
{
    eval_ceil_helper_impl<ret_type_ceil,M, typename M::value_type>::eval(ret, m);
};

template<class M>
void details::scalar_func_helper<M>::eval_round(matcl::Matrix& ret, const M& m)
{
    eval_round_helper_impl<ret_type_round,M, typename M::value_type>::eval(ret, m);
};

template<class M>
void details::scalar_func_helper<M>::eval_trunc(matcl::Matrix& ret, const M& m)
{
    eval_trunc_helper_impl<ret_type_trunc,M, typename M::value_type>::eval(ret, m);
};

template<class M>
void details::scalar_func_helper<M>::eval_ifloor(matcl::Matrix& ret, const M& m)
{
    eval_ifloor_helper_impl<ret_type_ifloor,M, typename M::value_type>::eval(ret, m);
};

template<class M>
void details::scalar_func_helper<M>::eval_iceil(matcl::Matrix& ret, const M& m)
{
    eval_iceil_helper_impl<ret_type_iceil,M, typename M::value_type>::eval(ret, m);
};

template<class M>
void details::scalar_func_helper<M>::eval_iround(matcl::Matrix& ret, const M& m)
{
    eval_iround_helper_impl<ret_type_iround,M, typename M::value_type>::eval(ret, m);
};

template<class M>
void details::scalar_func_helper<M>::eval_itrunc(matcl::Matrix& ret, const M& m)
{
    return eval_itrunc_helper_impl<ret_type_itrunc,M, typename M::value_type>::eval(ret, m);
};

template<class M>
void details::scalar_func_helper<M>::eval_sign(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_sign,M>::eval(ret, m,sign_helper<typename M::value_type>());
};

template<class M>
void details::scalar_func_helper<M>::eval_isign(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_isign,M>::eval(ret, m,isign_helper<typename M::value_type>());
};

template<class M>
void details::unary_helper_impl<M>::eval_minus(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_minus,M>::eval(ret, m,uminus_helper<typename M::value_type>());
    struct_flag sf = md::predefined_struct::uminus_cont(ret.get_struct(), m.get_struct());
    ret.set_struct(sf);
};

template<class M>
void details::unary_helper_impl<M>::eval_inv(matcl::Matrix& ret, const M& m)
{
    eval_functor_impl<ret_type_inv,M>::eval(ret, m,invs_helper<typename M::value_type>());
    struct_flag sf = md::predefined_struct::inv_cont(ret.get_struct(), m.get_struct());
    ret.set_struct(sf);
};

template<class M>
void details::unary_helper_impl<M>::eval_neg(matcl::Matrix& ret, const M& m)
{
    return eval_functor_impl<ret_type_neg,M>::eval(ret, m,op_neg_helper<typename M::value_type>());
};

template<class M>
void details::unary_helper_impl<M>::eval_is_true(matcl::Matrix& ret, const M& m)
{
    return eval_functor_impl<ret_type_is_true,M>::eval(ret, m,op_true_helper<typename M::value_type>());
};

template<class M>
void details::scalfunc_real_helper<M>::eval_real(matcl::Matrix& ret, const M& m)
{
    static const bool is_complex = md::is_complex<typename M::value_type>::value
                                    ||std::is_same<typename M::value_type,Object>::value;
    return eval_real_helper_impl<ret_type,M,is_complex>::eval(ret, m);
};

template<class M>
void details::scalfunc_real_helper<M>::eval_imag(matcl::Matrix& ret, const M& m)
{
    static const bool is_complex = md::is_complex<typename M::value_type>::value
                                    ||std::is_same<typename M::value_type,Object>::value;
    return eval_imag_helper_impl<ret_type_imag,M,is_complex>::eval(ret, m);
};

template<class M>
void details::scalfunc_real_helper<M>::eval_conj(matcl::Matrix& ret, const M& m)
{
    static const bool is_compl_or_obj = md::is_complex<typename M::value_type>::value
                                    ||std::is_same<typename M::value_type,Object>::value;
    return eval_conj_helper_impl<ret_type_conj,M,is_compl_or_obj>::eval(ret, m);
};

template<class M>
void details::scalfunc_real_helper<M>::eval_abs(matcl::Matrix& ret, const M& m)
{
    return eval_abs_helper_impl<ret_type,M, typename M::value_type>::eval(ret, m);
};

template<class M>
void details::scalfunc_real_helper<M>::eval_abs2(matcl::Matrix& ret, const M& m)
{
    return eval_abs2_helper_impl<ret_type,M, typename M::value_type>::eval(ret, m);
};

template<class M>
void details::scalfunc_real_helper<M>::eval_arg(matcl::Matrix& ret, const M& m)
{
    return eval_arg_helper_impl<ret_type_arg,M>::eval(ret, m);
};

template<class M>
void details::scalfunc_real_helper<M>::eval_eps(matcl::Matrix& ret, const M& m)
{
    return eval_eps_helper_impl<ret_type_arg,M>::eval(ret, m);
};

};};};

MACRO_INSTANTIATE_G_1(matcl::raw::details::scalfunc_isa_helper)
MACRO_INSTANTIATE_S_1(matcl::raw::details::scalfunc_isa_helper)

MACRO_INSTANTIATE_G_1(matcl::raw::details::scalar_func_helper)
MACRO_INSTANTIATE_S_1(matcl::raw::details::scalar_func_helper)

MACRO_INSTANTIATE_G_1(matcl::raw::details::unary_helper_impl)
MACRO_INSTANTIATE_S_1(matcl::raw::details::unary_helper_impl)

MACRO_INSTANTIATE_G_1(matcl::raw::details::scalfunc_real_helper)
MACRO_INSTANTIATE_S_1(matcl::raw::details::scalfunc_real_helper)
