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

#include "matcl-matrep/details/mpl.h"
#include "matcl-matrep/details/isa.h"
#include "matcl-matrep/details/utils.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"

#include "matcl-internals/func/raw_func_unary.h"

namespace matcl { namespace raw { namespace details
{

namespace md = matcl::details;

template<class MP, bool is_f_0>
struct func_ret_type
{
    using value_type    = typename MP::value_type;
    using struct_type   = typename MP::struct_type;

    using ret_val
        = typename matcl::details::select_if
        <
            std::is_same<value_type,Integer>::value,
            Real,
            value_type
        >::type;

    using ret_str
        = typename matcl::details::select_if
        <
            !is_f_0,
            struct_dense,
            struct_type
        >::type;

    using type = raw::Matrix<ret_val,ret_str>;
};

template<class M>
struct scalfunc_isa_helper
{
    using struct_type
        = typename matcl::details::select_if
        <
            std::is_same<typename M::value_type,Integer>::value,
            struct_sparse,
            typename M::struct_type
        >::type;

    using int_or_obj         
        = typename matcl::details::select_if
        <
            std::is_same<typename M::value_type,Object>::value,
            Object,
            Integer
        >::type;

    using ret_type_nan      = Matrix<int_or_obj,struct_type>;
    using ret_type_inf      = Matrix<int_or_obj,struct_type>;
    using ret_type_finite   = Matrix<int_or_obj,struct_dense>;

    // inplace is allowed; refcount must be increased for 
    // nontemporary objects; TODO
    static void eval_is_nan(matcl::Matrix& ret, const M& m);
    static void eval_is_inf(matcl::Matrix& ret, const M& m);
    static void eval_is_finite(matcl::Matrix& ret, const M& m);
    static void eval_is_regular(matcl::Matrix& ret, const M& m);
    static void eval_is_normal(matcl::Matrix& ret, const M& m);
    static void eval_is_int(matcl::Matrix& ret, const M& m);
    static void eval_is_real(matcl::Matrix& ret, const M& m);
};

template<class MP>
struct MATCL_MATREP_EXPORT unary_helper_impl
{	
    using str_type          = typename MP::struct_type;
    using val_type          = typename MP::value_type;
    using ret_type_minus    = MP;
    using ret_type_inv      = Matrix<typename md::unify_types<val_type, Float>::type, struct_dense>;    

    using str_neg
        = typename matcl::details::select_if
        < 
            std::is_same<str_type,struct_dense>::value,
            struct_sparse,
            struct_dense
        >::type;

    using int_or_obj
        = typename matcl::details::select_if
        < 
            std::is_same<val_type,Object>::value,
            Object,
            Integer
        >::type;

    using ret_type_neg      = Matrix<int_or_obj,str_neg>;
    using ret_type_is_true  = Matrix<int_or_obj,str_type>;

    // inplace is allowed; refcount must be increased for 
    // nontemporary objects; TODO
    static void eval_minus(matcl::Matrix& ret, const MP& m);
    static void eval_inv(matcl::Matrix& ret, const MP& m);
    static void eval_neg(matcl::Matrix& ret, const MP& m);
    static void eval_is_true(matcl::Matrix& ret, const MP& m);
};

template<class MP>
struct scalar_func_helper
{
    using struct_type       = typename MP::struct_type;
    using value_type        = typename MP::value_type;
    using integer_matrix    = typename md::select_if
                            <
                                md::is_object<value_type>::value,
                                Matrix<Object,struct_type>,
                                Matrix<Integer,struct_type>
                            > :: type;

    using ret_type_floor    = MP;
    using ret_type_ceil     = MP;
    using ret_type_round    = MP;
    using ret_type_trunc    = MP;
    using ret_type_sign     = MP;

    using ret_type_ifloor   = integer_matrix;
    using ret_type_iceil    = integer_matrix;
    using ret_type_iround   = integer_matrix;
    using ret_type_itrunc   = integer_matrix;
    using ret_type_isign    = integer_matrix;

    using ret_type_sqrt     = typename func_ret_type<MP,true>::type;
    using ret_type_log      = typename func_ret_type<MP,false>::type;
    using ret_type_log2     = typename func_ret_type<MP,false>::type;
    using ret_type_log10    = typename func_ret_type<MP,false>::type;
    using ret_type_exp      = typename func_ret_type<MP,false>::type;

    using ret_type_sin      = typename func_ret_type<MP,true>::type;
    using ret_type_cos      = typename func_ret_type<MP,false>::type;
    using ret_type_tan      = typename func_ret_type<MP,true>::type;
    using ret_type_cot      = typename func_ret_type<MP,false>::type;
    using ret_type_sec      = typename func_ret_type<MP,false>::type;
    using ret_type_csc      = typename func_ret_type<MP,false>::type;

    using ret_type_asin     = typename func_ret_type<MP,true>::type;
    using ret_type_acos     = typename func_ret_type<MP,false>::type;
    using ret_type_atan     = typename func_ret_type<MP,true>::type;
    using ret_type_acot     = typename func_ret_type<MP,false>::type;
    using ret_type_asec     = typename func_ret_type<MP,false>::type;
    using ret_type_acsc     = typename func_ret_type<MP,false>::type;

    using ret_type_sinh     = typename func_ret_type<MP,true>::type;
    using ret_type_cosh     = typename func_ret_type<MP,false>::type;
    using ret_type_tanh     = typename func_ret_type<MP,true>::type;
    using ret_type_coth     = typename func_ret_type<MP,false>::type;
    using ret_type_sech     = typename func_ret_type<MP,false>::type;
    using ret_type_csch     = typename func_ret_type<MP,false>::type;

    using ret_type_asinh    = typename func_ret_type<MP,true>::type;
    using ret_type_acosh    = typename func_ret_type<MP,false>::type;
    using ret_type_atanh    = typename func_ret_type<MP,true>::type;
    using ret_type_acoth    = typename func_ret_type<MP,false>::type;
    using ret_type_asech    = typename func_ret_type<MP,false>::type;
    using ret_type_acsch    = typename func_ret_type<MP,false>::type;

    // inplace is allowed; refcount must be increased for 
    // nontemporary objects; TODO

    // these functions can be called only on complex matrices    
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

    static void eval_sqrt(matcl::Matrix& ret, const MP& m);    
    static void eval_sqrt1pm1(matcl::Matrix& ret, const MP& m);    
    static void eval_log(matcl::Matrix& ret, const MP& m);    
    static void eval_log1p(matcl::Matrix& ret, const MP& m);   
    static void eval_log2(matcl::Matrix& ret, const MP& m);    
    static void eval_log10(matcl::Matrix& ret, const MP& m);    
    static void eval_exp(matcl::Matrix& ret, const MP& m);

    static void eval_floor(matcl::Matrix& ret, const MP& m);
    static void eval_ceil(matcl::Matrix& ret, const MP& m);
    static void eval_round(matcl::Matrix& ret, const MP& m);
    static void eval_trunc(matcl::Matrix& ret, const MP& m);
    static void eval_sign(matcl::Matrix& ret, const MP& m);		

    static void eval_ifloor(matcl::Matrix& ret, const MP& m);
    static void eval_iceil(matcl::Matrix& ret, const MP& m);
    static void eval_iround(matcl::Matrix& ret, const MP& m);
    static void eval_itrunc(matcl::Matrix& ret, const MP& m);
    static void eval_isign(matcl::Matrix& ret, const MP& m);

    static void eval_sin(matcl::Matrix& ret, const MP& m);
    static void eval_cos(matcl::Matrix& ret, const MP& m);
    static void eval_tan(matcl::Matrix& ret, const MP& m);
    static void eval_cot(matcl::Matrix& ret, const MP& m);
    static void eval_sec(matcl::Matrix& ret, const MP& m);
    static void eval_csc(matcl::Matrix& ret, const MP& m);

    static void eval_asin(matcl::Matrix& ret, const MP& m);    
    static void eval_acos(matcl::Matrix& ret, const MP& m);    
    static void eval_atan(matcl::Matrix& ret, const MP& m);
    static void eval_acot(matcl::Matrix& ret, const MP& m);
    static void eval_asec(matcl::Matrix& ret, const MP& m);    
    static void eval_acsc(matcl::Matrix& ret, const MP& m);    

    static void eval_sinh(matcl::Matrix& ret, const MP& m);
    static void eval_cosh(matcl::Matrix& ret, const MP& m);
    static void eval_tanh(matcl::Matrix& ret, const MP& m);
    static void eval_coth(matcl::Matrix& ret, const MP& m);
    static void eval_sech(matcl::Matrix& ret, const MP& m);
    static void eval_csch(matcl::Matrix& ret, const MP& m);

    static void eval_asinh(matcl::Matrix& ret, const MP& m);
    static void eval_acosh(matcl::Matrix& ret, const MP& m);    
    static void eval_atanh(matcl::Matrix& ret, const MP& m);    
    static void eval_acoth(matcl::Matrix& ret, const MP& m);    
    static void eval_asech(matcl::Matrix& ret, const MP& m);    
    static void eval_acsch(matcl::Matrix& ret, const MP& m);
};

};};};
