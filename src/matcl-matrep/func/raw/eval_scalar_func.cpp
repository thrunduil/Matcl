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

#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/base/utils.h"
#include "matcl-matrep/func/raw/eval_functor.h"

namespace matcl { namespace details
{

namespace mr = matcl::raw;

template<class struct_1,bool is_f_0>
struct ret_struct_type
{
    using type  = typename details::select_if
                <
                    !is_f_0,
                    struct_dense,
                    struct_1
                >::type;
};

template<class Ret, class Mat>
void function_eval_func<Ret,Mat>::eval(matcl::Matrix& result, const Mat& A, 
                                       const virtual_functor_base<Ret, In>& vf)
{
    using val_type      = typename Mat::value_type;
    using struct_type   = typename Mat::struct_type;
    using ti_type       = typename ti::get_ti_type<Mat>::type;

    val_type Z          = default_value<val_type>(ti::get_ti(A));        
    ti_type in_ti       = ti::get_ti(A);

    Ret ret             = vf.test_zero(Z);

    if (mrd::is_zero(ret) == true)
    {
        using struct_type_r = typename ret_struct_type<struct_type,true>::type;
        using ret_type      = raw::Matrix<Ret,struct_type_r>;
        using functor_type  = virtual_functor_base<Ret, In>;   

        raw::eval_functor_2_main<ret_type,Mat>::eval(result, ti::get_ti(ret), A, vf);
        return;
    }
    else
    {
        using struct_type_r = typename ret_struct_type<struct_type,false>::type;
        using ret_type      = raw::Matrix<Ret,struct_type_r>;
        using functor_type  = virtual_functor_base<Ret, In>;   

        raw::eval_functor_2_main<ret_type,Mat>::eval(result, ti::get_ti(ret), A, vf);
        return;
    };
};

#define MACRO_EXPAND_SCALAR(scal)                                   \
template struct function_eval_func<scal, mr::integer_dense>;        \
template struct function_eval_func<scal, mr::float_dense>;          \
template struct function_eval_func<scal, mr::real_dense>;           \
template struct function_eval_func<scal, mr::float_complex_dense>;  \
template struct function_eval_func<scal, mr::complex_dense>;        \
template struct function_eval_func<scal, mr::object_dense>;         \
template struct function_eval_func<scal, mr::integer_sparse>;       \
template struct function_eval_func<scal, mr::float_sparse>;         \
template struct function_eval_func<scal, mr::real_sparse>;          \
template struct function_eval_func<scal, mr::float_complex_sparse>; \
template struct function_eval_func<scal, mr::complex_sparse>;       \
template struct function_eval_func<scal, mr::object_sparse>;        \
template struct function_eval_func<scal, mr::integer_band>;         \
template struct function_eval_func<scal, mr::float_band>;           \
template struct function_eval_func<scal, mr::real_band>;            \
template struct function_eval_func<scal, mr::float_complex_band>;   \
template struct function_eval_func<scal, mr::complex_band>;         \
template struct function_eval_func<scal, mr::object_band>;

MACRO_EXPAND_SCALAR(Integer)
MACRO_EXPAND_SCALAR(Float)
MACRO_EXPAND_SCALAR(Real)
MACRO_EXPAND_SCALAR(Float_complex)
MACRO_EXPAND_SCALAR(Complex)
MACRO_EXPAND_SCALAR(Object)

};};
