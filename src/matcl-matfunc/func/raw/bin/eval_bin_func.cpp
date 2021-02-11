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

#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-matfunc/func/bin/op_info.h"
#include "matcl-matfunc/func/raw/bin/eval_bin_functor.h"
#include "matcl-matfunc/func/raw/bin/raw_func_helpers.h"

namespace matcl { namespace details
{

namespace mr = matcl::raw;

template<class Ret, class V1, class V2, bool Rev>
class virtual_functor2 : public virtual_functor_base<Ret, V1>
{
    private:
        using op_virtual = op_virtual_base<Ret, V1, V2>;

    private:
        const op_virtual&   m_op;
        const V2&           m_scal;

    public:
        virtual_functor2(const op_virtual& op, const V2& scal)
            :m_op(op), m_scal(scal)
        {};

        virtual Ret eval(const V1& v) const override        { return m_op.eval(v,m_scal); };        
        virtual Ret test_zero(const V1& v) const override   { return m_op.test(v,m_scal); };
};

template<class Ret, class V1, class V2>
class virtual_functor2<Ret, V1, V2, true> : public virtual_functor_base<Ret, V2>
{
    private:
        using op_virtual = op_virtual_base<Ret, V1, V2>;

    private:
        const op_virtual&   m_op;
        const V1&           m_scal;

    public:
        virtual_functor2(const op_virtual& op, const V1& scal)
            :m_op(op), m_scal(scal)
        {};

        virtual Ret eval(const V2& v) const override        { return m_op.eval(m_scal,v); };
        virtual Ret test_zero(const V2& v) const override   { return m_op.test(m_scal,v); };
};

template<class Ret, class M1, class M2>
void function_eval_bin_func<Ret, M1, M2>::eval(matcl::Matrix& ret, const M1& mat1, const M2& mat2, 
                                               const func& op)
{
    using value_type_1      = typename M1::value_type;
    using value_type_2      = typename M2::value_type;
    using UT                = typename md::unify_types<value_type_1,value_type_2>::type;

    typename ti::get_ti_type<M1>::type ti_1 = ti::get_ti(mat1);
    typename ti::get_ti_type<M2>::type ti_2 = ti::get_ti(mat2);

    if (mat1.rows() == 1 && mat1.cols() == 1)
    {			
        value_type_1 val(mat1(1,1));

        if (mat2.rows() == 1 && mat2.cols() == 1)
        {			
            value_type_2 val2(mat2(1,1));
            ret = op.eval(val, val2);
            return;
        };

        using functor   = function_eval_func<Ret, M2>;
        using virt_func = virtual_functor2<Ret, In1, In2, true>;

        virt_func vf(op,val);
        return functor::eval(ret, mat2, vf);
    };

    if (mat2.rows() == 1 && mat2.cols() == 1)
    {			
        value_type_2 val(mat2(1,1));

        using functor   = function_eval_func<Ret, M1>;
        using virt_func = virtual_functor2<Ret, In1, In2, false>;

        virt_func vf(op,val);
        return functor::eval(ret, mat1, vf);
    };

    Integer m1 = mat1.rows();
    Integer n1 = mat1.cols();
    Integer m2 = mat2.rows();
    Integer n2 = mat2.cols();
    error::check_eeop(m1,n1,m2,n2);

    using ret_type  = matcl::raw::Matrix<Ret,struct_dense>;
    using str_ret   = typename ret_type::struct_type;
    using str_1     = typename M1::struct_type;
    using str_2     = typename M2::struct_type;

    mrd::eval_bin_functor_impl<ret_type,func,M1,M2,false,false,false,false,
            str_ret, str_1, str_2>::eval(ret, mat1, mat2, op);

    return;
};

#define MACRO_EXPAND_SCALAR(scal)                                                           \
template struct function_eval_bin_func<scal, mr::integer_dense, mr::integer_dense>;         \
template struct function_eval_bin_func<scal, mr::integer_dense, mr::float_dense>;           \
template struct function_eval_bin_func<scal, mr::integer_dense, mr::real_dense>;            \
template struct function_eval_bin_func<scal, mr::integer_dense, mr::complex_dense>;         \
template struct function_eval_bin_func<scal, mr::integer_dense, mr::float_complex_dense>;   \
template struct function_eval_bin_func<scal, mr::float_dense, mr::integer_dense>;           \
template struct function_eval_bin_func<scal, mr::float_dense, mr::float_dense>;             \
template struct function_eval_bin_func<scal, mr::float_dense, mr::real_dense>;              \
template struct function_eval_bin_func<scal, mr::float_dense, mr::complex_dense>;           \
template struct function_eval_bin_func<scal, mr::float_dense, mr::float_complex_dense>;     \
template struct function_eval_bin_func<scal, mr::real_dense, mr::integer_dense>;            \
template struct function_eval_bin_func<scal, mr::real_dense, mr::float_dense>;              \
template struct function_eval_bin_func<scal, mr::real_dense, mr::real_dense>;               \
template struct function_eval_bin_func<scal, mr::real_dense, mr::complex_dense>;            \
template struct function_eval_bin_func<scal, mr::real_dense, mr::float_complex_dense>;      \
template struct function_eval_bin_func<scal, mr::complex_dense, mr::integer_dense>;         \
template struct function_eval_bin_func<scal, mr::complex_dense, mr::float_dense>;           \
template struct function_eval_bin_func<scal, mr::complex_dense, mr::real_dense>;            \
template struct function_eval_bin_func<scal, mr::complex_dense, mr::complex_dense>;         \
template struct function_eval_bin_func<scal, mr::complex_dense, mr::float_complex_dense>;   \
template struct function_eval_bin_func<scal, mr::float_complex_dense, mr::integer_dense>;   \
template struct function_eval_bin_func<scal, mr::float_complex_dense, mr::float_dense>;     \
template struct function_eval_bin_func<scal, mr::float_complex_dense, mr::real_dense>;      \
template struct function_eval_bin_func<scal, mr::float_complex_dense, mr::complex_dense>;   \
template struct function_eval_bin_func<scal, mr::float_complex_dense, mr::float_complex_dense>;\
template struct function_eval_bin_func<scal, mr::object_dense, mr::object_dense>;           \


MACRO_EXPAND_SCALAR(Integer)
MACRO_EXPAND_SCALAR(Float)
MACRO_EXPAND_SCALAR(Real)
MACRO_EXPAND_SCALAR(Float_complex)
MACRO_EXPAND_SCALAR(Complex)
MACRO_EXPAND_SCALAR(Object)

};};
