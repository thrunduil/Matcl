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

#include "matcl-matfunc/func/raw/bin/raw_func_atan2.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-matfunc/func/raw/bin/eval_op.h"
#include "matcl-internals/func/raw_func_op_helpers.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-matrep/details/matrix.inl"

namespace matcl { namespace raw { namespace details
{

namespace mdyf  = matcl::dynamic::functions;

template<class arg_type, bool rev>
class atan2_functor
{
    private:
        arg_type	arg;

    public:
        atan2_functor(const arg_type& arg) : arg(arg) {};

        bool is_special_case() const		{ return false; };

        template<class ret_type, class in_type>
        void eval_special_case(matcl::Matrix&, typename ti::get_ti_type<ret_type>::type,
                                   const in_type &) const
        {};

        template<class val_type>
        typename md::real_type<typename md::unify_types2<Float,val_type,arg_type>::type>::type
        eval(const val_type& val) const
        {
            using ret_type0 = typename md::unify_types2<Float,val_type,arg_type>::type;
            using ret_type  = typename md::real_type<ret_type0>::type;
            return eval_impl<ret_type,val_type,rev>::eval(val,arg);
        };

        template<class ret_type, class val_type, bool rev_>
        struct eval_impl
        {
            static ret_type eval(const val_type& val, const arg_type& arg)
            {
                return atan2_helper<val_type,arg_type>::eval(val,arg);
            };
        };

        template<class ret_type, class val_type>
        struct eval_impl<ret_type,val_type,true>
        {
            static ret_type eval(const val_type& val, const arg_type& arg)
            {
                return atan2_helper<arg_type,val_type>::eval(arg,val);
            };
        };
};

struct eval_atan2_functor
{
    static const bool ZZ = true;
    static const bool ZN = false;
    static const bool NZ = false;

    template<class ret, class T1, class T2>
    static ret eval(const T1& arg1, const T2& arg2)
    {
        return atan2_helper<T1,T2>::eval(arg1,arg2);
    };

    template<class ret, class T1, class TZ, bool zero_on_right>
    static ret eval_zero(ti::ti_type<ret> ti_ret, const T1& arg1, ti::ti_type<TZ> ti_z)
    {
        (void)ti_ret;
        TZ zero = md::default_value<TZ>(ti_z);
        return zero_on_right? atan2_helper<T1,TZ>::eval(arg1,zero)
                              :atan2_helper<TZ,T1>::eval(zero,arg1);
    };

    template<bool zero_on_right>
    using is_eval_zero_id   = std::false_type;

    template<class ret, bool is_inv, class TI1, class TI2>
    static ti::ti_type<ret> return_type(TI1 t1, TI2 t2)
    {
        if (md::runtime_value<bool, is_inv>::eval() == false)
            return ti::get_return_ti<ti::ti_type<ret>>(mdyf::atan2::eval(),t1,t2);
        else
            return ti::get_return_ti<ti::ti_type<ret>>(mdyf::atan2::eval(),t2,t1);
    };

    template<class T, bool is_rev>
    static atan2_functor<T,is_rev> get_scalar_functor(const T& val)
    {
        return atan2_functor<T,is_rev>(val);
    };

    static struct_flag op_struct(struct_flag, struct_flag, bool, bool, bool)
    {
        return struct_flag();
    };
};

template<class M1,class M2>
void atan2_helper_mat_mat_inpl<M1,M2>::eval(matcl::Matrix& ret, const M1& A, const M2& B)
{
    return eval_op_impl<ret_type_atan2,M1,M2,eval_atan2_functor>::eval(ret,A,B);
};

template<class M1,class M2>
void atan2_helper_scal_mat_inpl<M1,M2>::eval(matcl::Matrix& ret, const M1& A, const M2& B)
{
    return eval_scalar_impl<ret_type_atan2,M2,eval_atan2_functor,false>::eval(ret,B,true,A);
};

template<class M2,class M1>
void atan2_helper_mat_scal_inpl<M2,M1>::eval(matcl::Matrix& ret, const M2& B, const M1& A)
{
    return eval_scalar_impl<ret_type_atan2,M2,eval_atan2_functor,false>::eval(ret,B,false,A);
};

}}}

MACRO_INSTANTIATE_SG_DIAG_DENSE_2(matcl::raw::details::atan2_helper_scal_mat_inpl)
MACRO_INSTANTIATE_GS_DIAG_DENSE_2(matcl::raw::details::atan2_helper_mat_scal_inpl)
MACRO_INSTANTIATE_BIN_DIAG_DENSE(matcl::raw::details::atan2_helper_mat_mat_inpl)
