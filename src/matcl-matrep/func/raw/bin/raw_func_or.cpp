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

#include "matcl-matrep/func/raw/bin/raw_func_or.h"
#include "matcl-matrep/base/instantiate.h"
#include "matcl-matrep/func/raw/bin/eval_op.h"
#include "matcl-matrep/func/raw/bin/raw_func_op_helpers.h"
#include "matcl-scalar/details/matfunc_helpers.h"

namespace matcl { namespace raw { namespace details
{

template<class arg_type, bool rev>
class or_functor
{
    private:
        arg_type	arg;

    public:
        or_functor(const arg_type& arg) : arg(arg) {};

        bool is_special_case() const		
        { 
            return cast_bool_helper<arg_type>::eval(arg); 
        };

        template<class ret_type, class in_type>
        void eval_special_case(matcl::Matrix& ret, typename ti::get_ti_type<ret_type>::type ret_ti, 
                                   const in_type &x) const
        {
            ret = create_matrix<Integer>(ret_ti,1,x.rows(), x.cols());
        };

        template<class val_type>
        Integer eval(const val_type& val) const
        {
            return cast_bool_helper<val_type>::eval(val);
        }
};

template<bool rev>
class or_functor<Object,rev>
{
    private:
        using arg_type = Object;

    private:
        arg_type	arg;

    public:
        or_functor(const arg_type& arg) : arg(arg) {};

        bool is_special_case() const		
        { 
            return false;
        };

        template<class ret_type, class in_type>
        void eval_special_case(matcl::Matrix&, typename ti::get_ti_type<ret_type>::type, 
                                   const in_type &) const
        {};

        template<class val_type>
        arg_type eval(const val_type& val) const
        {
            return rev == false ? elem_or_helper<val_type,arg_type>::eval(val,arg)
                                : elem_or_helper<arg_type, val_type>::eval(arg,val);
        };
};

template<class T, bool zero_on_right>
struct eval_zero_or
{
    static T eval(T arg1, ti::ti_type<T> ti_z)
    {
        (void)ti_z;
        return cast_bool_helper<T>::eval(arg1);
    }
};

template<bool zero_on_right>
struct eval_zero_or<Object,zero_on_right>
{
    using T = Object;
    static T eval(T arg1, ti::ti_type<T> ti_z)
    {
        T zero = md::default_value<T>(ti_z);
        return zero_on_right? elem_or_helper<T,T>::eval(arg1,zero)
                              :elem_or_helper<T,T>::eval(zero,arg1);
    }
};

struct eval_or_functor
{
    static const bool ZZ = true;
    static const bool ZN = false;
    static const bool NZ = false;

    template<class ret, class T1, class T2>
    static ret eval(const T1& arg1, const T2& arg2)
    {
        return mrd::elem_or_helper<T1,T2>::eval(arg1, arg2);
    };

    template<class ret, class T1, class TZ, bool zero_on_right>
    static ret eval_zero(ti::ti_type<ret> ti_ret, const T1& arg1, ti::ti_type<TZ> ti_z)
    {
        (void)ti_ret;
        return eval_zero_or<ret,zero_on_right>::eval(arg1, ti_z);
    };

    template<bool zero_on_right>
    using is_eval_zero_id   = std::false_type;

    template<class ret, bool is_inv, class TI1, class TI2>
    static ti::ti_type<ret> return_type(TI1 t1, TI2 t2)
    {
        if (md::runtime_value<bool, is_inv>::eval() == false)
            return ti::get_return_ti<ti::ti_type<ret>>(mdyf::elem_and::eval(),t1,t2);
        else
            return ti::get_return_ti<ti::ti_type<ret>>(mdyf::elem_and::eval(),t2,t1);
    };

    template<class T, bool is_rev>
    static or_functor<T,is_rev> get_scalar_functor(const T& val)
    {
        return or_functor<T,is_rev>(val);
    };

    static struct_flag op_struct(struct_flag , struct_flag , bool, bool, bool)
    {
        return struct_flag();
    };    
};

template<class M1,class M2>
void or_helper_mat_mat_inpl<M1,M2>::eval(matcl::Matrix& ret, const M1& A, const M2& B)
{
    return eval_op_impl<ret_type_or,M1,M2,eval_or_functor>::eval(ret,A,B);
};

template<class M1,class M2>
void or_helper_scal_mat_inpl<M1,M2>::eval(matcl::Matrix& ret, const M1& A, const M2& B)
{
    return eval_scalar_impl<ret_type_or,M2,eval_or_functor,false>::eval(ret,B,true,A);
};

template<class M2,class M1>
void or_helper_mat_scal_inpl<M2,M1>::eval(matcl::Matrix& ret, const M2& B, const M1& A)
{
    return eval_scalar_impl<ret_type_or,M2,eval_or_functor,false>::eval(ret,B,false,A);
};

template struct or_helper_mat_mat_inpl<integer_dense, integer_dense>;
template struct or_helper_mat_mat_inpl<integer_dense, integer_sparse>;
template struct or_helper_mat_mat_inpl<integer_dense, integer_band>;
template struct or_helper_mat_mat_inpl<integer_sparse, integer_dense>;
template struct or_helper_mat_mat_inpl<integer_sparse, integer_sparse>;
template struct or_helper_mat_mat_inpl<integer_sparse, integer_band>;
template struct or_helper_mat_mat_inpl<integer_band, integer_dense>;
template struct or_helper_mat_mat_inpl<integer_band, integer_sparse>;
template struct or_helper_mat_mat_inpl<integer_band, integer_band>;

template struct or_helper_mat_scal_inpl<integer_dense, Integer>;
template struct or_helper_mat_scal_inpl<integer_band, Integer>;
template struct or_helper_mat_scal_inpl<integer_sparse, Integer>;

template struct or_helper_scal_mat_inpl<Integer, integer_dense>;
template struct or_helper_scal_mat_inpl<Integer, integer_sparse>;
template struct or_helper_scal_mat_inpl<Integer, integer_band>;

//
template struct or_helper_mat_mat_inpl<object_dense, object_dense>;
template struct or_helper_mat_mat_inpl<object_dense, object_sparse>;
template struct or_helper_mat_mat_inpl<object_dense, object_band>;
template struct or_helper_mat_mat_inpl<object_sparse, object_dense>;
template struct or_helper_mat_mat_inpl<object_sparse, object_sparse>;
template struct or_helper_mat_mat_inpl<object_sparse, object_band>;
template struct or_helper_mat_mat_inpl<object_band, object_dense>;
template struct or_helper_mat_mat_inpl<object_band, object_sparse>;
template struct or_helper_mat_mat_inpl<object_band, object_band>;

template struct or_helper_mat_scal_inpl<object_dense, Object>;
template struct or_helper_mat_scal_inpl<object_band, Object>;
template struct or_helper_mat_scal_inpl<object_sparse, Object>;

template struct or_helper_scal_mat_inpl<Object, object_dense>;
template struct or_helper_scal_mat_inpl<Object, object_sparse>;
template struct or_helper_scal_mat_inpl<Object, object_band>;

}}}
