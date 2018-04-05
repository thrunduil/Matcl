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

#include "matcl-matfunc/func/raw/bin/raw_func_pow.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-matfunc/func/raw/bin/eval_op.h"
#include "matcl-internals/func/raw_func_op_helpers.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-internals/func/converter.h"
#include "matcl-matrep/details/matrix.inl"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace raw { namespace details
{

template<class val_type, class arg_type, bool rev_, bool Is_compl>
struct pow_eval_impl
{
    using ret_type = typename pow_helper<val_type,arg_type>::ret;

    static ret_type eval(const val_type& val, const arg_type& arg)
    {
        return pow_helper<val_type,arg_type>::eval(val,arg);
    };
};

template<class val_type, class arg_type, bool rev_>
struct pow_eval_impl<val_type, arg_type, rev_, true>
{
    using ret_type = typename pow_c_helper<val_type,arg_type>::ret;

    static ret_type eval(const val_type& val, const arg_type& arg)
    {
        return pow_c_helper<val_type,arg_type>::eval(val,arg);
    };
};

template<class val_type, class arg_type>
struct pow_eval_impl<val_type,arg_type,true,false>
{
    using ret_type = typename pow_helper<arg_type,val_type>::ret;

    static ret_type eval(const val_type& val, const arg_type& arg)
    {
        return pow_helper<arg_type,val_type>::eval(arg,val);
    };
};

template<class val_type, class arg_type>
struct pow_eval_impl<val_type,arg_type,true,true>
{
    using ret_type = typename pow_c_helper<arg_type,val_type>::ret;

    static ret_type eval(const val_type& val, const arg_type& arg)
    {
        return pow_c_helper<arg_type,val_type>::eval(arg,val);
    };
};

template<class arg_type, bool rev, bool Is_complex>
class pow_functor
{
    private:
        arg_type	arg;

    public:
        pow_functor(const arg_type& arg) : arg(arg) {};

        bool is_special_case() const
        {
            return ( rev == false && (mrd::is_zero(arg) || mrd::is_one(arg))); 
        };

        template<class ret_type, class in_type>
        void eval_special_case(matcl::Matrix& ret, typename ti::get_ti_type<ret_type>::type ret_ti,
                                   const in_type &x) const
        {
            if (mrd::is_zero(arg))
            {
                using val_type  = typename ret_type::value_type;
                using val_real  = typename md::real_type<val_type>::type;
                val_real one    = md::one_value<val_real>(ret_ti);
                ret             = create_matrix<val_real>(ret_ti,one,x.rows(),x.cols());
                return;
            }
            else
            {
                using ret_val       = typename ret_type::value_type;
                using ret_val_real  = typename md::real_type<ret_val>::type;
                using in_val        = typename in_type::value_type;
            
                using struct_type   = typename in_type::struct_type;
                using value_type    = typename md::unify_types<in_val,ret_val_real>::type;
                using mat_type      = Matrix<value_type,struct_type>;

                ret = matcl::Matrix(converter<mat_type,in_type>::eval(x,ret_ti),true);
                return;
            };
        };

        template<class val_type>
        auto eval(const val_type& val) const 
            -> typename pow_eval_impl<val_type,arg_type,rev,Is_complex>::ret_type
        {            
            return pow_eval_impl<val_type,arg_type,rev,Is_complex>::eval(val,arg);
        };
};

template<bool rev, bool Is_complex>
class pow_functor<Object, rev, Is_complex>
{
    private:
        using arg_type  = Object;

    private:
        arg_type	arg;

    public:
        pow_functor(const arg_type& arg) : arg(arg) {};

        bool is_special_case() const
        {
            return false;
        };

        template<class ret_type, class in_type>
        void eval_special_case(matcl::Matrix&, typename ti::get_ti_type<ret_type>::type,
                                   const in_type &) const
        {};

        template<class val_type>
        Object eval(const val_type& val) const 
        {            
            return pow_eval_impl<val_type,arg_type,rev,Is_complex>::eval(val,arg);
        };
};

template<class Ret, class T, class TZ, bool zero_on_right>
struct eval_zero_pow
{
    static Ret eval(T arg1, ti::ti_type<TZ> ti_z)
    {
        TZ zero = md::default_value<TZ>(ti_z);

        if (zero_on_right)
        {
            auto tmp = pow_helper<T,TZ>::eval(arg1,zero);
            return mr::converter<Ret,decltype(tmp)>::eval(tmp);
        }
        else
        {
            auto tmp = pow_helper<TZ,T>::eval(zero,arg1);
            return mr::converter<Ret,decltype(tmp)>::eval(tmp);
        };
    }
};

template<bool zero_on_right>
struct eval_zero_pow<Object,Object,Object,zero_on_right>
{
    using T = Object;
    static T eval(const T& arg1, ti::ti_type<T> ti_z)
    {
        T zero = md::default_value<T>(ti_z);
        return zero_on_right? pow_helper<T,T>::eval(arg1,zero)
                              :pow_helper<T,T>::eval(zero,arg1);
    }
};

struct eval_pow_functor
{
    static const bool ZZ = false;
    static const bool ZN = false;
    static const bool NZ = false;

    template<class ret, class T1, class T2>
    static ret eval(const T1& arg1, const T2& arg2)
    {
        return pow_helper<T1,T2>::eval(arg1,arg2);
    };

    template<class ret, class T1, class TZ, bool zero_on_right>
    static ret eval_zero(ti::ti_type<ret> ti_ret, const T1& arg1, ti::ti_type<TZ> ti_z)
    {
        (void)ti_ret;
        return eval_zero_pow<ret,T1,TZ,zero_on_right>::eval(arg1, ti_z);
    };

    template<bool zero_on_right>
    using is_eval_zero_id   = std::false_type;

    template<class ret, bool is_inv, class TI1, class TI2>
    static ti::ti_type<ret> return_type(TI1 t1, TI2 t2)
    {
        if (md::runtime_value<bool, is_inv>::eval() == false)
            return ti::get_return_ti<ti::ti_type<ret>>(mdyf::pow::eval(),t1,t2);
        else
            return ti::get_return_ti<ti::ti_type<ret>>(mdyf::pow::eval(),t2,t1);
    };

    template<class T, bool is_rev>
    static pow_functor<T,is_rev,false> get_scalar_functor(const T& val)
    {
        return pow_functor<T,is_rev,false>(val);
    };

    static struct_flag op_struct(struct_flag, struct_flag, bool, bool, bool)
    {
        return struct_flag();
    };
};

template<class Ret, class T, class TZ, bool zero_on_right>
struct eval_zero_pow_c
{
    static Ret eval(T arg1, ti::ti_type<TZ> ti_z)
    {
        TZ zero = md::default_value<TZ>(ti_z);

        if (zero_on_right)
        {
            auto tmp = pow_helper<T,TZ>::eval(arg1,zero);
            return mr::converter<Ret,decltype(tmp)>::eval(tmp);
        }
        else
        {
            auto tmp = pow_helper<TZ,T>::eval(zero,arg1);
            return mr::converter<Ret,decltype(tmp)>::eval(tmp);
        };
    }
};

template<bool zero_on_right>
struct eval_zero_pow_c<Object,Object,Object,zero_on_right>
{
    using T = Object;
    static T eval(const T& arg1, ti::ti_type<T> ti_z)
    {
        T zero = md::default_value<T>(ti_z);
        return zero_on_right? pow_c_helper<T,T>::eval(arg1,zero)
                              :pow_c_helper<T,T>::eval(zero,arg1);
    }
};

struct eval_pow_c_functor
{
    static const bool ZZ = false;
    static const bool ZN = false;
    static const bool NZ = false;

    template<class ret, class T1, class T2>
    static ret eval(const T1& arg1, const T2& arg2)
    {
        return pow_c_helper<T1,T2>::eval(arg1,arg2);
    };

    template<class ret, class T1, class TZ, bool zero_on_right>
    static ret eval_zero(ti::ti_type<ret> ti_ret, const T1& arg1, ti::ti_type<TZ> ti_z)
    {
        (void)ti_ret;
        return eval_zero_pow_c<ret,T1,TZ,zero_on_right>::eval(arg1, ti_z);
    };

    template<bool zero_on_right>
    using is_eval_zero_id   = std::false_type;

    template<class ret, bool is_inv, class TI1, class TI2>
    static ti::ti_type<ret> return_type(TI1 t1, TI2 t2)
    {
        if (md::runtime_value<bool, is_inv>::eval() == false)
            return ti::get_return_ti<ti::ti_type<ret>>(mdyf::pow_c::eval(),t1,t2);
        else
            return ti::get_return_ti<ti::ti_type<ret>>(mdyf::pow_c::eval(),t2,t1);
    };

    template<class T, bool is_rev>
    static pow_functor<T,is_rev,true> get_scalar_functor(const T& val)
    {
        return pow_functor<T,is_rev,true>(val);
    };

    static struct_flag op_struct(struct_flag, struct_flag, bool, bool, bool)
    {
        return struct_flag();
    };
};

//-------------------------------------------------------------
//              MAT - MAT
//-------------------------------------------------------------
template<class M1,class M2>
void pow_helper_mat_mat_impl<M1,M2>::eval(matcl::Matrix& ret, const M1& A, const M2& B )
{
    return eval_op_impl<ret_type_pow,M1,M2,eval_pow_functor>::eval(ret,A,B);
};

template<class M1, class M2, bool is_obj>
struct pow_helper_mat_mat_impl_c
{
    static void eval(matcl::Matrix&, const M1&, const M2&) {};
};

template<class M1, class M2>
struct pow_helper_mat_mat_impl_c<M1, M2, true>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B )
    {
        using ret_type_pow = typename pow_helper_mat_mat_impl<M1,M2>::ret_type_pow;
        return eval_op_impl<ret_type_pow,M1,M2,eval_pow_c_functor>::eval(ret,A,B);
    };
};

template<class M1,class M2>
void pow_helper_mat_mat_impl<M1,M2>::eval_c(matcl::Matrix& ret, const M1& A, const M2& B )
{
    using VR                = typename ret_type_pow::value_type;
    static const bool iso   = md::is_object<VR>::value;
    return pow_helper_mat_mat_impl_c<M1, M2, iso>::eval(ret, A, B);
};

//-------------------------------------------------------------
//              SCAL - MAT
//-------------------------------------------------------------
template<class M1,class M2>
void pow_helper_scal_mat_impl<M1,M2>::eval(matcl::Matrix& ret, const M1& A, const M2& B)
{
    return eval_scalar_impl<ret_type_pow,M2,eval_pow_functor,false>::eval(ret,B,true,A);
};

template<class M1, class M2, bool is_obj>
struct pow_helper_scal_mat_impl_c
{
    static void eval(matcl::Matrix&, const M1&, const M2&) {};
};

template<class M1, class M2>
struct pow_helper_scal_mat_impl_c<M1, M2, true>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B )
    {
        using ret_type_pow = typename pow_helper_scal_mat_impl<M1,M2>::ret_type_pow;
        return eval_scalar_impl<ret_type_pow,M2,eval_pow_c_functor,false>::eval(ret,B,true,A);
    };
};

template<class M1,class M2>
void pow_helper_scal_mat_impl<M1,M2>::eval_c(matcl::Matrix& ret, const M1& A, const M2& B )
{
    using VR                = typename ret_type_pow::value_type;
    static const bool iso   = md::is_object<VR>::value;
    return pow_helper_scal_mat_impl_c<M1, M2, iso>::eval(ret, A, B);
};

//-------------------------------------------------------------
//              MAT - SCAL
//-------------------------------------------------------------
template<class M2,class M1>
void pow_helper_mat_scal_impl<M2,M1>::eval(matcl::Matrix& ret, const M2& B, const M1& A)
{
    return eval_scalar_impl<ret_type_pow,M2,eval_pow_functor,false>::eval(ret,B,false,A);
};

template<class M1, class M2, bool is_obj>
struct pow_helper_mat_scal_impl_c
{
    static void eval(matcl::Matrix&, const M1&, const M2&) {};
};

template<class M1, class M2>
struct pow_helper_mat_scal_impl_c<M1, M2, true>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B )
    {
        using ret_type_pow = typename pow_helper_scal_mat_impl<M1,M2>::ret_type_pow;
        return eval_scalar_impl<ret_type_pow,M2,eval_pow_c_functor,false>::eval(ret,B,false,A);
    };
};

template<class M2,class M1>
void pow_helper_mat_scal_impl<M2,M1>::eval_c(matcl::Matrix& ret, const M2& B, const M1& A)
{
    using VR                = typename ret_type_pow::value_type;
    static const bool iso   = md::is_object<VR>::value;
    return pow_helper_mat_scal_impl_c<M1, M2, iso>::eval(ret, A, B);
};

template struct mrd::pow_helper_mat_mat_impl<integer_band, complex_dense>;
template struct mrd::pow_helper_mat_mat_impl<integer_band, complex_sparse>;
template struct mrd::pow_helper_mat_mat_impl<integer_band, complex_band>;
template struct mrd::pow_helper_mat_mat_impl<integer_dense, complex_dense>;
template struct mrd::pow_helper_mat_mat_impl<integer_dense, complex_sparse>;
template struct mrd::pow_helper_mat_mat_impl<integer_dense, complex_band>;
template struct mrd::pow_helper_mat_mat_impl<integer_sparse, complex_dense>;
template struct mrd::pow_helper_mat_mat_impl<integer_sparse, complex_sparse>;
template struct mrd::pow_helper_mat_mat_impl<integer_sparse, complex_band>;

template struct mrd::pow_helper_mat_mat_impl<complex_dense , integer_band     >;
template struct mrd::pow_helper_mat_mat_impl<complex_sparse, integer_band     >;
template struct mrd::pow_helper_mat_mat_impl<complex_band  , integer_band     >;
template struct mrd::pow_helper_mat_mat_impl<complex_dense , integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<complex_sparse, integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<complex_band  , integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<complex_dense , integer_sparse   >;
template struct mrd::pow_helper_mat_mat_impl<complex_sparse, integer_sparse   >;
template struct mrd::pow_helper_mat_mat_impl<complex_band  , integer_sparse   >;

template struct mrd::pow_helper_mat_scal_impl<integer_band  , Complex   >;
template struct mrd::pow_helper_mat_scal_impl<integer_dense , Complex   >;
template struct mrd::pow_helper_mat_scal_impl<integer_sparse, Complex   >;

template struct mrd::pow_helper_scal_mat_impl<Complex, integer_band      >;
template struct mrd::pow_helper_scal_mat_impl<Complex, integer_dense     >;
template struct mrd::pow_helper_scal_mat_impl<Complex, integer_sparse    >;

template struct mrd::pow_helper_mat_mat_impl<real_dense,          integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<float_dense,         integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<complex_dense,       integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<float_complex_dense, integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<object_dense,        integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<real_sparse,         integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<float_sparse,        integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<complex_sparse,      integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<float_complex_sparse,integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<object_sparse,       integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<real_band,           integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<float_band,          integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<complex_band,        integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<float_complex_band,  integer_dense    >;
template struct mrd::pow_helper_mat_mat_impl<object_band,         integer_dense    >;

template struct mrd::pow_helper_mat_mat_impl<real_dense,          integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<float_dense,         integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<complex_dense,       integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<float_complex_dense, integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<object_dense,        integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<real_sparse,         integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<float_sparse,        integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<complex_sparse,      integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<float_complex_sparse,integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<object_sparse,       integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<real_band,           integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<float_band,          integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<complex_band,        integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<float_complex_band,  integer_sparse    >;
template struct mrd::pow_helper_mat_mat_impl<object_band,         integer_sparse    >;

template struct mrd::pow_helper_mat_mat_impl<real_dense,          integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<float_dense,         integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<complex_dense,       integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<float_complex_dense, integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<object_dense,        integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<real_sparse,         integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<float_sparse,        integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<complex_sparse,      integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<float_complex_sparse,integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<object_sparse,       integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<real_band,           integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<float_band,          integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<complex_band,        integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<float_complex_band,  integer_band    >;
template struct mrd::pow_helper_mat_mat_impl<object_band,         integer_band    >;

template struct mrd::pow_helper_mat_scal_impl<real_dense,          Integer    >;
template struct mrd::pow_helper_mat_scal_impl<float_dense,         Integer    >;
template struct mrd::pow_helper_mat_scal_impl<complex_dense,       Integer    >;
template struct mrd::pow_helper_mat_scal_impl<float_complex_dense, Integer    >;
template struct mrd::pow_helper_mat_scal_impl<object_dense,        Integer    >;
template struct mrd::pow_helper_mat_scal_impl<real_sparse,         Integer    >;
template struct mrd::pow_helper_mat_scal_impl<float_sparse,        Integer    >;
template struct mrd::pow_helper_mat_scal_impl<complex_sparse,      Integer    >;
template struct mrd::pow_helper_mat_scal_impl<float_complex_sparse,Integer    >;
template struct mrd::pow_helper_mat_scal_impl<object_sparse,       Integer    >;
template struct mrd::pow_helper_mat_scal_impl<real_band,           Integer    >;
template struct mrd::pow_helper_mat_scal_impl<float_band,          Integer    >;
template struct mrd::pow_helper_mat_scal_impl<complex_band,        Integer    >;
template struct mrd::pow_helper_mat_scal_impl<float_complex_band,  Integer    >;
template struct mrd::pow_helper_mat_scal_impl<object_band,         Integer    >;

template struct mrd::pow_helper_scal_mat_impl<Float, integer_dense >;
template struct mrd::pow_helper_scal_mat_impl<Float, integer_sparse>;
template struct mrd::pow_helper_scal_mat_impl<Float, integer_band  >;
template struct mrd::pow_helper_scal_mat_impl<Real, integer_dense >;
template struct mrd::pow_helper_scal_mat_impl<Real, integer_sparse>;
template struct mrd::pow_helper_scal_mat_impl<Real, integer_band  >;
template struct mrd::pow_helper_scal_mat_impl<Complex, integer_dense >;
template struct mrd::pow_helper_scal_mat_impl<Complex, integer_sparse>;
template struct mrd::pow_helper_scal_mat_impl<Complex, integer_band  >;
template struct mrd::pow_helper_scal_mat_impl<Float_complex, integer_dense >;
template struct mrd::pow_helper_scal_mat_impl<Float_complex, integer_sparse>;
template struct mrd::pow_helper_scal_mat_impl<Float_complex, integer_band  >;
template struct mrd::pow_helper_scal_mat_impl<Object, integer_dense >;
template struct mrd::pow_helper_scal_mat_impl<Object, integer_sparse>;
template struct mrd::pow_helper_scal_mat_impl<Object, integer_band  >;
}}}

MACRO_INSTANTIATE_SG_2(matcl::raw::details::pow_helper_scal_mat_impl)
MACRO_INSTANTIATE_SST_2(matcl::raw::details::pow_helper_scal_mat_impl)

MACRO_INSTANTIATE_GS_2(matcl::raw::details::pow_helper_mat_scal_impl)
MACRO_INSTANTIATE_STS_2(matcl::raw::details::pow_helper_mat_scal_impl)

MACRO_INSTANTIATE_BIN_ALL(matcl::raw::details::pow_helper_mat_mat_impl)

#pragma warning(pop)