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

#include "matcl-matfunc/func/raw/bin/eval_bin_functor.h"
#include "matcl-matrep/func/raw/eval_functor.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-internals/base/utils.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace raw { namespace details
{

template<class ret_type,class in_type,class functor,bool is_scal_str>
struct eval_scalar_impl
{
    template<class T>
    static void eval(matcl::Matrix& ret, const in_type& A,bool is_rev,const T& val)
    {
        using val_ret = typename ret_type::value_type;

        if (is_rev)
        {
            auto ti = functor::template return_type<val_ret,true>(ti::get_ti(A),ti::get_ti(val));
            return eval_functor_2_main<ret_type,in_type>
                    ::eval(ret, ti, A, functor::template get_scalar_functor<T,true>(val));
        }
        else
        {            
            auto ti = functor::template return_type<val_ret,false>(ti::get_ti(A),ti::get_ti(val));
            return eval_functor_2_main<ret_type,in_type>
                    ::eval(ret,ti,A,functor::template get_scalar_functor<T,false>(val));
        };
    };
};

template<class ret_type,class in_type,class functor>
struct eval_scalar_impl<ret_type,in_type,functor, true>
{
    template<class T>
    static void eval(matcl::Matrix& ret, const in_type& A,bool is_rev,const T& val)
    {
        using val_type      = typename ret_type::value_type;
        using FullMatrix    = Matrix<val_type,struct_dense>;

        eval_scalar_impl<FullMatrix,in_type,functor,false>::eval(ret, A,is_rev,val);
    };
};

template<class ret_type,class in_type,class functor,bool is_scal_str>
struct eval_scalar2_impl
{
    template<class T>
    static void eval(matcl::Matrix& ret, const in_type& A,bool is_rev,const T& val, const functor& func)
    {
        using val_ret       = typename ret_type::value_type;

        if (is_rev)
        {
            return eval_functor_2_main<ret_type,in_type>
                            ::eval(ret, functor::template return_type<val_ret,true>(ti::get_ti(A),ti::get_ti(val)),
                                   A,func.template get_scalar_functor<T,true>(val));
        }
        else
        {
            return eval_functor_2_main<ret_type,in_type>
                            ::eval(ret, functor::template return_type<val_ret,true>(ti::get_ti(A),ti::get_ti(val)),
                                   A,func.template get_scalar_functor<T,false>(val));
        };
    };
};

template<class ret_type,class in_type,class functor>
struct eval_scalar2_impl<ret_type,in_type,functor, true>
{
    template<class T>
    static void eval(matcl::Matrix& ret, const in_type& A,bool is_rev,const T& val, const functor& func)
    {
        using val_type      = typename ret_type::value_type;
        using FullMatrix    = Matrix<val_type,struct_dense>;

        eval_scalar2_impl<FullMatrix,in_type,functor,false>::eval(ret, A,is_rev,val,func);
    };
};

template<class ret_type,class M1,class M2,class functor>
struct eval_op_impl
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B)
    {
        using ret_val_type  = typename ret_type::value_type;
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;

        if (A.rows() == 1 && A.cols() == 1)
        {		
            static const bool is_str = !std::is_same<typename M1::struct_type,struct_dense>::value;

            val_type_1 s1(A(1,1));
            
            if (md::is_complex<val_type_1>::value 
                && mrd::is_zero(imag_helper<val_type_1>::eval(s1)))
            {
                using val_type = typename md::real_type<val_type_1>::type;	
                return eval_scalar_impl<ret_type,M2,functor,is_str>
                            ::eval(ret, B,true,real_helper<val_type_1>::eval(s1) );
            }
            else
            {
                return eval_scalar_impl<ret_type,M2,functor,is_str>::eval(ret,B,true,s1);
            };
        };

        if (B.rows() == 1 && B.cols() == 1)
        {
            static const bool is_str = !std::is_same<typename M2::struct_type,struct_dense>::value;
            val_type_2 s2(B(1,1));

            if (md::is_complex<val_type_2>::value 
                && mrd::is_zero(imag_helper<val_type_2>::eval(s2)))
            {
                using val_type = typename md::real_type<val_type_2>::type;
                return eval_scalar_impl<ret_type,M1,functor,is_str>
                            ::eval(ret, A,false,real_helper<val_type_2>::eval(s2));
            }
            else
            {
                return eval_scalar_impl<ret_type,M1,functor,is_str>::eval(ret, A,false,s2);
            };
        };

        static const bool ZZ = functor::ZZ;
        static const bool ZN = functor::ZN;
        static const bool NZ = functor::NZ;

        eval_bin_functor<ret_type,functor,M1,M2,ZZ,ZN,NZ>::eval(ret, A,B);
        
        bool is_square  = ret.is_square();
        struct_flag so  = functor::op_struct(A.get_struct(),B.get_struct(),
                                             is_real_matrix(A), is_real_matrix(B), is_square);
        ret.set_struct(so);

        return;
    };
};

}}}

#pragma warning( pop )