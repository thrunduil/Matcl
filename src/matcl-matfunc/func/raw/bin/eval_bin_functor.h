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

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant
#pragma warning(disable:4100)	// 'fun' : unreferenced formal parameter (strange warning)
#pragma warning(disable:4723)	// potential divide by 0

#include "matcl-internals/base/utils.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl { namespace raw
{

namespace md = matcl::details;
namespace mrd = matcl::raw::details;

template<class func,bool is_inv,class ret_type, class in1, class in2>
struct eval_F_impl
{
    static ret_type eval(const in1& a, const in2& b,const func& fun)
    {
        return fun.template eval<ret_type,in1,in2>(a,b);
    };
};

template<class func,class ret_type, class in1, class in2>
struct eval_F_impl<func,true,ret_type,in1,in2>
{
    static ret_type eval(const in1& a, const in2& b, const func& fun)
    {
        return fun.template eval<ret_type,in2,in1>(b,a);
    };
};

template<bool zero_on_right, class T1, class T2>
struct select_value_type 
{ 
    using type_1    = T1; 
    using type_2    = T2; 

    static ti::ti_type<T1> eval_ti_1(ti::ti_type<T1> t1,ti::ti_type<T2> t2)
    {
        return t1;
    };

    static ti::ti_type<T2> eval_ti_2(ti::ti_type<T1> t1,ti::ti_type<T2> t2)
    {
        return t2;
    };
};

template<class T1,class T2>
struct select_value_type<false,T1,T2> 
{ 
    using type_2    = T1; 
    using type_1    = T2; 

    static ti::ti_type<T2> eval_ti_1(ti::ti_type<T1> t1,ti::ti_type<T2> t2)
    {
        return t2;
    };

    static ti::ti_type<T1> eval_ti_2(ti::ti_type<T1> t1,ti::ti_type<T2> t2)
    {
        return t1;
    };
};

template<class func,bool is_inv,class ret_type, class in1, class in2>
struct eval_F
{
    const func&             f;
    ti::ti_type<ret_type>   ti_ret;
    ti::ti_type<in1>        ti_1;
    ti::ti_type<in2>        ti_2;

    eval_F (const func& f, ti::ti_type<ret_type> tr, ti::ti_type<in1> t1, ti::ti_type<in2> t2) 
        : f(f), ti_ret(tr), ti_1(t1), ti_2(t2)
    {};

    ret_type eval(const in1& a, const in2& b)
    {
        return eval_F_impl<func,is_inv,ret_type,in1,in2>::eval(a,b,f);
    };

    template<bool zero_on_right>
    using is_eval_zero_id   = typename func::template is_eval_zero_id<zero_on_right>;

    template<bool zero_on_right>
    ret_type eval_zero(const typename select_value_type<zero_on_right != is_inv,in1,in2>::type_1& val)
    {
        using T1    = typename select_value_type<zero_on_right != is_inv,in1,in2>::type_1;
        using T2    = typename select_value_type<zero_on_right != is_inv,in1,in2>::type_2;

        ti::ti_type<T2> ti_z = select_value_type<zero_on_right != is_inv,in1,in2>::eval_ti_2(ti_1, ti_2);
        return f.template eval_zero<ret_type,T1,T2,zero_on_right>(ti_ret,val,ti_z);
    };

    private:
        eval_F(const eval_F&) = delete;
        eval_F& operator=(const eval_F&) = delete;
};

namespace details
{

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,
            bool ZN,
            bool NZ,
            bool is_inv,
            class str_ret,
            class str_1,						
            class str_2
        >
struct eval_bin_functor_impl
{};

//--------------------------------------------------------------------
//				GM - GM
//--------------------------------------------------------------------

template<class M1, class M2, class Val_ret, class Val1, bool Is_inv, class Functor>
struct eval_DD_1_inpl
{
    static void eval(matcl::Matrix&, const M1&, const M2&, const Functor&)
    {
        //we should not be here
    };
};

template<class M1, class M2, class Val_ret, class Val2, bool Is_inv, class Functor>
struct eval_DD_2_inpl
{
    static void eval(matcl::Matrix&, const M1&, const M2&, const Functor&)
    {
        //we should not be here
    };
};

template<class M1, class M2, class Val, bool Is_inv, class Functor>
struct eval_DD_1_inpl<M1,M2,Val,Val,Is_inv,Functor>
{
    static void eval(matcl::Matrix& ret, M1& A, const M2& B, const Functor& func)
    {        
        using val_type_1    = Val;
        using val_type_2    = typename M2::value_type;
        using val_type_ret  = Val;

        Integer r = A.rows();
        Integer c = A.cols();

        error::check_eeop(r, c, B.rows(), B.cols()); 

        using eval_func     = eval_F<Functor,Is_inv,val_type_ret,val_type_1,val_type_2>;

        eval_func ef(func,ti::get_ti(A), ti::get_ti(A), ti::get_ti(B));        

        val_type_1* ptr_A       = A.ptr();
        const val_type_2* ptr_B = B.ptr();

        Integer ld_A    = A.ld();
        Integer ld_B    = B.ld();
        
        for (Integer j = 0; j < c; ++j)
        {
            for (Integer i = 0; i < r; ++i)
                ptr_A[i] = ef.eval(ptr_A[i],ptr_B[i]); 

            ptr_A       += ld_A;
            ptr_B       += ld_B;
        };
     
        ret = matcl::Matrix(A,false);
        ret.set_struct(struct_flag());
        return;
    };
};

template<class M1, class M2, class Val, bool Is_inv, class Functor>
struct eval_DD_2_inpl<M1,M2,Val,Val,Is_inv,Functor>
{
    static void eval(matcl::Matrix& ret, const M1& A, M2& B, const Functor& func)
    {        
        using val_type_1    = typename M1::value_type;
        using val_type_2    = Val;
        using val_type_ret  = Val;

        Integer r = A.rows();
        Integer c = A.cols();

        error::check_eeop(r, c, B.rows(), B.cols()); 

        using eval_func     = eval_F<Functor,Is_inv,val_type_ret,val_type_1,val_type_2>;

        eval_func ef(func,ti::get_ti(B), ti::get_ti(A), ti::get_ti(B));        

        const val_type_1* ptr_A = A.ptr();
        val_type_2* ptr_B       = B.ptr();

        Integer ld_A    = A.ld();
        Integer ld_B    = B.ld();
        
        for (Integer j = 0; j < c; ++j)
        {
            for (Integer i = 0; i < r; ++i)
                ptr_B[i] = ef.eval(ptr_A[i],ptr_B[i]); 

            ptr_A       += ld_A;
            ptr_B       += ld_B;
        };
     
        ret = matcl::Matrix(B,false);
        ret.set_struct(struct_flag());
        return;
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,	//0,1
            bool ZN,	//0,1
            bool NZ,	//0,1
            bool is_inv
            //class str_ret = struct_dense,
            //class str_1= struct_dense,
            //class str_2= struct_dense
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,ZZ,ZN,NZ,is_inv,
            struct_dense,struct_dense,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using val_type_ret  = typename ret_type::value_type;
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;

        ti::ti_type<val_type_ret> ret_ti 
                        = func.template return_type<val_type_ret,is_inv>(ti::get_ti(A),ti::get_ti(B));

        if (A.is_unique() && std::is_same<val_type_1, val_type_ret>::value && ret_ti == A.get_type())
        {
            M1 A2 = A.make_unique();
            return eval_DD_1_inpl<M1, M2, val_type_ret, val_type_1, is_inv, functor>
                        ::eval(ret, A2, B, func);
        };

        if (B.is_unique() && std::is_same<val_type_2, val_type_ret>::value && ret_ti == B.get_type())
        {
            M2 B2 = B.make_unique();
            return eval_DD_2_inpl<M1, M2, val_type_ret, val_type_2, is_inv, functor>
                        ::eval(ret, A, B2, func);
        };

        Integer r = A.rows();
        Integer c = A.cols();

        error::check_eeop(r, c, B.rows(), B.cols()); 

        ret_type res        = ret_type(ret_ti, r, c);      
        using eval_func     = eval_F<functor,is_inv,val_type_ret,val_type_1,val_type_2>;

        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        const val_type_1* ptr_A = A.ptr();
        const val_type_2* ptr_B = B.ptr();
        val_type_ret* ptr_res   = res.ptr();

        Integer ld_A    = A.ld();
        Integer ld_B    = B.ld();
        Integer ld_res  = res.ld();
        
        for (Integer j = 0; j < c; ++j)
        {
            for (Integer i = 0; i < r; ++i)
                ptr_res[i] = ef.eval(ptr_A[i],ptr_B[i]); 

            ptr_A       += ld_A;
            ptr_B       += ld_B;
            ptr_res     += ld_res;
        };
     
        ret = matcl::Matrix(res,true);
        return;
    };
};

//--------------------------------------------------------------------
//				GM - SM
//--------------------------------------------------------------------
template<class M1, class M2, class Val_ret, class Val1, bool Is_inv, class Functor>
struct eval_DS_1_inpl
{
    static void eval(matcl::Matrix&, M1&, const M2&, const Functor&)
    {
        //we should not be here
    };
};

template<class M1, class M2, class Val, bool Is_inv, class Functor>
struct eval_DS_1_inpl<M1,M2,Val,Val,Is_inv,Functor>
{
    static void eval(matcl::Matrix& ret, M1& A, const M2& B, const Functor& func)
    {
        using val_type_2    = typename M2::value_type;
        static const bool zero_on_right = !Is_inv;

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols()); 

        Val* ptr_A      = A.ptr();

        Integer r       = A.rows();
        Integer c       = A.cols();
        Integer A_ld    = A.ld();

        using eval_func         = eval_F<Functor,Is_inv,Val,Val,val_type_2>;
        static const bool is_id = eval_func::is_eval_zero_id<zero_on_right>::value;

        eval_func ef(func,ti::get_ti(A),ti::get_ti(A), ti::get_ti(B));        

        if (B.nnz() == 0)
        {
            if (is_id == false)
            {
                for (Integer j = 0; j < c; ++j)
                {
                    for (Integer i = 0; i < r; ++i)
                        ptr_A[i]    = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                    ptr_A   += A_ld;
                };
            };

            ret = matcl::Matrix(A,false);
            ret.set_struct(struct_flag());
            return;
        };
     
        const sparse_ccs<val_type_2>& Bd = B.rep();
        const Integer* Bd_c		    = Bd.ptr_c();
        const Integer* Bd_r		    = Bd.ptr_r();
        const val_type_2* Bd_x	    = Bd.ptr_x();

        if (is_id == false)
        {
            for (Integer j = 0; j < c; ++j)
            {
                Integer pos = 0;
                for (Integer k = Bd_c[j]; k < Bd_c[j + 1] ; ++k)
                {
                    Integer p = Bd_r[k];
                    while (pos < p)
                    {
                        ptr_A[pos]  = ef.template eval_zero<zero_on_right>(ptr_A[pos]);
                        ++pos;
                    }

                    ptr_A[pos] = ef.eval(ptr_A[pos],Bd_x[k]);
                    ++pos;
                };

                while (pos < r)
                {
                    ptr_A[pos] = ef.template eval_zero<zero_on_right>(ptr_A[pos]);
                    ++pos;
                };

                ptr_A   += A_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < c; ++j)
            {
                for (Integer k = Bd_c[j]; k < Bd_c[j + 1] ; ++k)
                {
                    Integer p   = Bd_r[k];
                    ptr_A[p]    = ef.eval(ptr_A[p],Bd_x[k]);
                };

                ptr_A   += A_ld;
            };
        };

        ret = matcl::Matrix(A,false);
        ret.set_struct(struct_flag());
        return;
    }
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,	//0,1
            bool ZN,	//0,1
            //bool NZ,	//0
            bool is_inv
            //class str_ret = struct_dense,
            //class str_1= struct_dense,
            //class str_2= struct_sparse
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,ZZ,ZN,false,is_inv,
                    struct_dense,struct_dense,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using val_type_ret  = typename ret_type::value_type;
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;

        static const bool zero_on_right = !is_inv;

        ti::ti_type<val_type_ret> ret_ti 
                        = func.template return_type<val_type_ret,is_inv>(ti::get_ti(A),ti::get_ti(B));

        if (A.is_unique() && std::is_same<val_type_1, val_type_ret>::value && ret_ti == A.get_type())
        {
            M1 A2 = A.make_unique();
            return eval_DS_1_inpl<M1, M2, val_type_ret, val_type_1, is_inv, functor>
                        ::eval(ret, A2, B, func);
        };

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols()); 

        ret_type res            = ret_type(ret_ti,A.rows(), A.cols()); 

        const val_type_1* ptr_A = A.ptr();
        val_type_ret* ptr_res   = res.ptr();

        Integer r       = A.rows();
        Integer c       = A.cols();
        Integer A_ld    = A.ld();
        Integer res_ld  = res.ld();

        using eval_func = eval_F<functor,is_inv,val_type_ret,val_type_1,val_type_2>;

        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        if (B.nnz() == 0)
        {
            for (Integer j = 0; j < c; ++j)
            {
                for (Integer i = 0; i < r; ++i)
                    ptr_res[i] = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                ptr_res += res_ld;
                ptr_A   += A_ld;
            };

            ret = matcl::Matrix(res,true);
            return;
        };
     
        const sparse_ccs<val_type_2>& Bd = B.rep();
        const Integer* Bd_c		    = Bd.ptr_c();
        const Integer* Bd_r		    = Bd.ptr_r();
        const val_type_2* Bd_x	    = Bd.ptr_x();

        for (Integer j = 0; j < c; ++j)
        {
            Integer pos = 0;
            for (Integer k = Bd_c[j]; k < Bd_c[j + 1] ; ++k)
            {
                Integer p = Bd_r[k];
                while (pos < p)
                {
                    ptr_res[pos] = ef.template eval_zero<zero_on_right>(ptr_A[pos]);
                    ++pos;
                }

                ptr_res[pos] = ef.eval(ptr_A[pos],Bd_x[k]);
                ++pos;
            };

            while (pos < r)
            {
                ptr_res[pos] = ef.template eval_zero<zero_on_right>(ptr_A[pos]);
                ++pos;
            };

            ptr_res += res_ld;
            ptr_A   += A_ld;
        };
 
        ret = matcl::Matrix(res,true);
        return;
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,	//0	
            bool ZN,	//0,1
            //bool NZ,	//1
            bool is_inv
            //class str_ret = struct_sparse,
            //class str_1= struct_dense,
            //class str_2= struct_sparse
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,false,ZN,true,is_inv,
                    struct_sparse,struct_dense,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using val_type_ret  = typename ret_type::value_type;
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;

        static const bool zero_on_right = !is_inv;

        ti::ti_type<val_type_ret> ret_ti 
                        = func.template return_type<val_type_ret,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols()); 

        ret_type res    = ret_type(ret_ti,A.rows(), A.cols(),B.nnz() + A.rows()); 

        if (A.rows() == 0 || A.cols() == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };

        Integer r       = A.rows();
        Integer c       = A.cols();
        Integer A_ld    = A.ld();
     
        const sparse_ccs<val_type_2>& Bd = B.rep();
        const Integer* Bd_c		    = Bd.ptr_c();
        const Integer* Bd_r		    = Bd.ptr_r();
        const val_type_2* Bd_x	    = Bd.ptr_x();

        sparse_ccs<val_type_ret>& d	= res.rep();
        Integer* d_c			= d.ptr_c();
        Integer* d_r			= d.ptr_r();
        val_type_ret* d_x		= d.ptr_x();

        Integer nz				= 0;

        using eval_func = eval_F<functor,is_inv,val_type_ret,val_type_1,val_type_2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        val_type_1 Z1           = matcl::details::default_value<val_type_1>(A.get_type());
        val_type_2 Z2           = matcl::details::default_value<val_type_2>(B.get_type());
        val_type_ret ret_zz		= ef.eval(Z1,Z2);

        const val_type_1* ptr_A = A.ptr();

        if (B.nnz() == 0)
        {
            for (Integer j = 0; j < c; ++j)
            {
                d_c[j]				= nz;

                if (nz + r > d.nzmax()) 
                {
                    d.add_memory( d.nzmax() + r);
                    d_r				= d.ptr_r();
                    d_x				= d.ptr_x();
                };

                for (Integer i = 0; i < r; ++i)
                {
                    if (mrd::is_zero(ptr_A[i]))
                    {
                        d_r[nz]	= i;
                        d_x[nz]	= ret_zz;
                        ++nz;
                    };
                };
                ptr_A += A_ld;
            };

            d_c[c] = nz;
            d.add_memory(-1);
     
            ret = matcl::Matrix(res,true);
            return;
        };

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j]				= nz;
            Integer pos	        = 0;

            if (nz + r > d.nzmax()) 
            {
                d.add_memory( d.nzmax() + r);
                d_r				= d.ptr_r();
                d_x				= d.ptr_x();
            };

            for (Integer k = Bd_c[j]; k < Bd_c[j + 1] ; ++k)
            {
                Integer p		= Bd_r[k];
                for (;pos < p; ++pos)
                {
                    if (mrd::is_zero(ptr_A[pos]))
                    {
                        d_r[nz]	= pos;
                        d_x[nz]	= ret_zz;
                        ++nz;
                    };
                };

                val_type_ret tmp = ef.eval(ptr_A[pos++],Bd_x[k]);
                if (!mrd::is_zero(tmp))
                {
                    d_r[nz]		= p;
                    d_x[nz]		= tmp;
                    ++nz;
                };
            };

            for (;pos < r; ++pos)
            {
                if (mrd::is_zero(ptr_A[pos]))
                {
                    d_r[nz]	= pos;
                    d_x[nz]	= ret_zz;
                    ++nz;
                };
            };

            ptr_A += A_ld;
        };

        d_c[c] = nz;
        d.add_memory(-1);
 
        ret = matcl::Matrix(res,true);
        return;
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,	//1
            bool ZN,	//0,1
            //bool NZ,	//1
            bool is_inv
            //class str_ret = struct_sparse,
            //class str_1= struct_dense,
            //class str_2= struct_sparse
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,true,ZN,true,is_inv,
                    struct_sparse,struct_dense,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using val_type_ret  = typename ret_type::value_type;
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;

        static const bool zero_on_right = !is_inv;

        ti::ti_type<val_type_ret> ret_ti 
                        = func.template return_type<val_type_ret,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols()); 

        ret_type res = ret_type(ret_ti,A.rows(), A.cols(),B.nnz() + A.rows()); 

        if (A.rows() == 0 || A.cols() == 0 || B.nnz() == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };

        Integer c       = A.cols();
        Integer A_ld    = A.ld();
     
        const sparse_ccs<val_type_2>& Bd = B.rep();
        const Integer* Bd_c		    = Bd.ptr_c();
        const Integer* Bd_r		    = Bd.ptr_r();
        const val_type_2* Bd_x	    = Bd.ptr_x();

        sparse_ccs<val_type_ret>& d	= res.rep();
        Integer* d_c			= d.ptr_c();
        Integer* d_r			= d.ptr_r();
        val_type_ret* d_x		= d.ptr_x();

        const val_type_1* ptr_A = A.ptr();

        Integer nz				= 0;

        using eval_func = eval_F<functor,is_inv,val_type_ret,val_type_1,val_type_2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j]				= nz;

            for (Integer k = Bd_c[j]; k < Bd_c[j + 1] ; ++k)
            {
                Integer p		= Bd_r[k];

                val_type_ret tmp = ef.eval(ptr_A[p],Bd_x[k]);
                if (!mrd::is_zero(tmp))
                {
                    d_r[nz]		= p;
                    d_x[nz]		= tmp;
                    ++nz;
                };
            };

            ptr_A += A_ld;
        };

        d_c[c] = nz;
        d.add_memory(-1);
 
        ret = matcl::Matrix(res,true);
        return;
    };
};

//--------------------------------------------------------------------
//				GM - B
//--------------------------------------------------------------------
template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,		//1
            bool ZN,		//0,1
            //bool NZ,		//1
            bool is_inv
            //class str_ret = struct_banded,
            //class str_1= struct_dense,
            //class str_2= struct_banded
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,true,ZN,true,is_inv,
                        struct_banded,struct_dense,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using VTR   = typename ret_type::value_type;
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;
        static const bool zero_on_right = !is_inv;

        ti::ti_type<VTR> ret_ti = func.template return_type<VTR,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols()); 
        
        Integer fdb     = B.first_diag();
        Integer ldb     = B.last_diag();

        ret_type res    = ret_type(ret_ti,B.rows(), B.cols(), fdb, ldb); 

        Integer c       = A.cols();
        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer res_ld  = res.ld();

        using eval_func = eval_F<functor,is_inv,VTR,VT1,VT2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        const VT1* ptr_A = A.ptr();

        if (fdb == ldb)
        {
            const VT2* ptr_B    = B.rep_ptr() + B.first_elem_diag(fdb);
            VTR* ptr_res        = res.rep_ptr() + res.first_elem_diag(fdb);
            Integer rc          = B.diag_length(fdb);

            for (Integer j = 0; j < rc; ++j)
            {
                ptr_res[0]  = ef.eval(ptr_A[0], ptr_B[0]);
                ptr_B       += B_ld;
                ptr_A       += A_ld+1;
                ptr_res     += res_ld;
            };
        }
        else
        {
            const VT2* ptr_B    = B.rep_ptr();
            VTR* ptr_res        = res.rep_ptr();

            for (Integer j = 0; j < c; ++j)
            {
                Integer first_row = B.first_row(j);
                Integer last_row  = B.last_row(j);
                Integer first_elem= B.first_elem_pos(j);

                for (Integer i = first_row, pos_B = first_elem; i <= last_row; ++i, ++pos_B)
                    ptr_res[pos_B] = ef.eval(ptr_A[i],ptr_B[pos_B]);

                ptr_A   += A_ld;
                ptr_B   += B_ld;
                ptr_res += res_ld;				
            };
        }; 

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class M1, class M2, class Val_ret, class Val1, bool Is_inv, class Functor>
struct eval_DB_1_inpl
{
    static void eval(matcl::Matrix&, const M1&, const M2&, const Functor&)
    {
        //we should not be here
    };
};

template<class M1, class M2, class Val, bool Is_inv, class Functor>
struct eval_DB_1_inpl<M1,M2,Val,Val,Is_inv,Functor>
{
    static void eval(matcl::Matrix& ret, M1& A, const M2& B, const Functor& func)
    {
        using VT2    = typename M2::value_type;
        static const bool zero_on_right = !Is_inv;

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols());         

        Integer r       = A.rows();
        Integer c       = A.cols();
        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();

        using eval_func         = eval_F<Functor,Is_inv,Val,Val,VT2>;
        static const bool is_id = eval_func::is_eval_zero_id<zero_on_right>::value;

        eval_func ef(func,ti::get_ti(A),ti::get_ti(A), ti::get_ti(B));

        Val* ptr_A          = A.ptr();                
        Integer fdb         = B.first_diag();
        Integer ldb         = B.last_diag();

        if (fdb == ldb)
        {
            Integer r0      = B.first_row_on_diag(fdb);
            Integer c0      = B.first_col_on_diag(fdb);
            Integer rc      = B.diag_length(fdb);

            const VT2* ptr_B = B.rep_ptr() + B.first_elem_diag(fdb);

            if (is_id == false)
            {
                for (Integer j = 0; j < c0; ++j)
                {
                    for (Integer i = 0; i < j; ++i)
                        ptr_A[i]    = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                    ptr_A   += A_ld;
                };

                for (Integer j = c0; j < c0 + rc; ++j, ++r0)
                {
                    for (Integer i = 0; i < r0; ++i)
                        ptr_A[i]    = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                    ptr_A[j]        = ef.eval(ptr_A[j],ptr_B[0]);

                    for (Integer i = r0 + 1; i < r; ++i)
                        ptr_A[i]    = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                    ptr_A   += A_ld;
                    ptr_B   += B_ld;
                };

                for (Integer j = c0 + rc; j < c; ++j)
                {
                    for (Integer i = 0; i < r; ++i)
                        ptr_A[i]    = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                    ptr_A   += A_ld;
                };
            }
            else
            {
                ptr_A           += c0 * A_ld;

                for (Integer j = 0; j < rc; ++j)
                {
                    ptr_A[0]    = ef.eval(ptr_A[0],ptr_B[0]);
                    ptr_A       += A_ld + 1;
                    ptr_B       += B_ld;
                };
            };
        }
        else
        {
            const VT2* ptr_B = B.rep_ptr();

            if (is_id == false)
            {
                for (Integer j = 0; j < c; ++j)
                {
                    Integer first_row = B.first_row(j);
                    Integer last_row  = B.last_row(j);
                    Integer first_elem= B.first_elem_pos(j);

                    for (Integer i = 0; i < first_row && i < r; ++i)
                        ptr_A[i]    = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                    for (Integer i = first_row, pos_B = first_elem; i <= last_row; ++i, ++pos_B)
                        ptr_A[i]    = ef.eval(ptr_A[i],ptr_B[pos_B]);

                    for (Integer i = last_row+1; i < r; ++i)
                        ptr_A[i]    = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                    ptr_A   += A_ld;
                    ptr_B   += B_ld;
                };
            }
            else
            {
                for (Integer j = 0; j < c; ++j)
                {
                    Integer first_row = B.first_row(j);
                    Integer last_row  = B.last_row(j);
                    Integer first_elem= B.first_elem_pos(j);

                    for (Integer i = first_row, pos_B = first_elem; i <= last_row; ++i, ++pos_B)
                        ptr_A[i]    = ef.eval(ptr_A[i],ptr_B[pos_B]);

                    ptr_A   += A_ld;
                    ptr_B   += B_ld;
                };
            };
        }; 

        ret = matcl::Matrix(A, false);
        ret.set_struct(struct_flag());
        return;
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,		//0,1
            bool ZN,		//0,1
            //bool NZ,		//0
            bool is_inv
            //class str_ret = struct_dense,
            //class str_1= struct_dense,
            //class str_2= struct_banded
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,ZZ,ZN,false,is_inv,
                        struct_dense,struct_dense,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using VTR   = typename ret_type::value_type;
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;
        static const bool zero_on_right = !is_inv;

        ti::ti_type<VTR> ret_ti = func.template return_type<VTR,is_inv>(ti::get_ti(A),ti::get_ti(B));

        if (A.is_unique() && std::is_same<VT1, VTR>::value && ret_ti == A.get_type())
        {
            M1 A2 = A.make_unique();
            return eval_DB_1_inpl<M1, M2, VTR, VT1, is_inv, functor>
                        ::eval(ret, A2, B, func);
        };

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols()); 
        
        ret_type res    = ret_type(ret_ti,A.rows(), A.cols()); 

        Integer r       = A.rows();
        Integer c       = A.cols();
        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer res_ld  = res.ld();

        using eval_func = eval_F<functor,is_inv,VTR,VT1,VT2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        Integer fdb     = B.first_diag();
        Integer ldb     = B.last_diag();

        const VT1* ptr_A    = A.ptr();
        VTR* ptr_res        = res.ptr();

        if (fdb == ldb)
        {
            const VT2* ptr_B    = B.rep_ptr() + B.first_elem_diag(fdb);
            Integer r0          = B.first_row_on_diag(fdb);
            Integer c0          = B.first_col_on_diag(fdb);            
            Integer rc          = B.diag_length(fdb);

            for (Integer j = 0; j < c0; ++j)
            {
                for (Integer i = 0; i < r; ++i)
                    ptr_res[i] = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                ptr_A       += A_ld;
                ptr_res     += res_ld;
            };            

            for (Integer j = c0; j < c0 + rc; ++j, ++r0)
            {
                for (Integer i = 0; i < r0; ++i)
                    ptr_res[i] = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                ptr_res[j] = ef.eval(ptr_A[j],ptr_B[0]);

                for (Integer i = r0+1; i < r; ++i)
                    ptr_res[i] = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                ptr_A   += A_ld;
                ptr_B   += B_ld;
                ptr_res += res_ld;
            };

            for (Integer j = c0 + rc; j < c; ++j)
            {
                for (Integer i = 0; i < r; ++i)
                    ptr_res[i] = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                ptr_A   += A_ld;
                ptr_res += res_ld;
            };
        }
        else
        {
            const VT2* ptr_B    = B.rep_ptr();

            for (Integer j = 0; j < c; ++j)
            {
                Integer first_row = B.first_row(j);
                Integer last_row  = B.last_row(j);
                Integer first_elem= B.first_elem_pos(j);

                for (Integer i = 0; i < first_row && i < r; ++i)
                    ptr_res[i] = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                for (Integer i = first_row, pos_B = first_elem; i <= last_row; ++i, ++pos_B)
                    ptr_res[i] = ef.eval(ptr_A[i],ptr_B[pos_B]);

                for (Integer i = last_row+1; i < r; ++i)
                    ptr_res[i] = ef.template eval_zero<zero_on_right>(ptr_A[i]);

                ptr_A   += A_ld;
                ptr_B   += B_ld;
                ptr_res += res_ld;
            };
        }; 

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,		//0
            bool ZN,		//0,1
            //bool NZ,		//1
            bool is_inv
            //class str_ret = struct_sparse,
            //class str_1= struct_dense,
            //class str_2= struct_banded
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,false,ZN,true,is_inv,
                        struct_sparse,struct_dense,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using VTR   = typename ret_type::value_type;
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;
        static const bool zero_on_right = !is_inv;

        ti::ti_type<VTR> ret_ti = func.template return_type<VTR,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols()); 
        
        Integer r   = A.rows();
        Integer c   = A.cols();

        ret_type res = ret_type(ret_ti,B.rows(), B.cols(), B.nnz() + r);

        if (B.rows() == 0 || B.cols() == 0 || B.nnz() == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };

        sparse_ccs<VTR>& d	= res.rep();

        VT1 zero_1		= md::default_value<VT1>(ti::get_ti(A));
        VT2 zero_2		= md::default_value<VT2>(ti::get_ti(B));

        using eval_func = eval_F<functor,is_inv,VTR,VT1,VT2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        VTR zz			= ef.eval(zero_1,zero_2);

        Integer nz		= 0;

        Integer* d_c	= d.ptr_c();
        Integer* d_r	= d.ptr_r();
        VTR* d_x		= d.ptr_x();

        const VT1* ptr_A = A.ptr();        
        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();

        Integer fdb     = B.first_diag();
        Integer ldb     = B.last_diag();        

        if (fdb == ldb)
        {
            const VT2* ptr_B = B.rep_ptr() + B.first_elem_diag(fdb);
            Integer rc      = B.diag_length(fdb);
            Integer r0      = B.first_row_on_diag(fdb);
            Integer c0      = B.first_col_on_diag(fdb);

            if (!mrd::is_zero(zz))
            {
                for (Integer j = 0; j < c0; ++j)
                {
                    d_c[j]	    = nz;
            
                    if (nz + r > d.nzmax()) 
                    {
                        d.add_memory( d.nzmax() + r);

                        d_r		= d.ptr_r();
                        d_x		= d.ptr_x();
                    };

                    for (Integer i = 0; i < r; ++i)
                    {
                        const VT1& tmp  = ptr_A[i];

                        if (mrd::is_zero(tmp))
                        {
                            d_r[nz]		= i;
                            d_x[nz]		= zz;
                            ++nz;
                        };
                    };

                    ptr_A   += A_ld;
                };                
            }
            else
            {
                for (Integer j = 0; j < c0; ++j)
                    d_c[j]	    = nz;

                ptr_A           += c0 * A_ld;
            };

            for (Integer j = c0; j < rc; ++j, ++r0)
            {
                d_c[j]			= nz;
                if (nz + r > d.nzmax()) 
                {
                    d.add_memory( d.nzmax() + r);

                    d_r			= d.ptr_r();
                    d_x			= d.ptr_x();
                };

                if (!mrd::is_zero(zz))
                {
                    for (Integer i = 0; i < r0; ++i)
                    {
                        const VT1& tmp  = ptr_A[i];
                        if (mrd::is_zero(tmp))
                        {
                            d_r[nz]		= i;
                            d_x[nz]		= zz;
                            ++nz;
                        };
                    };
                };

                VTR val = ef.eval(ptr_A[j],ptr_B[0]);
                if (!mrd::is_zero(val))
                {
                    d_r[nz]		= r0;
                    d_x[nz]		= val;
                    ++nz;
                };

                if (!mrd::is_zero(zz))
                {
                    for (Integer i = r0+1; i < r; ++i)
                    {
                        const VT1& tmp  = ptr_A[i];

                        if (mrd::is_zero(tmp))
                        {
                            d_r[nz]		= i;
                            d_x[nz]		= zz;
                            ++nz;
                        };
                    };
                };

                ptr_A += A_ld;
                ptr_B += B_ld;
            };

            if (!mrd::is_zero(zz))
            {
                for (Integer j = c0 + rc; j < c; ++j)
                {
                    d_c[j]			    = nz;

                    if (nz + r > d.nzmax()) 
                    {
                        d.add_memory( d.nzmax() + r);

                        d_r				= d.ptr_r();
                        d_x				= d.ptr_x();
                    };
                    
                    for (Integer i = 0; i < r; ++i)
                    {
                        const VT1& tmp  = ptr_A[i];
                        if (mrd::is_zero(tmp))
                        {
                            d_r[nz]		= i;
                            d_x[nz]		= zz;
                            ++nz;
                        };
                    };

                    ptr_A += A_ld;
                };
            }
            else
            {
                for (Integer j = rc; j < c; ++j)
                {
                    d_c[j]			    = nz;
                };
            };

            d_c[c] = nz;

            d.add_memory(-1);
            ret = matcl::Matrix(res,true);
            return;
        };

        const VT2* ptr_B = B.rep_ptr();

        for (Integer j = 0; j < c; ++j)
        {
            if (nz + r > d.nzmax()) 
            {
                d.add_memory( d.nzmax() + r);

                d_r				= d.ptr_r();
                d_x				= d.ptr_x();
            };

            d_c[j]			    = nz;

            Integer first_row = B.first_row(j);
            Integer last_row  = B.last_row(j);
            Integer first_elem= B.first_elem_pos(j);
            
            if (!mrd::is_zero(zz))
            {
                for (Integer i = 0; i < first_row && i < r; ++i)
                {
                    const VT1& tmp  = ptr_A[i];
                    if (mrd::is_zero(tmp))
                    {
                        d_r[nz]		= i;
                        d_x[nz]		= zz;
                        ++nz;
                    };
                };
            };

            for (Integer i = first_row, pos_B = first_elem; i <= last_row; ++i, ++pos_B)
            {
                const VT2& tmp   = ptr_B[pos_B];
                VTR val        = ef.eval(ptr_A[i],tmp);

                if (!mrd::is_zero(val))
                {
                    d_r[nz]		= i;
                    d_x[nz]		= val;
                    ++nz;
                };
            };

            if (!mrd::is_zero(zz))
            {
                for (Integer i = last_row+1; i < r; ++i)
                {
                    const VT1& tmp  = ptr_A[i];

                    if (mrd::is_zero(tmp))
                    {
                        d_r[nz]		= i;
                        d_x[nz]		= zz;
                        ++nz;
                    };
                };
            };

            ptr_A += A_ld;
            ptr_B += B_ld;
        };

        d_c[c] = nz;
        d.add_memory(-1);
 
        ret = matcl::Matrix(res,true);
        return;
    };
};

//--------------------------------------------------------------------
//				SM - GM
//--------------------------------------------------------------------
template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,			
            bool ZN,			
            bool NZ,			
            bool is_inv
            //class str_ret = struct_sparse,
            //class str_1= struct_sparse,
            //class str_2= struct_dense
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,ZZ,ZN,NZ,is_inv,
                        struct_sparse,struct_sparse,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        return eval_bin_functor_impl<ret_type,functor,M2,M1,ZZ,NZ,ZN,!is_inv,
                                        struct_sparse,struct_dense,struct_sparse>
                                        ::eval(ret,B,A,func);
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,
            bool ZN,
            bool NZ,
            bool is_inv
            //class str_ret = struct_dense,
            //class str_1= struct_sparse,
            //class str_2= struct_dense
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,ZZ,ZN,NZ,is_inv,
                        struct_dense,struct_sparse,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        return eval_bin_functor_impl<ret_type,functor,M2,M1,ZZ,NZ,ZN,!is_inv,
                                        struct_dense,struct_dense,struct_sparse>
                                        ::eval(ret,B,A,func);
    };
};

//--------------------------------------------------------------------
//				SM - SM
//--------------------------------------------------------------------
template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,		//1
            //bool ZN,		//0
            //bool NZ,		//0
            bool is_inv
            //class str_ret = struct_sparse,
            //class str_1= struct_sparse,
            //class str_2= struct_sparse
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,true,false,false,is_inv,
                        struct_sparse,struct_sparse,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using val_type_ret  = typename ret_type::value_type;
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;
        static const bool zero_on_right = !is_inv;

        ti::ti_type<val_type_ret> ret_ti 
                        = func.template return_type<val_type_ret,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols());

        Integer r = A.rows();
        Integer c = A.cols();

        if (r == 0 || c == 0)
        {
            ret = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };
        
        Integer nzret = B.nnz() + A.nnz();

        ret_type res(ret_ti, r, c, nzret);

        if (nzret == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };

        const sparse_ccs<val_type_1>&	Ad = A.rep();
        const sparse_ccs<val_type_2>&	Bd = B.rep();
        sparse_ccs<val_type_ret>&		d  = res.rep();

        Integer* d_c				= d.ptr_c();
        Integer* d_r				= d.ptr_r();
        val_type_ret* d_x			= d.ptr_x();

        Integer nz					= 0;

        using eval_func = eval_F<functor,is_inv,val_type_ret,val_type_1,val_type_2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        if (A.nnz() == 0)
        {
            const Integer* Bd_c			= Bd.ptr_c();
            const Integer* Bd_r			= Bd.ptr_r();
            const val_type_2* Bd_x		= Bd.ptr_x();

            for (Integer j = 0; j < c; ++j)
            {
                Integer kb				= Bd_c[j];
                d_c[j]					= nz;

                while (kb < Bd_c[j+1])
                {				
                    val_type_ret tmp = ef.template eval_zero<!zero_on_right>(Bd_x[kb]);
                    d_r[nz]			= Bd_r[kb]; 
                    d_x[nz]			= tmp;				
                    ++nz; 
                    ++kb; 
                };
            };

            d_c[c] = nz;
            d.add_memory(-1);
            ret = matcl::Matrix(res,true);
            return;
        };

        if (B.nnz() == 0)
        {
            const Integer* Ad_c			= Ad.ptr_c();
            const Integer* Ad_r			= Ad.ptr_r();
            const val_type_1* Ad_x		= Ad.ptr_x();

            for (Integer j = 0; j < c; ++j)
            {
                Integer ka				= Ad_c[j];
                d_c[j]					= nz;

                while (ka < Ad_c[j+1])
                {				
                    val_type_ret tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    d_r[nz]		= Ad_r[ka]; 
                    d_x[nz]		= tmp;
                    ++nz;
                    ++ka;
                };
            };
            d_c[c] = nz;
            d.add_memory(-1);
            ret = matcl::Matrix(res,true);
            return;
        };

        const Integer* Ad_c			= Ad.ptr_c();
        const Integer* Ad_r			= Ad.ptr_r();
        const val_type_1* Ad_x		= Ad.ptr_x();

        const Integer* Bd_c			= Bd.ptr_c();
        const Integer* Bd_r			= Bd.ptr_r();
        const val_type_2* Bd_x		= Bd.ptr_x();
         
        for (Integer j = 0; j < c; ++j)
        {
            Integer ka				= Ad_c[j]; 
            Integer kb				= Bd_c[j];

            d_c[j]					= nz;

            while ((ka < Ad_c[j+1]) && (kb < Bd_c[j+1]))
            {
                if (Ad_r[ka] < Bd_r[kb])
                {
                    val_type_ret tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    d_r[nz]		= Ad_r[ka]; 
                    d_x[nz]		= tmp;
                    ++nz; 
                    ++ka;
                }
                else if (Ad_r[ka] > Bd_r[kb])
                {
                    val_type_ret tmp = ef.template eval_zero<!zero_on_right>(Bd_x[kb]);
                    d_r[nz]		= Bd_r[kb]; 
                    d_x[nz]		= tmp;
                    ++nz; 
                    ++kb; 
                }
                else
                { 
                    val_type_ret tmp = ef.eval(Ad_x[ka],Bd_x[kb]);
                    if (!mrd::is_zero(tmp))
                    {
                        d_r[nz]		= Ad_r[ka]; 
                        d_x[nz]		= tmp;
                        ++nz; 
                    };
                    ++ka; 
                    ++kb;
                };
            };

            while (ka < Ad_c[j+1])
            {
                val_type_ret tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                d_r[nz]		= Ad_r[ka]; 
                d_x[nz]		= tmp;
                ++nz;
                ++ka;
            };

            while (kb < Bd_c[j+1])
            {				
                val_type_ret tmp = ef.template eval_zero<!zero_on_right>(Bd_x[kb]);
                d_r[nz]			= Bd_r[kb]; 
                d_x[nz]			= tmp;				
                ++nz; 
                ++kb; 
            };
        };

        d_c[c] = nz;
        d.add_memory(-1);
        ret = matcl::Matrix(res,true);
        return;
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,		//1
            //bool ZN,		//1
            //bool NZ,		//0
            bool is_inv
            //class str_ret = struct_sparse,
            //class str_1= struct_sparse,
            //class str_2= struct_sparse
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,true,true,false,is_inv,
                        struct_sparse,struct_sparse,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using val_type_ret  = typename ret_type::value_type;
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;
        static const bool zero_on_right = !is_inv;

        ti::ti_type<val_type_ret> ret_ti 
                        = func.template return_type<val_type_ret,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols());

        Integer r = A.rows();
        Integer c = A.cols();

        if (r == 0 || c == 0)
        {
            ret = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };
        
        Integer nzret = A.nnz();

        ret_type res(ret_ti, r, c, nzret);

        if (nzret == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };

        const sparse_ccs<val_type_1>&	Ad = A.rep();
        const sparse_ccs<val_type_2>&	Bd = B.rep();
        sparse_ccs<val_type_ret>&		d  = res.rep();

        Integer* d_c				= d.ptr_c();
        Integer* d_r				= d.ptr_r();
        val_type_ret* d_x			= d.ptr_x();

        Integer nz					= 0;

        using eval_func = eval_F<functor,is_inv,val_type_ret,val_type_1,val_type_2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        if (B.nnz() == 0)
        {
            const Integer* Ad_c			= Ad.ptr_c();
            const Integer* Ad_r			= Ad.ptr_r();
            const val_type_1* Ad_x		= Ad.ptr_x();

            for (Integer j = 0; j < c; ++j)
            {
                Integer ka				= Ad_c[j];
                d_c[j]					= nz;

                while (ka < Ad_c[j+1])
                {				
                    val_type_ret tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    d_r[nz]		= Ad_r[ka]; 
                    d_x[nz]		= tmp;
                    ++nz;
                    ++ka;
                };
            };

            d_c[c] = nz;
            d.add_memory(-1);
            ret = matcl::Matrix(res,true);
            return;
        };

        const Integer* Ad_c			= Ad.ptr_c();
        const Integer* Ad_r			= Ad.ptr_r();
        const val_type_1* Ad_x		= Ad.ptr_x();

        const Integer* Bd_c			= Bd.ptr_c();
        const Integer* Bd_r			= Bd.ptr_r();
        const val_type_2* Bd_x		= Bd.ptr_x();
         
        for (Integer j = 0; j < c; ++j)
        {
            Integer ka				= Ad_c[j]; 
            Integer kb				= Bd_c[j];

            d_c[j]					= nz;

            while ((ka < Ad_c[j+1]) && (kb < Bd_c[j+1]))
            {
                if (Ad_r[ka] < Bd_r[kb])
                {
                    val_type_ret tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    d_r[nz]		= Ad_r[ka]; 
                    d_x[nz]		= tmp;
                    ++nz; 
                    ++ka;
                }
                else if (Ad_r[ka] > Bd_r[kb])
                {
                    ++kb; 
                }
                else
                { 
                    val_type_ret tmp = ef.eval(Ad_x[ka],Bd_x[kb]);
                    if (!mrd::is_zero(tmp))
                    {
                        d_r[nz]		= Ad_r[ka]; 
                        d_x[nz]		= tmp;
                        ++nz; 
                    };
                    ++ka; 
                    ++kb;
                };
            };
            while (ka < Ad_c[j+1])
            {
                val_type_ret tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                d_r[nz]		= Ad_r[ka]; 
                d_x[nz]		= tmp;
                ++nz;
                ++ka;
            };
        };

        d_c[c] = nz;
        d.add_memory(-1);
        ret = matcl::Matrix(res,true);
        return;
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,		//1
            //bool ZN,		//1
            //bool NZ,		//1
            bool is_inv
            //class str_ret = struct_sparse,
            //class str_1= struct_sparse,
            //class str_2= struct_sparse
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,true,true,true,is_inv,
                        struct_sparse,struct_sparse,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using val_type_ret  = typename ret_type::value_type;
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;
        static const bool zero_on_right = !is_inv;

        ti::ti_type<val_type_ret> ret_ti 
                        = func.template return_type<val_type_ret,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols());

        Integer r = A.rows();
        Integer c = A.cols();

        if (r == 0 || c == 0)
        {
            ret = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };
        
        Integer nzret = std::min(A.nnz(),B.nnz());

        ret_type res(ret_ti, r, c, nzret);
        if (r == 0 || c == 0 || nzret == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };

        const sparse_ccs<val_type_1>&	Ad = A.rep();
        const sparse_ccs<val_type_2>&	Bd = B.rep();
        sparse_ccs<val_type_ret>&		d  = res.rep();

        const Integer* Ad_c			= Ad.ptr_c();
        const Integer* Ad_r			= Ad.ptr_r();
        const val_type_1* Ad_x		= Ad.ptr_x();

        const Integer* Bd_c			= Bd.ptr_c();
        const Integer* Bd_r			= Bd.ptr_r();
        const val_type_2* Bd_x		= Bd.ptr_x();

        Integer* d_c				= d.ptr_c();
        Integer* d_r				= d.ptr_r();
        val_type_ret* d_x			= d.ptr_x();

        Integer nz					= 0;

        using eval_func = eval_F<functor,is_inv,val_type_ret,val_type_1,val_type_2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));
         
        for (Integer j = 0; j < c; ++j)
        {
            Integer ka				= Ad_c[j]; 
            Integer kb				= Bd_c[j];

            d_c[j]					= nz;

            while ((ka < Ad_c[j+1]) && (kb < Bd_c[j+1]))
            {
                if (Ad_r[ka] < Bd_r[kb])
                {
                    ++ka;
                }
                else if (Ad_r[ka] > Bd_r[kb])
                {
                    ++kb; 
                }
                else
                { 
                    val_type_ret tmp = ef.eval(Ad_x[ka],Bd_x[kb]);
                    if (!mrd::is_zero(tmp))
                    {
                        d_r[nz]		= Ad_r[ka]; 
                        d_x[nz]		= tmp;
                        ++nz; 
                    };
                    ++ka; 
                    ++kb;
                };
            };
        };

        d_c[c] = nz;
        d.add_memory(-1);
        ret = matcl::Matrix(res,true);
        return;
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,		//0
            bool ZN,		//0,1
            bool NZ,		//0,1
            bool is_inv
            //class str_ret = struct_dense,
            //class str_1= struct_sparse,
            //class str_2= struct_sparse
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,false,ZN,NZ,is_inv,
                            struct_dense,struct_sparse,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using val_type_ret  = typename ret_type::value_type;
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;
        static const bool zero_on_right = !is_inv;

        ti::ti_type<val_type_ret> ret_ti 
                        = func.template return_type<val_type_ret,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols());

        Integer r = A.rows();
        Integer c = A.cols();

        if (r == 0 || c == 0)
        {
            ret = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };
        
        val_type_1 zero1    		= md::default_value<val_type_1>(A.get_type());
        val_type_2 zero2    		= md::default_value<val_type_2>(B.get_type());

        using eval_func = eval_F<functor,is_inv,val_type_ret,val_type_1,val_type_2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        val_type_ret zz				= ef.eval(zero1,zero2);        

        if (A.nnz() == 0 && B.nnz() == 0)
        {
            if (mrd::is_zero(zz))
            {
                using sparse_mat = Matrix<val_type_ret,struct_sparse>;

                sparse_mat res(ret_ti,r,c);
                ret = matcl::Matrix(res,false);
                return;
            }
            else
            {
                ret_type res(ret_ti,zz,r,c);
                ret = matcl::Matrix(res,false);
                return;
            };            
        };

        ret_type res(ret_ti,zz,r,c);

        const sparse_ccs<val_type_1>& Ad = A.rep();
        const sparse_ccs<val_type_2>& Bd = B.rep();

        val_type_ret* ptr_res   = res.ptr();
        Integer res_ld          = res.ld();

        if (A.nnz() == 0)
        {
            const Integer* Bd_c			= Bd.ptr_c();
            const Integer* Bd_r			= Bd.ptr_r();
            const val_type_2* Bd_x		= Bd.ptr_x();

            for (Integer j = 0; j < c; ++j)
            { 
                Integer kb				= Bd_c[j];

                while (kb < Bd_c[j+1])
                {
                    val_type_ret tmp    = ef.template eval_zero<!zero_on_right>(Bd_x[kb]);
                    ptr_res[Bd_r[kb]]   = tmp;
                    ++kb; 
                };

                ptr_res += res_ld;
            };

            ret = matcl::Matrix(res,true);
            return;
        };

        const Integer* Ad_c			= Ad.ptr_c();
        const Integer* Ad_r			= Ad.ptr_r();
        const val_type_1* Ad_x		= Ad.ptr_x();

        if (B.nnz() == 0)
        {
            for (Integer j = 0; j < c; ++j)
            {
                Integer ka			= Ad_c[j]; 

                while (ka < Ad_c[j+1])
                {
                    val_type_ret tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    ptr_res[Ad_r[ka]] = tmp;
                    ++ka; 
                };

                ptr_res += res_ld;
            };

            ret = matcl::Matrix(res,true);
            return;
        };

        const Integer* Bd_c			= Bd.ptr_c();
        const Integer* Bd_r			= Bd.ptr_r();
        const val_type_2* Bd_x		= Bd.ptr_x();
         
        for (Integer j = 0; j < c; ++j)
        {
            Integer ka				= Ad_c[j]; 
            Integer kb				= Bd_c[j];

            while (ka < Ad_c[j+1] && kb < Bd_c[j+1])
            {
                if (Ad_r[ka] < Bd_r[kb])
                {
                    val_type_ret tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    ptr_res[Ad_r[ka]] = tmp;
                    ++ka; 
                }
                else if (Ad_r[ka] > Bd_r[kb])
                {
                    val_type_ret tmp = ef.template eval_zero<!zero_on_right>(Bd_x[kb]);
                    ptr_res[Bd_r[kb]] = tmp;
                    ++kb; 
                }
                else
                { 
                    val_type_ret tmp = ef.eval(Ad_x[ka],Bd_x[kb]);
                    ptr_res[Bd_r[kb]] = tmp;
                    ++ka; 
                    ++kb;
                };
            };
            while (ka < Ad_c[j+1])
            {
                val_type_ret tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                ptr_res[Ad_r[ka]] = tmp;
                ++ka; 
            };
            while (kb < Bd_c[j+1])
            {
                val_type_ret tmp = ef.template eval_zero<!zero_on_right>(Bd_x[kb]);
                ptr_res[Bd_r[kb]] = tmp;
                ++kb; 
            };

            ptr_res += res_ld;
        };

        ret = matcl::Matrix(res,true);
        return;
    };
};

//--------------------------------------------------------------------
//				SM - B
//--------------------------------------------------------------------
template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,		//1
            //bool ZN,		//1
            //bool NZ,		//1
            bool is_inv
            //class str_ret = struct_banded,
            //class str_1= struct_sparse,
            //class str_2= struct_banded
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,true,true,true,is_inv,
                        struct_banded,struct_sparse,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using VTR   = typename ret_type::value_type;
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;
        static const bool zero_on_right = !is_inv;

        ti::ti_type<VTR> ret_ti = func.template return_type<VTR,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols());

        Integer r = A.rows();
        Integer c = A.cols();

        VTR zero_r = md::default_value<VTR>(ret_ti);

        if (r == 0 || c == 0)
        {
            ret_type out(ret_ti,r, c, 0, 0);
            ret = matcl::Matrix(out,false);
            return;
        };

        if (A.nnz() == 0)
        {
            ret_type res(ret_ti, zero_r, r,c, 0, 0);
            ret = matcl::Matrix(res,false);
            return;
        };
        
        Integer fdb = B.first_diag();
        Integer ldb = B.last_diag();

        ret_type res(ret_ti,r,c, fdb, ldb);

        const sparse_ccs<VT1>& Ad   = A.rep();
        const Integer* Ad_c			= Ad.ptr_c();
        const Integer* Ad_r			= Ad.ptr_r();
        const VT1* Ad_x             = Ad.ptr_x();

        using eval_func = eval_F<functor,is_inv,VTR,VT1,VT2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        Integer B_ld        = B.ld();
        Integer res_ld      = res.ld();

        if (fdb == ldb)
        {
            VTR* ptr_res        = res.rep_ptr() + res.first_elem_diag(fdb);
            const VT2* ptr_B    = B.rep_ptr() + B.first_elem_diag(fdb);

            Integer rc          = B.diag_length(fdb);

            for (Integer j = 0; j < rc; ++j)
            {
                Integer ka;
                bool exist = Ad.has_element(j,j,ka);

                if (!exist)
                    ptr_res[0]  = zero_r;
                else
                    ptr_res[0] = ef.eval(Ad_x[ka], ptr_B[0]);

                ptr_B   += B_ld;
                ptr_res += res_ld;
            };
        }
        else
        {
            VTR* ptr_res        = res.rep_ptr();
            const VT2* ptr_B    = B.rep_ptr();

            for (Integer j = 0; j < c; ++j, ptr_B += B_ld, ptr_res += res_ld)
            {
                Integer first_row   = B.first_row(j);
                Integer last_row    = B.last_row(j);
                Integer pos_B       = B.first_elem_pos(j);
                Integer kb          = first_row;

                if (first_row >= r)
                    continue;

                Integer ka;
                Ad.has_element(kb,j,ka);

                if (ka == Ad_c[j+1] || Ad_r[ka] > last_row)
                {
                    for (Integer i = first_row; i <= last_row; ++i, ++pos_B)
                        ptr_res[pos_B] = zero_r;

                    continue;
                };
                while ((ka < Ad_c[j+1]) && (kb <= last_row))
                {
                    if (Ad_r[ka] > kb)
                    {
                        ptr_res[pos_B] = zero_r;
                        ++kb; 
                        ++pos_B;
                    }
                    else
                    { 
                        VTR tmp = ef.eval(Ad_x[ka],ptr_B[pos_B]);
                        ptr_res[pos_B] = tmp;
                        ++ka; 
                        ++kb;
                        ++pos_B;
                    };
                };
                while (kb <= last_row)
                {
                    ptr_res[pos_B] = zero_r;
                    ++kb; 
                    ++pos_B;
                };
            };
        };

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,		//1
            //bool ZN,		//1
            //bool NZ,		//0
            bool is_inv
            //class str_ret
            //class str_1= struct_sparse,
            //class str_2= struct_banded
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,true,true,false,is_inv,
                        struct_sparse,struct_sparse,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using val_type_ret  = typename ret_type::value_type;
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;
        static const bool zero_on_right = !is_inv;

        ti::ti_type<val_type_ret> ret_ti 
                        = func.template return_type<val_type_ret,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols());

        Integer r = A.rows();
        Integer c = A.cols();

        if (r == 0 || c == 0)
        {
            ret = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };
        
        Integer nnze = A.nnz();

        ret_type res(ret_ti,r,c, nnze);
        if (nnze == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };

        const sparse_ccs<val_type_1>& Ad = A.rep();

        sparse_ccs<val_type_ret>& d		= res.rep();
        Integer* d_c				= d.ptr_c();
        Integer* d_r				= d.ptr_r();
        val_type_ret* d_x			= d.ptr_x();

        Integer nz = 0;

        const val_type_2* ptr_B     = B.rep_ptr();
        Integer B_ld                = B.ld();

        using eval_func = eval_F<functor,is_inv,val_type_ret,val_type_1,val_type_2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        const Integer* Ad_c			= Ad.ptr_c();
        const Integer* Ad_r			= Ad.ptr_r();
        const val_type_1* Ad_x		= Ad.ptr_x();

        for (Integer j = 0; j < c; ++j)
        {
            Integer ka				= Ad_c[j]; 
            Integer first_row		= B.first_row(j);
            Integer last_row		= std::min(B.last_row(j),r-1);
            Integer kb				= first_row;
            Integer pos_B			= B.first_elem_pos(j);

            d_c[j]					= nz;

            while ((ka < Ad_c[j+1]) && (kb <= last_row))
            {
                if (Ad_r[ka] < kb)
                {
                    val_type_ret tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    d_r[nz]			= Ad_r[ka];
                    d_x[nz]			= tmp;
                    ++nz;
                    ++ka; 
                }
                else if (Ad_r[ka] > kb)
                {
                    ++kb; 
                    ++pos_B;
                }
                else
                { 
                    val_type_ret tmp = ef.eval(Ad_x[ka],ptr_B[pos_B]);
                    if (!mrd::is_zero(tmp))
                    {
                        d_r[nz]		= kb;
                        d_x[nz]		= tmp;
                        ++nz;
                    };
                    ++ka; 
                    ++kb;
                    ++pos_B;
                };
            };
            while (ka < Ad_c[j+1])
            {
                val_type_ret tmp    = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                d_r[nz]				= Ad_r[ka];
                d_x[nz]				= tmp;
                ++nz;
                ++ka; 
            };

            ptr_B += B_ld;
        };

        d_c[c] = nz;

        d.add_memory(-1);

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,		//1
            //bool ZN,		//0
            //bool NZ,		//1
            bool is_inv
            //class str_ret
            //class str_1= struct_sparse,
            //class str_2= struct_banded
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,true,false,true,is_inv,
                        struct_banded,struct_sparse,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using VTR   = typename ret_type::value_type;
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;
        static const bool zero_on_right = !is_inv;

        ti::ti_type<VTR> ret_ti = func.template return_type<VTR,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols());

        Integer r = A.rows();
        Integer c = A.cols();

        if (r == 0 || c == 0)
        {
            ret_type out(ret_ti,r, c, 0, 0);
            ret = matcl::Matrix(out,false);
            return;
        };
        
        Integer fdb     = B.first_diag();
        Integer ldb     = B.last_diag();

        ret_type res(ret_ti,r,c, fdb, ldb);

        const sparse_ccs<VT1>& Ad = A.rep();

        const Integer* Ad_c	= Ad.ptr_c();
        const Integer* Ad_r	= Ad.ptr_r();
        const VT1* Ad_x		= Ad.ptr_x();

        using eval_func = eval_F<functor,is_inv,VTR,VT1,VT2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        Integer B_ld        = B.ld();
        Integer res_ld      = res.ld();

        if (A.nnz() == 0)
        {
            VTR* ptr_res        = res.rep_ptr();
            const VT2* ptr_B    = B.rep_ptr();

            for (Integer j = 0; j < c; ++j, ptr_B += B_ld, ptr_res += res_ld)
            {
                Integer first_row		= B.first_row(j);
                Integer last_row		= B.last_row(j);
                Integer pos_B			= B.first_elem_pos(j);

                for (Integer i = first_row; i <= last_row; ++i, ++pos_B)
                {
                    //op(0,x) != 0
                    VTR tmp = ef.template eval_zero<!zero_on_right>(ptr_B[pos_B]);
                    ptr_res[pos_B] = tmp;
                };
            };

            ret = matcl::Matrix(res,false);
            return;
        };

        if (fdb == ldb)
        {
            VTR* ptr_res        = res.rep_ptr() + res.first_elem_diag(fdb);
            const VT2* ptr_B    = B.rep_ptr() + B.first_elem_diag(fdb);
            Integer rc          = B.diag_length(fdb);

            for (Integer j = 0; j < rc; ++j)
            {
                Integer ka;
                bool exist = Ad.has_element(j,j,ka);

                //op(0,x) != 0

                if (!exist)
                {
                    VTR tmp     = ef.template eval_zero<!zero_on_right>(ptr_B[0]);
                    ptr_res[0]  = tmp;
                }
                else
                {
                    ptr_res[0]  = ef.eval(Ad_x[ka],ptr_B[0]);
                };

                ptr_B           += B_ld;
                ptr_res         += res_ld;
            };
        }
        else
        {
            VTR* ptr_res        = res.rep_ptr();
            const VT2* ptr_B    = B.rep_ptr();

            for (Integer j = 0; j < c; ++j, ptr_B += B_ld, ptr_res += res_ld)
            {
                Integer first_row		= B.first_row(j);
                Integer last_row		= B.last_row(j);
                Integer pos_B			= B.first_elem_pos(j);
                Integer kb				= first_row;

                if (first_row >= r)
                    continue;

                Integer ka;
                Ad.has_element(kb,j,ka);

                if (ka == Ad_c[j+1] || Ad_r[ka] > last_row)
                {
                    for (Integer i = first_row; i <= last_row; ++i, ++pos_B)
                    {
                        //op(0,x) != 0
                        VTR tmp = ef.template eval_zero<!zero_on_right>(ptr_B[pos_B]);
                        ptr_res[pos_B] = tmp;
                    };
                    continue;
                };

                while ((ka < Ad_c[j+1]) && (kb <= last_row))
                {
                    if (Ad_r[ka] > kb)
                    {
                        VTR tmp = ef.template eval_zero<!zero_on_right>(ptr_B[pos_B]);
                        ptr_res[pos_B] = tmp;
                        ++kb; 
                        ++pos_B;
                    }
                    else
                    { 
                        VTR tmp = ef.eval(Ad_x[ka],ptr_B[pos_B]);
                        ptr_res[pos_B] = tmp;
                        ++ka; 
                        ++kb;
                        ++pos_B;
                    };
                };
                while (kb <= last_row)
                {
                    VTR tmp = ef.template eval_zero<!zero_on_right>(ptr_B[pos_B]);
                    ptr_res[pos_B] = tmp;
                    ++kb; 
                    ++pos_B;
                };
            };
        };

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,		//1
            //bool ZN,		//0
            //bool NZ,		//0
            bool is_inv
            //class str_ret = struct_sparse,
            //class str_1= struct_sparse,
            //class str_2= struct_banded
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,true,false,false,is_inv,
                        struct_sparse,struct_sparse,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using VTR   = typename ret_type::value_type;
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;
        static const bool zero_on_right = !is_inv;

        ti::ti_type<VTR> ret_ti = func.template return_type<VTR,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols());

        Integer r   = A.rows();
        Integer c   = A.cols();

        if (r == 0 || c == 0)
        {
            ret = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };
        
        Integer nnze = A.nnz() + B.nnz();

        ret_type res(ret_ti,r,c, nnze);

        if (nnze == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };

        const sparse_ccs<VT1>& Ad = A.rep();

        sparse_ccs<VTR>& d	    = res.rep();
        Integer* d_c			= d.ptr_c();
        Integer* d_r			= d.ptr_r();
        VTR* d_x			    = d.ptr_x();

        Integer nz              = 0;
        Integer B_ld            = B.ld();        
        Integer fdb             = B.first_diag();
        Integer ldb             = B.last_diag();

        using eval_func = eval_F<functor,is_inv,VTR,VT1,VT2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        if (A.nnz() == 0)
        {            
            if (fdb == ldb)
            {
                Integer c0          = B.first_col_on_diag(fdb);
                const VT2* ptr_B    = B.rep_ptr() + B.first_elem_diag(fdb);
                Integer rc          = B.diag_length(fdb);

                for (Integer j = 0; j < c0; ++j)
                    d_c[j]			= nz;

                for (Integer j = c0; j < rc; ++j)
                {
                    d_c[j]			= nz;

                    const VT2& val2 = ptr_B[0];

                    if (!mrd::is_zero(val2))
                    {
                        VTR tmp     = ef.template eval_zero<!zero_on_right>(val2);
                        d_r[nz]	    = j;
                        d_x[nz]		= tmp;
                        ++nz;
                    };

                    ptr_B           += B_ld;
                };

                for (Integer j = c0 + rc; j < c; ++j)
                    d_c[j]	        = nz;

                d_c[c] = nz;
                d.add_memory(-1);
            }
            else
            {
                const VT2* ptr_B        = B.rep_ptr();

                for (Integer j = 0; j < c; ++j)
                {					
                    Integer first_row	= B.first_row(j);
                    Integer last_row	= std::min(B.last_row(j),r-1);
                    Integer kb			= first_row;
                    Integer pos_B		= B.first_elem_pos(j);

                    d_c[j]				= nz;

                    while (kb <= last_row)
                    {
                        const VT2& val2  = ptr_B[pos_B];
                        if (!mrd::is_zero(val2))
                        {
                            VTR tmp = ef.template eval_zero<!zero_on_right>(val2);

                            d_r[nz]		= kb;
                            d_x[nz]		= tmp;
                            ++nz;
                        };
                        ++kb; 
                        ++pos_B;
                    };

                    ptr_B += B_ld;
                };

                d_c[c] = nz;
                d.add_memory(-1);
            };

            ret = matcl::Matrix(res,true);
            return;
        };

        const Integer* Ad_c		= Ad.ptr_c();
        const Integer* Ad_r		= Ad.ptr_r();
        const VT1* Ad_x		    = Ad.ptr_x();

        if (fdb == ldb)
        {
            Integer c0          = B.first_col_on_diag(fdb);
            Integer r0          = B.first_row_on_diag(fdb);
            const VT2* ptr_B    = B.rep_ptr() + B.first_elem_diag(fdb);
            Integer rc          = B.diag_length(fdb);

            for (Integer j = 0; j < c0; ++j)
            {
                Integer ka	    = Ad_c[j]; 
                d_c[j]			= nz;

                while (ka < Ad_c[j+1])
                {
                    VTR tmp     = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    d_r[nz]		= Ad_r[ka];
                    d_x[nz]		= tmp;
                    ++nz;
                    ++ka; 
                };
            };

            for (Integer j = c0; j < c0+rc; ++j, ptr_B += B_ld, ++r0)
            {
                Integer ka	    = Ad_c[j]; 
                d_c[j]			= nz;

                while (ka < Ad_c[j+1] && Ad_r[ka] < r0)
                {
                    VTR tmp     = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    d_r[nz]		= Ad_r[ka];
                    d_x[nz]		= tmp;
                    ++nz;
                    ++ka; 
                };

                if (ka < Ad_c[j+1] && Ad_r[ka] == r0)
                {
                    VTR tmp = ef.eval(Ad_x[ka],ptr_B[0]);
                    if (!mrd::is_zero(tmp))
                    {
                        d_r[nz] = r0;
                        d_x[nz]	= tmp;
                        ++nz;
                    };
                    ++ka; 
                }
                else
                {
                    const VT2& val2 = ptr_B[0];
                    if (!mrd::is_zero(val2))
                    {
                        VTR tmp = ef.template eval_zero<!zero_on_right>(val2);

                        d_r[nz]	= r0;
                        d_x[nz]	= tmp;
                        ++nz;
                    };
                };

                while (ka < Ad_c[j+1])
                {
                    VTR tmp     = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    d_r[nz]		= Ad_r[ka];
                    d_x[nz]		= tmp;
                    ++nz;
                    ++ka; 
                };
            };

            for (Integer j = c0 + rc; j < c; ++j)
            {
                d_c[j]	        = nz;

                for (Integer ka = Ad_c[j]; ka < Ad_c[j+1]; ++ka)
                {
                    VTR tmp     = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    d_r[nz]		= Ad_r[ka];
                    d_x[nz]		= tmp;
                    ++nz;
                };
            };

            d_c[c] = nz;

            d.add_memory(-1);
        }
        else
        {
            const VT2* ptr_B    = B.rep_ptr();

            for (Integer j = 0; j < c; ++j)
            {
                Integer ka		    = Ad_c[j]; 
                Integer first_row	= B.first_row(j);
                Integer last_row	= std::min(B.last_row(j),r-1);
                Integer kb			= first_row;
                Integer pos_B		= B.first_elem_pos(j);

                d_c[j]				= nz;

                while ((ka < Ad_c[j+1]) && (kb <= last_row))
                {
                    if (Ad_r[ka] < kb)
                    {
                        VTR tmp     = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                        d_r[nz]		= Ad_r[ka];
                        d_x[nz]		= tmp;
                        ++nz;
                        ++ka; 
                    }
                    else if (Ad_r[ka] > kb)
                    {
                        const VT2& val2 = ptr_B[pos_B];
                        if (!mrd::is_zero(val2))
                        {
                            VTR tmp = ef.template eval_zero<!zero_on_right>(val2);

                            d_r[nz]	= kb;
                            d_x[nz]	= tmp;
                            ++nz;
                        };
                        ++kb; 
                        ++pos_B;
                    }
                    else
                    { 
                        VTR tmp = ef.eval(Ad_x[ka],ptr_B[pos_B]);
                        if (!mrd::is_zero(tmp))
                        {
                            d_r[nz]	= kb;
                            d_x[nz]	= tmp;
                            ++nz;
                        };
                        ++ka; 
                        ++kb;
                        ++pos_B;
                    };
                };
                while (ka < Ad_c[j+1])
                {
                    VTR tmp     = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    d_r[nz]		= Ad_r[ka];
                    d_x[nz]		= tmp;
                    ++nz;
                    ++ka; 
                };
                while (kb <= last_row)
                {
                    const VT2& val2 = ptr_B[pos_B];
                    if (!mrd::is_zero(val2))
                    {
                        VTR tmp = ef.template eval_zero<!zero_on_right>(val2);

                        d_r[nz]	= kb;
                        d_x[nz]	= tmp;
                        ++nz;
                    };
                    ++kb; 
                    ++pos_B;
                };

                ptr_B += B_ld;
            };

            d_c[c] = nz;

            d.add_memory(-1);
        };

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,		//0
            bool ZN,		//0,1
            bool NZ,		//0,1
            bool is_inv
            //class str_ret = struct_dense,
            //class str_1= struct_sparse,
            //class str_2= struct_banded
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,false,ZN,NZ,is_inv,
                        struct_dense,struct_sparse,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using VTR   = typename ret_type::value_type;
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;
        static const bool zero_on_right = !is_inv;

        ti::ti_type<VTR> ret_ti = func.template return_type<VTR,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols());

        Integer r = A.rows();
        Integer c = A.cols();

        if (r == 0 || c == 0)
        {
            ret = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };
        
        VT1 zero1   = md::default_value<VT1>(ti::get_ti(A));
        VT2 zero2   = md::default_value<VT2>(ti::get_ti(B));

        using eval_func = eval_F<functor,is_inv,VTR,VT1,VT2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        VTR zz		= ef.eval(zero1,zero2);

        ret_type res(ret_ti,zz,r,c);

        const sparse_ccs<VT1>& Ad = A.rep();
                
        Integer res_ld      = res.ld();
        Integer B_ld        = B.ld();
        Integer fdb         = B.first_diag();
        Integer ldb         = B.last_diag();

        if (A.nnz() == 0)
        {
            if (fdb == ldb)
            {
                const VT2* ptr_B    = B.rep_ptr() + B.first_elem_diag(fdb);
                VTR* ptr_res        = res.ptr() + B.first_row_on_diag(fdb)
                                    + B.first_col_on_diag(fdb) * res_ld;
                Integer rc          = B.diag_length(fdb);

                for (Integer j = 0; j < rc; ++j)
                {
                    VTR tmp     = ef.template eval_zero<!zero_on_right>(ptr_B[0]);
                    ptr_res[0]  = tmp;
                    ptr_res     += res_ld+1;
                    ptr_B       += B_ld;
                }; 
            }
            else
            {
                const VT2* ptr_B    = B.rep_ptr();
                VTR* ptr_res        = res.ptr();

                for (Integer j = 0; j < c; ++j)
                {
                    Integer first_row   = B.first_row(j);
                    Integer last_row    = B.last_row(j);
                    Integer kb		    = first_row;
                    Integer pos_B	    = B.first_elem_pos(j);

                    while (kb <= last_row)
                    {
                        VTR tmp     = ef.template eval_zero<!zero_on_right>(ptr_B[pos_B]);
                        ptr_res[kb] = tmp;
                        ++kb; 
                        ++pos_B;
                    };
                    ptr_B   += B_ld;
                    ptr_res += res_ld;
                };
            };

            ret = matcl::Matrix(res,true);
            return;
        };

        const Integer* Ad_c		= Ad.ptr_c();
        const Integer* Ad_r		= Ad.ptr_r();
        const VT1* Ad_x		    = Ad.ptr_x();
         
        if (fdb == ldb)
        {
            Integer c0          = B.first_col_on_diag(fdb);
            Integer r0          = B.first_row_on_diag(fdb);
            const VT2* ptr_B    = B.rep_ptr() + B.first_elem_diag(fdb);
            VTR* ptr_res        = res.ptr();
            Integer rc          = B.diag_length(fdb);            

            for (Integer j = 0; j < c0; ++j)
            {
                Integer ka      = Ad_c[j]; 

                while (ka < Ad_c[j+1])
                {
                    VTR tmp     = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    ptr_res[Ad_r[ka]] = tmp;
                    ++ka; 
                };

                ptr_res         += res_ld;
            };

            for (Integer j = c0; j < c; ++j, ++r0)
            {
                Integer ka      = Ad_c[j]; 

                while (ka < Ad_c[j+1] && Ad_r[ka] < r0)
                {
                    VTR tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    ptr_res[Ad_r[ka]] = tmp;
                    ++ka; 
                };
                if (j < rc && ka < Ad_c[j+1] && Ad_r[ka] == r0)
                {
                    VTR tmp     = ef.eval(Ad_x[ka],ptr_B[0]);
                    ptr_res[j]  = tmp;
                    ++ka; 
                }
                else if (j < rc)
                {
                    VTR tmp     = ef.template eval_zero<!zero_on_right>(ptr_B[0]);
                    ptr_res[j]  = tmp;
                };
                
                while (ka < Ad_c[j+1])
                {
                    VTR tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    ptr_res[Ad_r[ka]] = tmp;
                    ++ka; 
                };

                ptr_res += res_ld;
                ptr_B   += B_ld;
            };
        }
        else
        {
            const VT2* ptr_B    = B.rep_ptr();
            VTR* ptr_res        = res.ptr();

            for (Integer j = 0; j < c; ++j)
            {
                Integer ka		    = Ad_c[j]; 
                Integer first_row   = B.first_row(j);
                Integer last_row    = B.last_row(j);
                Integer kb          = first_row;
                Integer pos_B       = B.first_elem_pos(j);

                if (first_row >= r)
                {
                    while (ka < Ad_c[j+1])
                    {
                        VTR tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                        ptr_res[Ad_r[ka]] = tmp;
                        ++ka; 
                    };
                    ptr_B   += B_ld;
                    ptr_res += res_ld;
                    continue;
                };

                while ((ka < Ad_c[j+1]) && (kb <= last_row))
                {
                    if (Ad_r[ka] < kb)
                    {
                        VTR tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                        ptr_res[Ad_r[ka]] = tmp;
                        ++ka; 
                    }
                    else if (Ad_r[ka] > kb)
                    {
                        VTR tmp = ef.template eval_zero<!zero_on_right>(ptr_B[pos_B]);
                        ptr_res[kb] = tmp;
                        ++kb; 
                        ++pos_B;
                    }
                    else
                    { 
                        VTR tmp = ef.eval(Ad_x[ka],ptr_B[pos_B]);
                        ptr_res[kb] = tmp;
                        ++ka; 
                        ++kb;
                        ++pos_B;
                    };
                };
                while (ka < Ad_c[j+1])
                {
                    VTR tmp = ef.template eval_zero<zero_on_right>(Ad_x[ka]);
                    ptr_res[Ad_r[ka]] = tmp;
                    ++ka; 
                };
                while (kb <= last_row)
                {
                    VTR tmp = ef.template eval_zero<!zero_on_right>(ptr_B[pos_B]);
                    ptr_res[kb] = tmp;
                    ++kb; 
                    ++pos_B;
                };


                ptr_B   += B_ld;
                ptr_res += res_ld;
            };
        };

        ret = matcl::Matrix(res,true);
        return;
    };
};

//--------------------------------------------------------------------
//				B - GM
//--------------------------------------------------------------------
template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,
            bool ZN,
            bool NZ,
            bool is_inv
            //class str_ret = struct_banded,
            //class str_1= struct_banded,
            //class str_2= struct_dense
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,ZZ,ZN,NZ,is_inv,
                        struct_banded,struct_banded,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        return eval_bin_functor_impl<ret_type,functor,M2,M1,ZZ,NZ,ZN,!is_inv,
                                        struct_banded,struct_dense,struct_banded>
                                        ::eval(ret,B,A,func);
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,
            bool ZN,
            bool NZ,
            bool is_inv
            //class str_ret = struct_dense,
            //class str_1= struct_banded,
            //class str_2= struct_dense
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,ZZ,ZN,NZ,is_inv,
                        struct_dense,struct_banded,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        return eval_bin_functor_impl<ret_type,functor,M2,M1,ZZ,NZ,ZN,!is_inv,
                                        struct_dense,struct_dense,struct_banded>
                                        ::eval(ret,B,A,func);
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,
            bool ZN,
            bool NZ,
            bool is_inv
            //class str_ret = struct_sparse,
            //class str_1= struct_banded,
            //class str_2= struct_dense
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,ZZ,ZN,NZ,is_inv,
                            struct_sparse,struct_banded,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        return eval_bin_functor_impl<ret_type,functor,M2,M1,ZZ,NZ,ZN,!is_inv,
                                        struct_sparse,struct_dense,struct_banded>
                                        ::eval(ret,B,A,func);
    };
};

//--------------------------------------------------------------------
//				B - SM
//--------------------------------------------------------------------
template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,
            bool ZN,
            bool NZ,
            bool is_inv
            //class str_ret = struct_banded,
            //class str_1= struct_banded,
            //class str_2= struct_sparse
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,ZZ,ZN,NZ,is_inv,
                            struct_banded,struct_banded,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        return eval_bin_functor_impl<ret_type,functor,M2,M1,ZZ,NZ,ZN,!is_inv,
                                        struct_banded,struct_sparse,struct_banded>
                                        ::eval(ret,B,A,func);
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,
            bool ZN,
            bool NZ,
            bool is_inv
            //class str_ret = struct_dense,
            //class str_1= struct_banded,
            //class str_2= struct_sparse
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,ZZ,ZN,NZ,is_inv,
                            struct_dense,struct_banded,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        return eval_bin_functor_impl<ret_type,functor,M2,M1,ZZ,NZ,ZN,!is_inv,
                                        struct_dense,struct_sparse,struct_banded>
                                        ::eval(ret,B,A,func);
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,
            bool ZN,
            bool NZ,
            bool is_inv
            //class str_ret = struct_sparse,
            //class str_1= struct_banded,
            //class str_2= struct_sparse
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,ZZ,ZN,NZ,is_inv,
                            struct_sparse,struct_banded,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        return eval_bin_functor_impl<ret_type,functor,M2,M1,ZZ,NZ,ZN,!is_inv,
                                        struct_sparse,struct_sparse,struct_banded>
                                        ::eval(ret,B,A,func);
    };
};

//--------------------------------------------------------------------
//				B - B
//--------------------------------------------------------------------
template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,		//0
            bool ZN,		//0,1
            bool NZ,		//0,1
            bool is_inv
            //class str_ret = struct_dense,
            //class str_1= struct_banded,
            //class str_2= struct_banded
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,false,ZN,NZ,is_inv,
                            struct_dense,struct_banded,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using VTR   = typename ret_type::value_type;
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;
        static const bool zero_on_right = !is_inv;

        ti::ti_type<VTR> ret_ti = func.template return_type<VTR,is_inv>(ti::get_ti(A),ti::get_ti(B));

        error::check_eeop(A.rows(), A.cols(), B.rows(), B.cols());

        Integer r = A.rows();
        Integer c = A.cols();

        if (r == 0 || c == 0)
        {
            ret = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };
        
        VT1 zero1			= md::default_value<VT1>(ti::get_ti(A));
        VT2 zero2			= md::default_value<VT2>(ti::get_ti(B));

        using eval_func = eval_F<functor,is_inv,VTR,VT1,VT2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        VTR zz				= ef.eval(zero1,zero2);

        ret_type res(ret_ti,zz,r,c);

        Integer A_ld        = A.ld();
        Integer B_ld        = B.ld();
        Integer res_ld      = res.ld();        

        Integer fda         = A.first_diag();
        Integer lda         = A.last_diag();
        Integer fdb         = B.first_diag();
        Integer ldb         = B.last_diag();

        Integer fd  = std::min(fda, fdb);
        Integer ld  = std::max(lda, ldb);

        if (fd > ld)
        {
            ret = matcl::Matrix(res, false);
            return;
        };

        if (fd == ld)
        {
            const VT1* ptr_A    = A.rep_ptr() + A.first_elem_diag(fd);
            const VT2* ptr_B    = B.rep_ptr() + B.first_elem_diag(fd);
            VTR* ptr_res        = res.ptr() + A.first_row_on_diag(fd) 
                                + A.first_col_on_diag(fd) * res_ld;
            Integer rc          = A.diag_length(fd);

            for (Integer i = 0; i < rc; ++i)
            {
                ptr_res[i]  = ef.eval(ptr_A[0],ptr_B[0]);

                ptr_A       += A_ld;
                ptr_B       += B_ld;
                ptr_res     += res_ld;
            };

            ret = matcl::Matrix(res,false);
            return;
        };
         
        const VT1* ptr_A        = A.rep_ptr();
        const VT2* ptr_B        = B.rep_ptr();
        VTR* ptr_res            = res.ptr();

        for (Integer j = 0; j < c; ++j)
        {
            Integer fr1			= A.first_row(j);
            Integer fr2			= B.first_row(j);
            Integer lr1			= std::min(A.last_row(j),r-1);
            Integer lr2			= std::min(B.last_row(j),r-1);
            Integer pa			= A.first_elem_pos(j);
            Integer pb			= B.first_elem_pos(j);

            Integer ka			= fr1;
            Integer kb			= fr2;
            Integer pos_ret		= std::min(ka,kb);

            while ((ka <= lr1) && (kb <= lr2))
            {
                if (ka < kb)
                {
                    VTR tmp = ef.template eval_zero<zero_on_right>(ptr_A[pa]);
                    ptr_res[pos_ret] = tmp;
                    ++pos_ret;
                    ++ka; 
                    ++pa;
                }
                else if (ka > kb)
                {
                    VTR tmp = ef.template eval_zero<!zero_on_right>(ptr_B[pb]);
                    ptr_res[pos_ret] = tmp;
                    ++pos_ret;
                    ++kb; 
                    ++pb;
                }
                else
                { 
                    VTR tmp = ef.eval(ptr_A[pa],ptr_B[pb]);
                    ptr_res[pos_ret] = tmp;
                    ++ka; 
                    ++kb;
                    ++pa;
                    ++pb;
                    ++pos_ret;
                };
            };
            while (ka <= lr1)
            {
                VTR tmp = ef.template eval_zero<zero_on_right>(ptr_A[pa]);
                ptr_res[pos_ret] = tmp;
                ++ka; 
                ++pa;
                ++pos_ret;
            };
            while ( kb <= lr2 )
            {
                VTR tmp = ef.template eval_zero<!zero_on_right>(ptr_B[pb]);
                ptr_res[pos_ret] = tmp;
                ++kb; 
                ++pb;
                ++pos_ret;
            };

            ptr_A   += A_ld;
            ptr_B   += B_ld;
            ptr_res += res_ld;
        };

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            //bool ZZ,		//1
            bool ZN,		//0,1
            bool NZ,		//0,1
            bool is_inv
            //class str_ret = struct_banded,
            //class str_1= struct_banded,
            //class str_2= struct_banded
        >
struct eval_bin_functor_impl<ret_type,functor,M1,M2,true,ZN,NZ,is_inv,
                        struct_banded,struct_banded,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func)
    {
        using VTR   = typename ret_type::value_type;
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;

        static const bool zero_on_right = !is_inv;

        ti::ti_type<VTR> ret_ti = func.template return_type<VTR,is_inv>(ti::get_ti(A),ti::get_ti(B));

        Integer r   = A.rows();
        Integer c   = A.cols();
        Integer fda = A.first_diag();
        Integer lda = A.last_diag();
        Integer fdb = B.first_diag();
        Integer ldb = B.last_diag();
 
        error::check_eeop(r, c, B.rows(), B.cols());
     
        Integer fd  = std::min(fda, fdb);
        Integer ld  = std::max(lda, ldb);

        if (ZN == true)
        {
            fd		= std::max(fd, fda);
            ld		= std::min(ld, lda);
        };

        if (NZ == true)
        {
            fd		= std::max(fd, fdb);
            ld		= std::min(ld, ldb);
        };

        if (fd > ld)
        {
            using VR    = typename md::real_type<VTR>::type;
            using SM    = raw::Matrix<VR, struct_sparse>;
            ret = matcl::Matrix(SM(ret_ti,r, c),false);
            return;
        };

        ret_type res(ret_ti, r, c, fd, ld);

        using eval_func = eval_F<functor,is_inv,VTR,VT1,VT2>;
        eval_func ef(func,ret_ti,ti::get_ti(A), ti::get_ti(B));

        Integer A_ld            = A.ld();
        Integer B_ld            = B.ld();
        Integer res_ld          = res.ld();

        if (fd == ld)
        {
            Integer rc  = A.diag_length(fd);

            const VT1* ptr_A    = A.rep_ptr() + A.first_elem_diag(fd);
            const VT2* ptr_B    = B.rep_ptr() + B.first_elem_diag(fd);
            VTR* ptr_res        = res.rep_ptr() + res.first_elem_diag(fd);

            for (Integer i = 0; i < rc; ++i)
            {
                ptr_res[0]  = ef.eval(ptr_A[0], ptr_B[0]);
                ptr_A       += A_ld;
                ptr_B       += B_ld;
                ptr_res     += res_ld;
            };

            ret = matcl::Matrix(res,true);
            return;
        };

        for (Integer d = fd; d <= ld; ++d)
        {
            VTR* ptr_res        = res.rep_ptr() + res.first_elem_diag(d);

            bool has_A          = A.has_diag(d);
            bool has_B          = B.has_diag(d);

            if (has_A && has_B)
            {
                const VT1* ptr_A= A.rep_ptr() + A.first_elem_diag(d);
                const VT2* ptr_B= B.rep_ptr() + B.first_elem_diag(d);
                Integer s       = A.diag_length(d);

                for (Integer i = 0; i < s; ++i)
                {
                    VTR tmp     = ef.eval(ptr_A[0], ptr_B[0]);
                    ptr_res[0]  = tmp;

                    ptr_A       += A_ld;
                    ptr_B       += B_ld;
                    ptr_res     += res_ld;
                };
            }
            else if (has_A)
            {
                const VT1* ptr_A= A.rep_ptr() + A.first_elem_diag(d);
                Integer s       = A.diag_length(d);

                for (Integer i = 0; i < s; ++i)
                {
                    VTR tmp     = ef.template eval_zero<zero_on_right>(ptr_A[0]);
                    ptr_res[0]  = tmp;

                    ptr_A       += A_ld;
                    ptr_res     += res_ld;
                };
            }
            else if (has_B)
            {
                const VT2* ptr_B= B.rep_ptr() + B.first_elem_diag(d);
                Integer s       = B.diag_length(d);

                for (Integer i = 0; i < s; ++i)
                {
                    VTR tmp     = ef.template eval_zero<!zero_on_right>(ptr_B[0]);
                    ptr_res[0]  = tmp;

                    ptr_B       += B_ld;
                    ptr_res     += res_ld;
                };
            }
        };
     
        ret = matcl::Matrix(res,true);
        return;
    };
};

} //namespace details

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,
            bool ZN,
            bool NZ
        >
struct eval_bin_functor
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B) 
    {
        using str_ret   = typename ret_type::struct_type;
        using str_1     = typename M1::struct_type;
        using str_2     = typename M2::struct_type;

        return details::eval_bin_functor_impl
                <ret_type,functor,M1,M2,ZZ,ZN,NZ,false,str_ret, str_1, str_2>
                ::eval(ret,A,B,functor());
    };
};

template<	class ret_type, 
            class functor,
            class M1, 
            class M2,
            bool ZZ,
            bool ZN,
            bool NZ
        >
struct eval_bin_functor2
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, const functor& func) 
    {
        using str_ret   = typename ret_type::struct_type;
        using str_1     = typename M1::struct_type;
        using str_2     = typename M2::struct_type;

        return details::eval_bin_functor_impl
                    <ret_type,functor,M1,M2,ZZ,ZN,NZ,false,str_ret,str_1,str_2>
                    ::eval(ret,A,B,func);
    };
};

}}

#pragma warning( pop )
