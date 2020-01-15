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

#include "matcl-internals/base/utils.h"
#include "matcl-matrep/details/struct_flag_predefined.h"

#pragma warning( push )
#pragma warning(disable:4100)	// unreferenced formal parameter (strange warning)
#pragma warning(disable:4723)	// potential divide by 0

namespace matcl { namespace raw 
{

namespace mr    = matcl::raw;
namespace mrd   = matcl::raw::details;
namespace md    = matcl::details;

namespace details
{

// inplace is allowed; refcount must be increased for 
// nontemporary objects; TODO
template<	class ret_type, 
            class in_type, 
            class ret_str,
            class in_str
        >
struct eval_functor_2_inpl {};

#pragma warning(push)
#pragma warning(disable:4244) //possible loss of data

template<class VR, class VI, class Func>
struct eval_inplace_dense
{
    using Mat   = raw::Matrix<VI, struct_dense>;

    static void eval_diag(matcl::Matrix& ret, const Mat &x, const Func& f)
    {
        //this case should aready be removed
        (void)ret;
        (void)x;
        (void)f;
        return;
    };

    static void eval(matcl::Matrix& ret, const Mat &x, const Func& f)
    {
        //this case should aready be removed
        (void)ret;
        (void)x;
        (void)f;
        return;
    };
};

template<class V, class Func>
struct eval_inplace_dense<V,V,Func>
{
    using Mat   = raw::Matrix<V, struct_dense>;

    static void eval_diag(matcl::Matrix& ret, const Mat &x, const Func& f)
    {
        Mat res             = x.make_unique();
        V* ptr_ret          = res.ptr();

        Integer r           = x.rows();
        Integer c           = x.cols();
        Integer rc          = std::min(r,c);
        Integer res_step    = res.ld() + 1;

        for (Integer j = 0; j < rc; ++j) 
        {
            ptr_ret[0]  = V(f.eval(ptr_ret[0]));
            ptr_ret     += res_step;
        };

        res.set_struct(struct_flag());
        ret = matcl::Matrix(res,false);
        return;
    };

    static void eval(matcl::Matrix& ret, const Mat &x, const Func& f)
    {
        Mat res             = x.make_unique();
        V* ptr_res          = res.ptr();

        Integer res_ld      = res.ld();

        for(Integer j = 0; j < res.cols(); ++j)    
        {
            for(Integer i = 0; i < res.rows(); ++i)
            {
                ptr_res[i] = V(f.eval(ptr_res[i])); 
            };

            ptr_res += res_ld;
        };

        res.set_struct(struct_flag());

        ret = matcl::Matrix(res,false);
        return;
    };
};

template<	class ret_type, 
            class in_type
        >
struct eval_functor_2_inpl<ret_type,in_type,struct_dense,struct_dense> 
{
    using ti_type = typename ti::get_ti_type<ret_type>::type;    

    template<class Func>
    static void eval(matcl::Matrix& ret, ti_type ret_ti, const in_type &x, const Func& f)
    {		
        using val_type  = typename ret_type::value_type;
        using in_val    = typename in_type::value_type;
        
        if (x.get_struct().is_diag() == true)
        {            
            in_val Z    =  md::default_value<in_val>(x.get_type());
            val_type fZ(f.eval(Z));

            if (mrd::is_zero(fZ) == false)
            {
                ret_type res(ret_ti, fZ, x.rows(), x.cols());

                Integer r   = x.rows();
                Integer c   = x.cols();
                Integer k   = std::min(r,c);

                Integer res_step    = res.ld() + 1;
                Integer in_step     = x.ld() + 1;
                val_type* ptr_ret   = res.ptr();
                const in_val* ptr_x = x.ptr();

                for (Integer i = 0; i < k; ++i)
                {
                    ptr_ret[0]      = val_type(f.eval(ptr_x[0]));
                    ptr_ret         += res_step;
                    ptr_x           += in_step;
                };

                ret = matcl::Matrix(res,false);
                return;
            }
            else
            {
                if (x.is_unique() == true 
                    && std::is_same<val_type, in_val>::value 
                    && std::is_same<val_type, Object>::value == false
                    )
                {
                    //inplace version
                    return eval_inplace_dense<val_type,in_val,Func>::eval_diag(ret, x, f);
                }

                using ret_band_type = Matrix<val_type, struct_banded>;
                ret_band_type res(ret_ti, x.rows(), x.cols(), 0, 0);

                val_type* ptr_ret   = res.rep_ptr();
                const in_val* ptr_x = x.ptr();

                Integer r           = x.rows();
                Integer c           = x.cols();
                Integer rc          = std::min(r,c);
                Integer in_step     = x.ld() + 1;

                for (Integer j = 0; j < rc; ++j) 
                {
                    ptr_ret[j]  = val_type(f.eval(ptr_x[0]));
                    ptr_x       += in_step;
                };

                ret = matcl::Matrix(res,false);
                return;
            };
        }
        else
        {
            if (x.is_unique() == true 
                && std::is_same<val_type, in_val>::value 
                && std::is_same<val_type, Object>::value == false
                )
            {
                //inplace version
                return eval_inplace_dense<val_type,in_val,Func>::eval(ret, x, f);
            }

            ret_type res(ret_ti, x.rows(), x.cols()); 

            val_type* ptr_res   = res.ptr();
            const in_val* ptr_x = x.ptr();

            Integer res_ld      = res.ld();
            Integer in_ld       = x.ld();

            for(Integer j = 0; j < x.cols(); ++j)    
            {
                for(Integer i = 0; i < x.rows(); ++i)
                {
                    ptr_res[i] = val_type(f.eval(ptr_x[i])); 
                };

                ptr_x   += in_ld;
                ptr_res += res_ld;
            };

            ret = matcl::Matrix(res,false);
            return;
        };	    
    };
};

template<	class ret_type, 
            class in_type
        >
struct eval_functor_2_inpl<ret_type,in_type,struct_dense,struct_sparse> 
{
    using ti_type = typename ti::get_ti_type<ret_type>::type;

    template<class Func>
    static void eval(matcl::Matrix& ret, ti_type ret_ti,const in_type &mat, const Func& f)
    {
        using value_type        = typename in_type::value_type;
        using value_type_ret    = typename ret_type::value_type;

        Integer Ar = mat.rows();
        Integer Ac = mat.cols();

        if (Ar == 0 || Ac == 0)
        {
            ret_type out(ret_ti, Ar, Ac);
            ret = matcl::Matrix(out,false);
            return;
        };

        value_type Z_in     = md::default_value<value_type>(ti::get_ti(mat));
        value_type_ret tmp(f.eval(Z_in));

        ret_type res(ret_ti, tmp, mat.rows(), mat.cols());

        const sparse_ccs<value_type>& Ad = mat.rep();
        const Integer* Ad_c		    = Ad.ptr_c();
        const Integer* Ad_r		    = Ad.ptr_r();
        const value_type* Ad_x	    = Ad.ptr_x();
        value_type_ret* ptr_res     = res.ptr();

        Integer res_ld              = res.ld();

        for (Integer j = 0; j < Ac; ++j)
        {
            for (Integer k = Ad_c[j] ; k < Ad_c[j + 1] ; ++k)
            {
                value_type_ret tmp2 = value_type_ret(f.eval(Ad_x[k]));
                ptr_res[Ad_r[k]]    = tmp2;
            };

            ptr_res += res_ld;
        };
 
        ret = matcl::Matrix(res,false);
        return;
    };
};

template<	class ret_type, 
            class in_type
        >
struct eval_functor_2_inpl<ret_type,in_type,struct_dense,struct_banded> 
{
    using ti_type = typename ti::get_ti_type<ret_type>::type;

    template<class Func>
    static void eval(matcl::Matrix& ret, ti_type ret_ti, const in_type& A, const Func& f)
    {
        using VT    = typename in_type::value_type;
        using VTR   = typename ret_type::value_type;

        Integer r   = A.rows();
        Integer c   = A.cols();

        if (r == 0 || c == 0)
        {
            ret_type out(ret_ti, r, c);
            ret = matcl::Matrix(out,false);
            return;
        };

        VT Z_in     = md::default_value<VT>(ti::get_ti(A));
        VTR res_0(f.eval(Z_in));

        ret_type res(ret_ti,res_0, A.rows(),A.cols());
    
        Integer fda = A.first_diag();
        Integer lda = A.last_diag();

        if (fda == lda)
        {
            Integer A_ld    = A.ld();
            Integer res_ld  = res.ld();

            const VT* ptr_A = A.rep_ptr() + A.first_elem_diag(fda);
            VTR* ptr_res    = res.ptr() + A.first_row_on_diag(fda) 
                            + A.first_col_on_diag(fda) * res_ld;
            Integer rc      = A.diag_length(fda);

            for (Integer j = 0; j < rc; ++j) 
            {
                ptr_res[j]  = VTR(f.eval(ptr_A[0]));
                ptr_A       += A_ld;
                ptr_res     += res_ld;
            };

            ret = matcl::Matrix(res,false);
            return;
        };

        const VT* ptr_A     = A.rep_ptr();
        VTR* ptr_res        = res.ptr();

        Integer A_ld        = A.ld();
        Integer res_ld      = res.ld();

        for (Integer j = 0; j < c; ++j) 
        {
            Integer row_f   = A.first_row(j);
            Integer row_l   = A.last_row(j);
            Integer row_p   = A.first_elem_pos(j);

            for (Integer i = row_f; i <= row_l; ++i, ++row_p)
            {
                ptr_res[i] = VTR(f.eval(ptr_A[row_p]));
            };	

            ptr_res += res_ld;
            ptr_A   += A_ld;
        };

        ret = matcl::Matrix(res,false);
        return;
    };
};

#pragma warning(pop)

template<	class ret_type, 
            class in_type
        >
struct eval_functor_2_inpl<ret_type,in_type,struct_sparse,struct_dense> 
{
    using ti_type = typename ti::get_ti_type<ret_type>::type;

    template<class Func>
    static void eval(matcl::Matrix& ret, ti_type ret_ti, const in_type &mat, const Func& f)
    {
        using value_type    = typename in_type::value_type;
        using val_type_ret  = typename ret_type::value_type;		        

        Integer Ar = mat.rows();
        Integer Ac = mat.cols();
        Integer nz = 0;

        if (Ar == 0 || Ac == 0)
        {
            ret_type out(ret_ti, Ar, Ac);
            ret = matcl::Matrix(out,false);
            return;
        };

        value_type Z_in     = md::default_value<value_type>(ti::get_ti(mat));
        val_type_ret res_0  = f.eval(Z_in);

        if (mrd::is_zero(res_0) == false)
        {
            using full_ret  = Matrix<val_type_ret,struct_dense>;
            eval_functor_2_inpl<full_ret, in_type, struct_dense, struct_dense>
                                            ::eval(ret, ret_ti,mat,f);

            return;
        };

        ret_type res(ret_ti,mat.rows(), mat.cols(), Ar); 

        sparse_ccs<val_type_ret>& d = res.rep();;
        Integer * d_c		        = d.ptr_c();
        Integer * d_r		        = d.ptr_r();
        val_type_ret * d_x	        = d.ptr_x();
        const value_type* ptr_mat   = mat.ptr();
        Integer mat_ld              = mat.ld();

        for (Integer j = 0; j < Ac; ++j)
        {
            if (nz + Ar > d.nzmax()) 
            {
                d.add_memory( d.nzmax() + Ar);

                d_r				= d.ptr_r();
                d_x				= d.ptr_x();
            };

            d_c[j]  = nz;

            for (Integer k = 0 ; k < Ar ; ++k)
            {
                val_type_ret val = f.eval(ptr_mat[k]);

                if (mrd::is_zero(val))
                    continue;

                d_r[nz] = k;
                d_x[nz] = val;

                ++nz;
            };

            ptr_mat += mat_ld;
        };

        d_c[Ac] = nz;

        d.add_memory(-1);
        ret = matcl::Matrix(res,false);
        return;
    };
};

template<class VR, class VI, class Func>
struct eval_inplace_sparse
{
    using Mat   = raw::Matrix<VI, struct_sparse>;

    static void eval(matcl::Matrix& ret, const Mat &x, const Func& f)
    {
        (void)ret;
        (void)x;
        (void)f;
        //this case should aready be removed
    };
};

template<class V, class Func>
struct eval_inplace_sparse<V,V,Func>
{
    using Mat   = raw::Matrix<V, struct_sparse>;

    static void eval(matcl::Matrix& ret, const Mat &x, const Func& f)
    {
        Mat res                 = x.make_unique();
        sparse_ccs<V>& d        = res.rep();;
        Integer* d_c			= d.ptr_c();
        V* d_x		            = d.ptr_x();
        Integer N               = res.cols();

        Integer kf              = d_c[0];
        Integer kl              = d_c[N];

        for (Integer k = kf; k < kl; ++k)
            d_x[k]  = V(f.eval(d_x[k]));

        res.set_struct(struct_flag());

        ret = matcl::Matrix(res,false);
        return;
    };
};

template<	class ret_type, 
            class in_type
        >
struct eval_functor_2_inpl<ret_type,in_type,struct_sparse,struct_sparse> 
{
    using ti_type = typename ti::get_ti_type<ret_type>::type;
    
    template<class Func>
    static void eval(matcl::Matrix& ret, ti_type ret_ti, const in_type & mat, const Func& f)  
    { 
        using value_type    = typename in_type::value_type;
        using val_type_ret  = typename ret_type::value_type;		

        Integer Ar = mat.rows();
        Integer Ac = mat.cols();
        Integer nz = 0;

        if (Ar == 0 || Ac == 0)
        {
            ret_type out(ret_ti, Ar, Ac);
            ret = matcl::Matrix(out,false);
            return;
        };

        value_type Z_in = md::default_value<value_type>(ti::get_ti(mat));

        val_type_ret res_0(f.eval(Z_in));

        if (mrd::is_zero(res_0) == false)
        {
            using full_ret  = Matrix<val_type_ret,struct_dense>;
            eval_functor_2_inpl<full_ret,in_type,struct_dense,struct_sparse>
                                        ::eval(ret, ret_ti,mat,f);
            return;
        };

        if (mat.is_unique() == true 
            && std::is_same<value_type, val_type_ret>::value 
            && std::is_same<value_type, Object>::value == false
            )
        {
            //inplace version
            return eval_inplace_sparse<val_type_ret,value_type,Func>::eval(ret, mat, f);
        }

        ret_type res(ret_ti, mat.rows(), mat.cols(), mat.nnz()); 		

        const sparse_ccs<value_type>& Ad = mat.rep();
        sparse_ccs<val_type_ret>& d      = res.rep();;
        const Integer* Ad_c		    = Ad.ptr_c();
        const Integer* Ad_r		    = Ad.ptr_r();
        const value_type* Ad_x	    = Ad.ptr_x();

        Integer* d_c			    = d.ptr_c();
        Integer* d_r			    = d.ptr_r();
        val_type_ret* d_x		    = d.ptr_x();

        for (Integer j = 0; j < Ac; ++j)
        {
            d_c[j] = nz;

            for (Integer k = Ad_c[j] ; k < Ad_c[j + 1] ; ++k)
            {
                val_type_ret val(f.eval(Ad_x[k]));

                d_r[nz] = Ad_r[k];
                d_x[nz] = val;

                ++nz;
            };
        }

        d_c[Ac] = nz;

        d.add_memory(-1);
        ret = matcl::Matrix(res,false);
        return;
    };
};

template<	class ret_type, 
            class in_type
        >
struct eval_functor_2_inpl<ret_type,in_type,struct_sparse,struct_banded> 
{
    using ti_type = typename ti::get_ti_type<ret_type>::type;

    template<class Func>
    static void eval(matcl::Matrix& ret, ti_type ret_ti,const in_type& A, const Func& f)
    {
        using value_type    = typename in_type::value_type;
        using val_type_ret  = typename ret_type::value_type;		

        if (A.rows() == 0 || A.cols() == 0)
        {
            ret_type out(ret_ti, A.rows(), A.cols());
            ret = matcl::Matrix(out,false);
            return;
        }

        value_type Z_in = md::default_value<value_type>(ti::get_ti(A));
        val_type_ret res_0(f.eval(Z_in));

        if (mrd::is_zero(res_0))
        {
            using band_matrix   = Matrix<val_type_ret,struct_banded>;
            eval_functor_2_inpl<band_matrix,in_type,struct_banded,struct_banded>
                                            ::eval(ret, ret_ti,A,f);
            return;
        }
        else
        {
            using full_ret  = Matrix<val_type_ret,struct_dense>;
            eval_functor_2_inpl<full_ret, in_type, struct_dense, struct_banded>
                                            ::eval(ret, ret_ti,A,f);
            return;
        };
    };
};

template<	class ret_type, 
            class in_type
        >
struct eval_functor_2_inpl<ret_type,in_type,struct_banded,struct_dense> 
{};

template<	class ret_type, 
            class in_type
        >
struct eval_functor_2_inpl<ret_type,in_type,struct_banded,struct_sparse> 
{};

template<class VR, class VI, class Func>
struct eval_inplace_band
{
    using Mat   = raw::Matrix<VI,struct_banded>;

    static void eval(matcl::Matrix& ret, const Mat& A, const Func& f)
    {
        //this case should already be removed
        (void)ret;
        (void)A;
        (void)f;
        return;
    };
};

template<class V, class Func>
struct eval_inplace_band<V,V,Func>
{
    using Mat   = raw::Matrix<V,struct_banded>;

    static void eval(matcl::Matrix& ret, const Mat& A, const Func& f)
    {
        Mat res         = A.make_unique();

        Integer c       = res.cols();
        V* ptr_ret      = res.rep_ptr();
        Integer res_ld  = res.ld();

        for (Integer j = 0; j < c; ++j) 
        {
            Integer row_f = res.first_row(j);
            Integer row_l = res.last_row(j);
            Integer row_p = res.first_elem_pos(j);

            for (Integer i = row_f; i <= row_l; ++i, ++row_p)
                ptr_ret[row_p] = V(f.eval(ptr_ret[row_p]));

            ptr_ret += res_ld;
        };

        res.set_struct(struct_flag());

        ret = matcl::Matrix(res,false);
        return;
    };
};

template<	class ret_type, 
            class in_type
        >
struct eval_functor_2_inpl<ret_type,in_type,struct_banded,struct_banded> 
{
    using ti_type = typename ti::get_ti_type<ret_type>::type;
    
    template<class Func>
    static void eval(matcl::Matrix& ret, ti_type ret_ti, const in_type& A, const Func& f)
    {
        using VT    = typename in_type::value_type;
        using VTR   = typename ret_type::value_type;		

        if (A.rows() == 0 || A.cols() == 0)
        {
            ret_type out(ret_ti, A.rows(), A.cols(),0,0);
            ret = matcl::Matrix(out,false);
            return;
        };

        VT Z_in = md::default_value<VT>(ti::get_ti(A));
        VTR res_0(f.eval(Z_in));

        if (mrd::is_zero(res_0) == false)
        {
            using full_ret  = Matrix<VTR,struct_dense>;
            eval_functor_2_inpl<full_ret,in_type,struct_dense,struct_banded>
                                        ::eval(ret, ret_ti, A, f);

            return;
        };

        if (A.is_unique() == true && std::is_same<VTR, VT>::value 
                && std::is_same<VTR, Object>::value == false)
        {
            //inplace version
            return eval_inplace_band<VTR,VT,Func>::eval(ret, A, f);
        }

        Integer fda         = A.first_diag();
        Integer lda         = A.last_diag();
        ret_type res(ret_ti, A.rows(), A.cols(), fda, lda);

        Integer c           = A.cols();
        const VT* ptr_A     = A.rep_ptr();
        VTR* ptr_ret        = res.rep_ptr();

        Integer A_ld        = A.ld();
        Integer res_ld      = res.ld();

        for (Integer j = 0; j < c; ++j) 
        {
            Integer row_f   = A.first_row(j);
            Integer row_l   = A.last_row(j);
            Integer row_p   = A.first_elem_pos(j);

            for (Integer i = row_f; i <= row_l; ++i, ++row_p)
            {
                ptr_ret[row_p] = VTR(f.eval(ptr_A[row_p]));
            };	

            ptr_A   += A_ld;
            ptr_ret += res_ld;
        };

        ret = matcl::Matrix(res,false);
        return;
    };
};

};

template<	class ret_type, 
            class in_type
        >
struct eval_functor_impl
{ 
    template<class Func>
    static void eval(matcl::Matrix& ret, const in_type &x,const Func& f) 
    {
        using ret_str   = typename ret_type::struct_type;
        using in_str    = typename in_type::struct_type;
        using in_val    = typename in_type::value_type;
        using ret_ti_t  = typename ti::get_ti_type<ret_type>::type;

        ret_ti_t ret_ti = ti::get_return_ti<ret_ti_t>(Func::name(), x.get_type());

        details::eval_functor_2_inpl<ret_type,in_type,ret_str,in_str>::eval(ret, ret_ti, x, f);

        in_val Z        = md::default_value<in_val>(x.get_type());
        bool is_zero_id = mrd::is_zero(f.eval(Z))? true:false;
        struct_flag so  = md::predefined_struct::eval_struct(x.get_struct(),is_zero_id);

        so.add(ret.get_struct());
        ret.set_struct(so);

        return;
    };
};

template<	class ret_type, 
            class in_type
        >
struct eval_functor_2_main
{ 
    using ti_type = typename ti::get_ti_type<ret_type>::type;

    template<class Func>
    static void eval(matcl::Matrix& ret, ti_type ret_ti,const in_type &x,const Func& f) 
    {
        using ret_str   = typename ret_type::struct_type;
        using ret_val   = typename ret_type::value_type;
        using in_str    = typename in_type::struct_type;
        using in_val    = typename in_type::value_type;

        if (f.is_special_case())
        {
            f.template eval_special_case<ret_type,in_type>(ret, ret_ti,x);

            in_val Z        = md::default_value<in_val>(x.get_type());
            bool is_zero_id = mrd::is_zero(f.eval(Z))? true:false;
            struct_flag so  = md::predefined_struct::eval_struct(x.get_struct(),is_zero_id);

            so.add(ret.get_struct());
            ret.set_struct(so);

            return;
        };

        details::eval_functor_2_inpl<ret_type,in_type,ret_str,in_str>::eval(ret, ret_ti,x,f);

        in_val Z            = md::default_value<in_val>(x.get_type());
        bool is_zero_id     = mrd::is_zero(f.eval(Z))? true:false;
        struct_flag so      = md::predefined_struct::eval_struct(x.get_struct(),is_zero_id);

        so.add(ret.get_struct());
        ret.set_struct(so);

        return;
    };
};

};};

#pragma warning( pop )
