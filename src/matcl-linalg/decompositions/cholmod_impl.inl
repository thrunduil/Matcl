/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "matcl-internals/func/converter.h"
#include "matcl-internals/base/pv_constructor.h"
#include "matcl-core/general/fwd_decls.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-internals/func/inplace.h"

namespace matcl { namespace details
{

template<class VR>
struct cholmod_diag
{
    static inline void eval(cholmod_return_type& ret, const matcl::Matrix& handle, 
                                correction_type::type corr, Real tol_in)
    {
        matcl::Matrix A = real(handle);

        matcl::Matrix U,P;
        tie(U,P) = sort2(handle,1,false);
    
        U = full(U);

        Integer N   = (Integer)handle.numel();    
        VR* ptr     = U.get_array_unique<VR>();

        if (N == 0)
        {
            ret = cholmod_return_type(A, details::pv_constructor::make(P), 0, 0.0);
            return;
        };

        VR tol      = VR();

        if (tol_in < 0)
        {
            VR max_diag_abs = max(abs(ptr[0]), abs(ptr[N-1]));

            if (corr == correction_type::TYPE_I)
	        {
                tol = N * constants::eps<VR>() * max_diag_abs;
	        }
	        else
	        {
                VR tau = std::pow(constants::eps<VR>(), VR(2./3.0));
		        tol = tau * max_diag_abs;
	        };	
        }
        else
        {
            tol     = VR(tol_in);
        };


        Integer rank    = N;
        Real norm       = 0;
        VR sqrt_tol     = sqrt(tol);
    
        for(Integer i = 0; i < N; ++i)
        {
            if (ptr[i] > tol)
            {
                ptr[i] = sqrt(ptr[i]);
            }
            else
            {
                if (corr == correction_type::TYPE_I)
                {
                    VR d    = max(abs(ptr[i]), tol);
                    norm    = max(norm,d-ptr[i]);
                    ptr[i]  = sqrt(d);
                }
                else
                {
                    VR d    = tol - ptr[i];
                    norm    = max(norm,d);
                    ptr[i]  = sqrt_tol;
                };
                --rank;
            };
        };

        ret = cholmod_return_type(spdiag(U),details::pv_constructor::make(P),rank,norm);
        return;
    };
};

template<class T>
struct cholmod_struct
{
    static void eval(cholmod_return_type& ret, const matcl::Matrix& handle, const T& mat, bool upper, 
                        correction_alg::type alg, correction_type::type corr, Real tol)
    {
        using V     = typename T::value_type;
        using VR    = typename md::real_type<V>::type;

        Integer N = mat.rows();
        Integer M = mat.cols();

        if (N != M)
            throw error::square_matrix_required(N, M);

        if (mat.get_struct().is_id())
        {
            permvec P   = permvec::identity(N);
            ret         = cholmod_return_type(Matrix(mat,false),P,N,0);
            return;
        }
        else if (mat.get_struct().is_diag())
        {
            return cholmod_diag<VR>::eval(ret, get_diag(handle), corr, tol);
        }

        if ( mat.nnz() == 0)
        {
            permvec P = permvec::identity(N);
            if (tol < 0)
                tol = 0.0;

            if (tol == 0.0)
            {
                ret = cholmod_return_type(Matrix(mat,false),P,0,tol);
                return;
            };

            using full_matrix = raw::Matrix<VR,struct_dense>;
            full_matrix d = full_matrix(mat.get_type(), VR(tol), M, 1);

            Matrix dm = spdiag(Matrix(d,false));
            ret = cholmod_return_type(dm, P, 0, tol);
            return;
        }

        using FullMatrix = raw::Matrix<V,struct_dense>;
        FullMatrix A = raw::converter<FullMatrix,T>::eval(mat);
        return cholmod_struct<FullMatrix>::eval(ret,handle,A,upper,alg,corr,tol);
    };
};

template<class V>
struct cholmod_struct<raw::Matrix<V,struct_dense>>
{
    using Mat = raw::Matrix<V,struct_dense>;

    static void eval(cholmod_return_type& ret, const matcl::Matrix& handle, const Mat& mat, bool upper, 
                        correction_alg::type alg, correction_type::type corr, Real tol)
    {
        using VR    = typename md::real_type<V>::type;

        Integer N = mat.rows();
        Integer M = mat.cols();

        if (N != M)
            throw error::square_matrix_required(N, M);

        if (mat.get_struct().is_id())
        {
            permvec P   = permvec::identity(N);
            ret         = cholmod_return_type(Matrix(mat,false),P,N,0);
            return;
        }
        else if (mat.get_struct().is_diag())
        {
            return cholmod_diag<VR>::eval(ret, get_diag(handle), corr, tol);
        }

        matcl::details::cholmod_options opts;
        switch(alg)
        {
            case correction_alg::NONE:  opts.corr_alg(cholmod_options::E_NONE); break;
            case correction_alg::GMW:   opts.corr_alg(cholmod_options::E_GMW); break;
            case correction_alg::SE:    opts.corr_alg(cholmod_options::E_SE); break;
        };

        switch(corr)
        {
            case correction_type::TYPE_I:    opts.corr_type(cholmod_options::TYPE_I); break;
            case correction_type::TYPE_II:   opts.corr_type(cholmod_options::TYPE_II); break;
        };

        opts.tol(tol);

        Mat Ac(mat.get_type());
        if (mat.get_refstr()->is_unique() == true)
            Ac.assign_to_fresh(mat);
        else
            Ac.assign_to_fresh(mat.copy());

        Ac.get_struct().reset();

        Integer lda     = Ac.ld();        
        V* ptr          = Ac.ptr();

        md::cholmod_matrix<V> raw_mat(ptr,N,N,lda,
                   upper ? md::triang_part::UPPER: md::triang_part::LOWER);

        md::cholmod<V> chol(raw_mat,opts);

        Integer* piv    = chol.permutation().get();

        matcl::Matrix P = matcl::izeros(1,N);
        Integer step_I  = chol.rank();
        Real norm       = chol.norm_correction();

        Integer* ptr_P = P.get_array_unique<Integer>();
        for(Integer i = 0; i < N; ++i)
        {
            ptr_P[i] = piv[i]+1;
        };

        Matrix U;

        if (upper)
            matcl::raw::inplace::make_triu<Mat>::eval(U,Ac,0);
        else
            matcl::raw::inplace::make_tril<Mat>::eval(U,Ac,0);

        ret = cholmod_return_type(U,details::pv_constructor::make(P),step_I,norm);
        return;
    };
};

template<class T>
struct cholmod_impl
{
    static void eval(cholmod_return_type& ret, const matcl::Matrix& handle, const T& mat, bool upper, 
                        correction_alg::type alg, correction_type::type corr, Real tol)
    {
        return cholmod_struct<T>::eval(ret, handle, mat, upper, alg, corr, tol);
    };
};

template<class S>
struct cholmod_impl<raw::Matrix<Integer,S>>
{
    static void eval(cholmod_return_type& ret, const matcl::Matrix& handle, const raw::Matrix<Integer,S>& mat, 
                        bool upper, correction_alg::type alg, correction_type::type corr, Real tol)
    {
        using matrix_type   = raw::Matrix<Integer,S>;
        using FullMatrix    = raw::Matrix<Real,struct_dense>;

        FullMatrix A = raw::converter<FullMatrix,matrix_type>::eval(mat);
        return cholmod_impl<FullMatrix>::eval(ret,handle,A,upper,alg,corr,tol);
    };
};

template<class S>
struct cholmod_impl<raw::Matrix<Object,S>>
{
    static void eval(cholmod_return_type&, const matcl::Matrix& , const raw::Matrix<Object,S>&, bool,
                        correction_alg::type , correction_type::type , Real )
    {
        throw error::object_value_type_not_allowed("cholmod");
    };
};

};};