/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include "matcl-linalg/decompositions/lu/lu_impl.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-linalg/decompositions/lu/lu_struct.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-internals/base/pv_constructor.h"
#include "matcl-core/utils/workspace.h"

namespace matcl { namespace details
{

template<class V> 
void lu_dense<V>::eval(lu_return_type& ret, const raw::Matrix<V,struct_dense>& A, const options& opts)
{
    using VR            = typename md::real_type<V>::type;
    using Mat_D         = raw::Matrix<V,struct_dense>;
    using pivot_type    = opt::linsolve::pivot_type;

    Integer M           = A.rows();
    Integer N           = A.cols();
    struct_flag flag    = A.get_struct();    

    //-------------------------------------------------------------------
    //                  fast exit
    //-------------------------------------------------------------------
    if ( M == 0 || N == 0)
    {
        ti::ti_empty ti;
        Mat_D L(ti);
        Mat_D U(ti);

        if (M >= N) 
        { 
            // N = 0
            L.reset_unique(M, N); 
            U.reset_unique(N, N); 
        }
        else 
        { 
            // M = 0
            L.reset_unique(M, M); 
            U.reset_unique(M, N); 
        };

        permvec p = permvec::identity(A.rows());
        permvec q = permvec::identity(A.cols());

        ret = lu_return_type(Matrix(L,false), Matrix(U,false), p, q);
        return;
    }
     
    value_code vc           = matrix_traits::value_code<V>::value;

    if (flag.is_id())
    {
        Integer K           = std::min(M, N);
        matcl::Matrix L     = speye(M,K,vc);
        matcl::Matrix U     = speye(K,N,vc);

        permvec p = permvec::identity(A.rows());
        permvec q = permvec::identity(A.cols());

        ret = lu_return_type(L,U,p,q);
        return;
    }

    Integer A_ldiags    = get_ld(Matrix(A,false),0);
    Integer A_udiags    = get_ud(Matrix(A,false),0);    

    Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
    pivot_type piv      = (pivot_type)piv_int;
    VR tol              = (VR)opts.get_option<Real>(opt::linsolve::tol());
    bool isp            = (piv == pivot_type::partial);
    bool allow_nonunit  = opts.get_option<bool>(opt::linsolve::allow_nonunit_L());

    if (A_ldiags == 0 && A_udiags == 0)
        return lu_diag<V,struct_dense>::eval(ret, A, piv, tol);

    if (A_udiags == 0)
    {
        if (isp && allow_nonunit)
            return lu_tril_nonunit(ret, Matrix(A,false));
    }
    else if (A_ldiags == 0)
    {
        if (isp)
            return lu_triu(ret, Matrix(A,false));
    }
    
    //-------------------------------------------------------------------
    //                  factor
    //-------------------------------------------------------------------
    Mat_D A2        = A.make_unique();
    V *ptr_A        = A2.ptr();

    lapack::i_type info;
    raw::integer_dense raw_piv(ti::ti_empty(), 1, M);    

    Matrix P        = irange(1,M);
    Matrix Q        = irange(1,N);

    Integer *ptr_Q  = Q.get_array_unique<Integer>();
    Integer *ptr_P  = P.get_array_unique<Integer>();    
    Integer rank;

    switch(piv)
    {
        case pivot_type::partial:
        {
            lapack::getrf_rec(M, N, lap(ptr_A), A2.ld(), lap(raw_piv.ptr()), &info);
            if (info < 0)
            {
                matcl_assert(0,"invalid argument passed to getrf_rec");
                throw error::error_general("invalid argument");
            };
            rank    = std::min(M, N);
            break;
        }
        case pivot_type::rook:
        {
            VR TOLC = (VR)opts.get_option<Real>(opt::linsolve::tol_c());
            VR TOLR = (VR)opts.get_option<Real>(opt::linsolve::tol_r());
            VR TOLV = (VR)tol;
            
            V work_query;
            lapack::getrfr(M, N, lap(ptr_A), A2.ld(), lap(raw_piv.ptr()), lap(ptr_Q), TOLC, TOLR, TOLV, 
                           lap(&work_query), -1, &info);

            Integer lwork   = (Integer)real(work_query);

            using VTR       = pod_type<V>;
            using workspace = matcl::pod_workspace<VTR>;
            workspace WORK  = workspace(lwork);
            V* ptr_WORK     = reinterpret_cast<V*>(WORK.ptr());

            lapack::getrfr(M, N, lap(ptr_A), A2.ld(), lap(raw_piv.ptr()), lap(ptr_Q), TOLC, TOLR, TOLV, 
                           lap(ptr_WORK), lwork, &info);

            if (info < 0)
                throw error::error_general("invalid argument passed to getrfr");

            rank            = info;
            break;
        };
        case pivot_type::complete:
        {
            VR TOLV = (VR)tol;

            lapack::getrfc(M, N, lap(ptr_A), A2.ld(), lap(raw_piv.ptr()), lap(ptr_Q), TOLV, &info);
            
            if (info < 0)
                throw error::error_general("invalid argument passed to getrfc");

            rank            = info;
            break;
        }
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };    

    //-------------------------------------------------------------------
    //                  extract
    //-------------------------------------------------------------------
    Integer K           = std::min(M, N);
    Integer* ptr_piv    = raw_piv.ptr();    

    for (Integer i = 0; i < K; i++)
    {
        Integer sw      = ptr_P[i];
        Integer pos     = ptr_piv[i]-1;
        ptr_P[i]        = ptr_P[pos];
        ptr_P[pos]      = sw;
    }

    if (M > N)
    {
        ti::ti_empty ti;

        Mat_D U(ti, K, N);

        V *ptr_U        = U.ptr(); 
        ptr_A           = A2.ptr();

        for (Integer j = 0; j < N; ++j)
        {
            Integer jk  = std::min(rank, j + 1);

            memcpy(ptr_U, ptr_A, jk * sizeof(V));
            memset(ptr_U+j+1, 0, (K - jk) * sizeof(V));

            ptr_U       += U.ld();
            ptr_A       += A2.ld();
        };

        ptr_A           = A2.ptr();

        for (Integer j = 0; j < rank; ++j)
        {
            memset(ptr_A, 0, j * sizeof(V));
            ptr_A[j]    = 1.;

            ptr_A       += A2.ld();
        };
        for (Integer j = rank; j < K; ++j)
        {
            memset(ptr_A, 0, M * sizeof(V));
            ptr_A[j]    = 1.;

            ptr_A       += A2.ld();
        };

        Matrix L        = Matrix(A2,true);
        L               = L(colon(), colon(1,K));

        L.get_struct().set(predefined_struct_type::tril);
        U.get_struct().set(predefined_struct_type::triu);

        ret = lu_return_type(L, Matrix(U,true), details::pv_constructor::make(P), 
                            details::pv_constructor::make(Q));
        return;
    }
    else
    {
        ti::ti_empty ti;

        Mat_D L(ti, M, K);

        V *ptr_L = L.ptr();

        for (Integer j = 0; j < rank; ++j)
        {
            memset(ptr_L, 0, j * sizeof(V));
            ptr_L[j] = 1.;
            memcpy(ptr_L+j+1, ptr_A+j+1, (M-j-1)*sizeof(V));

            ptr_L   += L.ld();
            ptr_A   += A2.ld();
        };
        for (Integer j = rank; j < K; ++j)
        {
            memset(ptr_L, 0, M * sizeof(V));
            ptr_L[j] = 1.;
            ptr_L   += L.ld();
        };

        ptr_A       = A2.ptr();

        for (Integer j = 1; j < K; ++j)
        {
            memset(ptr_A+j, 0, (K - j) * sizeof(V));

            ptr_A   += A2.ld();
        }

        Matrix U        = Matrix(A2,true);
        U               = U(colon(1,K), colon());

        U.get_struct().set(predefined_struct_type::triu);
        L.get_struct().set(predefined_struct_type::tril);        

        ret = lu_return_type(Matrix(L,true), U, details::pv_constructor::make(P), 
                        details::pv_constructor::make(Q));
        return;
    };
}; 

template struct lu_dense<matcl::Real>;
template struct lu_dense<matcl::Float>;
template struct lu_dense<matcl::Complex>;
template struct lu_dense<matcl::Float_complex>;

};};
