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
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/decompositions/lu/lu_struct.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-internals/base/optim_params.h"
#include "matcl-internals/base/pv_constructor.h"
#include "matcl-internals/container/mat_d.h"

namespace matcl { namespace details
{

template<class V> 
void lu_band<V>::eval(lu_return_type& ret, const raw::Matrix<V,struct_banded>& A, const options& opts)
{
    using VR            = typename md::real_type<V>::type;
    using Mat_D         = raw::Matrix<V,struct_dense>;
    using Mat_B         = raw::Matrix<V,struct_banded>;
    using Mat_S         = raw::Matrix<V,struct_sparse>;
    using pivot_type    = opt::linsolve::pivot_type;

    Integer M           = A.rows();
    Integer N           = A.cols();

    Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
    pivot_type piv      = (pivot_type)piv_int;
    VR tol              = (VR)opts.get_option<Real>(opt::linsolve::tol());
    bool isp            = (piv == pivot_type::partial);
    bool allow_nonunit  = opts.get_option<bool>(opt::linsolve::allow_nonunit_L());
    value_code vc       = matrix_traits::value_code<V>::value;

    //-------------------------------------------------------------------
    //                  fast exit
    //-------------------------------------------------------------------
    if ( M == 0 || N == 0)
    {
        ti::ti_empty ti;
        Mat_D L(ti, M, 0);
        Mat_D U(ti, 0, N);

        permvec P = permvec::identity(M);
        permvec Q = permvec::identity(N);

        ret = lu_return_type(Matrix(L,false), Matrix(U,false), P, Q);
        return;
    }

    //struct switch
    struct_flag flag        = A.get_struct();
 
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

    if (A_ldiags == 0 && A_udiags == 0)
        return lu_diag<V,struct_banded>::eval(ret, A, piv, tol);

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

    if (!isp)
    {
        Real density = Real(A.nnz())/(M+1.)/(N+1.);

        if (density < optim_params::max_sparse_density_min)
        {            
            Mat_S smat = raw::converter<Mat_S,Mat_B>::eval(A);

            return lu_sparse<V>::eval(ret, smat, opts);
        }
        else
        {
            Mat_D smat = raw::converter<Mat_D, Mat_B>::eval(A);
            return lu_dense<V>::eval(ret, smat,opts);
        };
    };

    Integer ld  = A.number_subdiagonals();
    Integer ud  = A.number_superdiagonals();
    Integer N2  = max(N,ld + ud + 1);

    Mat_B AFM(ti::ti_empty(), M, N2, -ld, ud+ld);

    using VL    = typename matcl::details::lapack_value_type<V>::type;

    const VL* ptr_A = lap(A.rep_ptr());
    VL* ptr_AF      = lap(AFM.rep_ptr());

    for (Integer J = 1; J <= N; ++J)
    {
        Integer J1  = max( J-ud, 1 );
        Integer J2  = min( J+ld, M );

        matcl::lapack::copy(J2-J1+1, ptr_A+ud-J+J1, 1, ptr_AF+ld+ud-J+J1, 1);

        ptr_A   += A.ld();
        ptr_AF  += AFM.ld();
    };

    ptr_AF      = lap(AFM.rep_ptr());

    raw::integer_dense raw_piv(ti::ti_empty(), 1, M);

    lapack::i_type info;
    lapack::gbtrf(M, N, ld, ud, ptr_AF, AFM.ld(), lap(raw_piv.ptr()), &info);

    Integer K           = std::min(M, N);
    Integer* ptr_piv    = raw_piv.ptr(); 

    Matrix P            = matcl::irange(1, M);
    Integer* ptr_P      = P.get_array_unique<Integer>(); 

    for (Integer i = 0; i < K; i++)
    {        
        Integer pos     = ptr_piv[i]-1;

        if (pos != i)
        {
            Integer sw  = ptr_P[i];
            ptr_P[i]    = ptr_P[pos];
            ptr_P[pos]  = sw;
        };
    }

    Mat_S L(ti::ti_empty(), M, K, K*(ld+1));
    {
        raw::details::sparse_ccs<V>& L_rep = L.rep();

        Integer* ptr_c  = L_rep.ptr_c();
        Integer* ptr_r  = L_rep.ptr_r();
        V* ptr_x        = L_rep.ptr_x();
        V* ptr_A2       = AFM.rep_ptr();
        Integer nz      = 0;

        for (Integer i = 0; i < K; ++i)
        {
            ptr_c[i]    = nz;

            ptr_r[nz]   = i;
            ptr_x[nz]   = 1.;
            ++nz;

            Integer fr  = AFM.first_row(i);
            Integer lr  = AFM.last_row(i);
            Integer ii  = AFM.first_elem_pos(i) + (i+1-fr);
            fr          = i+1;

            for(Integer k = fr; k<= lr; ++k,++ii)
            {
                ptr_r[nz] = k;
                ptr_x[nz] = ptr_A2[ii];
                ++nz;
            };

            ptr_A2      += AFM.ld();
        };

        ptr_c[K]        = nz;

        //apply delayed permutations
        matcl::Matrix P_inv = irange(0, M - 1);

        Integer* ptr_Pi = P_inv.get_array_unique<Integer>();

        for (Integer i = K - 1; i >= 0; --i)
        {
            for (Integer k = ptr_c[i]+1; k < ptr_c[i+1]; ++k)
            {
                ptr_r[k] = ptr_Pi[ptr_r[k]];
            };

            Integer pos = ptr_piv[i]-1;
            if (pos != i)
            {
                Integer sw      = ptr_Pi[i];
                ptr_Pi[i]       = ptr_Pi[pos];
                ptr_Pi[pos]     = sw;
            };
        };

        L_rep.sort();
    };

    Mat_B U(ti::ti_empty(), K, N, 0, ud+ld);    
    {
        V* ptr_U        = U.rep_ptr();
        V* ptr_A2       = AFM.rep_ptr();

        for (Integer i = 0; i < N; ++i)
        {
            Integer fr  = AFM.first_row(i);
            Integer lr  = min(AFM.last_row(i),i);
            Integer ii  = AFM.first_elem_pos(i);
            Integer jj  = U.first_elem_pos(i);

            for(Integer k = fr; k<= lr; ++k,++ii, ++jj)
                ptr_U[jj] = ptr_A2[ii];

            ptr_A2  += AFM.ld();
            ptr_U   += U.ld();
        };
    };

    L.get_struct().set(predefined_struct_type::tril);
    U.get_struct().set(predefined_struct_type::triu);

    permvec Q = permvec::identity(N);
    ret = lu_return_type(Matrix(L,true),Matrix(U,true), details::pv_constructor::make(P) , Q);
    return;
}; 

};};

template struct matcl::details::lu_band<matcl::Real>;
template struct matcl::details::lu_band<matcl::Float>;
template struct matcl::details::lu_band<matcl::Complex>;
template struct matcl::details::lu_band<matcl::Float_complex>;
