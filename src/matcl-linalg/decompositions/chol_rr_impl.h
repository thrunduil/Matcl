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

#pragma once

#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-linalg/lusol/lusol.h"
#include "matcl-linalg/utils/utils.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-internals/base/optim_params.h"
#include "matcl-linalg/norms_error/norm.h"

namespace matcl { namespace details
{

template<class V, class S>
struct chol_rr_diag
{
    using Mat           = raw::Matrix<V,S>;
    using DenseMatrix   = raw::Matrix<V,struct_dense>;
    using VR            = typename md::real_type<V>::type;

    static void eval(chol_rr_return_type& ret, const Mat& A, VR tol)
    {        
        Integer K = A.rows();

        if (K == 0)
        {
            permvec p = details::pv_constructor::make(matcl::izeros(1,0));
            ret = chol_rr_return_type(Matrix(A,false), p, 0);
            return;
        };

        Matrix p;
        Matrix D        = Matrix(A.get_diag(),false);        
        tie(D,p)        = matcl::sort2(D,1,false);
        
        Integer rank    = 0;
        V* ptr          = D.get_array_unique<V>();        

        if (tol < 0)
            tol         = VR(10.0) * matcl::constants::eps<VR>() * matcl::abs(ptr[0]);

        for(Integer i = 0; i < K; ++i)
        {
            VR val      = matcl::real(ptr[i]);
            ptr[i]      = raw::details::sqrt_helper<V>::eval(std::max<VR>(val, 0.));

            if (val > tol)
                ++rank;
        }

        Matrix S_ = matcl::bdiag(D, 0);

        ret = chol_rr_return_type(S_,details::pv_constructor::make(p),rank);
    }
};

template<class V, class S> 
struct chol_rr_impl{};

//--------------------------------------------------------------------------------
//                                  DENSE
//--------------------------------------------------------------------------------
template<class V>  
struct chol_rr_impl<V,struct_dense> 
{
    using Mat       = raw::Matrix<V,struct_dense>;
    using real_type = typename md::real_type<V>::type;

    static void eval(chol_rr_return_type& ret, const Mat& A, bool upper, real_type tol)
    {
        if (A.get_struct().is_id())
        {
            permvec p   = permvec::identity(A.rows());
            Matrix S    = Matrix(A,false);

            if (tol < real_type(0.0))
                tol     = real_type(10.0) * constants::eps<real_type>();

            Integer rank= (1. <= tol) ? 0 : A.rows();
            ret = chol_rr_return_type(S,p,rank);
            return;
        }

        if(is_diag(Matrix(A,false)))
            return chol_rr_diag<V,struct_dense>::eval(ret, A, tol);

        Integer N = A.rows();
        if (N == 0)
        {
            permvec p   = details::pv_constructor::make(matcl::izeros(1,0));
            Matrix S    = Matrix(A,false);

            ret = chol_rr_return_type(S,p,0);
            return;
        };

        using workspace     = pod_workspace<real_type>;

        Mat Ac              = A.make_unique();
        Ac.get_struct().reset();

        Integer* ptr_PIV;
        matcl::Matrix PIV   = matcl::make_integer_dense_noinit(N,1,ptr_PIV);
        workspace WORK      = matcl::pod_workspace<real_type>(2*N,1);
        real_type* ptr_WORK = WORK.ptr();
        Integer rank        = 0;

        lapack::i_type info;

        char uplo           = upper ? 'U' : 'L';

        lapack::potfp3(&uplo, N, lap(Ac.ptr()), Ac.ld(), lap(ptr_PIV), *lap(&rank), 
                       tol, lap(ptr_WORK), info);

        Matrix S;
        if (upper)
        {
            Matrix U;
            matcl::raw::inplace::make_triu<Mat>::eval(U,Ac,0);
            S = U(colon(1,rank), colon());
        }
        else
        {
            Matrix U;
            matcl::raw::inplace::make_tril<Mat>::eval(U,Ac,0);
            S = U(colon(), colon(1,rank));
        }

        ret = chol_rr_return_type(S,details::pv_constructor::make(PIV),rank);
        return;
    };
};

//--------------------------------------------------------------------------------
//                                  SPARSE
//--------------------------------------------------------------------------------
template<class V> 
struct chol_rr_impl<V,struct_sparse> 
{
    using Mat       = raw::Matrix<V,struct_sparse>;
    using real_type = typename md::real_type<V>::type;

    static void eval(chol_rr_return_type& ret, const Mat& A, bool upper, real_type tol)
    {
        Integer N           = A.rows();
        struct_flag flag    = A.get_struct();

        if (flag.is_id())
        {
            permvec p       = permvec::identity(A.rows());

            if (tol < real_type(0.0))
                tol         = real_type(10.0) * constants::eps<real_type>();

            Matrix S        = Matrix(A,false);            
            Integer rank    = (1. <= tol) ? 0 : A.rows();
            
            ret = chol_rr_return_type(S,p,rank);
            return;
        }

        if(is_diag(Matrix(A,false)))
            return chol_rr_diag<V,struct_sparse>::eval(ret, A, tol);

        if (A.nnz() == 0)
        {
            Matrix S    = Matrix(A,false);
            permvec p   = details::pv_constructor::make(matcl::izeros(1,0));
            ret         = chol_rr_return_type(S,p,0);
            return;
        };

        //-------------------------------------------------------------------
        //                  factor and extract
        //-------------------------------------------------------------------
        int pivot  = lusol::LUSOL_PIVOT_TSP;

        using VT = typename lusol_value_type<V>::type;

        lusol::LUSOLrec<VT> *LUSOL = lusol::LUSOLrec<VT>::LUSOL_create(0, 0, pivot, 0);

        if (LUSOL == nullptr)
            throw error::alloc();

        real_type small = get_tolerance(tol, A);
        real_type Amax  = 0.; 

        try
        {
            LUSOL->luparm[lusol::LUSOL_IP_PRINTLEVEL]       = lusol::LUSOL_MSG_NONE;
            LUSOL->parmlu[lusol::LUSOL_RP_ZEROTOLERANCE]    = small;

            if (Amax >= 1) 
            {
                LUSOL->parmlu[lusol::LUSOL_RP_FACTORMAX_Lij] = Amax;
            }

            {
                const raw::details::sparse_ccs<V>& d = A.rep();
                const Integer* Ad_r = d.ptr_r();
                const Integer* Ad_c = d.ptr_c();
                const VT* Ad_x = reinterpret_cast<const VT*>(d.ptr_x());

                bool status = LUSOL->LUSOL_assign(Ad_r, Ad_c, Ad_x, N, N, N);

                if (status == false) 
                    throw error::alloc();
            };

            // Factor  A = L U.
            lapack::i_type inform;
            LUSOL->LU1FAC(&inform);

            if (inform > lusol::LUSOL_INFORM_SERIOUS) 
            {
                const char* msg = LUSOL->LUSOL_informstr(inform);
                throw error::error_lu(msg);
            }

            // Extract vector P
            Integer nrank = LUSOL->luparm[lusol::LUSOL_IP_RANK_U];         

            raw::Matrix<Integer,struct_dense> Q(ti::ti_empty(),0,1,N);
            {
                Integer* ptr = Q.ptr();

                Integer j;
                for (Integer k = 0; k < N; k++) 
                {
                    j       = LUSOL->iq[k+1]-1;
                    ptr[k]  = j+1;
                }
            };

            permvec q    = details::pv_constructor::make(Matrix(Q,false));
            permvec qi   = invperm(q);

            const Integer* ptr_qi   = qi.to_array();

            // Extract U(p,:)^T (i.e. transpose of U)
            Integer Unz = LUSOL->luparm[lusol::LUSOL_IP_NONZEROS_U];
            raw::Matrix<V,struct_sparse> Upt(ti::ti_empty(),N,N,Unz);

            if (Unz > 0)
            {
                raw::details::sparse_ccs<V>& Ud = Upt.rep();
                Integer* Ui = Ud.ptr_r();
                Integer* Up = Ud.ptr_c();
                V* Ux       = Ud.ptr_x();

                Integer k, pos, j;
                for (j = 0, k = 0; j < nrank; j++) 
                {
                    Up[j]   = k;
                    k       += LUSOL->lenr[LUSOL->ip[j+1]];
                }
                for (j = nrank; j <= N; j++) 
                    Up[j]   = k;

                for (k = 0, pos = 1; k < Up[nrank]; k++, pos++) 
                {
                    Integer rr = ptr_qi[LUSOL->indr[pos]-1]-1;
                    Ui[k]   = rr;
                    Ux[k]   = reinterpret_cast<V&>(LUSOL->a[pos]);
                }

                Ud.sort();

                //convert to cholesky factor
                for (j = 0; j < nrank; j++) 
                {
                    Integer diag_pos    = Up[j];
                    V diag              = Ux[diag_pos];
                    real_type diag_r    = matcl::real(diag);
                    real_type scal      = real_type(1.0) /matcl::sqrt(std::max(diag_r, real_type()));

                    if (upper == true)
                    {
                        for (Integer k2 = diag_pos; k2 < Up[j+1]; ++k2)
                        {
                            Ux[k2]      = Ux[k2] * scal;
                        };
                    }
                    else
                    {
                        for (Integer k2 = diag_pos; k2 < Up[j+1]; ++k2)
                        {
                            Ux[k2]      = conj(Ux[k2] * scal);
                        };
                    }
                }
            };

            permvec p = q;

            matcl::Matrix out_U = Matrix(Upt,false);
            if (upper == true)
            {
                out_U               = trans(out_U);
                out_U.add_struct(struct_flag(predefined_struct_type::triu));
            }
            else
            {
                out_U.add_struct(struct_flag(predefined_struct_type::tril));
            };

            using SparseMatrix = raw::Matrix<V,struct_sparse>;

            LUSOL->LUSOL_free();            

            ret = chol_rr_return_type(out_U, p, nrank);
            return;
        }
        catch(...)
        {
            LUSOL->LUSOL_free();
            throw;
        };    
    };

    static real_type get_tolerance(real_type tol, const Mat& A)
    {
        if (tol >= 0.0)
            return tol;

        tol                 = real_type(10.0) * constants::eps<real_type>();

        Integer MN          = std::min(A.rows(), A.cols());
        MN                  = std::max(1, MN);
        real_type norm_F    = real_type(norm(Matrix(A,false), -1));
        norm_F              = norm_F / std::sqrt(real_type(MN));
        tol                 = norm_F * tol;
        
        return tol;
    };
};

//--------------------------------------------------------------------------------
//                                  BAND
//--------------------------------------------------------------------------------
template<class V> 
struct chol_rr_impl<V,struct_banded> 
{
    using Mat       = raw::Matrix<V,struct_banded>;
    using real_type = typename md::real_type<V>::type;

    static void eval(chol_rr_return_type& ret, const Mat& A, bool upper, real_type tol)
    {
        struct_flag flag = A.get_struct();

        if (flag.is_id())
        {
            permvec p       = permvec::identity(A.rows());

            Matrix S        = Matrix(A,false);            

            if (tol < real_type(0.0))
                tol         = real_type(10.0) * constants::eps<real_type>();

            Integer rank    = (1. <= tol) ? 0 : A.rows();
            ret = chol_rr_return_type(S,p,rank);
            return;
        }

        if(is_diag(Matrix(A,false)))
            return chol_rr_diag<V,struct_banded>::eval(ret, A, tol);

        if (A.rows() == 0)
        {
            permvec p   = details::pv_constructor::make(matcl::izeros(1,0));
            Matrix S    = Matrix(A,false);

            ret = chol_rr_return_type(S,p,0);
            return;
        };

         Real density = Real(A.nnz())/(A.rows()+1.)/(A.cols()+1.);

         if (density < optim_params::max_sparse_density_min)
         {
             using SparseMatrix = raw::Matrix<V,struct_sparse>;
             SparseMatrix smat  = raw::converter<SparseMatrix,raw::Matrix<V,struct_banded>>::eval(A);
             return chol_rr_impl<V, struct_sparse>::eval(ret, smat,upper,tol);
         }
         else
         {
             using DenseMatrix  = raw::Matrix<V,struct_dense>;
             using VR           = md::real_type<V>::type;

             DenseMatrix smat = raw::converter<DenseMatrix,raw::Matrix<V,struct_banded>>::eval(A);
             return chol_rr_impl<V, struct_dense>::eval(ret, smat,upper,VR(tol));
         };
    };
};

//--------------------------------------------------------------------------------
//                                  VALUE TYPE
//--------------------------------------------------------------------------------
template<class V, class S>
struct chol_rr
{
    using M         = raw::Matrix<V,S>;
    using real_type = typename md::real_type<V>::type;

    static void eval(chol_rr_return_type& ret, const M& A, bool upper, real_type tol)
    {
        return chol_rr_impl<V,S>::eval(ret, A, upper, tol);
    };
};

template<class S>
struct chol_rr<Integer,S>
{
    using M = raw::Matrix<Integer,S>;

    static void eval(chol_rr_return_type& ret, const M& A, bool upper, Real tol)
    {
        using MC    = raw::Matrix<Real,S>;
        MC AC       = raw::converter<MC,M>::eval(A);
        return chol_rr<Real,S>::eval(ret, AC, upper, tol);
    };
};

};};
