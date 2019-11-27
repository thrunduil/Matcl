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

#include "matcl-linalg/decompositions/lu/lu_impl.h"
#include "matcl-linalg/decompositions/lu/lu_struct.h"
#include "matcl-linalg/lusol/lusol.h"
#include "matcl-linalg/utils/utils.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-internals/base/pv_constructor.h"
#include "matcl-linalg/norms_error/norm.h"
#include "lu_superlu.h"
#include "matcl-linalg/special_matrices/struct_flag_linalg.h"

#undef small

namespace matcl { namespace details
{

using matcl::operator-;

template<class V> 
void lu_sparse<V>::eval(lu_return_type& ret, const raw::Matrix<V,struct_sparse>& A, const options& opts)
{ 
    using VR            = typename md::real_type<V>::type;
    using Mat_S         = raw::Matrix<V,struct_sparse>;
    using Mat_D         = raw::Matrix<V,struct_dense>;
    using pivot_type    = opt::linsolve::pivot_type;
    using lu_solver_type= opt::linsolve::lu_solver_type;

    Integer M           = A.rows();
    Integer N           = A.cols();
    Integer K           = std::min(M, N);

    Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
    pivot_type piv      = (pivot_type)piv_int;
    VR tol              = (VR)opts.get_option<Real>(opt::linsolve::tol());
    bool isp            = (piv == pivot_type::partial);
    bool allow_nonunit  = opts.get_option<bool>(opt::linsolve::allow_nonunit_L());
    value_code vc       = matrix_traits::value_code<V>::value;

    //-------------------------------------------------------------------
    //                  fast exit
    //-------------------------------------------------------------------
    if ( M == 0 || N == 0 || A.nnz() == 0)
    {
        value_code vcr  = matrix_traits::value_code<VR>::value;

        Matrix L    = speye(M, K, vcr);
        Matrix U    = spzeros(K, N, 0, vcr);
        permvec P   = permvec::identity(M);
        permvec Q   = permvec::identity(N);

        ret = lu_return_type(L, U, P, Q);
        return;
    }

    struct_flag flag        = A.get_struct();

    if (flag.is_id())
    {
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
        return lu_diag<V,struct_sparse>::eval(ret, A, piv, tol);

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

    if (isp == false)
        return eval_lusol(ret, A, opts);

    lu_solver_type solver   = (lu_solver_type)opts.get_option<Integer>(opt::linsolve::lu_solver());

    switch(solver)
    {
        case lu_solver_type::lusol:
            return eval_lusol(ret, A, opts);
        case lu_solver_type::superlu:
        default:
        {
            //superlu is broken for M < N
            if (M >= N)
                return eval_superlu(ret, A, opts);
            else
                return eval_lusol(ret, A, opts);
        }
    };
}

template<class V> 
void lu_sparse<V>::eval_superlu(lu_return_type& ret, const raw::Matrix<V,struct_sparse>& A, const options& opts)
{
    return superlu_wrap<V>::eval_lu(A, ret, opts);
};

template<class V> 
void lu_sparse<V>::eval_lusol(lu_return_type& ret, const raw::Matrix<V,struct_sparse>& A, const options& opts)
{
    using VR            = typename md::real_type<V>::type;
    using Mat_S         = raw::Matrix<V,struct_sparse>;
    using Mat_D         = raw::Matrix<V,struct_dense>;
    using pivot_type    = opt::linsolve::pivot_type;

    VR tol_r            = (VR)opts.get_option<Real>(opt::linsolve::tol_r());
    VR tol              = (VR)opts.get_option<Real>(opt::linsolve::tol());
    Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
    pivot_type piv      = (pivot_type)piv_int;

    Integer M           = A.rows();
    Integer N           = A.cols();

    //-------------------------------------------------------------------
    //                  factor and extract
    //-------------------------------------------------------------------
    VR TOL          = tol;
    VR NORM_F       = 1.0;
    Integer K       = std::min(M, N);

    if (TOL < VR(0.0))
    {
        NORM_F      = (VR)norm(Matrix(A,false), -1) / sqrt(VR(K));
        tol         = VR(10.0) * constants::eps<VR>() * NORM_F;
    };      

    if ( NORM_F == VR(0.0) )
    {
        value_code vcr  = matrix_traits::value_code<VR>::value;

        Matrix L    = speye(M, K, vcr);
        Matrix U    = spzeros(K, N, 0, vcr);
        permvec P   = permvec::identity(M);
        permvec Q   = permvec::identity(N);
        ret         = lu_return_type(L, U, P, Q);
        return;
    }

    VR small    = TOL;
    VR Amax     = VR(1.)/tol_r; 
    int pivot   = -1;

    switch(piv)
    {
        case pivot_type::partial:   pivot = lusol::LUSOL_PIVOT_TPP; break;
        case pivot_type::rook:      pivot = lusol::LUSOL_PIVOT_TRP; break;
        case pivot_type::complete:  pivot = lusol::LUSOL_PIVOT_TCP; break;
        default:
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
    };         

    lusol::LUSOLrec<VL> *LUSOL = lusol::LUSOLrec<VL>::LUSOL_create(0, 0, pivot, 0);

    if (LUSOL == nullptr)
        throw error::alloc();

    try
    {
        LUSOL->luparm[lusol::LUSOL_IP_PRINTLEVEL]       = lusol::LUSOL_MSG_NONE;
        LUSOL->parmlu[lusol::LUSOL_RP_ZEROTOLERANCE]    = small;
        
        if (Amax >= 1) 
            LUSOL->parmlu[lusol::LUSOL_RP_FACTORMAX_Lij] = Amax;

        {
            const raw::details::sparse_ccs<V>& d = A.rep();
            const Integer* Ad_r = d.ptr_r();
            const Integer* Ad_c = d.ptr_c();
            const VL* Ad_x = reinterpret_cast<const VL*>(d.ptr_x());

            bool status = LUSOL->LUSOL_assign(Ad_r, Ad_c, Ad_x, M, N, N);

            if (status == false) 
                throw error::alloc();
        };

        // Factor  A = L U.
        lusol::INT inform;
        LUSOL->LU1FAC(&inform);

        if (inform > lusol::LUSOL_INFORM_SERIOUS) 
        {
            const char* msg = LUSOL->LUSOL_informstr(inform);
            throw error::error_lu(msg);
        }

        // Extract vectors P and Q        
        raw::Matrix<Integer,struct_dense> P(ti::ti_empty(), 0, 1, M);   
        raw::Matrix<Integer,struct_dense> Q(ti::ti_empty(), 0, 1, N);    

        {
            Integer* ptr = P.ptr();

            Integer j;
            for (Integer k = 0; k < M; k++) 
            {
                j = LUSOL->ip[k+1]-1;
                ptr[k] = j+1;
            }
        };

        {
            Integer *ptr = Q.ptr();
            Integer j;

            for (Integer k = 0; k < N; k++) 
            {
                j = LUSOL->iq[k+1]-1;
                ptr[k] = j+1;
            };
        };

        permvec p = details::pv_constructor::make(Matrix(P,false));
        permvec q = details::pv_constructor::make(Matrix(Q,false));

        Matrix out_L, out_U;

        extract_factors(LUSOL, out_L, out_U, p, q);

        LUSOL->LUSOL_free();

        ret = lu_return_type(out_L,out_U,p,q);
        return;
    }
    catch(std::bad_alloc&)
    {
        LUSOL->LUSOL_free();
        throw matcl::error::alloc(0);
    }
    catch(...)
    {
        LUSOL->LUSOL_free();
        throw;
    };    
}; 

template<class V> 
void lu_sparse<V>::extract_factors(lusol::LUSOLrec<VL>* LUSOL, Matrix& out_L, Matrix& out_U,
                                   const permvec& p, const permvec& q)
{
    permvec pi = invperm(p);
    permvec qi = invperm(q);

    const Integer* ptr_p    = p.is_id() ? nullptr : p.to_array();
    const Integer* ptr_pi   = pi.is_id() ? nullptr : pi.to_array();
    const Integer* ptr_qi   = qi.is_id() ? nullptr : qi.to_array();

    Integer M       = p.length();
    Integer N       = q.length();
    Integer K       = std::min(M, N);

    // Extract L    
    Integer Lnz     = LUSOL->luparm[lusol::LUSOL_IP_NONZEROS_L];        
    Integer nrank   = LUSOL->luparm[lusol::LUSOL_IP_RANK_U];

    Integer Lm      = M;
    Integer Ln      = K;
    Integer Um      = K;
    Integer Un      = N;        

    raw::Matrix<V,struct_sparse> L(ti::ti_empty(),Lm,Ln,Lnz+Ln);
    {
        //speye(Am,Am)-out_L(colon(P),colon(P));

        raw::details::sparse_ccs<V>& Ld = L.rep();
        Integer* Li = Ld.ptr_r();
        Integer* Lp = Ld.ptr_c();
        V* Lx       = Ld.ptr_x();

        Integer numL0 = LUSOL->luparm[lusol::LUSOL_IP_COLCOUNT_L0]; 

        std::vector<Integer> Lstart(M + 1, 0);
        std::vector<Integer> Llen(M + 1, 0);

        Integer i, pos, j, k, len;
        for (i = 0, pos = LUSOL->lena; i < numL0; i++, pos -= len) 
        {
            j           = LUSOL->indr[pos]-1;
            len         = LUSOL->lenc[i+1];
            Lstart[j]   = pos;
            Llen[j]     = len;
        }

        for (j = 0, k = 0; j < nrank; j++) 
        {
            Lp[j]       = k;
            Li[k]       = j;
            Lx[k]       = 1.;

            Integer cc = ptr_p? ptr_p[j] - 1 : j; 

            ++k;

            if (ptr_pi)
            {
                for (Integer i2 = 0, pos2 = Lstart[cc]; i2 < Llen[cc]; ++i2, ++k, --pos2) 
                {
                    Integer rr  = ptr_pi[LUSOL->indc[pos2]-1]-1;
                    Li[k]       = rr;
                    Lx[k]       = -reinterpret_cast<V&>(LUSOL->a[pos2]);
                }
            }
            else
            {
                for (Integer i2 = 0, pos2 = Lstart[cc]; i2 < Llen[cc]; ++i2, ++k, --pos2) 
                {
                    Integer rr  = LUSOL->indc[pos2]-1;
                    Li[k]       = rr;
                    Lx[k]       = -reinterpret_cast<V&>(LUSOL->a[pos2]);
                }
            }
        };
        for (j = nrank; j < Ln; ++j, ++k) 
        {
            Lp[j]       = k;
            Li[k]       = j;
            Lx[k]       = 1.;
        };

        Lp[Ln]          = k;

        Ld.sort();
    };    

    using Mat_S         = raw::Matrix<V,struct_sparse>;

    // Extract U(p,:)^T (i.e. transpose of U)
    Integer Unz = LUSOL->luparm[lusol::LUSOL_IP_NONZEROS_U];
    Mat_S Upt(ti::ti_empty(), Un, Um, Unz);

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
        for (j = nrank; j < K; j++) 
            Up[j]   = k;

        Up[K]       = k;

        if (ptr_qi)
        {
            for (k = 0, pos = 1; k < Unz; k++, pos++) 
            {
                Integer rr  = ptr_qi[LUSOL->indr[pos]-1]-1;
                Ui[k]       = rr;
                Ux[k]       = reinterpret_cast<V&>(LUSOL->a[pos]);
            }
        }
        else
        {
            for (k = 0, pos = 1; k < Unz; k++, pos++) 
            {
                Integer rr  = LUSOL->indr[pos]-1;
                Ui[k]       = rr;
                Ux[k]       = reinterpret_cast<V&>(LUSOL->a[pos]);
            }
        }

        Ud.sort();
    };

    out_L = Matrix(L,true);
    out_U = Matrix(Upt,true);
    out_U = trans(out_U);

    out_L.add_struct(struct_flag(predefined_struct_type::tril));
    out_U.add_struct(struct_flag(predefined_struct_type::triu));
};

};};

template struct matcl::details::lu_sparse<matcl::Real>;
template struct matcl::details::lu_sparse<matcl::Float>;
template struct matcl::details::lu_sparse<matcl::Complex>;
template struct matcl::details::lu_sparse<matcl::Float_complex>;
