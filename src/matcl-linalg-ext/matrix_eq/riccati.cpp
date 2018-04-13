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

#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-linalg/matrix_eq/riccati.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/func/test_inf_nan.h"

#include "matcl-linalg/decompositions/cholmod.h"
#include "matcl-linalg/matrix_eq/lyapunov.h"
#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/norms_error/norm.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

//TODO
#if 0
//#include "extern_lib_src/slicot/include/slicot.h"

namespace matcl { namespace details
{

enum class riccati_kind
{
    DARE, CARE
};

template<class V>
mat_tup_3 solve_riccati_val(const V& A, const V& B, const V& Q, const V& R, 
                            const V& L, const riccati_kind kind)
{
    using constants::nan;

    typedef typename V::value_type  vt;
    typedef typename V::struct_type st;

    bool isv =    A.all_finite() && B.all_finite() && Q.all_finite() && R.all_finite() && L.all_finite();

    if (isv == false)
    {
        Matrix cond = constants::nan();
        Matrix eig  = repmat(constants::nan(), A.rows(), 1);
        return mat_tup_3(repmat(constants::nan(), A.rows(), A.cols()), cond, eig);
    }

    Integer M = R.rows();       
    Integer N = A.cols();
    
    V   x(ti::ti_real(),N,N);
    V   alfar(ti::ti_real(),2*N,1);
    V   alfai(ti::ti_real(),2*N,1);
    V   beta(ti::ti_real(),2*N,1);
    V   s(ti::ti_real(),max(1,2*N + M),2*N + M);
    V   t(ti::ti_real(),max(1,2*N + M),2*N);
    V   u(ti::ti_real(),max(1,2*N),2*N);

    x.set_struct(struct_flag());

    Integer info = 0;

    char dico   = kind == riccati_kind::DARE ? 'D' : 'C' ; // discrete or continuous time
    char jobb   = 'B'; // Matrix B provided
    char fact   = 'N'; // R, Q, L provided non-factorized
    char uplo   = 'U'; 
    char jobl   = 'N';
    char sort   = 'S'; // stable spectrum solution, i.e. we need a non negative definite solution

    lapack::i_type lda  = A.ld();
    lapack::i_type ldb  = B.ld();
    lapack::i_type ldq  = Q.ld();
    lapack::i_type ldr  = R.ld();
    lapack::i_type ldl  = L.ld();
    lapack::i_type ldx  = x.ld();
    lapack::i_type lds  = s.ld();
    lapack::i_type ldt  = t.ld();
    lapack::i_type ldu  = u.ld();

    Real rcond      = -1;
    Real tol        = 1e-50;

    Integer ldwork  = max(max(max(7*(2*N+1)+16,16*N),2*N+M),3*M);
    
    raw::integer_dense  iwork(ti::ti_int(), max(max(1,M),2*N),1);
    V                   dwork(A.get_type(),ldwork,1);
    raw::integer_dense  bwork(ti::ti_int(), 2*N, 1);

    Integer             trash = 0;

    slicot::sb02od_(&dico, &jobb, &fact, &uplo, &jobl, &sort, 
                    lap(&N), lap(&M), lap(&trash), 
                    lap(A.ptr()), &lda, lap(B.ptr()), &ldb, lap(Q.ptr()), &ldq, 
                    lap(R.ptr()), &ldr, lap(L.ptr()), &ldl, lap(&rcond), 
                    lap(x.ptr()), &ldx, lap(alfar.ptr()), lap(alfai.ptr()), lap(beta.ptr()), 
                    lap(s.ptr()), &lds, lap(t.ptr()), &ldt, lap(u.ptr()), &ldu, lap(&tol), 
                    lap(iwork.ptr()), lap(dwork.ptr()), lap(&ldwork), lap(bwork.ptr()), lap(&info));

    if (info != 0 )
    {
        throw error::error_ricc(info);
    }

    return mat_tup_3(Matrix(x,true), rcond, 
                     div((Matrix(alfar,true) + constants::i() * Matrix(alfai,true)), Matrix(beta,true)));
}

mat_tup_3 solve_riccati(const Matrix& A, const Matrix& B, const Matrix& Q, const Matrix& R, 
                        const Matrix& L, const riccati_kind kind)
{
    if (!A.is_square() || !Q.is_square() || !R.is_square() || 
        A.rows() != B.rows() || B.rows() != Q.rows() || Q.rows() != L.rows() ||
        B.cols() != R.cols() || R.cols() != L.cols())
    {
        throw error::error_size_ricc();
    }

    matcl::value_code vA    = A.get_value_code();
    matcl::value_code vB    = B.get_value_code();
    matcl::value_code vQ    = Q.get_value_code();
    matcl::value_code vR    = R.get_value_code();
    matcl::value_code vL    = L.get_value_code();
    matcl::value_code v     = (matcl::value_code) max(max(Integer(vA),Integer(vB)),
                                                      max(max(Integer(vQ),Integer(vR)),Integer(vL)));

    switch (v)
    {
        case value_code::v_integer:
        case value_code::v_float:
        {
            //TODO: impl float
        }
        case value_code::v_real:
        {
            Matrix Ac   = convert(A,mat_code::real_dense);
            Matrix Bc   = convert(B,mat_code::real_dense);
            Matrix Qc   = convert(Q,mat_code::real_dense);
            Matrix Rc   = convert(R,mat_code::real_dense);
            Matrix Lc   = convert(L,mat_code::real_dense);

            typedef matcl::raw::Matrix<Real,struct_dense> DM;

            return details::solve_riccati_val<DM>(Ac.get_impl<DM>(),Bc.get_impl<DM>(),
                                                  Qc.get_impl<DM>(),Rc.get_impl<DM>(),
                                                  Lc.get_impl<DM>(), kind);
        }
        case value_code::v_float_complex:
        {
            //TODO: impl float
        }
        case value_code::v_complex:
        {
            throw error::error_complex_value_type_not_allowed();
        }
        case value_code::v_object:
        {
            throw error::object_value_type_not_allowed("riccati");
        }
    };

    throw error::error_general("impossible type case in solve_dare/care");
}

void check_value_types_ARE(const std::vector<Matrix>& to_check)
{
    for (std::vector<Matrix>::const_iterator it = to_check.begin(); it < to_check.end(); it++)
    {
        switch (it->get_value_code())
        {
            case value_code::v_integer:
            case value_code::v_float:
            {
                //TODO: impl float
            }
            case value_code::v_real:
                break; // OK
            case value_code::v_float_complex:
            {
                //TODO: impl float
            }
            case value_code::v_complex:
            {
                throw error::error_complex_value_type_not_allowed();
            }
            case value_code::v_object:
            {
                throw error::object_value_type_not_allowed("riccati");
            }
            default:
                throw error::error_general("impossible type case in solve_dare/care");
        }
    }
}
bool check_finite_ARE(const std::vector<Matrix>& to_check)
{
    for (std::vector<Matrix>::const_iterator it = to_check.begin(); it < to_check.end(); it++)
    {
        if (any(any(is_nan(*it), 1), 2) || any(any(is_inf(*it), 1), 2))
            return false;
    }

    return true;
}

Matrix sparse_solve_driccati(const Matrix& A, const Matrix& B, const Matrix& C, const Matrix& Q,
                             const Matrix& R, const Matrix& S, const Real lyapunov_tol, const Integer kmax, 
                             const Integer l_zero_cap)
{
    if (!A.is_square() || A.rows() != B.rows() || A.rows() != C.cols() || C.rows() != Q.rows()
        || !Q.is_square() || R.rows() != B.cols() || !R.is_square() || 
        S.rows() != B.cols() || S.cols() != Q.rows())
    {
        throw error::error_size_ricc();
    }

    if (A.numel() == 0) 
        return A;

    std::vector<Matrix> to_check;

    to_check.push_back(A);
    to_check.push_back(B);
    to_check.push_back(C);
    to_check.push_back(Q);
    to_check.push_back(R);
    to_check.push_back(S);

    check_value_types_ARE(to_check);

    Integer n = A.rows();

    if (!check_finite_ARE(to_check))
    {
        return  repmat(constants::nan(), n, n);
    }

    Matrix      Q_wave;
    permvec     Q_wave_p;
    Integer     Q_rank;
    Real        Q_e;

    std::tie(Q_wave, Q_wave_p, Q_rank, Q_e) = cholmod(Q,true);

    if (Q_e > 1e-8) 
        throw error::error_nonposdef(false);

    //TODO: change to linsolve_rev if possible
    Matrix to_factor = R - 
                       S * (linsolve(Q_wave, linsolve(trans(Q_wave),trans(S)(Q_wave_p, colon())))(Q_wave_p.invperm(), colon()));

    Matrix      W; 
    permvec     W_p;
    Integer     W_rank;
    Real        W_e;

    std::tie(W,W_p, W_rank, W_e) = cholmod(to_factor,true);

    if (W_e > 1e-8) 
        throw error::error_nonposdef(false);

    W           = trans(W(colon(), W_p.invperm()));
    Matrix L    = trans(Q_wave(colon(), Q_wave_p.invperm()));
    Matrix K    = spzeros(n, R.rows());

    // equivalent to X_0 = 0, does not always work
    // TODO: investigate, which is better
    //Matrix K = trans(C) * trans(linsolve(trans(R), S));    

    Matrix lastU, lastF;

    Matrix sp_B     = sparse(B);
    Matrix Z        = spzeros(n, 0);

    for (int k = 0 ; k < kmax; k++)
    {
        Matrix F    = A - sp_B * sparse(trans(K));
        Matrix U    = (mat_row(), trans(C) * L - K * trans(linsolve(L, trans(S))), - K * W);

        if (k > 0 && (     norm(lastF - F, -2) / norm(lastF, -2) < lyapunov_tol 
                        && norm(lastU - U, -2) / norm(lastU, -2) < lyapunov_tol)    ) 
            return Z;

        //TODO: try to do something with these transposes
        Z       = sparse_dale(F, U, l_zero_cap);
        K       = real(trans(linsolve(trans(R + trans(sp_B) * (Z * (ctrans(Z) * sp_B))),
                                       trans(trans(A) * (Z * (ctrans(Z) * sp_B)) + trans(C) * trans(S)))));
        lastF   = F;
        lastU   = U;
    }

    return Z;
}
Matrix sparse_solve_criccati(const Matrix& A, const Matrix& B, const Matrix& C, const Matrix& Q,
                             const Matrix& R, const Real lyapunov_tol, const Integer kmax, 
                             const Integer l_zero_cap, const Real V_update_tol)
{
    if (!A.is_square() || A.rows() != B.rows() || A.rows() != C.cols() || C.rows() != Q.rows()
        || !Q.is_square() || R.rows() != B.cols() || !R.is_square())
    {
        throw error::error_size_ricc();
    }

    if (A.numel() == 0) 
        return A;

    std::vector<Matrix> to_check;

    to_check.push_back(A);
    to_check.push_back(B);
    to_check.push_back(C);
    to_check.push_back(Q);
    to_check.push_back(R);

    check_value_types_ARE(to_check);
    
    Integer n = A.rows();

    if (! check_finite_ARE(to_check))
    {
        return  repmat(constants::nan(), n, n);
    }

    Matrix      Q_wave;
    Matrix      R_wave;
    permvec     Q_wave_p;
    permvec     R_wave_p;
    Integer     Q_rank;
    Integer     R_rank;
    Real        Q_e;
    Real        R_e;

    std::tie(Q_wave, Q_wave_p, Q_rank, Q_e) = cholmod(Q,true);
    std::tie(R_wave, R_wave_p, R_rank, R_e) = cholmod(R,true);

    if (Q_e > 1e-8)
        throw error::error_nonposdef(false);

    if (R_e > 1e-8)
        throw error::error_nonposdef(false);

    Q_wave      = trans(Q_wave(colon(), Q_wave_p.invperm()));
    R_wave      = trans(R_wave(colon(), R_wave_p.invperm()));

    Matrix K    = spzeros(n, R.rows());
    Matrix Z    = spzeros(n, 0);

    Matrix lastG, lastF;

    Matrix sp_B = sparse(B);

    for (int k = 0 ; k < kmax; k++)
    {
        Matrix F = A - sp_B * sparse(trans(K));
        Matrix G = horzcat(trans(C) * Q_wave, K * R_wave);

        if (k > 0 && (     norm(lastF - F, -2) / norm(lastF, -2) < lyapunov_tol 
                        && norm(lastG - G, -2) / norm(lastG, -2) < lyapunov_tol)    ) 
            return Z;

        Z       = sparse_cale(F, G, l_zero_cap, V_update_tol);
        K       = real(Z * trans (linsolve (trans(R), trans(ctrans(Z) * B))));
        lastF   = F;
        lastG   = G;
    }
    
    return Z;
}

} // end namespace details


mat_tup_3 solve_dare(const Matrix& A, const Matrix& B, const Matrix& Q,
                     const Matrix& R, const Matrix& L)
{
    return details::solve_riccati(A, B, Q, R, L, details::riccati_kind::DARE);
}

mat_tup_3 solve_care(const Matrix& A, const Matrix& B, const Matrix& Q,
                     const Matrix& R, const Matrix& L)
{
    return details::solve_riccati(A, B, Q, R, L, details::riccati_kind::CARE);
}

Matrix sparse_dare(const Matrix& A, const Matrix& B, const Matrix& C, const Matrix& Q,
                   const Matrix& R, const Matrix& S, const Real lyapunov_tol, const Integer kmax, 
                   const Integer l_zero_cap)
{
    return details::sparse_solve_driccati(A, B, C, Q, R, S, lyapunov_tol, kmax, l_zero_cap);
}

Matrix sparse_care(const Matrix& A, const Matrix& B, const Matrix& C, const Matrix& Q,
                   const Matrix& R, const Real lyapunov_tol, const Integer kmax, 
                   const Integer l_zero_cap, const Real V_update_tol)
{
    return details::sparse_solve_criccati(A, B, C, Q, R, lyapunov_tol, kmax, l_zero_cap, V_update_tol);
}

}

#endif