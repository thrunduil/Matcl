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
#pragma once

#include "matcl-matrep/matrix/matrix.h"
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-linalg/special_matrices/matrix_functors.h"
#include "matcl-core/options/matcl_options.h"
#include "matcl-linalg/options/options_speigs.h"
#include "matcl-linalg/decompositions/schur.h"

#pragma warning(push)
#pragma warning(disable:4251)  //needs to have dll-interface to be used by clients of class

namespace matcl
{

// which eigenvalues should be found
enum class cluster_type
{
    LM,     // want the k eigenvalues of largest absolute value
    SM,     // want the k eigenvalues of smallest absolute value

    LR,     // want the k eigenvalues of largest real part
            // not available for symmetric problems
    SR,     // want the k eigenvalues of smallest real part
            // not available for symmetric problems
    
    LI,     // want the k eigenvalues of largest imaginary part
            // not available for symmetric problems
    SI,     // want the k eigenvalues of smallest imaginary part
            // not available for symmetric problems

    LA,     // compute the k largest eigenvalues (taking signs into account)
            // not available for complex and nonsymmetric problems
    SA,     // compute the k smallest eigenvalues (taking signs into account)
            // not available for complex and nonsymmetric problems
    
    BE,     // compute k eigenvalues, half from each end of the spectrum;
            // when k is odd, compute one more from the high end than from the low end
            // not available for complex and nonsymmetric problems
};

// perform partial Schur decomposition of a linear operator A in the form
//     A * U = U * T,  U' * U = I      (1)
// where T is upper quasi triangular or diagonal if A is real and hermitian;
// this allows to find a few eigenvalues and eigenvectors of the operator A;
// this decomposition is found using Implicitly Restarted Arnoldi Method; 
// this class is a wrapper around Arpack;
// Arpack is not reentrant, operator A cannot call any function of this class
class MATCL_LINALG_EXPORT pschur_decomposition : protected schur_decomposition
{
    protected:
        using impl_type = std::shared_ptr<details::pschur_impl>;
        using base_type = schur_decomposition;

        friend details::pschur_impl;

    protected:
        impl_type               m_impl;
        bool                    m_initialized;

    public:
        pschur_decomposition();
        ~pschur_decomposition();

        //------------------------------------------------------------------------
        //                      FACTORIZATION
        //------------------------------------------------------------------------

        // find k eigenvalues of a linear operator A of size N x N; specialized 
        // algorithm is used if A is hermitian and real; 1 <= k <= N - 3;
        // see namespace opt::speigs for available options
        pschur_decomposition(const linear_operator& A, Integer k, cluster_type ec = cluster_type::LM,
                const options& opts = options());

        // find k eigenvalues of a linear operator A of size N x N; initial vector
        // of the iteration is additionally supplied
        pschur_decomposition(const linear_operator& A, Integer k, const Matrix& x0, 
                cluster_type ec = cluster_type::LM, const options& opts = options());

        // compute decomposition of another matrix; see constructors for details
        pschur_decomposition&  operator()(const linear_operator& A, Integer k, 
                                    cluster_type ec = cluster_type::LM, 
                                           const options& opts = options());
        pschur_decomposition&  operator()(const linear_operator& A, Integer k, const Matrix& x0, 
                                    cluster_type ec = cluster_type::LM, 
                                    const options& opts = options());

        //------------------------------------------------------------------------
        //                 SHIFT-INVERSE FACTORIZATION
        //------------------------------------------------------------------------

        // find k eigenvalues of a matrix A of size N x N; if sing = false, 
        // then As = A - sig*I must be nonsinguar, otherwise can be singular; 
        //
        // cluster type ec defines which eigenvalues of inv(As) should be found, 
        // these eigenvalues are given by mu_i = 1/(lam_i - sig) with lam - eigenvalue
        // of A, mu - eigenvalue of inv(As); if ec = LM, then largest eigenvalues of
        // As will be found, i.e. smallest eigenvalues of A;
        //
        // specialized algorithm is used if A is hermitian and real; 1 <= k <= N - 3
        // see namespace opt::speigs for available options; 
        // see also linsolve_shift_inverse
        pschur_decomposition(const Matrix& A, Real sig, bool sing, Integer k, 
                cluster_type ec = cluster_type::LM, const options& opts = options());

        // find k eigenvalues of a matrix A of size N x N; initial vector of
        // the iteration is additionally supplied;
        pschur_decomposition(const Matrix& A, Real sig, bool sing, Integer k, 
                const Matrix& x0, cluster_type ec = cluster_type::LM, const options& opts = options());

        // compute decomposition of another matrix; see constructors for details
        pschur_decomposition&  operator()(const Matrix& A, Real sig, bool sing, 
                                    Integer k, cluster_type ec = cluster_type::LM, 
                                    const options& opts = options());
        pschur_decomposition&  operator()(const Matrix& A, Real sig, bool sing, 
                                    Integer k, const Matrix& x0, cluster_type ec = cluster_type::LM, 
                                    const options& opts = options());

        //------------------------------------------------------------------------
        //              FUNCTIONS FROM SCHUR FACTORIZATION
        //------------------------------------------------------------------------
        // see schur_factorization for details

        // return the unitary factor such that (1) holds
        using base_type::U;

        // return the quasi-uppertriangular factor of A; TA is real diagonal if
        // A is real hermitian
        using base_type::TA;

        // return the eigenvalues of A (and also of T in (1) )
        using base_type::eig;

        // reorder eigenvalues; eigenvalus at positions i with ind(i) == 1 will be moved to 
        // upper block of T
        using base_type::select;

        // compute left eigenvectors
        using base_type::left_eigenvectors;

        // compute right eigenvectors
        using base_type::right_eigenvectors;

        // compute left and right eigenvectors; equivalent to 
        // tuple<Matrix,Matrix>(left_eigenvectors(), right_eigenvectors())
        using base_type::eigenvectors;

        // estimate reciprocal condition number for eigenvalues;
        using base_type::rcond_eig;

        // estimate reciprocal condition number for eigenvectors
        using base_type::rcond_vec;

        // return an approximate bound on the average absolute error of the selected 
        // cluster of first M eigenvalues
        using base_type::rcond_projector;

        // return estimated separation
        using base_type::separation;

        // return a lower bound on the smallest Frobenius norm of perturbation E for which
        // the first M eigenvalues may move and coalesce with the last N-M eigenvalues for the 
        // perturbed matrix (A + E); 
        using base_type::pert_bound;

        //------------------------------------------------------------------------
        //                      DIAGONOSTICS
        //------------------------------------------------------------------------

        // return true if all required eigenvalues converged
        bool                    converged() const;

        // return number of converged eigenvalues
        Integer                 number_converged_eigenvalues() const;

        // profiling statistics
        Integer                 number_Arnoldi_iterations() const;
        Integer                 number_operations_Ax() const;
        Integer                 number_reorthogonalizations() const;

    protected:
        void                    initialize(const linear_operator& A, const options& opts);
        void                    compute(Integer k, bool return_nonconverged);
        void                    transform(const Matrix& A);
        void                    set_x0(const Matrix& x0);
        void                    set_cluster(cluster_type ec);
        void                    check() const;
};

// perform partial Schur decomposition of a linear operator A in the form
//     A * U = U * T,  U' * B * U = I                          (1)
// where T is upper quasi triangular and B is hermitian and semi positive 
// definite; from (1) we have that U spans invariant subspace of the operator
// A and eigenvalues of T are also eigenvalues of A. The partial Schur
// decomposition can be obtained by taking qr decompostion: Q*R = U, then
//     A * Q = Q * (R * T * R^-1), Q' * Q = I
//
// if the operator A is hermitian with respect to inner product given by the
// matrix B, then T is diagonal;  A is hermitian with respect to inner product
// given by B iff
//     B * A = A' * B, or equivalently < x,Ay > = < Ax,y >     (2)
// where <z,w> = z'Bw.
//
// representation (1) may help finding eigenvalues of the matrix pair (X,Y)
// using shift and invert method by taking:
//     A = inv[X - sigma*Y] * Y, B = Y or
//     A = inv[Y] * X, B = Y
// where sigma is a shift, such that X - sigma*Y is invertible, (such shift
// always exists if the matrix pair (X, Y) is regular), especially if Y is
// nearly singular or singular; however if Y can be factored into a Cholesky
// factorization Y = LL', then the pschur decomposition should be used for
// A = inv[L] * X * inv[L']
//
// pbschur_decomposition allows to find a few eigenvalues and eigenvectors 
// of the operator A; this decomposition is found using Implicitly Restarted
// Arnoldi Method;  this class is a wrapper around Arpack;
// Arpack is not reentrant, operators A, B cannot call any function from this
// class
class MATCL_LINALG_EXPORT pbschur_decomposition : public pschur_decomposition
{
    private:
        using base_type = pschur_decomposition;

    public:
        pbschur_decomposition();
        ~pbschur_decomposition();

        //------------------------------------------------------------------------
        //                      FACTORIZATION
        //------------------------------------------------------------------------

        // find k eigenvalues of a linear operator A of size N x N; specialized 
        // algorithm is used if A is real and hermitian with respect to semi-inner 
        // product given by semi positive definite hermitian and real operator B;
        // this case is assumed if hermitian = true; 1 <= k <= N - 3
        pbschur_decomposition(const linear_operator& A, const linear_operator& B, bool hermitian,
                              Integer k, cluster_type ec = cluster_type::LM, 
                              const options& opts = options());

        // find k eigenvalues of a linear operator A of size N x N; specialized 
        // algorithm is used if A is real and hermitian with respect to semi-inner 
        // product given by semi positive definite hermitian and real operator B;
        // this case is assumed if hermitian = true; 1 <= k <= N - 3; 
        // initial vector of the iteration is supplied;
        pbschur_decomposition(const linear_operator& A, const linear_operator& B, bool hermitian, 
                              Integer k, const Matrix& x0, cluster_type ec = cluster_type::LM,
                              const options& opts = options());

        // compute decomposition of another matrix; see constructors for details
        pbschur_decomposition&  operator()(const linear_operator& A, const linear_operator& B, 
                                    bool hermitian, Integer k, cluster_type ec = cluster_type::LM, 
                                    const options& opts = options());
        pbschur_decomposition&  operator()(const linear_operator& A, const linear_operator& B,
                                    bool hermitian, Integer k, const Matrix& x0, 
                                    cluster_type ec = cluster_type::LM, const options& opts = options());

        // return the U factor such that (1) holds (U is not unitary)
        Matrix                  U() const;

        //------------------------------------------------------------------------
        //                      DIAGONOSTICS
        //------------------------------------------------------------------------

        // profiling statistics
        Integer                 number_operations_Bx() const;

    private:
        void                    initialize_G(const linear_operator& A, const linear_operator& B, 
                                           bool hermitian, const options& opts);

        friend details::pschur_impl;
};

// construct linsolve_obj for a square matrix As = A - sig * I
// this function calls make_linsolve_obj for the matrix As
//     lo  = linsolve_shift_invert(A, sig, sing, opts)
// A       - a square matrix;
// sig     - a shift
// sing    = true:   As is possibly singular, in this case rank revealing
//           factorization is used and small perturbation is added if As
//           is detected to be singular
//         = false:  As is nonsingular, exception is thrown if As is 
//           detected to be singular
// opts    - options used by factorization routines, see namespace 
//           opt::linsolve for details; see also make_linsolve_obj
// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_shift_invert(const Matrix& A, Real sig, bool sing,
                                                       const options& opts = options());

};

#pragma warning(pop)