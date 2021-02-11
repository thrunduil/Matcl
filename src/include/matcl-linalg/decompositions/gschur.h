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

#include "matcl-matrep/matrix/matrix.h"
#include "matcl-linalg/general/config_linalg.h"

namespace matcl
{

// perform generalized schur decomposition of the matrix pair (A,B) 
// in the form
//
//          A * Z = Q * TA,    Q' * Q = I
//          B * Z = Q * TB,    Z' * Z = I
//
//  where Q, Z are unitary matrices, TA is quasi-uppertriangular, 
//  and TB is uppertriangular
class MATCL_LINALG_EXPORT gschur_decomposition
{
    public:
        gschur_decomposition();        
        ~gschur_decomposition();

        //------------------------------------------------------------------------
        //                      FACTORIZATION
        //------------------------------------------------------------------------

        // compute decomposition of the matrix pair A, B; if with_QZ is false, then
        // unitary matrices Q, Z are not computed
        // not available for sparse and band matrices
        gschur_decomposition(const Matrix &A, const Matrix &B, bool with_QZ = true);
        gschur_decomposition(Matrix &&A, const Matrix& B, bool with_QZ = true);
        gschur_decomposition(const Matrix &A, Matrix&& B, bool with_QZ = true);
        gschur_decomposition(Matrix&& A, Matrix&& B, bool with_QZ = true);
        
        // create from existing decomposition, Q, Z are unitary matrices,
        // TA is quasi-uppertriangular and TB is uppertriangular
        // Q,Z need not be square
        gschur_decomposition(const Matrix &Q, const Matrix &Z, const Matrix &TA, const Matrix &TB);

        // compute decomposition of another matrix pair; see constructor for details
        gschur_decomposition&   operator()(const Matrix &A, const Matrix &B, bool with_QZ = true);
        gschur_decomposition&   operator()(Matrix&& A, const Matrix &B, bool with_QZ = true);
        gschur_decomposition&   operator()(const Matrix &A, Matrix&& B, bool with_QZ = true);
        gschur_decomposition&   operator()(Matrix&& A, Matrix&& B, bool with_QZ = true);
        
        // change decomposition; Q,Z must be unitary matrices, TA quasi-uppertriangular, TB is
        // triangular; Q,Z need not be square
        gschur_decomposition&   set_factors(const Matrix &Q, const Matrix &Z, const Matrix &TA, 
                                        const Matrix &TB);        

        // create from existing decomposition; TA quasi-uppertriangular, TB is triangular;
        // unitary matrix will not be available
        gschur_decomposition&   set_factors(const Matrix &TA, const Matrix &TB);        

        // return the left unitary factor
        Matrix                  Q() const;
        
        // return the right unitary factor.
        Matrix                  Z() const;
        
        // return the quasi-uppertriangular factor of A.
        Matrix                  TA() const;
                
        // return the uppertriangular factor of B.
        Matrix                  TB() const;
        
        // return the generalized eigenvalues.
        Matrix                  eig() const;
        
        // return alpha's for which alpha A - beta B is singular.
        Matrix                  alpha() const;
        
        // return beta's for which alpha A - beta B is singular.
        Matrix                  beta() const;
        
        // reorder eigenvalues; eigenvalus at positions i with ind(i) == 1 will be moved to 
        // upper block of TA and TB
        void                    select(const Matrix& ind);

        //------------------------------------------------------------------------
        //                      EIGENVECTORS
        //------------------------------------------------------------------------

        // compute all left eigenvectors of the input matrix pair (A,B); if unitary 
        // matrices are not stored, then it is assumed, that Q, Z are the identity matrix
        // (therefore eigenvectors of TA and TB are returned); left eigenvector y for
        // i-th eigenvalue l_i is defined as:
        //         ctrans(y) * TA = l_i * ctrans(y) * TB
        // Matrix of eigenvectors Y satisfies
        //         Y' * TA  = E * Y' * TB, where E = diag(eig())
        Matrix                  left_eigenvectors() const;

        // compute some left eigenvectors of the input matrix pair (A,B) specified by 
        // the matrix ind; if ind(i) == 1, then compute eignevector for i-th generalized 
        // eigenvalue; return matrix of size N x M, where N is size of the matrix A, and
        // M is number of selected eigenvalues; matrix of eigenvectors Y satisfies
        //         Y' * TA  = E * Y' * TB, where E = diag(eig()(J)) with J = find(ind)
        Matrix                  left_eigenvectors(const Matrix& ind) const;
        Matrix                  left_eigenvectors(Matrix&& ind) const;

        // compute all right eigenvectors of the input matrix pair (A,B); if unitary 
        // matrice Q, Z are not stored, then it is assumed, that Q,Z are the identity 
        // matrice (therefore eigenvectors of TA and TB are returned); right eigenvector
        // y for i-th generalized eigenvalue l_i is defined as:
        //         TA * y = l_i * TB * y
        // Matrix of eigenvectors Y satisfies
        //         TA * Y = TB * Y * E, where E = diag(eig())
        Matrix                  right_eigenvectors() const;

        // compute some right eigenvectors of the input matrix pair (A,B) specified by the 
        // matrix ind; if ind(i) == 1, then compute eignevector for i-th generalized eigenvalue;
        // return matrix of size N x M, where N is size of the matrix A, and M is number of 
        // selected eigenvalues; matrix of eigenvectors Y satisfies
        //         TA * Y = TB* Y * E, where E = diag(eig()(J)) with J = find(ind)
        Matrix                  right_eigenvectors(const Matrix& ind) const;
        Matrix                  right_eigenvectors(Matrix&& ind) const;

        // compute left and right eigenvectors; equivalent to 
        // tuple<Matrix,Matrix>(left_eigenvectors(), right_eigenvectors())
        tuple<Matrix,Matrix>    eigenvectors() const;

        // compute left and right eigenvectors; equivalent to 
        // tuple<Matrix,Matrix>(left_eigenvectors(ind), right_eigenvectors(ind))
        tuple<Matrix,Matrix>    eigenvectors(const Matrix& ind) const;
        tuple<Matrix,Matrix>    eigenvectors(Matrix&& ind) const;

        //------------------------------------------------------------------------
        //                      CONDITION NUMBERS
        //------------------------------------------------------------------------

        // estimate reciprocal condition number for all eigenvalues; return matrix of size N 
        // (size of the matrix A); function requires left eigenvectors VL and right eigenvectors
        // VR as returned by left_eigenvectors() and right_eigenvectors()
        //
        // The reciprocal of the condition number of an eigenvalue lambda is defined as in Lapack's
        // function tgsna:
        //         S(lambda) = sqrt(|u'*A*v|^2 + |u'*B*V|^2) / (norm(u)*norm(v))
        // where u and v are the right and left eigenvectors of (A,B) corresponding to lambda;
        // v' denotes the conjugate transpose of v, and norm(u) denotes the Euclidean norm.        
        // An approximate error bound on the chordal distance between the i-th computed generalized
        // eigenvalue lambda w and the exact eigenvalue lambda is
        //          chord(w, lambda) <=  EPS * norm([A,B]) / S(i)
        // where EPS is the machine precision, and norm denotes the second norm of a matrix
        Matrix                  rcond_eig(const Matrix& VL, const Matrix& VR) const;

        // estimate reciprocal condition number for selected eigenvalues specified
        // by the matrix ind; if ind(i) == 1, then i-th eigenvalue is selected; return matrix of
        // size M, where M is number of selected eigenvalues; function requires left eigenvectors VL
        // and right eigenvectors VR as returned by left_eigenvectors(ind) and right_eigenvectors(ind)
        Matrix                  rcond_eig(const Matrix& VL, const Matrix& VR, const Matrix& ind) const;

        // estimate reciprocal condition number for all eigenvalues; return matrix of size N 
        // (size of the matrix A); this function internally computes left and right eigenvectors
        Matrix                  rcond_eig() const;

        // estimate reciprocal condition number for selected eigenvalues specified
        // by the matrix ind; if ind(i) == 1, then i-th eigenvalue is selected; return matrix of
        // size M, where M is number of selected eigenvalues; this function internally computes
        // left and right eigenvectors
        Matrix                  rcond_eig(const Matrix& ind) const;

        // estimate reciprocal condition number for all eigenvectors; return matrix of size N 
        // (size of the matrix A)
        //
        // The reciprocal of the condition number of the right eigenvector u and left eigenvector
        // v corresponding to the generalized eigenvalue lambda is defined as in Lapack's function
        // tgsna
        //
        //  An approximate error bound for a computed right and left eigenvector VL(i), R(i) 
        // is given by
        //             EPS * norm([A,B]) / DIF(i)
        // where EPS is the machine precision, norm denotes the second norm and DIF is estimated
        // condition number
        Matrix                  rcond_vec() const;

        // estimate reciprocal condition number for selected eigenvectors specified
        // by the matrix ind; if ind(i) == 1, then eigenvector associated with i-th eigenvalue
        // is selected; return matrix of size M, where M is number of selected eigenvalues
        Matrix                  rcond_vec(const Matrix& ind) const;
        
        // return the reciprocal norm of the projectors on the left (PL) and right (PR) eigenspaces 
        // associated with (A11, B11), where A11 and B11 are upper left blocks of TA and TB of size
        // M (or M+1 is the M-th eigenvalue is complex); this gives an approximate (asymptotic) bound
        // on the average absolute error S of the selected cluster of first M eigenvalues given by
        //               EPS * norm([A,B]) / PL
        // where EPS is the machine precision and norm is the second norm; the reciprocal norms PL, PR
        // satisfy 0 <= PL, PR <= 1; they are defined as in the Lapack's function tgsen
        // usage: [PL, PR] = rcond_eig_cluster(M);
        tuple<Real,Real>        rcond_projector(Integer M) const;

        // return upper bounds on Dif_u and Dif_l based on Frobenius norm estimate, where Dif_u is the
        // separation between matrix pairs (A11, B11) and (A22, B22), where A11 and B11 are upper left
        // blocks of TA and TB of size M (or M+1 is the M-th eigenvalue is complex) and A22, B22 are 
        // lower right blocks of TA, TB of size N-M (or N-M-1):
        //     Dif_u[(A11,B11), (A22,B22)] = inf_{||[L,R]||_F=1} ||[A11*R-L*A22, B11*R-L*B22]||_F
        //     Dif_l[(A11,B11),(A22,B22)]  = Dif_u[(A22,B22),(A11,B11)]
        //
        // when Dif_l is small, small changes in (A, B) can cause large changes in the deflating subspace.
        // An approximate (asymptotic) bound on the maximum angular error in the computed deflating 
        // subspaces is
        //                 EPS * norm([A, B]) / Dif_l
        // where EPS is the machine precision and norm is the second norm
        // usage: [Dif_u, Dif_l] = separation_est(M);
        tuple<Real,Real>        separation_est(Integer M) const;

        // estimate Dif_u and Dif_l based on 1-norm estimate; see separation_est for details
        tuple<Real,Real>        separation(Integer M) const;

        // return a lower bound on the smallest Frobenius norm of perturbation (E,F) for which
        // the first M eigenvalues may move and coalesce with the last N-M eigenvalues for the 
        // perturbed matrix pair (A + E, B + F); this bound depends on reciprocal norm of the
        // projectors and separations Dif_l, Dif_u; if use_est then separations are estimated
        // using separation_est function, otherwise are calculated using separation function 
        Real                    pert_bound(Integer M, bool use_est = true) const;

    private:
        void                    set_alpha_beta_eig();

        void                    set(const Matrix &Q, const Matrix &Z, const Matrix &TA, const Matrix &TB, 
                                        bool with_QZ);
        void                    compute(const Matrix &A, const Matrix &B, bool with_QZ);
        template<class V> void  compute_val(const V &A, const V &B);
        template<class V> void  select_val(const Matrix& ind, const bool no_fast_exit = false);

        Matrix                  comp_eig(const Matrix& alpha, const Matrix& beta) const;
        bool                    test_factors();
        void                    clear();
        Matrix                  convert_value(const Matrix& mat, value_code vc) const;
        void                    comp_eigenvectors(Matrix& ind, Matrix& XL, Matrix& XR, 
                                    bool comp_left, bool comp_right) const;
        template<class Val>
        Matrix                  calc_eig_inf(const Matrix& alpha, const Matrix& I) const;
        void                    rcond_est(Integer M, Integer type, Real& ret1, Real& ret2) const;

    private:
        Matrix                  m_Q_factor;
        Matrix                  m_Z_factor;
        Matrix                  m_TA;
        Matrix                  m_TB;
        Matrix                  m_eig;
        Matrix                  m_alpha;
        Matrix                  m_beta;
        bool                    m_has_QZ;
        bool                    m_is_nan;

        template<class Val>
        friend struct details::gschur_str;
};

};