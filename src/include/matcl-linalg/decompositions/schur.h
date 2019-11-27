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
#include "matcl-linalg/linear_eq/linsolve_object.h"
#include "matcl-linalg/special_matrices/unitary_matrix.h"
#include "matcl-matrep/details/enablers.h"

namespace matcl
{

namespace md = matcl::details;

/// algorithm used for symmetric eigenvalue problem
enum class schur_sym_alg
{
    qr,     /// qr algorithm
    dc,     /// divide-and-conquer
    rrr     /// relatively robust representation, this algorithm is faster but may yield 
            /// non-unitary (non-orthogonal) U; unitary flag will not be set for unitary matrix U
            /// in the schur decomposition; rrr algorithm is not available for band matrices, in
            /// this case dc algorithm is used
};

/// perform Schur decomposition of the matrix A in the form
///
///          A  * U = U * T,    U' * U = I
///
///  where U is unitary matrix, T is quasi-uppertriangular
class MATCL_LINALG_EXPORT schur_decomposition
{
    public:
        schur_decomposition();
        ~schur_decomposition();

        //------------------------------------------------------------------------
        //                      FACTORIZATION
        //------------------------------------------------------------------------

        /// compute decomposition of the matrix A; uses specialized algoritms for
        /// symmetric/hermitian matrices if symetric/hermitian structure is assigned
        /// to A;
        /// alg describe algorithm used for symmetric problems
        /// if with_U is false, then unitary matrix is not computed
        /// not available for sparse matrices and band general matrices
        schur_decomposition(const Matrix &A, schur_sym_alg alg = schur_sym_alg::dc, bool with_U = true);
        schur_decomposition(Matrix &&A, schur_sym_alg alg = schur_sym_alg::dc, bool with_U = true);

        /// create from existing decomposition, U is unitary matrix,
        /// T is quasi-uppertriangular; U need not be square
        schur_decomposition(const Matrix &U, const Matrix &T);
        
        /// compute decomposition of another matrix; see constructor for details
        schur_decomposition&    operator()(const Matrix &A, schur_sym_alg  alg = schur_sym_alg::dc,
                                           bool with_U = true);
        schur_decomposition&    operator()(Matrix &&A, schur_sym_alg alg = schur_sym_alg::dc,
                                           bool with_U = true);
        
        /// change decomposition; U must be unitary but not necessary square, and 
        /// T quasi-uppertriangular
        schur_decomposition&    set_factors(const Matrix &T, const Matrix &U);

        /// create from existing decomposition; T is quasi-uppertriangular;
        /// unitary matrix will not be available
        schur_decomposition&    set_factors(const Matrix &T);

        /// return the unitary factor.
        Matrix                  U() const;

        /// return the quasi-uppertriangular factor of A.
        Matrix                  TA() const;

        /// return the eigenvalues.
        Matrix                  eig() const;

        /// reorder eigenvalues; eigenvalus at positions i with ind(i) == 1 will be moved to 
        /// upper block of T
        void                    select(const Matrix& ind);

        //------------------------------------------------------------------------
        //                      EIGENVECTORS
        //------------------------------------------------------------------------

        /// compute all left eigenvectors of the input matrix A; if unitary matrix U
        /// is not stored, then it is assumed, that U is the identity matrix (therefore 
        /// eigenvectors of TA are returned); left eigenvector y for i-th eigenvalue l_i 
        /// is defined as:
        ///         ctrans(y) * A = l_i * ctrans(y)
        /// Matrix of eigenvectors Y satisfies
        ///         Y' * A  = E * Y', where E = diag(eig())
        Matrix                  left_eigenvectors() const;

        /// compute some left eigenvectors of the input matrix A specified by the matrix ind;
        /// if ind(i) == 1, then compute eignevector for i-th eigenvalue; return matrix of size
        /// N x M, where N is size of the matrix A, and M is number of selected eigenvalues
        /// Matrix of eigenvectors Y satisfies
        ///         Y' * A  = E * Y', where E = diag(eig()(J)) with J = find(ind)
        Matrix                  left_eigenvectors(const Matrix& ind) const;
        Matrix                  left_eigenvectors(Matrix&& ind) const;

        /// compute all right eigenvectors of the input matrix A; if unitary matrix U
        /// is not stored, then it is assumed, that U is the identity matrix (therefore 
        /// eigenvectors of TA are returned); right eigenvector y for i-th eigenvalue l_i 
        /// is defined as:
        ///         A * y = l_i * y
        /// Matrix of eigenvectors Y satisfies
        ///         A * Y = Y * E, where E = diag(eig())
        Matrix                  right_eigenvectors() const;

        /// compute some right eigenvectors of the input matrix A specified by the matrix ind;
        /// if ind(i) == 1, then compute eignevector for i-th eigenvalue; return matrix of size
        /// N x M, where N is size of the matrix A, and M is number of selected eigenvalues
        /// Matrix of eigenvectors Y satisfies
        ///         A * Y = Y * E, where E = diag(eig()(J)) with J = find(ind)
        Matrix                  right_eigenvectors(const Matrix& ind) const;
        Matrix                  right_eigenvectors(Matrix&& ind) const;

        /// compute left and right eigenvectors; equivalent to 
        /// tuple<Matrix,Matrix>(left_eigenvectors(), right_eigenvectors())
        tuple<Matrix,Matrix>    eigenvectors() const;

        /// compute left and right eigenvectors; equivalent to 
        /// tuple<Matrix,Matrix>(left_eigenvectors(ind), right_eigenvectors(ind))
        tuple<Matrix,Matrix>    eigenvectors(const Matrix& ind) const;
        tuple<Matrix,Matrix>    eigenvectors(Matrix&& ind) const;

        //------------------------------------------------------------------------
        //                      CONDITION NUMBERS
        //------------------------------------------------------------------------

        /// estimate reciprocal condition number for all eigenvalues; return matrix of size N 
        /// (size of the matrix A); function requires left eigenvectors VL and right eigenvectors
        /// VR as returned by left_eigenvectors() and right_eigenvectors()
        ///
        /// The reciprocal of the condition number of an eigenvalue lambda is defined as in Lapack's
        /// function trsna:
        ///         S(lambda) = |v'*u| / (norm(u)*norm(v))
        /// where u and v are the right and left eigenvectors of T corresponding to lambda;
        /// v' denotes the conjugate transpose of v, and norm(u) denotes the Euclidean norm. 
        /// These reciprocal condition numbers always lie between zero (very badly conditioned) 
        /// and one (very well conditioned). If n = 1, S(lambda) is defined to be 1.
        /// An approximate error bound for a computed eigenvalue W(i) is given by
        ///             EPS * norm(TA) / S(i)
        /// where EPS is the machine precision
        Matrix                  rcond_eig(const Matrix& VL, const Matrix& VR) const;

        /// estimate reciprocal condition number for selected eigenvalues specified
        /// by the matrix ind; if ind(i) == 1, then i-th eigenvalue is selected; return matrix of
        /// size M, where M is number of selected eigenvalues; function requires left eigenvectors VL
        /// and right eigenvectors VR as returned by left_eigenvectors(ind) and right_eigenvectors(ind)
        Matrix                  rcond_eig(const Matrix& VL, const Matrix& VR, const Matrix& ind) const;

        /// estimate reciprocal condition number for all eigenvalues; return matrix of size N 
        /// (size of the matrix A); this function internally computes left and right eigenvectors
        Matrix                  rcond_eig() const;

        /// estimate reciprocal condition number for selected eigenvalues specified
        /// by the matrix ind; if ind(i) == 1, then i-th eigenvalue is selected; return matrix of
        /// size M, where M is number of selected eigenvalues; this function internally computes
        /// left and right eigenvectors
        Matrix                  rcond_eig(const Matrix& ind) const;

        /// estimate reciprocal condition number for all eigenvectors; return matrix of size N 
        /// (size of the matrix A)
        ///
        /// The reciprocal of the condition number of the right eigenvector u corresponding to 
        /// lambda is defined as in Lapack's function trsna:
        ///     SEP( lambda, T22 ) = sigma-min( T22 - lambda*I )
        /// where sigma-min denotes the smallest singular value. T22 is submatrix of TAO obtained
        /// by reordering given eigenvalue to the first position and T22 = TAO(2:end,2:end). 
        /// If n = 1, SEP(1) is defined to be abs(TA(1,1)). 
        ///
        ///  An approximate error bound for a computed right eigenvector VR(i) is given by
        ///             EPS * norm(TA) / SEP(i)
        /// where EPS is the machine precision
        Matrix                  rcond_vec() const;

        /// estimate reciprocal condition number for selected eigenvectors specified
        /// by the matrix ind; if ind(i) == 1, then eigenvector associated with i-th eigenvalue
        /// is selected; return matrix of size M, where M is number of selected eigenvalues
        Matrix                  rcond_vec(const Matrix& ind) const;

        /// return the reciprocal norm of the projector (S) on the eigenspace associated with A11,
        /// where A11 is upper left blocks of TA of size M (or M+1 is the M-th eigenvalue is 
        /// complex); this gives an approximate bound on the average absolute error of the 
        /// selected cluster of first M eigenvalues given by
        ///               EPS * norm(TA) / S
        /// where EPS is the machine precision and norm is the second norm; the reciprocal norms S
        /// satisfies 0 <= S <= 1; it is defined as in the Lapack's function trsen
        Real                    rcond_projector(Integer M) const;

        /// return estimated separation between matrices A11 and A22, where A11 is upper left
        /// blocks of TA of size M (or M+1 is the M-th eigenvalue is complex) and A22 is lower
        /// right blocks of TA of size N-M (or N-M-1):
        ///             Sep(A11, A22) = inf_{||X||_F=1} ||A11*X-X*A22||_F
        ///
        /// when Sep is small, small changes in TA can cause large changes in the invariant subspace;
        /// an approximate bound on the maximum angular error in the computed right invariant subspace
        /// is       
        ///                 EPS * norm(TA) / SEP
        /// where EPS is the machine precision
        Real                    separation(Integer M) const;

        /// return a lower bound on the smallest Frobenius norm of perturbation E for which
        /// the first M eigenvalues may move and coalesce with the last N-M eigenvalues for the 
        /// perturbed matrix (A + E); 
        Real                    pert_bound(Integer M) const;

    private:
        void                    set_eig();
        void                    set(const Matrix &T, const Matrix &U, bool with_U);
        void                    compute(const Matrix &A, schur_sym_alg alg, bool with_U);        
        void                    comp_eigenvectors(Matrix& ind, Matrix& XL, Matrix& XR, 
                                    bool comp_left, bool comp_right) const;
        void                    select_val(const Matrix& ind, bool no_fast_exit = false);
        void                    clear();
        bool                    test_factors();

    private:
        Matrix                  m_U_factor;
        mutable Matrix          m_TA;
        Matrix                  m_eig;
        bool                    m_is_nan;
        bool                    m_has_U;

        template<class V, class S>
        friend struct details::schur_str;
        friend struct details::reorder_diag;
};

/// compute the eigendecomposition of a 2-by-2 hermitian matrix
///   [  A   B  ]
///   [  B'  C  ].
/// on return, eig_1 is the eigenvalue of larger absolute value, eig_2 is the
/// eigenvalue of smaller absolute value, and (cos,sin) is the unit right
/// eigenvector for eig_1, giving the decomposition
///
///   [ cos  sin ] [  A   B  ] [ cos -sin ]  =  [ eig_1  0     ]
///   [-sin' cos ] [  B'  C  ] [ sin' cos ]     [  0     eig_2 ]
///
///     schur_22_sym(A, B, C, A_22, sin, cos, eig_1, eig_2)
///
/// A, B, C         - (input) elements of the 2 x 2 hermitian matrix
/// cos, sin        - (output) cosine and sine of plane rotation
/// eig_1, eig_2    - (output) defines larger and smaller eigenvalue
///
template<class V, class Enable = typename md::enable_if_float<V,void>::type, 
        class VR = typename md::real_type<V>::type>
MATCL_LINALG_EXPORT void schur_22_sym(const VR& A, const V& B, const VR& C, 
                                    VR& cos, V& sin, VR& eig_1, VR& eig_2);

/// compute the Schur factorization of a 2-by-2 nonsymmetric matrix in 
/// standard form:
///     [ CS   SN ] [ A  B ] [ CS  -SN ] =  [ AA  BB ] 
///     [ -SN' CS ] [ C  D ] [ SN'  CS ]    [ CC  DD ] 
///
/// where either
///     1) CC = 0 so that AA and DD are eigenvalues of the matrix, or
///     2) the matrix is real and AA = DD and BB * CC < 0, so that 
///         AA +- sqrt(BB*CC) are complex conjugate eigenvalues.
///
///     schur_22(A, B, C, D, CS, SN, eig_1, eig_2)
///
/// A, B, C, D      - (input/output) elements of the 2 x 2 matrix
/// CS, SN          - (output) cosine and sine of plane rotation
/// eig_1, eig_2    - (output) complex eigenvalues
///
template<class V, class Enable = typename md::enable_if_float<V,void>::type, 
        class VR = typename md::real_type<V>::type,
        class VC = typename md::complex_type<V>::type>
MATCL_LINALG_EXPORT void schur_22(V& A, V& B, V& C, V& D, VR& CS, V& SN,
                                VC& eig_1, VC& eig_2);

/// compute eigenvalues of a 2-by-2 hermitian matrix
///   [  A   B  ]
///   [  B'  C  ].
/// on return, eig_1 is the eigenvalue of larger absolute value, eig_2 is the
/// eigenvalue of smaller absolute value
///
///     eig_22_sym(A, B, C, A_22, eig_1, eig_2)
///
/// A, B, C         - (input) elements of the 2 x 2 hermitian matrix
/// eig_1, eig_2    - (output) defines larger and smaller eigenvalue
///
template<class V, class Enable = typename md::enable_if_float<V,void>::type, 
        class VR = typename md::real_type<V>::type>
MATCL_LINALG_EXPORT void eig_22_sym(const VR& A, const V& B, const VR& C, VR& eig_1, VR& eig_2);

/// compute the eigenvalues of the Schur factorization of a 2-by-2 nonsymmetric matrix
///     [ A  B ]
///     [ C  D ]
///
///     eig_22(A, B, C, D, eig_1, eig_2)
///
/// A, B, C, D      - (input) elements of the 2 x 2 matrix
/// eig_1, eig_2    - (output) complex eigenvalues
///
template<class V, class Enable = typename md::enable_if_float<V,void>::type, 
        class VR = typename md::real_type<V>::type,
        class VC = typename md::complex_type<V>::type>
MATCL_LINALG_EXPORT void eig_22(const V& A, const V& B, const V& C, const V& D, VC& eig_1, VC& eig_2);

/// construct linsolve_obj from Schur factors such that
///     A  * U = U * T,    U' * U = I
///
///     lo  = linsolve_schur(A, U, T, opts)
/// A       - the factored matrix
/// U       - unitary matrix
/// T       - square quasi upper triangular matrix
/// opts    - options, see namespace opt::linsolve for details
/// lo      - a linsolve_obj
/// 
/// if T is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_schur(const Matrix& A, const Matrix& U, const Matrix& T,
                                                const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_schur(const Matrix& A, const unitary_matrix& U, const Matrix& T,
                                                const options& opts = options());

};