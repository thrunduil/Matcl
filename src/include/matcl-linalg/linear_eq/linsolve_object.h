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
#include "matcl-matrep/matrix/permvec.h"
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-linalg/options/options_linsolve.h"
#include "matcl-linalg/details/linalg_fwd.h"
#include "matcl-linalg/norms_error/norm.h"

#pragma warning (push)
#pragma warning (disable: 4251) //needs to have dll-interface to be used by

namespace matcl
{

/// internal representation of linsolve_obj
class MATCL_LINALG_EXPORT linsolve_obj_data
        : public std::enable_shared_from_this<linsolve_obj_data>
{
    public:
        /// shared pointer type
        using data_ptr  = std::shared_ptr<linsolve_obj_data>;

    public:
        virtual ~linsolve_obj_data(){};

        /// get number of rows
        virtual Integer         rows() const = 0;

        /// get number of columns
        virtual Integer         cols() const = 0;

        /// code of stored elements type
        virtual value_code      get_value_code() const = 0;

        /// return type_info of stored elements
        virtual ti::ti_object   get_type() const;

        /// check if all elements are finite
        virtual bool            all_finite() const = 0;

        /// should return true if the matrix is real symmetric or complex hermitian
        virtual bool            is_hermitian() const
                                    { return false; };

        /// should return true if the matrix is positive definite hermitian matrix
        virtual bool            is_posdef() const
                                    { return false; };

        /// return true if small perturbations was added to make factors nonsingular
        virtual bool            is_modified() const = 0;

        /// return true if this object represents a direct solver and false if this
        /// object represents an iterative solver
        virtual bool            is_direct() const = 0;

        /// this function is called when expected value type is different 
        /// than value type of stored elements and conversion is expected;
        /// given linsolve_obj may ignore this request or perform required
        /// conversion, in this case this function must be overriden;
        virtual data_ptr        convert(value_code new_val_code) const = 0;

        /// create inverse matrix; default implementation solves A * Ainv = I
        /// for the inverse Ainv
        virtual matcl::Matrix   inv() const;

        /// create the matrix A, default implementation calls mmul_right(I),
        /// where I is the identity matrix;
        virtual matcl::Matrix   base_matrix() const;

        /// solve op(A) * Y = X, where A is the matrix represented by this object
        virtual matcl::Matrix   solve(const Matrix& X, trans_type tA) const = 0;
        virtual matcl::Matrix   solve(Matrix&& X, trans_type tA) const = 0;

        /// solve Y * op(A) = X, where A is the matrix represented by this object;
        /// defalt implementation solves transposed problem
        virtual matcl::Matrix   solve_rev(const Matrix& X, trans_type tA) const;
        virtual matcl::Matrix   solve_rev(Matrix&& X, trans_type tA) const;

        /// return log |det(A)|, i.e. logarithm of absolute value of the determinant
        /// of the matrix A
        virtual Real            log_det() const = 0;

        /// evaluate mmul(A, X, t), i.e. trans(A,t) * X
        virtual Matrix          mmul_right(const Matrix& X, trans_type t) const = 0;
        virtual Matrix          mmul_right(Matrix&& X, trans_type t) const;

        /// evaluate mmul(X, A, t), i.e. X * trans(A,t)
        virtual Matrix          mmul_left(const Matrix& X, trans_type t) const = 0;
        virtual Matrix          mmul_left(Matrix&& X, trans_type t) const;

        /// estimate norm-1 of the operator inv(A), default implementation uses
        /// matcl::normest_1; derived class may provide different algorithm
        virtual Real            normest_1() const;

        /// estimate norm-2 of the operator inv(A), default implementation uses
        /// matcl::normest_2; derived class may provide different algorithm
        virtual Real            normest_2() const;

        /// estimate infinity norm of the operator inv(A), default implementation uses
        /// matcl::normest_inf; derived class may provide different algorithm
        virtual Real            normest_inf() const;

        /// estimate norm-1 of the operator A, default implementation uses
        /// matcl::normest_1; derived class may provide different algorithm
        virtual Real            mat_normest_1() const;

        /// estimate norm-2 of the operator A, default implementation uses
        /// matcl::normest_2; derived class may provide different algorithm
        virtual Real            mat_normest_2() const;

        /// estimate infinity norm of the operator A, default implementation uses
        /// matcl::normest_inf; derived class may provide different algorithm
        virtual Real            mat_normest_inf() const;
};

/// class containing factorization of a square matrix A allowing to
/// solve linear equations
/// scalars are always treated as 1x1 matrix, i.e. A.solve(x), where A is created from 
/// a scalar and x is not a 1xN matrix, will produced error
class MATCL_LINALG_EXPORT linsolve_obj
{
    public: 
        /// type of internal representation
        using linsolve_data_ptr = std::shared_ptr<linsolve_obj_data>;

    private:
        linsolve_data_ptr   m_impl;

    public:
        /// create linsolve functor from real scalar 1.0
        linsolve_obj();

        /// create linsolve_obj from representation rep; rep cannot be empty
        explicit linsolve_obj(const linsolve_data_ptr& rep);
        explicit linsolve_obj(linsolve_data_ptr&& rep);

        /// standard copy and move constructor
        linsolve_obj(const linsolve_obj& mat);
        linsolve_obj(linsolve_obj&& mat);

        /// standard assignment and move assignment operator
        linsolve_obj&           operator=(const linsolve_obj&) &;
        linsolve_obj&           operator=(linsolve_obj&&) &;

        /// standard destructor
        ~linsolve_obj();

        //--------------------------------------------------------------------
        //          member functions specific to linsolve_obj
        //--------------------------------------------------------------------
        //op(A, t) denotes matcl function trans(A, t)

        /// this function is called when expected value type is different 
        /// than value type of stored elements and conversion is expected;
        /// however given linsolve_obj may ignore this request
        linsolve_obj            convert(value_code new_val_code) const;

        /// should return true if the matrix is real symmetric or complex hermitian
        bool                    is_hermitian() const;

        /// should return true if the matrix is real symmetric or complex hermitian
        /// and positive definite
        bool                    is_posdef() const;

        /// return true if small perturbations was added to make factors nonsingular;
        /// only modifications made in matcl are reported
        bool                    is_modified() const;

        /// return true if this object represents a direct solver and false if this
        /// object represents an iterative solver
        bool                    is_direct() const;

        /// return log |det(A)|, i.e. logarithm of absolute value of the determinant
        /// of the matrix A;
        /// log_det is not available for linsolve_obj implementing iterative methods,
        /// in such case exception is thrown
        Real                    log_det() const;

        /// create inverse of the matrix represented by this object
        matcl::Matrix           inv() const;

        /// create the matrix A, notice that A need not be available, in this 
        /// case this is equivalent to mmul_right(eye(N));
        matcl::Matrix           base_matrix() const;

        /// solve op(A) * Y = X, where A is the matrix represented by this object
        matcl::Matrix           solve(const Matrix& X, trans_type tA = trans_type::no_trans) const;
        matcl::Matrix           solve(Matrix&& X, trans_type tA = trans_type::no_trans) const;

        /// solve Y * op(A) = X, where A is the matrix represented by this object
        matcl::Matrix           solve_rev(const Matrix& X, trans_type tA = trans_type::no_trans) const;
        matcl::Matrix           solve_rev(Matrix&& X, trans_type tA = trans_type::no_trans) const;

        /// improves the computed solution Y to a system of linear equations
        /// op(A) * Y = X; see iterative refinement options in opt::linsolve
        /// for available options
        matcl::Matrix           iterative_refinement(const Matrix& Y, const Matrix& X, 
                                        trans_type tA = trans_type::no_trans, 
                                        const options& opts = options()) const;
        matcl::Matrix           iterative_refinement(Matrix&& Y, const Matrix& X, 
                                        trans_type tA = trans_type::no_trans,
                                        const options& opts = options()) const;

        /// improves the computed solution Y to a system of linear equations
        /// Y * op(A) = X; see iterative refinement options in opt::linsolve
        /// for available options
        matcl::Matrix           iterative_refinement_rev(const Matrix& Y, const Matrix& X, 
                                        trans_type tA = trans_type::no_trans,
                                        const options& opts = options()) const;
        matcl::Matrix           iterative_refinement_rev(Matrix&& Y, const Matrix& X, 
                                        trans_type tA = trans_type::no_trans,
                                        const options& opts = options()) const;

        /// improve the computed inverse matrix I using one step of Schulz iteration
        matcl::Matrix           improve_inv(const Matrix& I) const;

        /// evaluate mmul(A, X, t), i.e. trans(A,t) * X
        Matrix                  mmul_right(const Matrix& X, trans_type tA = trans_type::no_trans) const;
        Matrix                  mmul_right(Matrix&& X, trans_type tA = trans_type::no_trans) const;

        /// evaluate mmul(X, A, t), i.e. X * trans(A,t)
        Matrix                  mmul_left(const Matrix& X, trans_type tA = trans_type::no_trans) const;
        Matrix                  mmul_left(Matrix&& X, trans_type tA = trans_type::no_trans) const;

        /// convert to linear_operator that represents the matrix A, multiplication
        /// is implemented using funcions linsolve_obj::mmul_right
        linear_operator         linear_operator_mat() const;

        /// convert to linear_operator that represents the matrix inv(A), multiplication
        /// is implemented using funcions linsolve_obj::solve; note that this
        /// is equivalent to the operator obtained from implicit conversion of 
        /// linsolve_obj to linear_operator
        linear_operator         linear_operator_inv() const;

        //--------------------------------------------------------------------
        //                      error analysis
        //--------------------------------------------------------------------

        /// estimation of the norm of inv(op(A,t))
        Real                    normest(basic_vector_norm p, trans_type t = trans_type::no_trans) const;

        /// estimation of the norm of op(A,t)
        Real                    mat_normest(basic_vector_norm p, trans_type t = trans_type::no_trans) const;

        /// reciprocal condition number estimate in p-norm defined as
        /// ki = 1.0/(|A|_p x |inv(A)|_p), where p is the first norm, the second
        /// norm or the infinity norm; if A is well conditioned, ki is near 1.0;
        /// if A is badly conditioned, ki is near 0.0; this condition number is
        /// obtained based on norm estimates returned by normest and mat_normest.
        Real                    rcond_est(basic_vector_norm p) const;

        /// return Skeel condition number of the matrix A defined as
        ///     Cond(op(A,t)) = || abs(inv(op(A,t))) * abs(op(A,t)) ||_p
        /// the matrix A and the inverse matrix inv(A) are created
        Real                    skeel_cond(basic_vector_norm p, trans_type t = trans_type::no_trans) const;

        /// return Skeel condition number defined for each column i of X as
        ///     Cond(op(A,t), X_i) = |abs(inv(op(A,t))) * abs(op(A,t)) * abs(X_i)|_inf
        ///                        / |X_i|_inf
        /// the matrix A is created but not the inverse matrix; return matrix of
        /// size 1xN, where N is the number of columns of X
        Matrix                  skeel_vec_cond(const Matrix& X, trans_type t = trans_type::no_trans) const;

        /// return Skeel condition number of the matrix A and a matrix E defined as
        ///     Cond_E(op(A,t)) = || abs(inv(op(A,t))) * abs(op(E,t)) ||_p
        /// the inverse matrix inv(A) is created
        Real                    skeel_gen_cond(const Matrix& E, basic_vector_norm p, 
                                        trans_type t = trans_type::no_trans) const;

        /// return Skeel condition defined for each column i of X as
        ///     Cond_E(op(A,t), X_i) = | abs(inv(op(A,t))) * abs(X_i) |_inf
        ///                         / |X_i|_inf
        /// the inverse matrix inv(A) is not created
        Matrix                  skeel_gen_vec_cond(const Matrix& X, trans_type t = trans_type::no_trans) const;

        /// return residuals r = op(A, t) * Y - X, where Y is returned by solve(X, t)
        Matrix                  resid(const Matrix& Y, const Matrix& X, 
                                        trans_type t = trans_type::no_trans) const;

        /// return residuals r = Y * op(A, t) - X, where Y is returned by solve_rev(X, t)
        Matrix                  resid_rev(const Matrix& Y, const Matrix& X,
                                        trans_type t = trans_type::no_trans) const;

        /// return |op(A, t) * Y - X|_p, where Y is returned by solve(X, t) and p defines
        /// vector norm, see norm_col for details
        Matrix                  resid_norm(const Matrix& Y, const Matrix& X, basic_vector_norm p,
                                        trans_type t = trans_type::no_trans) const;

        /// return |Y * op(A, t) - X|_p, where Y is returned by solve_rev(X, t) and p defines
        /// vector norm, see norm_row for details
        Matrix                  resid_rev_norm(const Matrix& Y, const Matrix& X, basic_vector_norm p,
                                        trans_type t = trans_type::no_trans) const;

        /// return normwise backward error bN and forward error fN of the solution 
        /// Y = solve(X, t), i.e.
        ///       bN = min{ bN : (op(A,t) + DA) * Y = X + DX, where
        ///                    ||DA||_p <= bN * ||op(A,t)||_p, |DX|_p <= bN * |X|_p }
        ///       fN >= |Y* - Y|_p / |Y*|_p
        /// where Y* is true solution, |*|_p is vector norm, and ||*||_p is the 
        /// corresponding subordinate matrix norm estimate (see norm_col for details); 
        /// bN and fN are given by
        ///         bN = |op(A,t) * Y - X|_p / ( ||op(A,t)||_p * |Y|_p + |X|_p ) 
        ///         fN = 2 * kA / (1 - kA), kA = bN / rcond( op(A,t), p)
        /// fN can take negative value, in this case given component of Y is meaningless
        /// [bN, fN] = norm_error(X, Y, p, t)
        mat_tup_2               norm_error(const Matrix& Y, const Matrix& X, basic_vector_norm p, 
                                        trans_type t = trans_type::no_trans) const;

        /// return normwise backward error bN and forward error fN of the solution 
        /// Y = solve_rev(X, t), i.e.
        ///       bN =  min{ bN : Y * (op(A,t) + DA) = X + DX, where
        ///                    ||DA'||_p <= bN * ||op(A,t)'||_p, |DX|_p <= bN * |X|_p }
        ///       fN >= |Y* - Y|_p / |Y*|_p
        /// where Y* is true solution, |*|_p is vector norm, and ||*||_p is the corresponding
        /// subordinate matrix norm estimate (see norm_row for details); bN is given by
        ///         bN = |Y * op(A,t) - X|_p / ( ||op(A,t)'||_p * |Y|_p + |X|_p ) 
        ///         fN = 2 * kA / (1 - kA), kA = bN / rcond( op(A,t)', p)
        mat_tup_2               norm_error_rev(const Matrix& Y, const Matrix& X, basic_vector_norm p, 
                                        trans_type t = trans_type::no_trans) const;

        /// return componentwise backward error of the solution Y = solve(X, t) with respect
        /// to matrix perturbation E >= 0 and vector perturbation F >= 0  whose elements
        /// provide relative tolerance againts which the components DA and DX
        /// are measured, i.e.
        ///         min{ bC : (op(A,t) + DA) * Y = X + DX, where
        ///                   |DA| <= bC * op(E,t), |DX| <= bN * F }
        /// where |P| < Q means |p_ij| < q_ij for all elements p_ij, q_ij of matrices P, Q
        /// bC is given by
        ///         bC = max_i |op(A,t) * Y - X|_i / ( op(E,t) * abs(Y) + F )_i 
        /// usually E = abs(A) and F = abs(X).
        Matrix                  comp_bacward_error(const Matrix& Y, const Matrix& X, const Matrix& E, 
                                        const Matrix& F, trans_type t = trans_type::no_trans) const;

        /// return componentwise backward error of the solution Y = solve_rev(X, t) with
        /// respect to matrix perturbation E >= 0 and vector perturbation F >= 0 whose 
        /// elements provide relative tolerance againts which the components DA
        /// and DX are measured, i.e.
        ///         min{ bC : Y * (op(A,t) + DA) = X + DX, where
        ///                   |DA| <= bC * op(E,t), |DX| <= bN * F }
        /// where |P| < Q means |p_ij| < q_ij for all elements p_ij, q_ij of matrices P, Q
        /// bC is given by
        ///         bC = max_i |Y * op(A,t) - X|_i / ( abs(Y) * op(E,t) + F )_i 
        /// usually E = abs(A) and F = abs(X).
        Matrix                  comp_bacward_error_rev(const Matrix& Y, const Matrix& X, const Matrix& E, 
                                        const Matrix& F, trans_type t = trans_type::no_trans) const;

        /// return componentwise backward bC and forward fC error of the solution 
        /// Y = solve(X, t) with respect to matrix perturbation E = abs(A) and vector 
        /// perturbation F = abs(X) in the infinity norm; the backward error is calculated
        /// as in comp_bacward_error; the forward error is defined as
        ///     |Y - Y*|_inf / |Y|_inf <= fC
        /// where Y* is the true solution and is given by
        ///     fC_i = | abs(inv(op(A,t))) * (abs(r_i) + g * t_i)) |_inf / |Y_i|_inf
        ///     r_i  = op(A,t) * Y_i - X_i
        ///     t_i  = abs(X_i) + abs(op(A,t)) * abs(Y_i)
        ///     g    = NZ * EPS
        /// where EPS is the machine precision for value type of A, NZ is the maximum number
        /// of nonzeros in any row of A, plus 1 (this is the forward error estimator used in
        /// Lapack); the matrix inv(op(A)) is not constructed explicitly
        /// [bC, fC] = comp_error(X, Y, t)
        mat_tup_2               comp_error(const Matrix& Y, const Matrix& X,
                                        trans_type t = trans_type::no_trans) const;

        /// return componentwise backward bN and forward fN error of the solution 
        /// Y = solve_rev(X, t) with respect to matrix perturbation E = A and vector 
        /// perturbation F = X in the infinity norm; the backward error is calculated
        /// as in comp_bacward_error; the forward error is defined as
        ///     |Y - Y*|_inf / |Y|_inf <= fN
        /// where Y* is the true solution and is given by
        ///     fN_i = | (abs(r_i) + g * t_i) * abs( inv(op(A,t)) ) ) |_inf / |Y_i|_inf
        ///     r    = Y_i * op(A,t) - X_i
        ///     t    = abs(X_i) + abs(Y_i) * abs(op(A,t))
        ///     g    = NZ * EPS
        /// where EPS is the machine precision for value type of A, NZ is the maximum number
        /// of nonzeros in any row of A, plus 1 (this is the forward error estimator used in
        /// Lapack); the matrix inv(op(A)) is not constructed explicitly
        /// [bC, fC] = comp_error_rev(X, Y, t)
        mat_tup_2               comp_error_rev(const Matrix& Y, const Matrix& X,
                                        trans_type t = trans_type::no_trans) const;

        //--------------------------------------------------------------------
        //          const functions defined for all matrix types
        //--------------------------------------------------------------------
        /// get number of rows
        Integer                 rows() const;

        /// get number of columns
        Integer                 cols() const;
        
        /// larger or rows() and cols(); return zero for empty matrices
        Integer                 length() const;
                
        /// number of rows times number of columns
        Real                    numel() const;

        /// check if all elements are finite
        bool                    all_finite() const;

        /// true if this is 0xn or mx0 matrix
        bool                    is_empty() const;

        /// true is this is a scalar or 1x1 matrix
        bool                    is_scalar() const;

        /// true if rows() == cols()
        bool                    is_square() const;

        /// true if rows() == 1 or cols() == 1
        bool                    is_vector() const;

        /// code of stored elements type
        value_code              get_value_code() const;

        /// return type_info of stored elements
        ti::ti_object           get_type() const;
};

//--------------------------------------------------------------------
//              operations defined on linsolve_obj
//--------------------------------------------------------------------
/// makes transpose of given mat
MATCL_LINALG_EXPORT linsolve_obj trans(const linsolve_obj& mat);

/// makes conjugate transpose of given mat
MATCL_LINALG_EXPORT linsolve_obj ctrans(const linsolve_obj& mat);

/// make transposition of mat given of type by t
MATCL_LINALG_EXPORT linsolve_obj trans(const linsolve_obj& mat, trans_type t);

/// make transposition of mat given of type by t
MATCL_LINALG_EXPORT linsolve_obj trans(const linsolve_obj& mat, trans_type_ext t);

/// makes conjugate of given mat
MATCL_LINALG_EXPORT linsolve_obj conj(const linsolve_obj& mat);

/// return linsolve_obj representing -A
MATCL_LINALG_EXPORT linsolve_obj operator-(const linsolve_obj& A);

/// return linsolve_obj representing -A
MATCL_LINALG_EXPORT linsolve_obj uminus(const linsolve_obj& A);

/// return linsolve_obj representing A * B
MATCL_LINALG_EXPORT linsolve_obj operator*(const linsolve_obj& A, const linsolve_obj& B);

/// return linsolve_obj representing A * B
MATCL_LINALG_EXPORT linsolve_obj mmul(const linsolve_obj& A, const linsolve_obj& B,
                                         trans_type tA = trans_type::no_trans,
                                         trans_type tB = trans_type::no_trans);

/// return linsolve_obj representing alpha * A, where alpha is a scalar or 1x1 matrix
/// alpha * A is interpreted as element-by-element multiplication by alpha
/// notice that operator*(alpha,A) and mmul(alpha,A) will interpret alpha as a 1x1 matrix
/// and error will be thrown if A does not represent 1xN matrix
MATCL_LINALG_EXPORT linsolve_obj scale(const Matrix& alpha, const linsolve_obj& A);

template<class S, class Enable = typename md::enable_if_scalar<S,void>::type>
linsolve_obj                     scale(const S& alpha, const linsolve_obj& A) 
                                            { return scale(Matrix(alpha), A); };

//--------------------------------------------------------------------
//              rank-k  updates
//--------------------------------------------------------------------
/// form linsolve_obj representing A + U*V'
/// Q = I + V'*inv(A)*U must be invertible, opts are used to form 
/// linsolve_obj of Q
/// if U, V has size N x 1, then this gives the Sherman-Morrison formula
MATCL_LINALG_EXPORT linsolve_obj update(const linsolve_obj& A, const Matrix& U, const Matrix& V,
                                           const options& opt = options());

/// form linsolve_obj representing A + U*U', 
/// Q = I + U'*inv(A)*U must be invertible, opts are used to form 
/// linsolve_obj of Q
/// if U has size N x 1, then this gives the Sherman-Morrison formula
MATCL_LINALG_EXPORT linsolve_obj update(const linsolve_obj& A, const Matrix& U,
                                           const options& opt = options());

/// form linsolve_obj representing A + U*C*V' for invertible C,
/// Q = inv(C) + V'*inv(A)*U must be invertible, opts are used to form 
/// linsolve_obj of Q
/// this is the Woodbury formula
MATCL_LINALG_EXPORT linsolve_obj update(const linsolve_obj& A, const Matrix& U, 
                                           const linsolve_obj& C, const Matrix& V,
                                           const options& opt = options());

/// form linsolve_obj representing A + U*C*U' for invertible C,
/// Q = inv(C) + U'*inv(A)*U must be invertible, opts are used to form 
/// linsolve_obj of Q
/// this is the Woodbury formula
MATCL_LINALG_EXPORT linsolve_obj update(const linsolve_obj& A, const Matrix& U, 
                                           const linsolve_obj& C,
                                           const options& opt = options());

/// form linsolve_obj representing A + U*C*V' for singular (or non-square) C,
/// Q = I + C * V' * inv(A) * U must be invertible, opts are used to form 
/// linsolve_obj of Q
/// this is the binomial inverse formula
MATCL_LINALG_EXPORT linsolve_obj update_BI(const linsolve_obj& A, const Matrix& U, 
                                           const Matrix& C, const Matrix& V,
                                           const options& opt = options());

/// form linsolve_obj representing A + U*C*U' for singular (or non-square) C,
/// Q = I + C * U' * inv(A) * U must be invertible, opts are used to form 
/// linsolve_obj of Q
/// this is the binomial inverse formula
MATCL_LINALG_EXPORT linsolve_obj update_BI(const linsolve_obj& A, const Matrix& U, 
                                           const Matrix& C, const options& opt = options());

/// form linsolve_obj representing a block matrix
///         [A  B]
///         [C  D]
/// with A, D square and A invertible, the matrix Q = D - C*inv(A)*B must also be
/// invertible; opts are used to form linsolve_obj of Q
MATCL_LINALG_EXPORT linsolve_obj block_invert(const linsolve_obj& A, const Matrix& B, 
                                           const Matrix& C, const Matrix& D,
                                           const options& opt = options());

/// form linsolve_obj representing a block matrix
///         [A  B]
///         [B' D]
/// with A, D square and A invertible, the matrix Q = D - B'*inv(A)*B must also be
/// invertible; opts are used to form linsolve_obj of Q
MATCL_LINALG_EXPORT linsolve_obj block_invert(const linsolve_obj& A, const Matrix& B, 
                                           const Matrix& D, const options& opt = options());

/// form linsolve_obj representing a block matrix
///         [A  B]
///         [C  D]
/// with A, D square and D invertible, the matrix Q = A - B*inv(D)*C must also be
/// invertible; opts are used to form linsolve_obj of Q
MATCL_LINALG_EXPORT linsolve_obj block_invert(const Matrix& A, const Matrix& B, 
                                           const Matrix& C, const linsolve_obj& D,
                                           const options& opt = options());

/// form linsolve_obj representing a block matrix
///         [A  B]
///         [B' D]
/// with A, D square and D invertible, the matrix Q = A - B*inv(D)*B' must also be
/// invertible; opts are used to form linsolve_obj of Q
MATCL_LINALG_EXPORT linsolve_obj block_invert(const Matrix& A, const Matrix& B, 
                                           const linsolve_obj& D,
                                           const options& opt = options());

};

#pragma warning (pop)