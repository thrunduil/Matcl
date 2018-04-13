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
#include "matcl-linalg/special_matrices/unitary_matrix.h"
#include "matcl-linalg/details/linalg_fwd.h"

#pragma warning(push)
#pragma warning(disable:4251)   //needs to have dll-interface to be used by clients of class

namespace matcl
{

/// abstract class representing a matrix like object A with defined multiplication 
/// operator; derived classes must implement matrix multiplications:
/// op(A) * X, where op(A) = A or op(A) = ctrans(A) and X is a vector
class MATCL_LINALG_EXPORT linear_operator_data 
            : private std::enable_shared_from_this<linear_operator_data>
{
    public:
        /// shared pointer type
        using data_ptr  = std::shared_ptr<linear_operator_data>;

    public:
        virtual ~linear_operator_data() {};

        /// value code of elements stored in matrix A
        virtual value_code  get_value_code() const = 0;

        /// this function is called when expected value type is different 
        /// than value type of stored elements and conversion is expected;
        /// given linear operator may ignore this request or perform required
        /// conversion, in this case this function must be overriden;
        virtual data_ptr    convert(value_code new_val_code) const;

        /// number of rows of the matrix A
        virtual Integer     rows() const = 0;

        /// number of columns of the matrix A
        virtual Integer     cols() const = 0;

        /// should return true if the matrix is real symmetric or complex hermitian
        virtual bool        is_hermitian() const
                                { return false; };

        /// evaluate mmul(A, X, t), i.e. trans(A,t) * X
        virtual Matrix      mmul_right(const Matrix& X, trans_type t) const = 0;

        /// derived class can implement this version of mmul_right if inplace
        /// multiplication is possible
        virtual Matrix      mmul_right(Matrix&& X, trans_type t) const 
                                { return mmul_right(X, t); };

        /// evaluate y(:) = mmul(A, X, t), i.e. y(:) = trans(A,t) * X; results
        /// must be stored in array held by a dense unique matrix y
        virtual void        mmul_right(const Matrix& X, trans_type t, Matrix& y) const;
};

/// shared pointer like class representing a linear operator A, i.e. matrix-like
/// object with defined multiplication by a matrix; in order to define linear
/// operator one must derive from linear_operator_data class
///
/// scalars are always treated as 1x1 matrix, i.e. A * x, where A is created from 
/// a scalar and x is not a 1xN matrix, will produced error
class MATCL_LINALG_EXPORT linear_operator
{
    public: 
        /// type of internal representation
        using data_ptr  = std::shared_ptr<linear_operator_data>;

    private:
        data_ptr        m_impl;
        trans_type_ext  m_trans;

    public:
        /// create linear_operator representing a real scalar 0.0
        linear_operator();

        /// create linear_operator from a type convertible to matrix
        template<class T, class Enable = typename details::enable_convertible_to_matrix<T,void>::type>
        linear_operator(const T& mat)               { from_matrix(Matrix(mat));};

        /// create linear_operator from a Matrix
        linear_operator(const Matrix& mat)          { from_matrix(Matrix(mat));};

        /// create linear_operator from a unitary matrix
        linear_operator(const unitary_matrix& mat)  { from_unitary(mat);};

        /// create linear_operator from a linsolve object
        /// A * X is defined as mat.solve(X)
        linear_operator(const linsolve_obj& mat)    { from_linsolve(mat);};

        /// create linear_operator from representation ret; rep cannot be empty
        explicit linear_operator(const data_ptr& rep, trans_type_ext t = trans_type_ext::no_trans);
        explicit linear_operator(data_ptr&& rep, trans_type_ext t = trans_type_ext::no_trans);

        /// standard copy and move constructor
        linear_operator(const linear_operator& mat);
        linear_operator(linear_operator&& mat);

        /// standard assignment and move assignment operator
        linear_operator&    operator=(const linear_operator&) &;
        linear_operator&    operator=(linear_operator&&) &;

        /// standard destructor
        ~linear_operator();

    public:
        /// value code of elements stored in matrix A
        value_code          get_value_code() const;

        /// this function is called when expected value type is different 
        /// than value type of stored elements and conversion is expected;
        /// however given linear operator may ignore this request
        linear_operator     convert(value_code new_val_code) const;

        /// number of rows of the matrix A
        Integer             rows() const;

        /// number of columns of the matrix A
        Integer             cols() const;

        /// return true if the matrix is marked as real symmetric or complex hermitian
        bool                is_hermitian() const;

        /// evaluate mmul(A, X, t), i.e. trans(A,t) * X
        Matrix              mmul_right(const Matrix& X, trans_type t) const;
        Matrix              mmul_right(Matrix&& X, trans_type t) const;

        /// evaluate y(:) = mmul(A, X, t), i.e. y(:) = trans(A,t) * X; results
        /// must be stored in array held by a dense unique matrix y
        void                mmul_right(const Matrix& X, trans_type t, Matrix& y) const;

    private:
        void                from_matrix(const Matrix& mat);
        void                from_unitary(const unitary_matrix& mat);
        void                from_linsolve(const linsolve_obj& mat);

    //internal use
    public:
        trans_type_ext      get_trans() const           { return m_trans; };
        void                set_trans(trans_type_ext t) { m_trans = t; };
};

//--------------------------------------------------------------------
//              operations defined on linear_operator
//--------------------------------------------------------------------
/// makes transpose of given mat
MATCL_LINALG_EXPORT linear_operator linop_trans(const linear_operator& mat);

/// makes conjugate transpose of given mat
MATCL_LINALG_EXPORT linear_operator linop_ctrans(const linear_operator& mat);

/// make transposition of mat given of type by t
MATCL_LINALG_EXPORT linear_operator linop_trans(const linear_operator& mat, trans_type t);

/// makes conjugate of given mat
MATCL_LINALG_EXPORT linear_operator linop_conj(const linear_operator& mat);

/// return linear operator representing A + B
MATCL_LINALG_EXPORT linear_operator linop_plus(const linear_operator& A, const linear_operator& B);

/// return linear operator representing A - B
MATCL_LINALG_EXPORT linear_operator linop_minus(const linear_operator& A, const linear_operator& B);

/// return linear operator representing -A
MATCL_LINALG_EXPORT linear_operator linop_uminus(const linear_operator& A);

/// return linear operator representing A * B
MATCL_LINALG_EXPORT linear_operator linop_mmul(const linear_operator& A, const linear_operator& B,
                                         trans_type tA = trans_type::no_trans,
                                         trans_type tB = trans_type::no_trans);

/// return linear operator representing alpha * A, where alpha is a scalar or 1x1 matrix
/// alpha * A is interpreted as element-by-element multiplication by alpha
/// notice that operator*(alpha,A) and mmul(alpha,A) will interpret alpha as a 1x1 matrix
/// and error will be thrown if A does not represent 1xN matrix
MATCL_LINALG_EXPORT linear_operator linop_scale(const Matrix& alpha, const linear_operator& A);

template<class S, class Enable = typename md::enable_if_scalar<S,void>::type>
linear_operator                     linop_scale(const S& alpha, const linear_operator& A) 
                                            { return scale(Matrix(alpha), A); };

/// return linear operator representing A * A' (trans = false) or A' * A (trans = true)
MATCL_LINALG_EXPORT linear_operator linop_symprod(const linear_operator& A, bool trans = true);

/// return linear operator representing A + A'
MATCL_LINALG_EXPORT linear_operator linop_symsum(const linear_operator& A);

/// return linear operator representing [0 A'; A 0] (if trans = false) or [0 A; A' 0]
/// (if trans = true)
MATCL_LINALG_EXPORT linear_operator linop_symcat(const linear_operator& A, bool trans = true);

};

#pragma warning(pop)