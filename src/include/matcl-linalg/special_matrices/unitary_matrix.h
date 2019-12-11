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
#include "matcl-core/matrix/enums.h"
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-matrep/IO/serialization_helper.h"

#pragma warning (push)
#pragma warning (disable: 4251) //needs to have dll-interface to be used by

namespace matcl
{

// internal representation of unitary matrix
class MATCL_LINALG_EXPORT unitary_matrix_data
{
    public:
        using data_ptr  = std::shared_ptr<unitary_matrix_data>;

    public:
        virtual ~unitary_matrix_data(){};

    public:
        // get number of rows
        virtual Integer         rows() const = 0;

        // get number of columns
        virtual Integer         cols() const = 0;

        // check if all elements are finite
        virtual bool            all_finite() const = 0;

        // return true if unitary_matrix was constructed from a matrix
        virtual bool            is_matrix() const   { return false; };

        // return true if unitary_matrix represents an identity matrix
        virtual bool            is_id() const       { return false; };

        // code of stored elements type
        virtual value_code      get_value_code() const = 0;

        // return type_info of stored elements
        virtual ti::ti_object   get_type() const;

        // convert to matcl matrix
        virtual void            to_matrix(matcl::Matrix& ret) const = 0;

        // this function is called when expected value type is different 
        // than value type of stored elements and conversion is expected;
        // given linear operator may ignore this request or perform required
        // conversion, in this case this function must be overriden;
        virtual data_ptr        convert(value_code new_val_code) const = 0;

        // form op(this) * mat; if mat is unique, then can be modified inplace
        virtual void            mult_right(matcl::Matrix& ret, const matcl::Matrix& mat,
                                        trans_type t_unitary) const = 0;

        // form mat * op(this); if mat is unique, then can be modified inplace
        virtual void            mult_left(matcl::Matrix& ret, const matcl::Matrix& mat,
                                        trans_type t_unitary) const = 0;

        // derived class must return 
        //     serialization_helper<unitary_matrix_data>::get<derived>(unique_id);
        // where unique_id is selected string uniquely representing derived class;
        virtual serialization_helper<unitary_matrix_data>*
                                get_serialization_helper() const = 0;

        // serialize unitary_matrix representation using
        // boost::serialization library.
        virtual void            save(oarchive& os) const = 0;

        // output stream operator
        virtual void            save(std::ostream& os) const = 0;

        // derived classes must also implement
        // static unitary_matrix_data* load(std::istream& is);
        // static unitary_matrix_data* load(iarchive& ar);
};

// class representing an unitary matrix, not necessary square, 
// is matrix is not square, then it is submatrix of some square unitary 
// matrix of size K x K obtained by taking first M rows and N columns, 
// where K = max(M,N)
class MATCL_LINALG_EXPORT unitary_matrix
{
    //internal use
    public:
        using impl_type = unitary_matrix_data;
        using impl_ptr  = std::shared_ptr<impl_type>;

    public: 
        // type of internal representation
        using unitary_matrix_data_ptr = std::shared_ptr<unitary_matrix_data>;

    private:
        trans_type_ext  m_trans;
        impl_ptr        m_impl;

    public:
        // create unitary matrix from real scalar 1.0
        unitary_matrix();

        // create unitary matrix from matcl Matrix mat representing
        // a valid unitary matrix; matrix mat should contain only finite
        // values; if test_finite is true, then validity test is performed
        // and if mat fails this test, then special unitary matrix representing
        // matrix with NaN values is created
        explicit unitary_matrix(const Matrix& mat, bool test_finite);
        explicit unitary_matrix(Matrix&& mat, bool test_finite);        

        // create unitary_matrix from representation rep; rep cannot be empty
        explicit unitary_matrix(const unitary_matrix_data_ptr& rep);
        explicit unitary_matrix(unitary_matrix_data_ptr&& rep);

        // standard copy and move constructor
        unitary_matrix(const unitary_matrix& mat);
        unitary_matrix(unitary_matrix&& mat);

        // standard assignment and move assignment operator
        unitary_matrix&             operator=(const unitary_matrix&) &;
        unitary_matrix&             operator=(unitary_matrix&&) &;

        // standard destructor
        ~unitary_matrix();

    public:
        //--------------------------------------------------------------------
        //          member functions specific to unitary_matrix
        //--------------------------------------------------------------------
        
        // convert to matcl matrix
        matcl::Matrix               to_matrix() const;

        // this function is called when expected value type is different 
        // than value type of stored elements and conversion is expected;
        // however given linear operator may ignore this request
        unitary_matrix              convert(value_code new_val_code) const;

        // return unitary matrix representing matrix with NaN values
        static unitary_matrix       from_nan(Integer M, Integer N, value_code vc);

        //--------------------------------------------------------------------
        //          const functions defined for all matrix types
        //--------------------------------------------------------------------
        // get number of rows
        Integer                     rows() const;

        // get number of columns
        Integer                     cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer                     length() const;
                
        // number of rows times number of columns
        Real                        numel() const;

        // check if all elements are finite
        bool                        all_finite() const;

        // true if this is 0xn or mx0 matrix
        bool                        is_empty() const;

        // true is this is a scalar or 1x1 matrix
        bool                        is_scalar() const;

        // true if rows() == cols()
        bool                        is_square() const;

        // true if rows() == 1 or cols() == 1
        bool                        is_vector() const;

        // code of stored elements type
        value_code                  get_value_code() const;

        // return type_info of stored elements
        ti::ti_object               get_type() const;

    //internal use
    public:
        trans_type_ext              get_trans() const           { return m_trans; };
        void                        set_trans(trans_type_ext t) { m_trans = t; };
        const impl_ptr&             get_impl() const            { return m_impl; };
        void                        set_impl(const impl_ptr& p) { m_impl = p; };
};

//--------------------------------------------------------------------
//              operations defined on unitary_matrix
//--------------------------------------------------------------------
// makes transpose of given mat
MATCL_LINALG_EXPORT unitary_matrix  trans(const unitary_matrix& mat);

// makes conjugate transpose of given mat
MATCL_LINALG_EXPORT unitary_matrix  ctrans(const unitary_matrix& mat);

// makes conjugate of given mat
MATCL_LINALG_EXPORT unitary_matrix  conj(const unitary_matrix& mat);

// make transposition of mat given of type by t
MATCL_LINALG_EXPORT unitary_matrix  trans(const unitary_matrix& mat, trans_type t);

// make transposition of mat given of type by t
MATCL_LINALG_EXPORT unitary_matrix  trans(const unitary_matrix& mat, trans_type_ext t);

// matrix multiply; not implemented if B is sparse or band matrix
MATCL_LINALG_EXPORT matcl::Matrix   operator*(const unitary_matrix& A, const matcl::Matrix& B);
MATCL_LINALG_EXPORT matcl::Matrix   operator*(const unitary_matrix& A, matcl::Matrix&& B);
MATCL_LINALG_EXPORT matcl::Matrix   operator*(const matcl::Matrix& A, const unitary_matrix& B);
MATCL_LINALG_EXPORT matcl::Matrix   operator*(matcl::Matrix&& A, const unitary_matrix& B);
MATCL_LINALG_EXPORT unitary_matrix  operator*(const unitary_matrix& A, const unitary_matrix& B);

// matrix multiply, op(A) * B, where 
//     op(A) = A           if t_A = trans_type::no_trans
//     op(A) = trans(A)    if t_A = trans_type::trans
//     op(A) = ctrans(A)   if t_A = trans_type::conj_trans
//
// not implemented for sparse and band matrix or t_B != no_trans
MATCL_LINALG_EXPORT matcl::Matrix   mmul(const unitary_matrix& A, const matcl::Matrix& B, 
                                        trans_type tA = trans_type::no_trans,
                                         trans_type tB = trans_type::no_trans);
MATCL_LINALG_EXPORT matcl::Matrix   mmul(const unitary_matrix& A, matcl::Matrix&& B, 
                                        trans_type tA = trans_type::no_trans,
                                        trans_type tB = trans_type::no_trans);

// matrix multiply, A * op(B), where 
//     op(B) = B           if t_B = trans_type::no_trans
//     op(B) = trans(B)    if t_B = trans_type::trans
//     op(B) = ctrans(B)   if t_B = trans_type::conj_trans
//
// not implemented for sparse and band matrix or t_A != no_trans
MATCL_LINALG_EXPORT matcl::Matrix   mmul(const matcl::Matrix& A, const unitary_matrix& B, 
                                         trans_type tA = trans_type::no_trans,
                                         trans_type tB = trans_type::no_trans);
MATCL_LINALG_EXPORT matcl::Matrix   mmul(matcl::Matrix&& A, const unitary_matrix& B, 
                                         trans_type tA = trans_type::no_trans,
                                         trans_type tB = trans_type::no_trans);

// multiply two unitary matrices
MATCL_LINALG_EXPORT unitary_matrix  mmul(const unitary_matrix& A, const unitary_matrix& B, 
                                         trans_type tA = trans_type::no_trans,
                                         trans_type tB = trans_type::no_trans);

// serialize unitary_matrix using boost::serialization library.
MATCL_LINALG_EXPORT void            save(oarchive& ar,const unitary_matrix& mat);

// deserialize unitary_matrix using boost::serialization library.
MATCL_LINALG_EXPORT void            load(iarchive& ar,unitary_matrix& mat);

// output stream operator for an unitary_matrix
MATCL_LINALG_EXPORT std::ostream&   operator<<(std::ostream& os, const unitary_matrix& mat);
    
// input stream operator for an unitary_matrix
MATCL_LINALG_EXPORT std::istream&   operator>>(std::istream& is, unitary_matrix& mat);

// generate random dense unitary matrix of size K x K with real values
MATCL_LINALG_EXPORT unitary_matrix  rand_unitary(Integer K);

// generate random dense unitary matrix of size K x K with float values
MATCL_LINALG_EXPORT unitary_matrix  frand_unitary(Integer K);

// generate random dense unitary matrix of size K x v with complex values
MATCL_LINALG_EXPORT unitary_matrix  crand_unitary(Integer K);

// generate random dense unitary matrix of size K x K with float complex values
MATCL_LINALG_EXPORT unitary_matrix  fcrand_unitary(Integer K);

// generate random dense unitary matrix of size K x K with real values 
// given by value code vc
MATCL_LINALG_EXPORT unitary_matrix  rand_unitary(Integer K, value_code vc);
};

#pragma warning (pop)