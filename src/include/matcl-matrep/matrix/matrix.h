/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-matrep/general/config.h"
#include "matcl-matrep/details/matrix_details.h"
#include "matcl-matrep/details/isa.h"
#include "matcl-matrep/matrix/struct_flag.h"
#include "matcl-matrep/general/tuple.h"
#include "matcl-matrep/details/enablers.h"

namespace matcl
{

// predefined tuples
using mat_tup_2 = tuple<Matrix,Matrix>;
using mat_tup_3 = tuple<Matrix,Matrix,Matrix>;
using mat_tup_4 = tuple<Matrix,Matrix,Matrix,Matrix>;
using mat_tup_5 = tuple<Matrix,Matrix,Matrix,Matrix,Matrix>;
using mat_tup_6 = tuple<Matrix,Matrix,Matrix,Matrix,Matrix,Matrix>;
using int_tup_2 = tuple<Integer,Integer>;

// Matrix is a container storing scalars or dense, sparse, or band matrices of
// Integer, Float, Real, Complex, Float_complex, or Object numbers. Matrix 
// representation can change, when nonconst methods are called. Matrix bahaves 
// according to copy-on-write scheme. Matrix is thread-safe, but given instance 
// cannot be shared between threads.
//
// MATCL preserves structure of sparse matrices, i.e. zeroes are not implicitly 
// dropped except binary operations (all functions declared in func_binary, kron,
// sh_prod); all unary, manip, concat, assign functions preserve sparsity (for 
// example assigning zero to sparse matrix can add new element to the structure); 
// conversion of band matrix, assignment band matrix to sparse matrix or concat
// with band matrix will add all structurally nonzero elements of band matrix to
// the structure of sparse matrix; in conversion, assign, concat of dense matrix
// test for zeroes are always performed; getting submatrix of band and sparse matrices
// also preserves structure; notice however that sparse matrix can be converted to
// dense matrix if is too dense
class MATCL_MATREP_EXPORT Matrix : private details::matrix_base
{
    private:
        using mat_ptr   = details::matrix_container_base;
        using base_type = details::matrix_base;

        const mat_ptr*  get_mat_ptr_prv() const;
        
        friend matcl::sub_matrix;
        friend matcl::sub_matrix_1;
        friend matcl::sub_matrix_2;
        friend details::matrix_data_accesser;        

    public:
        //--------------------------------------------------------------------
        //      constructors, destructor, assignments
        //--------------------------------------------------------------------

        // default constructor create Real zero scalar
        Matrix();

        // conversion from bool to Matrix is explicit
        explicit Matrix(bool val);

        // conversion from scalar types
        template<class T>            
        Matrix(T&& val, typename details::enable_if_is_conv_to_mat<T, true>::type = 0);        
        
        // explicit conversion from typed object scalars
        template<class T>
        explicit Matrix(const object_type<T>& obj);

        template<class T>
        explicit Matrix(object_type<T>&& obj);

        // conversion from internal types convertible to Matrix
        template<class T>            
        explicit Matrix(T&& val,bool allow_conversions,
                typename details::enable_if_is_conv_to_mat<T, false>::type = 0);

        // conversion from dense matrix
        template<class T, bool Is_safe>
        explicit Matrix(const dense_matrix<T, Is_safe>& mat);

        template<class T, bool Is_safe>
        explicit Matrix(dense_matrix<T, Is_safe>&& mat);

        // conversion from sparse matrix
        template<class T, bool Is_const>
        explicit Matrix(const sparse_matrix<T, Is_const>& mat);
        
        template<class T, bool Is_const>
        explicit Matrix(sparse_matrix<T, Is_const>&& mat);

        // conversion from band matrix
        template<class T, bool Is_const>
        explicit Matrix(const band_matrix<T, Is_const>& mat);

        template<class T, bool Is_const>
        explicit Matrix(band_matrix<T, Is_const>&& mat);

        // conversion from object matrix
        template<class T>
        explicit Matrix(const object_matrix<T>& mat);
        
        template<class T>
        explicit Matrix(object_matrix<T>&& mat);

        // copy constructor, data are not copied, only refcount is increased
        Matrix(const Matrix&);
        
        // move constructor
        Matrix(Matrix&&);
        
        // conversion from submatrix, a new matrix is usually created
        Matrix(const matcl::sub_matrix&);
        Matrix(const matcl::sub_matrix_1&);
        Matrix(const matcl::sub_matrix_2&);

        // conversion from unique matrix
        Matrix(unique_matrix&&);

        // conversion from mat_row (i.e. concatenation helper)
        Matrix(const mat_row&);
        
        // conversion from mat_col (i.e. concatenation helper)
        Matrix(const mat_col&);

        // conversions from typed concatenation helpers and submatrices
        template<class T>
        explicit Matrix(const sub_dense_matrix<T>&);
        
        template<class T>
        explicit Matrix(const dense_row<T>&);
        
        template<class T>
        explicit Matrix(const dense_col<T>&);

        // conversions from typed concatenation helpers and submatrices
        template<class T>
        explicit Matrix(const sub_sparse_matrix<T>&);
        
        template<class T>
        explicit Matrix(const sub_sparse_matrix_1<T>&);
        
        template<class T>
        explicit Matrix(const sub_sparse_matrix_2<T>&);
        
        template<class T>
        explicit Matrix(const sparse_row<T>&);
        
        template<class T>
        explicit Matrix(const sparse_col<T>&);

        // conversions from typed concatenation helpers and submatrices
        template<class T>
        explicit Matrix(const sub_band_matrix<T>&);
        
        template<class T>
        explicit Matrix(const sub_band_matrix_1<T>&);
        
        template<class T>
        explicit Matrix(const sub_band_matrix_2<T>&);

        // conversions from typed concatenation helpers and submatrices
        template<class T>
        explicit Matrix(const object_row<T>&);
        
        template<class T>
        explicit Matrix(const object_col<T>&);
        
        template<class T>
        explicit Matrix(const sub_object_matrix<T>&);
        
        template<class T>
        explicit Matrix(const sub_object_matrix_1<T>&);        
        
        template<class T>
        explicit Matrix(const sub_object_matrix_2<T>&);

        // destructor, memory is freed if refcount drops to zero
        ~Matrix();

        // assignments; pointer to representation is changed to pointer
        // in rhs
        Matrix&                     operator=(const Matrix&) &;
        Matrix&                     operator=(Matrix&&) &;

        // assignments from submatrices; submatrices are converted to matrices
        // and then standard assignment take place
        Matrix&                     operator=(const matcl::sub_matrix&) &;
        Matrix&                     operator=(const matcl::sub_matrix_1&) &;
        Matrix&                     operator=(const matcl::sub_matrix_2&) &;

        //--------------------------------------------------------------------
        //      functions specific to Matrix
        //--------------------------------------------------------------------
        // get scalar of given type; matrix must store scalars or 1x1 matrices;
        // conversions are allowed
        template<class V, class Enable = typename details::enable_if_matcl_scalar<V,void>::type>
        V                           get_scalar() const;

        // get scalar of given type; matrix must store scalars or 1x1 matrices;
        // conversions are allowed; make_unique is called
        template<class V, class Enable = typename details::enable_if_matcl_scalar<V,void>::type>
        V&                          get_scalar_unique();

        // get array of given type; matrix is converted to dense explicit matrix
        // (i.e. if matrix is a noncontinuous view, then new matrix is created);
        // conversions are allowed
        template<class T, class Enable = typename details::enable_if_matcl_scalar<T,void>::type>
        const T*                    get_array() const;

        // get array of given type; matrix is converted to dense explicit and unique matrix
        // (i.e. if matrix is a noncontinuous view, then new matrix is created);
        // conversions are allowed; all assigned structures are removed
        template<class T, class Enable = typename details::enable_if_matcl_scalar<T,void>::type>
        T*                          get_array_unique();

        //--------------------------------------------------------------------
        //          const functions defined for all matrix types
        //--------------------------------------------------------------------
        // get submatrix, a new matrix is usually created, unless matrix views 
        // can be created
        const Matrix                operator()(Integer r) const;
        const Matrix                operator()(Integer r, Integer c) const;
        const Matrix                operator()(const colon& r) const;        
        const Matrix                operator()(const colon& r, const colon& c) const;
        
        // get diagonal
        const Matrix                diag(Integer d = 0) const;

        // delete rows specified by colon c.
        const Matrix                delrows(const colon&) const &;
        const Matrix                delrows(const colon&) const &&;

        // delete columns specified in colon c.
        const Matrix                delcols(const colon&) const &;
        const Matrix                delcols(const colon&) const &&;

        // delete rows specified by colon c1 and columns specified by colon c2.
        const Matrix                delrowscols(const colon& c1, const colon& c2) const &;
        const Matrix                delrowscols(const colon& c1, const colon& c2) const &&;

        // get number of rows; return 1 for scalars
        Integer                     rows() const;

        // get number of columns; return 1 for scalars
        Integer                     cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer                     length() const;
        
        // estimate of number of nonzeroes elements based on representation type
        // and struct flags (without testing of stored values)
        Integer                     structural_nnz() const;
        
        // estimate of number of subdiagonal based on representation type
        // and struct flags if use_flags = true, but without testing of stored values
        Integer                     structural_ldiags(bool use_flags = true) const;
        
        // estimate of number of superdiagonals based on representation type
        // and struct flags if use_flags = true, but without testing of stored values
        Integer                     structural_udiags(bool use_flags = true) const;
        
        // number of rows times number of columns
        Real                        numel() const;

        // check if all elements are finite
        bool                        all_finite() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit                    operator bool() const;

        // true if this is 0xn or mx0 matrix
        bool                        is_empty() const;

        // true is this is a scalar or 1x1 matrix
        bool                        is_scalar() const;

        // true if rows() == cols()
        bool                        is_square() const;

        // true if rows() == 1 or cols() == 1 or is_empty() == true
        bool                        is_vector() const;

        // true if given matrix is represented as dense, sparse or band matrix
        bool                        is_matrix_type() const;

        // true if given matrix is represented as scalar
        bool                        is_scalar_type() const;

        // true if reference count is 1
        bool                        is_unique() const;

        // code of stored elements type
        value_code                  get_value_code() const;

        // code if matrix representation (scalar, band, dense, or sparse)
        struct_code                 get_struct_code() const;

        // matrix code (value code and struct code)
        mat_code                    get_matrix_code() const     { return m_type;}

        // return type_info of stored elements
        ti::ti_object               get_type() const;

        // return additional structures assigned to this matrix
        struct_flag&                get_struct() const;        

        // struct flag must be valid, data are not changed
        void                        set_struct(const struct_flag&) const;
        
        // struct flag must be valid, data are not changed
        void                        add_struct(const struct_flag&) const;        

        // if matrix is unique (i.e. refcount is 1) then do nothing; otherwise
        // make copy of this matrix, note that refcount of this instance will
        // also drop to one; TODO
        const Matrix&               make_unique() const;

        // make independent copy
        const Matrix                clone() const;        

        //--------------------------------------------------------------------
        //          non const functions defined for all matrix types
        //--------------------------------------------------------------------
        // get submatrix, change submatrix, or modify structure of sparse
        // and band matrices
        matcl::sub_matrix_1         operator()(Integer r);
        matcl::sub_matrix_2         operator()(Integer r, Integer c);
        matcl::sub_matrix           operator()(const colon& r);
        matcl::sub_matrix           operator()(const colon& r, const colon& c);

        // get or change diagonal
        matcl::sub_matrix           diag(Integer d = 0);

        // change matrix size; if number of rows and columns is less than capacity,
        // then no copy is done
        void                        resize(Integer r, Integer c);

        // increase rows and columns capacity
        void                        reserve(Integer r, Integer c);

        // change matrix size including first and last diagonal index; 
        // if stored matrix is not a band matrix, then equivalent to resize(r,c);
        void                        resize_band(Integer r, Integer c, Integer fd, Integer ld);

        // increase rows columns and first and last diagonal capacity; 
        // if stored matrix is not a band matrix, then equivalent to reserve(r,c) 
        void                        reserve_band(Integer r, Integer c, Integer fd, Integer ld);                

    //internal use
    public:
        // internal representation of this matrix must be M, matrix must 
        // be temporary and unique
        template<class M> 
        details::rvalue_holder<M>   move_impl();
        
        // internal representation of this matrix must be M
        template<class M> const M&  get_impl() const;
        
        // internal representation of this matrix must be M,
        // makes matrix unique; structure flags are not cleared
        template<class M> M&        get_impl_unique() &;
                
        // mark or unmark the matrix as being effectively unique
        void                        mark_unique(bool unique);
};

};
