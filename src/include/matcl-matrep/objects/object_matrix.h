/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/matrix/struct_flag.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-matrep/matrix/colon.h"
#include "matcl-scalar/matcl_scalar.h"
#include "matcl-matrep/objects/details/object_matrix_subs.h"
#include "matcl-matrep/matrix/matrix_concat.h"
#include "matcl-matrep/lib_functions/eval_functors.h"
#include "matcl-core/options/matcl_options.h"
#include "matcl-matrep/matrix/matrix_rep_dense.h"
#include "matcl-dynamic/matcl_dynamic.h"

namespace matcl
{

// matrix horizontal concatenation helper, i.e. constructor of [A, B, ...]
// instances of this class cannot be shared between threads
template<class T>
class object_row : public mat_row
{
    private:
        using base_type = mat_row;

    public:
        // construct empty matrix of Integer type
        object_row();        

        // standard destructor
        ~object_row();        

        // standard copy constructor and move constructor
        object_row(const object_row&);
        object_row(object_row&&);

        // standard assignment and move assignment;
        object_row&             operator=(const object_row&);
        object_row&             operator=(object_row&&);        

        // add scalar to the list of matrices
        object_row&             add(const T&);
        object_row&             add(T&&);
        object_row&             add(const object_type<T>&);
        object_row&             add(object_type<T>&&);
        object_row&             operator,(const T& v)               { return add(v); };
        object_row&             operator,(T&& v)                    { return add(std::move(v)); };
        object_row&             operator,(const object_type<T>& v)  { return add(v); };
        object_row&             operator,(object_type<T>&& v)       { return add(std::move(v)); };

        // add matrix to the list of matrices
        object_row&             add(const object_matrix<T>&);
        object_row&             add(object_matrix<T>&&);
        object_row&             operator,(const object_matrix<T>& v){ return add(v); };
        object_row&             operator,(object_matrix<T>&& v)     { return add(std::move(v)); };

        // add concatenation helpers to the list of matrices
        object_row&             add(const object_row&);
        object_row&             add(object_row&&);
        object_row&             add(const object_col<T>&);
        object_row&             add(object_col<T>&&);
        object_row&             operator,(const object_row& v)      { return add(v); };
        object_row&             operator,(object_row&& v)           { return add(std::move(v)); };
        object_row&             operator,(const object_col<T>& v)   { return add(v); };
        object_row&             operator,(object_col<T>&& v)        { return add(std::move(v)); };

        // convert to dense_matrix
        const object_matrix<T>  to_matrix() const;
};

// matrix vertical concatenation helper, i.e. constructor of [A; B; ...]
// instances of this class cannot be shared between threads
template<class T>
class object_col : public mat_col
{
    private:
        using base_type = mat_col;

    public:
        // construct empty matrix of Integer type
        object_col();        

        // standard copy constructor and move constructor
        object_col(const object_col&);
        object_col(object_col&&);

        // standard assignment and move assignment;
        object_col&             operator=(const object_col&);
        object_col&             operator=(object_col&&);

        // standard destructor
        ~object_col();

        // add scalar to the list of matrices
        object_col&             add(const T&);
        object_col&             add(T&&);
        object_col&             add(const object_type<T>&);
        object_col&             add(object_type<T>&&);
        object_col&             operator,(const T& v)               { return add(v); };
        object_col&             operator,(T&& v)                    { return add(std::move(v)); };
        object_col&             operator,(const object_type<T>& v)  { return add(v); };
        object_col&             operator,(object_type<T>&& v)       { return add(std::move(v)); };

        // add matrix to the list of matrices
        object_col&             add(const object_matrix<T>&);
        object_col&             add(object_matrix<T>&&);
        object_col&             operator,(const object_matrix<T>& v){ return add(v); };
        object_col&             operator,(object_matrix<T>&& v)     { return add(std::move(v)); };

        // add concatenation helpers to the list of matrices
        object_col&             add(const object_row<T>&);
        object_col&             add(object_row<T>&&);
        object_col&             add(const object_col&);
        object_col&             add(object_col&&);
        object_col&             operator,(const object_row<T>& v)   { return add(v); };
        object_col&             operator,(object_row<T>&& v)        { return add(std::move(v)); };
        object_col&             operator,(const object_col& v)      { return add(v); };
        object_col&             operator,(object_col&& v)           { return add(std::move(v)); };

        // convert to dense_matrix
        const object_matrix<T>  to_matrix() const;
};

// object_matrix is a container which behaves as general matrix but stores only numbers
// of type object_type<T>. Matrix representation can change, when nonconst methods are
// called. One may need to define properties of the type T by partially specializing
// template matcl::dynamic::object_type_traits.
// Matrix is thread-safe, but given instance cannot be shared between threads unless 
// all accesses are read only.
template<class T>
class object_matrix
{
    private:
        using sub_object_matrix     = sub_object_matrix<T>;
        using sub_object_matrix_1   = sub_object_matrix_1<T>;
        using sub_object_matrix_2   = sub_object_matrix_2<T>;
        using object_type           = object_type<T>;

    public:
        using value_type            = T;

    private:
        matcl::Matrix       m_matrix;

        friend matcl::sub_object_matrix<T>;
        friend matcl::sub_object_matrix_1<T>;
        friend matcl::sub_object_matrix_2<T>;

    public:    
        //--------------------------------------------------------------------
        //      constructors, destructor, assignments
        //--------------------------------------------------------------------

        // default constructor create 1x1 matrix with zero scalar
        object_matrix();

        // conversion from scalar to object_matrix is explicit
        explicit object_matrix(const T& val);

        // conversion from typed object
        object_matrix(const object_type&);
        object_matrix(object_type&&);

        // convert general matrix to object_matrix; if mat is not dense, then
        // conversions take place
        explicit object_matrix(const matcl::Matrix& mat);
        explicit object_matrix(matcl::Matrix&& mat);

        // copy constructor and move constructor
        object_matrix(const object_matrix& mat);
        object_matrix(object_matrix&& mat);

        // conversion from submatrix, a new matrix is usually created
        object_matrix(const sub_object_matrix&);
        object_matrix(const sub_object_matrix_1&);
        object_matrix(const sub_object_matrix_2&);        

        // conversion from object_row (i.e. concatenation helper)
        object_matrix(const object_row<T>&);

        // conversion from object_col (i.e. concatenation helper)
        object_matrix(const object_col<T>&);

        // assignments 
        object_matrix&      operator=(const object_matrix&) &;
        object_matrix&      operator=(object_matrix&&) &;

        // destructor, memory is freed if refcount drops to zero
        ~object_matrix();

        //--------------------------------------------------------------------
        //      functions specific to object_matrix
        //--------------------------------------------------------------------

        // get scalar of given type; matrix must store scalars or 1x1 matrices
        object_type         get_scalar() const;

        // get scalar of given type; matrix must store scalars or 1x1 matrices;
        // make_unique is called
        object_type&        get_scalar_unique();

        // get array of given type; matrix is converted to dense explicit matrix
        // (i.e. if matrix is a noncontinuous view, then new matrix is created);
        const object_type*  get_array() const;

        // get array of given type; matrix is converted to dense explicit and unique matrix
        // (i.e. if matrix is a noncontinuous view, then new matrix is created);
        object_type*        get_array_unique();

        // convert to general matrix
        const Matrix&       to_matrix() const &;
        Matrix&&            to_matrix() &&;

        //--------------------------------------------------------------------
        //          const functions defined for all matrix types
        //--------------------------------------------------------------------
        
        // number of rows
        Integer             rows() const;
        
        // number of columns
        Integer             cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer             length() const;
        
        // estimate of number of nonzeroes elements based on representation type
        // and struct flags (without testing of stored values)
        Integer             structural_nnz() const;
        
        // estimate of number of subdiagonal based on representation type
        // and struct flags if use_flags = true, but without testing of stored values
        Integer             structural_ldiags(bool use_flags = true) const;
        
        // estimate of number of superdiagonals based on representation type
        // and struct flags if use_flags = true, but without testing of stored values
        Integer             structural_udiags(bool use_flags = true) const;
        
        // number of rows times number of columns
        Real                numel() const;        

        // check if all elements are finite
        bool                all_finite() const;

        // return additional structures assigned to this matrix
        const struct_flag   get_struct() const;
        
        // struct flag must be valid, data are not changed
        void                set_struct(const struct_flag&) const;
        
        // struct flag must be valid, data are not changed
        void                add_struct(const struct_flag&) const;

        // return type_info of stored elements
        ti::ti_object       get_type() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit            operator bool() const;

        // true if this is 0xn or mx0 matrix
        bool                is_empty() const;

        // true is this is a scalar or 1x1 matrix
        bool                is_scalar() const;

        // true if rows() == cols()
        bool                is_square() const;

        // true if rows() == 1 or cols() == 1 or is_empty() == true
        bool                is_vector() const;

        // true if given matrix is represented as dense, sparse or band matrix
        bool                is_matrix_type() const;

        // true if given matrix is represented as scalar
        bool                is_scalar_type() const;

        // true if reference count is 1
        bool                is_unique() const;

        // code of stored elements type
        value_code          get_value_code() const;

        // code if matrix representation (scalar, band, dense, or sparse)
        struct_code         get_struct_code() const;

        // matrix code (value code and struct code)
        mat_code            get_matrix_code() const;

        // delete rows specified by colon c.
        const object_matrix delrows(const colon&) const &;
        const object_matrix delrows(const colon&) const &&;

        // delete columns specified in colon c.
        const object_matrix delcols(const colon&) const &;
        const object_matrix delcols(const colon&) const &&;

        // delete rows specified by colon c1 and columns specified by colon c2.
        const object_matrix delrowscols(const colon& c1, const colon& c2) const &;
        const object_matrix delrowscols(const colon& c1, const colon& c2) const &&;

        // get element at position pos; checks are performed
        object_type         operator()(Integer pos) const;
        
        // get element at given row r and column c; checks are performed
        object_type         operator()(Integer r, Integer c) const;
        
        // get submatrix, a new matrix is usually created, unless matrix views 
        // can be created
        const object_matrix operator()(const colon& r) const;        
        const object_matrix operator()(const colon& r, const colon& c) const;
        
        // get diagonal
        const object_matrix diag(Integer d = 0) const;

        // if matrix is unique (i.e. refcount is 1) then do nothing; otherwise
        // make copy of this matrix, note that refcount of this instance will
        // also drop to one
        const object_matrix& make_unique() const;

        // make independent copy
        const object_matrix clone() const;

        //--------------------------------------------------------------------
        //          non const functions defined for all matrix types
        //--------------------------------------------------------------------

        // get or change submatrix
        sub_object_matrix_1 operator()(Integer r);
        sub_object_matrix_2 operator()(Integer r, Integer c);
        sub_object_matrix   operator()(const colon& r);
        sub_object_matrix   operator()(const colon& r, const colon& c);

        // get or change diagonal
        sub_object_matrix   diag(Integer d = 0);

        // change matrix size; if number of rows and columns is less than capacity,
        // then no copy is done
        void                resize(Integer r, Integer c);

        // increase rows and columns capacity
        void                reserve(Integer r, Integer c);

        // change matrix size including first and last diagonal index; 
        // if stored matrix is not a band matrix, then equivalent to resize(r,c);
        void                resize_band(Integer r, Integer c, Integer fd, Integer ld);

        // increase rows columns and first and last diagonal capacity; 
        // if stored matrix is not a band matrix, then equivalent to reserve(r,c) 
        void                reserve_band(Integer r, Integer c, Integer fd, Integer ld);

        //--------------------------------------------------------------------
        //          static functions
        //--------------------------------------------------------------------
        // create dense matrix with zeroes of size rxc
        static object_matrix    zeros(Integer r, Integer c);

        // create sparse matrix with zeroes of size rxc, with allocated space for nnz 
        // elements
        static object_matrix    spzeros(Integer r, Integer c,Integer nnz = 0);

        // create band matrix with Real zeroes of size rxc, with first diag fd and last diag ld
        static object_matrix    bzeros(Integer r, Integer c,Integer fd, Integer ld);

        // create dense matrix with ones of size rxc
        static object_matrix    ones(Integer r, Integer c);

        // create sparse matrix with ones of size rxc
        static object_matrix    spones(Integer r, Integer c);

        // create band matrix with Real ones of size rxc
        static object_matrix    bones(Integer r, Integer c);

        // create dense matrix with ones on diagonal of size rxc
        static object_matrix    eye(Integer r, Integer c);

        // create dense matrix with ones on diagonal of size rxr
        static object_matrix    eye(Integer r);

        // create sparse matrix with ones on diagonal of size rxc
        static object_matrix    speye(Integer r, Integer c);

        // create sparse matrix with ones on diagonal of size rxr
        static object_matrix    speye(Integer r);

        // create band matrix with ones on diagonal of size rxc with first diag fd
        // and last diag ld
        static object_matrix    beye(Integer r, Integer c, Integer fd, Integer ld);

        // create band matrix with ones on diagonal of size rxr with first diag fd
        // and last diag ld
        static object_matrix    beye(Integer r, Integer fd, Integer ld);

        // create dense matrix with elements on diagonal d given by v of size kxk, where 
        // k = length(v) + abs(d)
        static object_matrix    diag(const object_matrix& v, Integer d = 0);

        // create sparse matrix with elements on diagonal d given by v of size kxk, where 
        // k = length(v) + abs(d)
        static object_matrix    spdiag(const object_matrix &v, Integer d = 0);

        // create band matrix with elements on diagonal d given by v of size kxk, where 
        // k = length(v) + abs(d)
        static object_matrix    bdiag(const object_matrix &v, Integer d = 0);

        // create rxc dense matrix from the columns of A and places them along the diagonals 
        // specified by d 
        static object_matrix    diags(const object_matrix& A, const Matrix &d, Integer r, Integer c);

        // create rxc sparse matrix from the columns of A and places them along the diagonals 
        // specified by d 
        static object_matrix    spdiags(const object_matrix &A, const Matrix &d, Integer r, Integer c);

        // create rxc band matrix from the columns of A and places them along the diagonals 
        // specified by d 
        static object_matrix    bdiags(const object_matrix &A, const Matrix &d, Integer r, Integer c);

        // create dense matrix of size rows x cols with zeroes
        static object_matrix    make_dense(Integer rows,Integer cols);

        // create dense matrix of size rows x cols with elements copied from array
        // arr (column oriented storage assumed); optional argument set leading dimension 
        // (distance between consecutive columns)
        static object_matrix    make_dense(Integer rows,Integer cols, const Object *arr);
        static object_matrix    make_dense(Integer rows,Integer cols, const Object *arr, Integer ld);

        // create dense matrix of size rows x cols with elements stored in the array arr 
        // (column oriented storage assumed) with leading dimension ld (distance between
        // consecutive columns) array arr is not internally copied;
        //
        // this is unsafe function; one should ensure that this matrix and all matrices 
        // created from this matrix directly or indirectly are destroyed before destroying
        // the array arr
        static object_matrix    make_dense_foreign(Integer rows,Integer cols, Object* arr, 
                                        Integer ld);

        // create dense matrix of size rows x cols with elements val
        static object_matrix    make_dense(const T& val, Integer rows, Integer cols);
        static object_matrix    make_dense(const object_type& val, Integer rows, Integer cols);

        // create band matrix of size rows x cols with first diag fd and last diag ld,
        // with zeroes
        static object_matrix    make_band(Integer rows,Integer cols, Integer fd, Integer ld);

        // create band matrix of size rows x cols with first diag fd and last diag ld,
        // with elements val
        static object_matrix    make_band(const T& val,Integer rows,Integer cols, Integer fd, Integer ld);
        static object_matrix    make_band(const object_type& val,Integer rows,Integer cols, 
                                    Integer fd, Integer ld);

        // create sparse matrix of size rows x cols with nnz = 0, and capacity nzmax, 
        // with zero elements
        static object_matrix    make_sparse(Integer rows,Integer cols, Integer nzmax = 0);

        // create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_x
        // are arrays of length at least nnz; allocate space for at least nzmax elements
        static object_matrix    make_sparse(const Integer *trip_r, const Integer *trip_c, 
                                    const Object* trip_x, Integer r, Integer c, Integer nnz, 
                                    Integer nzmax = 0);

        // create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_x
        // are vectors of equal length nnz; allocate space for at least nzmax elements
        static object_matrix    make_sparse(const Matrix& trip_r, const Matrix& trip_c, 
                                    const Matrix& trip_x, Integer r, Integer c, Integer nzmax = 0);

        // create dense matrix of size rows x cols; elements are not initialized; 
        // on return ptr_data is the pointer to stored data
        static object_matrix    make_dense_noinit(Integer rows,Integer cols, Object*& ptr_data);

        // create band matrix of size rows x cols with first diag fd and last diag ld; 
        // elements are not initialized
        static object_matrix    make_band_noinit(Integer rows,Integer cols, Integer fd, Integer ld);

        // create sparse matrix of size rows x cols with capacity for nzmax nonzero elements; 
        // elements are not initialized;
        static object_matrix    make_sparse_noinit(Integer rows,Integer cols, Integer nzmax = 0);
};

};

#include "matcl-matrep/objects/object_matrix_functions.h"
#include "matcl-matrep/objects/details/object_matrix.inl"
