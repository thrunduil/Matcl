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
#include "matcl-internals/container/mat_base.h"

namespace matcl { namespace raw 
{

template<class value_type_>
class MATCL_MATREP_EXPORT Matrix<value_type_, struct_banded> : public dense_matrix_base<value_type_>
{
    public:
        using base_type     = dense_matrix_base<value_type_>;
        using value_type    = value_type_;
        using struct_type   = struct_banded;
        using refcount_str  = matcl::details::refcount_str_default;

        using this_type     = Matrix<value_type_, struct_banded>;

        static const matcl::mat_code matrix_code
                            = md::type_to_code<this_type>::value;

    private:
        using DenseMatrix   = Matrix<value_type,struct_dense>;
        using tinfo         = ti::ti_type<value_type>;

        using base_type::ptr;
        using base_type::size;

    public:
        Matrix(tinfo ti);

        Matrix(tinfo ti, Integer r, Integer c, Integer fd, Integer ld);
        Matrix(tinfo ti,const value_type &val, Integer r, Integer c, Integer fd, Integer ld);

        Matrix(const Matrix&&) = delete;   
        Matrix(Matrix&&);   
        Matrix(const Matrix&);        

        Matrix              copy(bool keep_bufor = false) const;
        Matrix              clone(bool keep_bufor = false) const;
        Matrix              make_unique(bool keep_bufor = false) const;
        refcount_str*       get_refstr() const      { return base_type::get_refstr(); };
        void                destroy_data();
        bool                is_same_matrix(const Matrix& other) const;
        
        inline
        const value_type&   operator()(Integer i, Integer j) const; //1 base

        //only for raw matrix not stored in a container, for example allocated on stack
        void                assign_to_fresh(const Matrix&);
        void                assign_to_fresh(Matrix&&);

        inline Integer      rows() const;
        inline Integer      cols() const;
        inline Integer      length() const;
        inline Integer      impl_size() const       { return base_type::size(); };
        inline Integer      size() const            { return m_rows*m_cols; };

        // index of first nonzero diagonal (negative values for subdiagonals
        // zero for main diagonal, positive values for superdiagonals)
        inline Integer      first_diag() const;

        // index of last nonzero diagonal (negative values for subdiagonals
        // zero for main diagonal, positive values for superdiagonals)
        // last_diag >= first_diag
        inline Integer      last_diag() const;

        // equivalent to max(0, last_diag())
        inline Integer      number_superdiagonals() const;

        // equivalent to max(0, -first_diag())
        inline Integer      number_subdiagonals() const;
        
        Integer             nnz() const;

        // Get the leading dimmension of the representation of a banded matrix
        inline Integer      ld() const;

        inline Integer      max_udiags() const;
        inline Integer      max_ldiags() const;

        inline Integer      first_row(Integer col) const;
        inline Integer      last_row(Integer col) const;
        inline Integer      first_col(Integer row) const;
        inline Integer      last_col(Integer row) const;
        inline Integer      first_elem_pos(Integer col) const;
        inline Integer      first_elem_pos_row(Integer row) const;

        inline Integer      first_nonzero_column() const;
        inline Integer      last_nonzero_column() const;
        inline Integer      first_nonzero_row() const;
        inline Integer      last_nonzero_row() const;

        // offset of first element on the diagonal d from the root of
        // pointer to representation (d = 0 : main diagonal, d > 0: 
        // upper diagonals, d < 0: lower diagonals)
        inline Integer      first_elem_diag(Integer d) const;

        // check if given diagonal d is nonzero
        inline bool         has_diag(Integer d) const;

        // length of given diagonal d
        inline Integer      diag_length(Integer d) const;

        // row index of the first element on diagonal d (0-based)
        inline Integer      first_row_on_diag(Integer d) const;

        // columns index of the first element on diagonal d (0-based)
        inline Integer      first_col_on_diag(Integer d) const;

        // offset of given element from the root that contains element 
        // on row r and column c if exists (r, c are 0-based)
        inline Integer      element_pos(Integer r, Integer c) const;

        // Get the pointer to the representation of a banded matrix
        inline 
        const value_type*   rep_ptr() const;
        inline value_type*  rep_ptr();

        DenseMatrix         get_diag(Integer = 0) const;
        matcl::Matrix       fast_optim() const;
        bool                all_finite() const;

        Matrix              make_view(Integer rcs, Integer re, Integer ce) const;
        Matrix              make_view(Integer rcs, Integer re, Integer ce, Integer fd, Integer ld) const;

        // Returns a copy of a banded matrix with memory reserved for a matrix of given size
        // r - number of rows, c - number of columns, with unchanged number of sub- and 
        // superdiagonals
        Matrix              reserve(Integer r, Integer c) const;

        // Returns a copy of a banded matrix with memory reserved for a matrix of given size
        // r - number of rows, c - number of columns, fl - first diagonal, ld - last diagonal
        Matrix              reserve(Integer r, Integer c, Integer fd, Integer ld) const;

        // Returns a resized copy of a banded matrix, r - number of rows, c - number of columns
        Matrix              resize(Integer r, Integer c) const;
        
        // Returns a resized copy of a banded matrix, r - number of rows, c - number of columns
        Matrix              resize(Integer r, Integer c);

        // new version of resize
        // returns a resized copy of a banded matrix, r - number of rows, c - number of columns,
        // fd - first diagonal, ld - last diagonal
        Matrix              resize(Integer r, Integer c, Integer fd, Integer ld) const;

        // Returns a resized copy of a banded matrix, r - number of rows, c - number of columns,
        // fd - first diagonal, ld - last diagonal
        Matrix              resize(Integer r, Integer c, Integer fd, Integer ld);

        void                serialize(oarchive_impl & ar, const unsigned int version) const;
        void                serialize(iarchive_impl & ar, const unsigned int version); 

        // update struct flag based on number of sub- and super-diagonals
        void                update_struct() const;

    protected:
        Integer             m_ldiags;
        Integer             m_udiags;
        Integer             m_rows;
        Integer             m_cols;      

    private:
        Matrix&             operator=(const Matrix&) = delete;
        Matrix              get_diag_band() const;

        Matrix&             reset_unique();
        Matrix&             reset_unique(Integer r, Integer c, Integer fd, Integer ld);
};

};};

#include "matcl-internals/container/mat_b.inl"
