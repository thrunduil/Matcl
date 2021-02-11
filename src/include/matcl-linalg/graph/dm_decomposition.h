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

#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-linalg/details/linalg_fwd.h"

namespace matcl
{

/// finds the Dulmage-Mendelsohn decomposition of A; compute row and column
/// permutations p, q such that A(p,q) has block upper triangular form; A(p,q)
/// is decomposed into underdetermined, exactly determined, and overdetermined
/// parts:
///
///     A(p,q) = [ A_11 A_12  A_13]
///              [ 0    A_22  A_23]
///              [ 0    0     A_33]
///
/// underdetermined part (A_11) is rectangular matrix of size HR x HC, HC > HR
/// with maximum structural row rank; exactly determined part (A_22) is a square
/// matrix with full structural rank and overdetermined part (A_33) is rectangular
/// matrix of size VR x VC, VR > VC, with full structural column rank; some parts
/// may be empty; 
///
/// matrices A_11, A_33, are further decomposed into block diagonal matrices with
/// rectangular blocks of size m_i x n_i; in A_11 part m_i < n_i and i-th block 
/// has structure [A_i, D_i], where D_i is square matrix with zero free diagonal, 
/// m_i can be zero; in A_33 part m_i > n_i and i-th block has structure
/// [D_i; A_i], where D_i is square matrix with zero free diagonal, n_i can be zero; 
/// A_22 matrix is decomposed into block upper triangula matrices with square 
/// diagonal blocks, each diagonal block has zero free diagonal
///
/// algorithm: A. Pothen and C.J. Fan, "Computing the Block Triangular form of a Sparse 
//             Matrix", ACM Trans. on Mathematical Software, 16 (4), pp 303-324, 1990.
class MATCL_LINALG_EXPORT dm_decomp
{
    public:
        /// tuple of four integers
        using int_tup_4 = tuple<Integer,Integer,Integer,Integer>;

    private:
        permvec     m_rows_perm;
        permvec     m_cols_perm;
        Integer     m_comp_under;
        Integer     m_comp_square;
        Integer     m_comp_over;
        Matrix      m_split_rows;
        Matrix      m_split_cols;

        template<class V, class S>
        friend struct details::dmperm_impl;

    public:
        /// create uninitialized object
        dm_decomp();

        /// performe Dulmage-Mendelsohn decomposition of a matrix
        dm_decomp(const Matrix& A);

        /// standard destructor
        ~dm_decomp();

        /// get row permutations p
        permvec     row_perms() const;

        /// get column permutations q
        permvec     col_perms() const;

        /// get number of connected components in underdetermined part A_11
        Integer     components_under() const;

        /// get number of connected components in exactly determined part A_22
        Integer     components_square() const;

        /// get number of connected components in exactly overdetermined part A_33
        Integer     components_over() const;

        /// get total number of connected components
        Integer     components_total() const;

        /// get indices of first row in diagonal blocks; integer matrix of
        /// size (components_total() + 1) x 1
        Matrix      row_split() const;

        /// get indices of first column in diagonal blocks; integer matrix of
        /// size (components_total() + 1) x 1
        Matrix      col_split() const;

        /// get number of rows of underdetermined part A_11
        Integer     under_rows() const;

        /// get number of columns of underdetermined part A_11
        Integer     under_cols() const;

        /// get number of rows and columns of exactly determined part A_22
        Integer     square_size() const;

        /// get number of rows of overdetermined part A_33
        Integer     over_rows() const;

        /// get number of columns of overdetermined part A_33
        Integer     over_cols() const;

        /// get number of rows of the matrix A
        Integer     rows() const;

        /// get number of columns of the matrix A
        Integer     cols() const;

        /// return true if the matrix A has underdetermined part A_11
        bool        has_under() const;

        /// return true if the matrix A has exactly determined part A_22
        bool        has_square() const;

        /// return true if the matrix A has overdetermined part A_33
        bool        has_over() const;

        /// return structural rank of the matrix A (see sprank function
        /// for details)
        Integer     sprank() const;

        /// return [first_row, last_row, first_col, last_col] boundary of
        /// component number 1 <= comp <= components_under()
        int_tup_4   get_component_under(Integer comp) const;

        /// return [first_row, last_row, first_col, last_col] boundary of
        /// component number 1 <= comp <= components_square()
        int_tup_4   get_component_square(Integer comp) const;

        /// return [first_row, last_row, first_col, last_col] boundary of
        /// component number 1 <= comp <= components_over()
        int_tup_4   get_component_over(Integer comp) const;

        /// return [first_row, last_row, first_col, last_col] boundary of
        /// underdetermined part A_11
        int_tup_4   get_under_block() const;

        /// return [first_row, last_row, first_col, last_col] boundary of
        /// exactly determined part A_22
        int_tup_4   get_square_block() const;

        /// return [first_row, last_row, first_col, last_col] boundary of
        /// overdetermined part A_33
        int_tup_4   get_over_block() const;
};

};