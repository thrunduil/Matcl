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

/// type of gschur_sym decomposition
enum class gschur_sym_type
{
    A_B,        /// eigenproblem:   Ax  = lambda Bx
    AB,         /// eigenproblem:   ABx = lambda x
    BA          /// eigenproblem:   BAx = lambda x
};

/// perform generalized Schur decomposition of the matrix pair (A,B) 
/// where A and B are symmetric/hermitian and B is positive definite
/// in the form:
///
/// type:       eigen problem:      decomposition:
/// A_B:        Ax =  lambda Bx.   A * V = B * V * D    and V' * B * V = I;
/// AB:         ABx = lambda x.    A * B * V = V * D    and V' * B * V = I;
/// BA:         BAx = lambda x.    B * A * V = V * D    and V' * inv(B) * V = I
///
/// where V is square and D is diagonal matrix
/// this solver reduce this problem to the standard symmetric eigenvalue
/// problem implicitly inverting the matrix B, therefore B should be 
/// sufficiently positive definite,
///
/// not available for sparse matrices; available for band matrices if decomposition
/// type is A_B
class MATCL_LINALG_EXPORT gschur_sym_decomposition
{
    public:
        gschur_sym_decomposition();        
        ~gschur_sym_decomposition();

        //------------------------------------------------------------------------
        //                      FACTORIZATION
        //------------------------------------------------------------------------

        /// compute decomposition of the matrix pair A, B; if with_V is false, then
        /// unitary matricex V is not computed        
        gschur_sym_decomposition(const Matrix &A, const Matrix &B, gschur_sym_type type, bool with_V = true);
        gschur_sym_decomposition(Matrix &&A, const Matrix& B, gschur_sym_type type, bool with_V = true);
        gschur_sym_decomposition(const Matrix &A, Matrix&& B, gschur_sym_type type, bool with_V = true);
        gschur_sym_decomposition(Matrix&& A, Matrix&& B, gschur_sym_type type, bool with_V = true);

        /// compute decomposition of another matrix pair; see constructor for details
        gschur_sym_decomposition&   operator()(const Matrix &A, const Matrix &B, gschur_sym_type type, 
                                               bool with_V = true);
        gschur_sym_decomposition&   operator()(Matrix&& A, const Matrix &B, gschur_sym_type type, bool with_V = true);
        gschur_sym_decomposition&   operator()(const Matrix &A, Matrix&& B, gschur_sym_type type, bool with_V = true);
        gschur_sym_decomposition&   operator()(Matrix&& A, Matrix&& B, gschur_sym_type type, bool with_V = true);

        /// return the eigenvector matrix
        Matrix                  V() const;
                
        /// return the diagonal matrix D with eigenvalues on diagonal sorted in asceding
        /// order
        Matrix                  D() const;
                        
        /// return the generalized eigenvalues sorted in asceding order
        Matrix                  eig() const;
        
    private:
        void                    clear();
        void                    compute(const Matrix &A, const Matrix &B, gschur_sym_type type, bool with_V);
        bool                    test_factors();        

    private:
        bool                    m_has_Q;
        bool                    m_is_nan;
        Matrix                  m_V;
        Matrix                  m_D;
        Matrix                  m_eig;

        template<class Val, class Str>
        friend struct details::gschur_sym_str;
};

};