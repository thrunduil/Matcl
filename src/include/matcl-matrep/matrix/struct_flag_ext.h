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

#include "matcl-matrep/general/config.h"
#include "matcl-matrep/matrix/struct_flag.h"
#include "matcl-matrep/details/struct_flag_predefined.h"
#include "matcl-core/matrix/enums.h"

namespace matcl
{

// set of binary operations on struct_flags
struct MATCL_MATREP_EXPORT predefined_struct_ext : details::predefined_struct
{
    public:
        // matrix with struct_flag sf was multiplied by scalar value of class vt
        static struct_flag  get_scal_mult(const struct_flag& sf, value_struct_class vt);

        // rows or columns of a matrix with struct_flag sf was scaled
        static struct_flag  get_scal_rowcol(const struct_flag& sf);

        // kronecker product of matrices A and B; matrix A has struct f1, is real if is_re_1=true,
        // and is square if is_square_1 = true; matrix B has struct f2, is real if is_re_2 = true,
        // and is square if is_square_2 = true; if is_square_ret = true, then resulting
        // matrix is square
        static struct_flag  kron_struct(struct_flag f1, struct_flag f2, bool is_re_1, bool is_re_2, 
                                bool is_square_1, bool is_square_2, bool is_square_ret);

        // eval mmul(X, X, t1, t2); matrix X has struct f, value code vc, and is square if 
        // is_square = true; if is_square_ret is true, then resulting matrix is a square matrix
        static struct_flag  seft_mult_struct(struct_flag f, trans_type t1, trans_type t2, value_code vc, 
                                bool is_square, bool is_square_ret);

        // eval mmul(A, B, t1, t2); matrix A has struct f1, is real if is_re_1=true,
        // and is square if is_square_1 = true; matrix B has struct f2, is real if is_re_2 = true,
        // and is square if is_square_2 = true; if is_square_ret = true, then resulting
        // matrix is square
        static struct_flag  mult_struct(struct_flag f1, struct_flag f2, trans_type t1, trans_type t2,
                                bool is_re_1, bool is_re_2, bool is_square_1, bool is_square_2, 
                                bool is_square_ret);

        // evelement by element product of matrices A and B; matrix A has struct f1 and is real 
        // if is_re_1=true; matrix B has struct f2 and is real if is_re_2 = true;if is_square = true,
        // then resulting matrix is a square matrix
        static struct_flag  dmult_struct(struct_flag f1, struct_flag f2, bool is_re_1, bool is_re_2,
                                bool is_square);

        // eval A + B; matrix A has struct f1 and is real if is_re_1=true;
        // matrix B has struct f2 and is real if is_re_2 = true; if is_square is true, then resulting
        // matrix is a square matrix
        static struct_flag  plus_struct(struct_flag f1, struct_flag f2, bool is_re_1, bool is_re_2,
                                bool is_square);

        // eval A - B; matrix A has struct f1 and is real if is_re_1=true;
        // matrix B has struct f2 and is real if is_re_2 = true; if is_square is true, then resulting
        // matrix is a square matrix
        static struct_flag  minus_struct(struct_flag f1, struct_flag f2, bool is_re_1, bool is_re_2,
                                bool is_square);

    private:   
        enum op_type { op_mult, op_plus, op_minus, op_dmult, op_kron };

        //is_square_2 is used only if op = op_kron
        static struct_flag  op_struct(struct_flag f1, struct_flag f2, op_type op, 
                                bool is_real_1, bool is_real_2, bool is_square_1, bool is_square_2, 
                                bool is_square_ret);
        //is_square_2 is used only if op = op_kron
        static struct_flag  op_diags(struct_flag f1, struct_flag f2, op_type op, bool is_square_2);
        static struct_flag  op_symher(struct_flag f1, struct_flag f2, op_type op, 
                                      bool is_square_1, bool is_square_2, bool is_re_1, bool is_re_2);
};

// mark real or complex matrix as unitary, i.e. A*ctrans(A) = ctrans(A)*A = I
struct MATCL_MATMULT_EXPORT unitary_flag : user_flag_config, user_flag
{
    public:
        unitary_flag();

    private:
        virtual bool        test(const matcl::Matrix& mat) const override;
        virtual std::string tag() const override;

        virtual user_flag   conj(struct_flag sf) const override;
        virtual user_flag   trans(struct_flag sf) const override;
        virtual user_flag   ctrans(struct_flag sf) const override;
        virtual user_flag   uminus(struct_flag sf) const override;

        virtual user_flag   mult_both(struct_flag sf_left, struct_flag sf_right, bool is_square) const override;
        virtual user_flag   kron_both(struct_flag sf_left, struct_flag sf_right) const override;
        virtual user_flag   kron_1(struct_flag sf_left, struct_flag sf_right) const;
        virtual user_flag   kron_2(struct_flag sf_left, struct_flag sf_right) const;
};

// check if unitary flag was set
bool MATCL_MATMULT_EXPORT is_unitary(const struct_flag& sf);

};