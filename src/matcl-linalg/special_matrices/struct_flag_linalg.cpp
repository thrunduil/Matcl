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

#include "matcl-linalg/special_matrices/struct_flag_linalg.h"
#include "matcl-linalg/decompositions/chol.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/lib_functions/func_matrix.h"
#include "matcl-matrep/matrix/colon.h"
#include "matcl-linalg/norms_error/norm.h"

namespace matcl 
{

namespace md = matcl::details;

//---------------------------------------------------------------------
//                      posdef_flag
//---------------------------------------------------------------------
bool posdef_flag::test(const matcl::Matrix& mat) const
{
    Matrix S;
    permvec p;
    Integer rank;
    tie(S,p,rank) = chol_rr(mat, false, 0.0);

    Integer N = mat.rows();

    if (rank < N)
        return false;
    else
        return true;
};

std::string posdef_flag::tag() const
{
    return "possitive_definite";
}

user_flag posdef_flag::real(struct_flag) const
{
    return posdef_flag();
};
user_flag posdef_flag::conj(struct_flag) const
{
    return posdef_flag();
};

user_flag posdef_flag::abs(struct_flag sf, bool is_square) const
{
    if (is_square && sf.is_diag())
        return posdef_flag();
    else
        return user_flag();
};
user_flag posdef_flag::trans(struct_flag) const
{
    return posdef_flag();
};
user_flag posdef_flag::ctrans(struct_flag) const
{
    return posdef_flag();
};

user_flag posdef_flag::plus_both(struct_flag, struct_flag) const
{
    return posdef_flag();
};

user_flag posdef_flag::plus_1(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_left;
    if (is_semi_posdef(sf_right) || sf_right.is_id())
        return posdef_flag();
    else
        return user_flag();
};

user_flag posdef_flag::plus_2(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_right;

    if (is_semi_posdef(sf_left) || sf_left.is_id())
        return posdef_flag();
    else
        return user_flag();
};

user_flag posdef_flag::scal(struct_flag, value_struct_class vc) const
{
    if (vc == value_struct_class::vc_one || vc == value_struct_class::vc_pos_real)
        return posdef_flag();
    else
        return user_flag();
};

user_flag posdef_flag::kron_both(struct_flag, struct_flag) const
{
    return posdef_flag();
};
user_flag posdef_flag::kron_1(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_left;

    if (sf_right.is_id())
        return posdef_flag();
    if (is_semi_posdef(sf_right))
        return semi_posdef_flag();
    else
        return user_flag();
};
user_flag posdef_flag::kron_2(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_right;

    if (sf_left.is_id())
        return posdef_flag();
    if (is_semi_posdef(sf_left))
        return semi_posdef_flag();
    else
        return user_flag();
};

user_flag posdef_flag::precision_increased(struct_flag) const
{
    return posdef_flag();
};

user_flag posdef_flag::mult_both(struct_flag sf_left, struct_flag sf_right, bool is_square) const
{
    if (is_square && sf_left.is_diag() && sf_right.is_diag())
        return posdef_flag();
    else
        return user_flag();
};
user_flag posdef_flag::mult_1(struct_flag sf_left, struct_flag sf_right, bool is_square) const
{
    if (is_square && sf_left.is_diag() && sf_right.is_diag())
    {
        if (is_semi_posdef(sf_right))
            return semi_posdef_flag();
    };

    return user_flag();
};
user_flag posdef_flag::mult_2(struct_flag sf_left, struct_flag sf_right, bool is_square) const
{
    if (is_square && sf_left.is_diag() && sf_right.is_diag())
    {
        if (is_semi_posdef(sf_left))
            return semi_posdef_flag();
    };

    return user_flag();
};

user_flag posdef_flag::seft_mult(struct_flag sf, trans_type t1, trans_type t2, value_code vc, bool is_square) const
{
    bool is_real    = matcl::matrix_traits::is_float_real(vc);

    if (t1 == trans_type::no_trans)
    {
        if (t2 == trans_type::conj_trans || t2 == trans_type::trans && is_real)
            return posdef_flag();
    }
    else if (t2 == trans_type::no_trans)
    {
        if (t1 == trans_type::conj_trans || t1 == trans_type::trans && is_real)
            return posdef_flag();
    };

    if (is_square && sf.is_diag())
        return posdef_flag();

    return user_flag();
};

bool posdef_flag::register_visit_scal() const
{
    return true;
};
bool posdef_flag::visit_scal(struct_flag mat_flags, value_struct_class vc) const
{
    if (vc == value_struct_class::vc_pos_real && mat_flags.is_id())
        return true;
    else 
        return false;
};

bool posdef_flag::register_visit_plus() const
{
    return true;
};
bool posdef_flag::visit_plus(struct_flag sf_X, struct_flag sf_Y) const
{
    if (sf_X.is_id() && sf_Y.is_id())
        return true;
    else
        return false;
};

bool matcl::is_posdef(const struct_flag& sf)
{
    return sf.get_user(posdef_flag());
};

//---------------------------------------------------------------------
//                      semi_posdef_flag
//---------------------------------------------------------------------
bool semi_posdef_flag::test(const matcl::Matrix& mat) const
{
    Matrix S;
    permvec p;
    Integer rank;
    tie(S,p,rank) = chol_rr(mat, false, 0.0);

    Integer N = mat.rows();

    if (rank == N)
        return true;

    value_code vc   = matcl::matrix_traits::real_value_type(mat.get_value_code());
    vc              = matcl::matrix_traits::unify_value_types(vc, value_code::v_float);
    Matrix dif      = S*matcl::ctrans(S) - mat(p,p);
    Real n_dif      = norm(dif, 1);
    Real n_A        = norm(mat, 1);
    // backward error from Higham, Accuracy and stability of numerical methods, (eq. 10.7)
    Real tol        = constants::eps(vc) * 4 * N * (3 * N + 1) * std::sqrt(N) * n_A;

    if (n_dif > tol * 10.0)
        return false;
    else
        return true;
};

std::string semi_posdef_flag::tag() const
{
    return "semi_possitive_definite";
}

user_flag semi_posdef_flag::abs(struct_flag sf, bool is_square) const
{
    if (is_square && sf.is_diag())
        return semi_posdef_flag();
    else
        return user_flag();
};

user_flag semi_posdef_flag::real(struct_flag) const
{
    return semi_posdef_flag();
};
user_flag semi_posdef_flag::conj(struct_flag) const
{
    return semi_posdef_flag();
};

user_flag semi_posdef_flag::trans(struct_flag) const
{
    return semi_posdef_flag();
};
user_flag semi_posdef_flag::ctrans(struct_flag) const
{
    return semi_posdef_flag();
};

user_flag semi_posdef_flag::plus_both(struct_flag, struct_flag) const
{
    return semi_posdef_flag();
};

user_flag semi_posdef_flag::plus_1(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_left;
    if (is_posdef(sf_right) || sf_right.is_id())
        return posdef_flag();
    else
        return user_flag();
};

user_flag semi_posdef_flag::plus_2(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_right;

    if (is_posdef(sf_left) || sf_left.is_id())
        return posdef_flag();
    else
        return user_flag();
};

user_flag semi_posdef_flag::scal(struct_flag, value_struct_class vc) const
{
    if (vc == value_struct_class::vc_one || vc == value_struct_class::vc_pos_real)
        return semi_posdef_flag();
    else
        return user_flag();
};

user_flag semi_posdef_flag::kron_both(struct_flag, struct_flag) const
{
    return semi_posdef_flag();
};
user_flag semi_posdef_flag::kron_1(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_left;

    if (sf_right.is_id() || is_posdef(sf_right))
        return semi_posdef_flag();
    else
        return user_flag();
};
user_flag semi_posdef_flag::kron_2(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_right;

    if (sf_left.is_id() || is_posdef(sf_left))
        return semi_posdef_flag();
    else
        return user_flag();
};

user_flag semi_posdef_flag::seft_mult(struct_flag sf, trans_type t1, trans_type t2, value_code vc, 
                                      bool is_square) const
{
    bool is_real    = matcl::matrix_traits::is_float_real(vc);

    if (t1 == trans_type::no_trans)
    {
        if (t2 == trans_type::conj_trans || t2 == trans_type::trans && is_real)
            return semi_posdef_flag();
    }
    else if (t2 == trans_type::no_trans)
    {
        if (t1 == trans_type::conj_trans || t1 == trans_type::trans && is_real)
            return semi_posdef_flag();
    };

    if (is_square && sf.is_diag())
        return semi_posdef_flag();

    return user_flag();
};

user_flag semi_posdef_flag::mult_both(struct_flag sf_left, struct_flag sf_right, bool is_square) const
{
    if (is_square && sf_left.is_diag() && sf_right.is_diag())
        return semi_posdef_flag();
    else
        return user_flag();
};
user_flag semi_posdef_flag::mult_1(struct_flag sf_left, struct_flag sf_right, bool is_square) const
{
    if (is_square && sf_left.is_diag() && sf_right.is_diag())
    {
        if (is_posdef(sf_right))
            return semi_posdef_flag();
    };

    return user_flag();
};
user_flag semi_posdef_flag::mult_2(struct_flag sf_left, struct_flag sf_right, bool is_square) const
{
    if (is_square && sf_left.is_diag() && sf_right.is_diag())
    {
        if (is_semi_posdef(sf_left))
            return semi_posdef_flag();
    };

    return user_flag();
};

user_flag semi_posdef_flag::precision_increased(struct_flag) const
{
    return semi_posdef_flag();
};

bool semi_posdef_flag::register_visit_abs() const
{
    return true;
};
bool semi_posdef_flag::visit_abs(struct_flag mat_flags, bool is_square) const
{
    if (is_square && mat_flags.is_diag())
        return true;
    else
        return false;
};

bool semi_posdef_flag::register_visit_self_mult() const
{
    return true;
};
bool semi_posdef_flag::visit_self_mult(struct_flag mat_flags, trans_type t1, trans_type t2, value_code vc,
                                       bool is_square) const
{
    bool is_float   = matcl::matrix_traits::is_float(vc);

    // do not set this flag for integer or object matrices; in case of integer matrices
    // overlow may happened
    if (is_float == false)
        return false;

    bool is_real    = matcl::matrix_traits::is_float_real(vc);

    if (t1 == trans_type::no_trans)
    {
        if (t2 == trans_type::conj_trans || t2 == trans_type::trans && is_real)
            return true;
    }
    else if (t2 == trans_type::no_trans)
    {
        if (t1 == trans_type::conj_trans || t1 == trans_type::trans && is_real)
            return true;
    };

    if (is_square && mat_flags.is_diag())
    {
        if (is_real || t1 == trans_type::no_trans && t2 == trans_type::conj_trans
                    || t2 == trans_type::no_trans && t1 == trans_type::conj_trans)
        {
            return true;
        };
    };

    return false;
};

bool matcl::is_semi_posdef(const struct_flag& sf)
{
    return sf.get_user(semi_posdef_flag());
};

};
