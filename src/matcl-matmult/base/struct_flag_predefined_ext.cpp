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

#include "matcl-matrep/matrix/struct_flag_ext.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-matrep/lib_functions/func_matrix.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-matrep/lib_functions/vecfunc.h"

#include <algorithm>

namespace matcl 
{

namespace md = matcl::details;

unitary_flag::unitary_flag() 
    : user_flag(user_flag::get<unitary_flag>()) 
{};

bool unitary_flag::test(const matcl::Matrix& mat) const
{
    if (mat.rows() != mat.cols())
        return false;

    if (mat.rows() == 0)
        return true;

    value_code vc   = mat.get_value_code();
    Matrix tmp      = mat*matcl::ctrans(mat) - speye(mat.rows(), mat.cols(), vc);
    Matrix tmp2     = max_abs_vec(tmp);

    if (tmp2.get_scalar<Real>() <= 10.0 * mat.rows() * constants::eps(vc) )
        return true;
    else
        return false;
};

std::string unitary_flag::tag() const
{
    return "unitary";
}

user_flag unitary_flag::conj(struct_flag) const
{
    return unitary_flag();
}

user_flag unitary_flag::trans(struct_flag) const
{
    return unitary_flag();
}

user_flag unitary_flag::ctrans(struct_flag) const
{
    return unitary_flag();
}

user_flag unitary_flag::uminus(struct_flag) const
{
    return unitary_flag();
}

user_flag unitary_flag::mult_both(struct_flag sf_left, struct_flag sf_right, bool is_square) const
{
    (void)sf_left;
    (void)sf_right;
    (void)is_square;
    return unitary_flag();
}

user_flag unitary_flag::kron_both(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_left;
    (void)sf_right;
    return unitary_flag();
};

user_flag unitary_flag::kron_1(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_left;
    
    if (sf_right.is_id())
        return unitary_flag();
    else
        return user_flag();
};

user_flag unitary_flag::kron_2(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_right;
    
    if (sf_left.is_id())
        return unitary_flag();
    else
        return user_flag();
};

bool matcl::is_unitary(const struct_flag& sf)
{
    return sf.get_user(unitary_flag());
};

};
