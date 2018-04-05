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

#include "matcl-matrep/details/mpl.h"
#include "matcl-matrep/details/isa.h"
#include "matcl-matrep/details/utils.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"

namespace matcl { namespace raw { namespace details
{

namespace md = matcl::details;

template<class M1, class M2>
struct kron_helper
{
    using val_type_1    = typename M1::value_type;
    using str_type_1    = typename M1::struct_type;
    using val_type_2    = typename M2::value_type;
    using str_type_2    = typename M2::struct_type;
    using val_ret       = typename md::unify_types<val_type_1,val_type_2>::type;

    using str_ret       = typename md::select_if
                        <
                            std::is_same<str_type_2,struct_dense>::value,
                            typename md::select_if
                            <
                                std::is_same<str_type_1,struct_banded>::value,
                                struct_sparse,
                                str_type_1
                            >::type,
                            struct_sparse
                        >::type;

    using ret_type      = Matrix<val_ret,str_ret>;

    static void eval(matcl::Matrix& ret, const M1& A, const M2& B);
};

}}}
