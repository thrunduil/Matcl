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

#include "matcl-matrep/details/mpl.h"
#include "matcl-matrep/details/isa.h"
#include "matcl-matrep/details/utils.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"

namespace matcl { namespace raw { namespace details
{

namespace md = matcl::details;

template<class MP>
struct MATCL_MATREP_EXPORT scalfunc_real_helper
{
    using value_type            = typename MP::value_type;
    using struct_type           = typename MP::struct_type;
    using real_value_type       = typename md::real_type<value_type>::type;
    using real_int_value_type   = typename md::real_type_int_real<value_type>::type;
    using ret_type_conj         = MP;
    using ret_type_arg          = Matrix<real_int_value_type,struct_type>;

    using ret_type
        = typename matcl::details::select_if
        <
            md::is_complex<value_type>::value,
            Matrix<real_value_type,struct_type>,
            MP
        >::type;    

    using ret_type_imag
        = typename matcl::details::select_if
        <
            md::is_complex<value_type>::value,
            Matrix<real_value_type,struct_type>,
            Matrix<real_value_type,struct_sparse>
        >::type;    

    // inplace is allowed; refcount must be increased for 
    // nontemporary objects; TODO
    static void eval_real(matcl::Matrix& ret, const MP& m);
    static void eval_imag(matcl::Matrix& ret, const MP& m);
    static void eval_conj(matcl::Matrix& ret, const MP& m);
    static void eval_abs(matcl::Matrix& ret, const MP& m);
    static void eval_abs2(matcl::Matrix& ret, const MP& m);
    static void eval_arg(matcl::Matrix& ret, const MP& m);
    static void eval_eps(matcl::Matrix& ret, const MP& m);
};

};};};
