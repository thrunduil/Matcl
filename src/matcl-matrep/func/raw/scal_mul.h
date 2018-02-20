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

//TODO
/*
template<class Val_ret, class M1, class M2,class str_1, class str_2>
struct eval_mult { };
*/

template<class Val_ret, class M1, class M2, class str_1>
struct eval_mult_abs { };

}}};

namespace matcl { namespace raw { namespace details
{

template<class T1, class T2>
struct mult_helper_mat_scal
{
    static void eval(matcl::Matrix& ret, const T1& mat, const T2& scal, trans_type t_A, trans_type t_B);
};

template<class T1, class T2>
struct mult_helper_scal_mat
{
    static void eval(matcl::Matrix& ret, const T1& scal, const T2& mat, trans_type t_A, trans_type t_B);
};

}}}
