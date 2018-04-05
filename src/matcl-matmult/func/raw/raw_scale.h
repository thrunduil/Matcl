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
#include "matcl-matrep/details/utils.h"

namespace matcl { namespace raw { namespace details
{

namespace md = matcl::details;

template<class MT>
struct scaling_helper
{	
    using V     = typename MT::value_type;
    using S     = typename MT::struct_type;
    using Mat   = raw::Matrix<V,S>;
    using Mat_D = raw::Matrix<V,struct_dense>;

    static void eval_rows(matcl::Matrix& ret, const Mat& m, const Mat_D& Dr);
    static void eval_cols(matcl::Matrix& ret, const Mat& m, const Mat_D& Dr);
    static void eval_rowscols(matcl::Matrix& ret, const Mat& m, const Mat_D& Dr, const Mat_D& Dc);
};

template<class MT>
struct scaling_real_helper
{	
    using V     = typename MT::value_type;
    using S     = typename MT::struct_type;
    using VR    = typename md::real_type<V>::type;
    using Mat   = raw::Matrix<V,S>;
    using Mat_D = raw::Matrix<VR,struct_dense>;

    static void eval_rows(matcl::Matrix& ret, const Mat& m, const Mat_D& Dr);
    static void eval_cols(matcl::Matrix& ret, const Mat& m, const Mat_D& Dr);
    static void eval_rowscols(matcl::Matrix& ret, const Mat& m, const Mat_D& Dr, const Mat_D& Dc);
};

};};};