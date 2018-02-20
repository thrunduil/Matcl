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

#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"

namespace matcl { namespace raw { namespace details
{

namespace md = matcl::details;

template<class MP>
struct sort_helper
{    
    using val_type      = typename MP::value_type;
    using str_type      = typename MP::struct_type;
    using ret_str_type  = typename matcl::details::select_if
                        <
                            std::is_same<str_type,struct_banded>::value,
                            struct_sparse,
                            str_type
                        >::type;

    using ret_type_sort     = Matrix<val_type,ret_str_type>;
    using ret_type_sortrows = Matrix<val_type,ret_str_type>;
    using ret_type_sortcols = Matrix<val_type,ret_str_type>;
    using ret_type_issorted = Matrix<Integer,struct_dense>;

    using ret_type_sortrows_2   = std::pair<ret_type_sortrows,integer_dense>;
    using ret_type_sortcols_2   = std::pair<ret_type_sortcols,integer_dense>;


    static void     eval_sort(matcl::Matrix& i, const MP& m, Integer dim, bool asceding);
    static void     eval_sort_2(matcl::Matrix& i, matcl::Matrix& x, const MP& m, Integer dim,
                                bool asceding);

    static void     eval_sortrows(matcl::Matrix& i, const MP& m);
    static void     eval_sortrows_2(matcl::Matrix& i, matcl::Matrix& x, const MP& m);

    static void     eval_sortrows(matcl::Matrix& i, const MP& m, const integer_dense& dims);
    static void     eval_sortrows_2(matcl::Matrix& i, matcl::Matrix& x, const MP& m, 
                                const integer_dense& dims);

    static void     eval_sortcols(matcl::Matrix& i, const MP& m);
    static void     eval_sortcols_2(matcl::Matrix& i, matcl::Matrix& x, const MP& m);

    static void     eval_sortcols(matcl::Matrix& i, const MP& m, const integer_dense& dims);
    static void     eval_sortcols_2(matcl::Matrix& i, matcl::Matrix& x, const MP& m, 
                                const integer_dense& dims);

    static void     eval_issorted(matcl::Matrix& ret, const MP& m, Integer dim, bool asceding);
    static bool     eval_issorted_rows(const MP& m);
    static bool     eval_issorted_cols(const MP& m);
};

};};};
