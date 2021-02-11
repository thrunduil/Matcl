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

#include "matcl-linalg/graph/graph_manip.h"

namespace matcl { namespace details
{

struct wrappers_utils
{
    using Mat_I = raw::Matrix<Integer, struct_dense>;

    void    build_node_weights(Mat_I& node_weights, const Matrix& W, bool only_nonnegative,
                               Integer max_val);
    void    build_edge_weights(Mat_I& edge_weights, const Matrix& A, bool only_positive,
                               Integer max_val);
    void    get_structure(const Matrix& A, const Integer*& ptr_c, const Integer*& ptr_r);

    Matrix  convert_to_partitioning(const Mat_I& mat_res, Integer n_part, bool with_sep, 
                                    bool rem_empty, value_code vc);

    template<class V>
    void    build_node_weights_impl(const Matrix& W, Mat_I& node_weights, bool only_nonnegative,
                                    Integer max_val);

    template<class V>
    void    get_structure_impl(const Matrix& A, const Integer*& ptr_c, const Integer*& ptr_r);

    template<class V>
    Matrix  convert_to_partitioning_impl(const Mat_I& mat_res, Integer n_part, bool with_sep, 
                                    bool rem_empty);
};

}};