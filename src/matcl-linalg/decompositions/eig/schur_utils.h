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

#pragma once

#include "matcl-matrep/matrix/matrix.h"
#include "matcl-matrep/details/isa.h"

namespace matcl { namespace details
{
    namespace md = matcl::details;

    void schur_reorder_check(const Matrix& I, const Integer N, Integer& M);
    bool schur_is_trivial_reorder(const Matrix& I);

    template<class V, bool Is_real = md::is_float_real_scalar<V>::value> 
    struct make_complex_eigenvectors
    {};

    template<class V> 
    struct make_complex_eigenvectors<V,true>
    {
        using Mat = raw::Matrix<V,struct_dense>;

        static Integer calc_sel_size(const Mat& TA, const Integer* SELECT, Integer N);

        static void eval(const Mat& TA, const Matrix& ind, const Matrix& mat_VL, const Matrix& mat_VR, 
                         Matrix& XL, Matrix& XR, bool comp_left, bool comp_right);

        static void complex_vectors_to_real(const Mat& TA, const Matrix& WL, const Matrix& WR, 
                                            Mat& WL_R, Mat& WR_R, const Matrix& ind);
    };

    template<class V> 
    struct make_complex_eigenvectors<V,false>
    {
        using Mat = raw::Matrix<V,struct_dense>;

        static Integer calc_sel_size(const Mat& TA, const Integer* SELECT, Integer N);

        static void eval(const Mat& TA, const Matrix& ind, const Matrix& mat_VL, const Matrix& mat_VR, 
                         Matrix& XL, Matrix& XR, bool comp_left, bool comp_right);

        static void complex_vectors_to_real(const Mat& TA, const Matrix& WL, const Matrix& WR, 
                                            Mat& WL_R, Mat& WR_R, const Matrix& ind);
    };

    /// apply unitary transformations to make symmetric tridiagonal matrix real
    template<class V>
    Matrix make_tridiag_subdiag_real(const Matrix& D1, V* ptr_U, bool is_D1_compl);
}}