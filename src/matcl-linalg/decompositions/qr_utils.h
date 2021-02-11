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

#include "matcl-matrep/matrix/matrix.h"
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-linalg/special_matrices/unitary_matrix.h"

namespace matcl { namespace details
{

template<class VTR, class VM>
struct make_copy
{
    using ret_type  = raw::Matrix<VTR,struct_dense>;
    using in_type   = raw::Matrix<VM,struct_dense>;

    static ret_type eval(const in_type& X, Integer M, Integer N)
    {
        if (M == X.rows() && N == X.cols())
        {
            return raw::converter<ret_type,in_type>::eval(X);
        }
        else
        {
            ret_type out(X.get_type(), M, N);

            VTR* ptr_Y      = out.ptr();
            const VM* ptr_X = X.ptr();
            Integer Y_ld    = out.ld();
            Integer X_ld    = X.ld();

            Integer N0      = X.cols();
            Integer K       = X.rows();

            for (Integer j = 0; j < N0; ++j)
            {
                for (Integer i = 0; i < K; ++i)
                {
                    ptr_Y[i]    = VTR(ptr_X[i]);
                };
                for (Integer i = K; i < M; ++i)
                {
                    ptr_Y[i]    = VTR(0.0);
                };

                ptr_X   += X_ld;
                ptr_Y   += Y_ld;
            };

            for (Integer j = N0; j < N; ++j)
            {
                for (Integer i = 0; i < M; ++i)
                {
                    ptr_Y[i]    = VTR(0.0);
                };

                ptr_Y   += Y_ld;
            };

            return out;
        };
    };
};
template<class VTR>
struct make_copy<VTR,VTR>
{
    using ret_type  = raw::Matrix<VTR,struct_dense>;
    using in_type   = raw::Matrix<VTR,struct_dense>;

    static ret_type eval(const in_type& X, Integer M, Integer N)
    {
        if (X.is_unique() && M == X.rows() && N == X.cols())
            return ret_type(X, ret_type::copy_is_safe());
        else if (M == X.rows() && N == X.cols())
            return X.copy();
        else
            return X.make_unique().resize(M,N);
    };
};

}};