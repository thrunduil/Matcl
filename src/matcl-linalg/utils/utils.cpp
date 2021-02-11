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

#include "matcl-linalg/utils/utils.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/utils/linalg_utils.h"

namespace matcl
{

Matrix invperm(const Matrix& p)
{
    Matrix pinv     = izeros(p.rows(),p.cols());
    pinv(p)         = irange(1,p.length());

    return pinv;
};

template<class Val>
Matrix details::make_nan_matrix(Integer rows, Integer cols)
{
    using VR    = typename md::real_type_int_real<Val>::type;
    VR nan      = constants::nan<VR>();

    Integer k   = min(rows,cols);
    Matrix I    = repmat(nan,k,1);

    Matrix ret = spdiags(I, 0, rows, cols);
    return ret;
};

Matrix details::make_nan_matrix(Integer rows, Integer cols, value_code vt)
{
    Integer k   = min(rows,cols);

    Matrix I;
    switch (vt)
    {
        case value_code::v_float:
        case value_code::v_float_complex:
        {
            Float nan   = constants::nan<Float>();    
            I           = repmat(nan,k,1);
            break;
        }
        default:
        {
            Real nan    = constants::nan<Real>();    
            I           = repmat(nan,k,1);
            break;
        };
    };

    Matrix ret = spdiags(I, 0, rows, cols);
    return ret;
};

template Matrix matcl::details::make_nan_matrix<Integer>(Integer, Integer);
template Matrix matcl::details::make_nan_matrix<Real>(Integer, Integer);
template Matrix matcl::details::make_nan_matrix<Float>(Integer, Integer);
template Matrix matcl::details::make_nan_matrix<Complex>(Integer, Integer);
template Matrix matcl::details::make_nan_matrix<Float_complex>(Integer, Integer);

};
