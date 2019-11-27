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

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "utils.h"

#include <algorithm>

namespace matcl { namespace lapack
{

template<class V>
i_type lapack::check_hess_type(i_type N, const V* A, i_type LDA, const V* B, i_type LDB, i_type ILO, i_type IHI)
{    
    V ZERO = V(0.0);

    // check whether B if upper triangular
    B   = B + (ILO-1)*LDB;

    for (i_type I = ILO; I <= IHI; ++I)
    {
        for (i_type J = I + 1; J <= N; ++J)
        {
            if (B[J-1] != ZERO)
            {
                //QR and hess are required
                return 2;
            };
        };

        B   += LDB;
    };

    //check whether A is upper hessenberg
    A   = A + (ILO-1)*LDA;

    for (i_type I = ILO; I <= IHI; ++I)
    {
        for (i_type J = I + 2; J <= N; ++J)
        {
            if (A[J-1] != ZERO)
            {
                //QR is not required but hess is required
                return 1;
            };
        };

        A   += LDB;
    };

    //QR and hess are not required
    return 0;
};

template
i_type lapack::check_hess_type(i_type N, const d_type* A, i_type LDA, const d_type* B, i_type LDB, i_type ILO, i_type IHI);

template
i_type lapack::check_hess_type(i_type N, const s_type* A, i_type LDA, const s_type* B, i_type LDB, i_type ILO, i_type IHI);

template
i_type lapack::check_hess_type(i_type N, const c_type* A, i_type LDA, const c_type* B, i_type LDB, i_type ILO, i_type IHI);

template
i_type lapack::check_hess_type(i_type N, const z_type* A, i_type LDA, const z_type* B, i_type LDB, i_type ILO, i_type IHI);
}};
