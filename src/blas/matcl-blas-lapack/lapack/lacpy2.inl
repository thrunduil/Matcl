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
#include "matcl-blas-lapack/blas/details/config_blas_lib.h"
#include "matcl-simd/simd.h"
#include <algorithm>

namespace matcl { namespace lapack
{

template<class V>
struct lacpy2
{
    static void eval(const char *uplo, i_type m, i_type n, const V* a, 
                     i_type lda, V* b, i_type ldb)
    {
        if (uplo[0] == 'u' || uplo[0] == 'U')
        {
            for (i_type j = 0; j < n; ++j)
            {
                i_type maxi = std::min(m,j+1);
                for (i_type i = 0; i < maxi; ++i)
                    b[i]    = a[i];

                b           += ldb;
                a           += lda;
            };
        }
        else if (uplo[0] == 'l' || uplo[0] == 'L')
        {
            for (i_type j = 0; j < n; ++j)
            {
                for (i_type i = j; i < m; ++i)
                    b[i]    = a[i];

                b           += ldb;
                a           += lda;
            };
        }
        else
        {
            for (i_type j = 0; j < n; ++j)
            {
                for (i_type i = 0; i < m; ++i)
                    b[i]    = a[i];

                b           += ldb;
                a           += lda;
            };
        };
    }
};

};};