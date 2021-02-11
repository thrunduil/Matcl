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

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-blas-lapack/blas/config_blas.h"

namespace matcl { namespace lapack
{

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
qtriu2gschur(i_type n, i_type m, V* T, i_type ldt, V* S, i_type lds, bool wantq, bool wantz, 
             V* Q, i_type ldq, V* Z, i_type ldz)
{
    for (i_type i = 0; i < n-1; ++i)
    {
        bool is_compl   = (T[i+1 + i * ldt] != 0);

        if (is_compl == false)
            continue;

        V ALPHAR[2];
        V ALPHAI[2];
        V BETA[2];

        V SR, CR, SL, CL;

        //generalized schur decomposition of 2x2 block
        lapack::lagv2(T + i + i*ldt, ldt, S+ i + i*lds, lds, ALPHAR, ALPHAI, BETA, CL, SL, CR, SR);

        if (SL != 0 && n - i - 2 > 0)
            lapack::rot(n - i - 2, T + i + (i+2) * ldt, ldt, T + i+1 + (i+2) * ldt, ldt, CL, SL );

        if (SR != 0 && i > 0)
            lapack::rot(i, T + i*ldt, 1, T + (i+1)*ldt, 1, CR, SR);

        if (SL != 0 && n - i - 2 > 0)
            lapack::rot<V>(n - i - 2, S + i + (i+2)*lds, lds, S + i + 1 + (i+2)*lds, lds, CL, SL);            

        if (SR != 0 && i > 0)
            lapack::rot<V>(i, S + i * lds, 1, S + (i+1) * lds, 1, CR, SR);

        if(wantq && SL != 0)
            lapack::rot<V>(m, Q + i * ldq, 1, Q + (i+1)*ldq, 1, CL, SL );

        if (wantz &&  SR != 0)
            lapack::rot<V>(m, Z + i * ldz, 1, Z + (i+1)*ldz, 1, CR, SR );

        //ignore next row/columns
        ++i;
    };
};

//----------------------------------------------------------------
//         instantiations and specializations
//----------------------------------------------------------------
template<> 
BLAS_EXT_EXPORT void qtriu2gschur<c_type>(i_type n, i_type m, c_type* T, i_type ldt, c_type* S, i_type lds, 
             bool wantq, bool wantz, c_type* Q, i_type ldq, c_type* Z, i_type ldz)
{};
template<> 
BLAS_EXT_EXPORT void qtriu2gschur<z_type>(i_type n, i_type m, z_type* T, i_type ldt, z_type* S, i_type lds, 
             bool wantq, bool wantz, z_type* Q, i_type ldq, z_type* Z, i_type ldz)
{};

template void BLAS_EXT_EXPORT
qtriu2gschur<d_type>(i_type n, i_type m, d_type* T, i_type ldt, d_type* S, i_type lds, 
             bool wantq, bool wantz, d_type* Q, i_type ldq, d_type* Z, i_type ldz);

template void BLAS_EXT_EXPORT
qtriu2gschur<s_type>(i_type n, i_type m, s_type* T, i_type ldt, s_type* S, i_type lds, 
             bool wantq, bool wantz, s_type* Q, i_type ldq, s_type* Z, i_type ldz);

};};