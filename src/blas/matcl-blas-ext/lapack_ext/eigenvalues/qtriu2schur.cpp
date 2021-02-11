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
lapack::qtriu2schur(i_type n, V* T, i_type ldt, bool wantq, V* Q, i_type QR, i_type ldq)
{
    for (i_type i = 0; i < n-1; ++i)
    {
        bool is_compl   = (T[i+1 + i * ldt] != 0);

        if (is_compl == false)
            continue;
        
        //find givens rotation
        V* A    = T + i     + i*ldt;
        V* B    = T + i     + (i+1)*ldt;
        V* C    = T + i + 1 + i*ldt;
        V* D    = T + i + 1 + (i+1)*ldt;

        V RT1R, RT1I, RT2R, RT2I;
        V CS, SN;

        lapack::lanv2(*A,*B,*C,*D, RT1R, RT1I, RT2R, RT2I, CS, SN);

        if (SN != 0)
        {
            // apply the transformation to the rest of T.
            i_type m    = n - i - 2;
            if ( m > 0)
                lapack::rot<V>(m, T + i + (i+2)*ldt, ldt, T + i + 1 + (i+2)*ldt, ldt, CS, SN );            

            if (i > 0)
                lapack::rot<V>(i, T + i * ldt, 1, T + (i+1) * ldt, 1, CS, SN );

            //apply givens rotations to q matrix
            if (wantq == true)
                lapack::rot(QR, Q + i * ldq, 1, Q + (i+1)*ldq, 1, CS, SN );
        };

        //ignore next row/columns
        ++i;
    };
};

//----------------------------------------------------------------
//         instantiations and specializations
//----------------------------------------------------------------
template<> 
BLAS_EXT_EXPORT void qtriu2schur<c_type>(i_type n, c_type* T, i_type ldt, bool wantq, c_type* Q, i_type QR, i_type ldq)
{};
template<> 
BLAS_EXT_EXPORT void qtriu2schur<z_type>(i_type n, z_type* T, i_type ldt, bool wantq, z_type* Q, i_type QR, i_type ldq)
{};

template void BLAS_EXT_EXPORT
qtriu2schur<d_type>(i_type n, d_type* T, i_type ldt, bool wantq, d_type* Q, i_type QR, i_type ldq);

template void BLAS_EXT_EXPORT
qtriu2schur<s_type>(i_type n, s_type* T, i_type ldt, bool wantq, s_type* Q, i_type QR, i_type ldq);

};};