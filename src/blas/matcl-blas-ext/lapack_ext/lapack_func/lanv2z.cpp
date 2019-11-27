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
#include "matcl-blas-ext/lapack_ext/blas_ext.h"

#include <algorithm>

namespace matcl { namespace lapack
{

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid_complex<void,V>::type
lapack::lanv2z(V& A, V& B, V& C, V& D, V& RT1, V& RT2, typename details::real_type<V>::type& CS, V& SN)
{
    using VR    = typename details::real_type<V>::type;

    const VR RONE   = VR(1.0);
    const VR RZERO  = VR(0.0);
    const V ZERO    = V(0.0);
    const V ONE     = V(1.0);
    const VR HALF   = VR(0.5);

    //  Initialize CS and SN
    CS              = RONE;
    SN              = ZERO;

    if (C == ZERO )
        goto lab_10;

    if (B == ZERO )
    {
        // Swap rows and columns
        CS          = RZERO;
        SN          = ONE;
        V TEMP      = D;
        D           = A;
        A           = TEMP;
        B           = -C;
        C           = ZERO;

        goto lab_10;
    }
    else if ( (A-D) == ZERO )
    {
        V TEMP      = std::sqrt(B * C);
        A           = A + TEMP;
        D           = D - TEMP;

        if ( ( B+C ) == ZERO )
        {
            CS      = std::sqrt( HALF );
            SN      = V(RZERO, RONE) * CS;
        }
        else
        {
            TEMP    = std::sqrt( B+C );
            V TEMP2 = std::sqrt( B ) / TEMP;
            CS      = real( TEMP2 );
            SN      = std::sqrt( C ) / TEMP;
        };

        B           = B - C;
        C           = ZERO;
        goto lab_10;
    }
    else
    {
        // Compute eigenvalue closest to D
        V T         = D;
        V U         = B*C;
        V X         = HALF*( A-T );
        V Y         = std::sqrt( X*X+U );

        if ( real(X) * real(Y) + imag(X) * imag(Y) < RZERO )
            Y       = -Y;

        V TEMP      = U / ( X+Y );
        T           = T - TEMP;

        // Do one QR step with exact shift T - resulting 2 x 2 in
        // triangular form.
        V AA;
        lartg( A-T, C, &CS, &SN, &AA );

         D          = D - T;
         V BB       = CS*B + SN*D;
         V DD       = -conj( SN )*B + CS*D;

         A          = AA*CS + BB*conj( SN ) + T;
         B          = -AA*SN + BB*CS;
         C          = ZERO;
         D          = T;
    };

  lab_10:

    // Store eigenvalues in RT1 and RT2.
    RT1             = A;
    RT2             = D;
};

template BLAS_EXT_EXPORT void
lapack::lanv2z(z_type& A, z_type& B, z_type& C, z_type& D, z_type& RT1, z_type& RT2, d_type& CS, z_type& SN);

template BLAS_EXT_EXPORT void
lapack::lanv2z(c_type& A, c_type& B, c_type& C, c_type& D, c_type& RT1, c_type& RT2, s_type& CS, c_type& SN);


}};