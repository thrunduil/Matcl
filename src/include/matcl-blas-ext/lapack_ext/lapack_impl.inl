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
#include "matcl-blas/lapack_ext.h"
#include "matcl-blas/blas_ext.h"

#include <algorithm>

// inline versions of functions implemented in LAPACK

namespace matcl { namespace lapack
{

template<class V>
typename details::enable_if_valid_real<V,V>::type
lapack::lapy2(const V& X, const V& Y)
{
    V XABS  = std::abs(X);
    V YABS  = std::abs(Y);
    V W     = std::max(XABS, YABS);
    V Z     = std::min(XABS, YABS);

    if (Z == V(0.0) )
        return W;
    else
        return W * std::sqrt(V(1.0) + pow_2(Z/W) );
};

template<class V>
typename details::enable_if_valid_real<V,V>::type
lapack::lapy3(const V& X, const V& Y, const V& Z)
{
    V XABS  = std::abs(X);
    V YABS  = std::abs(Y);
    V ZABS  = std::abs(Z);
    V W     = std::max(std::max(XABS, YABS), ZABS);

    if (W == V(0.0) )
        return XABS + YABS + ZABS;
    else
        return W * std::sqrt( pow_2(XABS/W) + pow_2(YABS/W) + pow_2(ZABS/ W) );
};

template<class V>
typename details::enable_if_valid<void,V>::type
lae2(typename details::real_type<V>::type A, V B0, typename details::real_type<V>::type C, 
     typename details::real_type<V>::type& RT1, typename details::real_type<V>::type& RT2)
{
    using VR    = typename details::real_type<V>::type;

    const VR ONE    = VR(1.0);
    const VR TWO    = VR(2.0);
    const VR ZERO   = VR(0.0);
    const VR HALF   = VR(0.0);

    // Compute the eigenvalues
    VR B            = abs(B0);
    VR SM           = A + C;
    VR DF           = A - C;
    VR ADF          = abs(DF);
    VR TB           = B + B;
    VR AB           = abs( TB );

    VR ACMX, ACMN;

    if ( abs( A ) > abs( C ) )
    {
        ACMX    = A;
        ACMN    = C;
    }
    else
    {
        ACMX    = C;
        ACMN    = A;
    };

    VR RT;

    if ( ADF > AB )
         RT     = ADF * std::sqrt(ONE + pow_2(AB/ADF) );
    else if ( ADF < AB )
         RT     = AB * std::sqrt(ONE + pow_2(ADF/AB) );
    else
        // Includes case AB=ADF=0
         RT     = AB * std::sqrt( TWO );

    if (SM < ZERO )
    {
        RT1     = HALF * ( SM-RT );

        // Order of execution important.
        // To get fully accurate smaller eigenvalue,
        // next line needs to be executed in higher precision.
         RT2    = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B;
    }
    else if (SM > ZERO )
    {
        RT1    = HALF * ( SM+RT );

        // Order of execution important.
        // To get fully accurate smaller eigenvalue,
        // next line needs to be executed in higher precision.
         RT2    = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B;
    }
    else
    {
        // Includes case RT1 = RT2 = 0
         RT1    = HALF*RT;
         RT2    = -HALF*RT;
    };
};

template<class V>
typename details::enable_if_valid_real<void,V>::type
lapack::las2(V F, V G, V H, V& SSMIN, V& SSMAX )
{
    const V ZERO    = V(0.0);
    const V ONE     = V(1.0);
    const V TWO     = V(2.0);

    V FA            = abs(F);
    V GA            = abs(G);
    V HA            = abs(H);
    V FHMN          = std::min( FA, HA );
    V FHMX          = std::max( FA, HA );
      
    if( FHMN == ZERO )
    {
        SSMIN       = ZERO;

        if( FHMX == ZERO )
            SSMAX   = GA;
        else
            SSMAX   = std::max(FHMX, GA) * std::sqrt(ONE + pow_2(std::min(FHMX, GA) / std::max(FHMX, GA) ) );
    }
    else
    {
        if (GA < FHMX )
        {
            V AS    = ONE + FHMN / FHMX;
            V AT    = ( FHMX-FHMN ) / FHMX;
            V AU    = pow_2( GA / FHMX );
            V C     = TWO / ( std::sqrt( AS*AS+AU )+std::sqrt( AT*AT+AU ) );
            SSMIN   = FHMN*C;
            SSMAX   = FHMX / C;
        }
        else
        {
            V AU    = FHMX / GA;

            if ( AU == ZERO )
            {
                // Avoid possible harmful underflow if exponent range
                // asymmetric (true SSMIN may not underflow even if AU underflows)
                SSMIN   = ( FHMN*FHMX ) / GA;
                SSMAX   = GA;
            }
            else
            {
                V AS    = ONE + FHMN / FHMX;
                V AT    = ( FHMX-FHMN ) / FHMX;
                V C     = ONE/( std::sqrt(ONE + pow_2(AS*AU)) + std::sqrt(ONE + pow_2(AT*AU) ) );
                SSMIN   = ( FHMN*C )*AU;
                SSMIN   = SSMIN + SSMIN;
                SSMAX   = GA / ( C+C );
            };
        };
    };
};

}};