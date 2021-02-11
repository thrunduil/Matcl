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
#include "matcl-blas-ext/lapack_ext/blas_ext.h"

#include <algorithm>

namespace matcl { namespace lapack
{

/*
*  Purpose
*  =======
*
*  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
*  eigenvector for RT1, giving the decomposition
*
*     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
*     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  B       (input) DOUBLE PRECISION
*          The (1,2) element and the conjugate of the (2,1) element of
*          the 2-by-2 matrix.
*
*  C       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  RT1     (output) DOUBLE PRECISION
*          The eigenvalue of larger absolute value.
*
*  RT2     (output) DOUBLE PRECISION
*          The eigenvalue of smaller absolute value.
*
*  CS1     (output) DOUBLE PRECISION
*  SN1     (output) DOUBLE PRECISION
*          The vector (CS1, SN1) is a unit right eigenvector for RT1.
*
* this is C++ version of Lapack function dlaev2
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid_real<void,V>::type
lapack::laev2(const V& A, const V& B, const V& C, V& RT1, V& RT2, V& CS1, V& SN1 )
{
    const V ZERO    = V(0.0);
    const V ONE     = V(1.0);
    const V TWO     = V(2.0);
    const V HALF    = V(0.5);

    // Compute the eigenvalues
    V SM        = A + C;
    V DF        = A - C;
    V ADF       = abs( DF );
    V TB        = B + B;
    V AB        = abs( TB );

    V  ACMX, ACMN;

    if( abs(A) > abs(C) )
    {
        ACMX    = A;
        ACMN    = C;
    }
    else
    {
        ACMX    = C;
        ACMN    = A;
    };

    V RT;

    if (ADF > AB )
         RT     = ADF * std::sqrt(ONE + matcl::lapack::pow_2( AB / ADF ) );
    else if (ADF < AB )
         RT     = AB * std::sqrt(ONE +  matcl::lapack::pow_2( ADF / AB ) );
    else
    {
        // Includes case AB=ADF=0
        RT      = AB * std::sqrt(TWO);
    };

    i_type SGN1;

    if (SM < ZERO )
    {
        RT1     = HALF * ( SM-RT );
        SGN1    = -1;

        // Order of execution important.
        // To get fully accurate smaller eigenvalue,
        // next line needs to be executed in higher precision.

        RT2     = ( ACMX / RT1 ) * ACMN - ( B / RT1 )*B;
    }
    else if( SM > ZERO )
    {
        RT1     = HALF * ( SM+RT );
        SGN1    = 1;

        // Order of execution important.
        // To get fully accurate smaller eigenvalue,
        // next line needs to be executed in higher precision.
        RT2     = ( ACMX / RT1 ) * ACMN - ( B / RT1 )*B;
    }
    else
    {
        // Includes case RT1 = RT2 = 0
        RT1     = HALF * RT;
        RT2     = -HALF * RT;
        SGN1    = 1;
    };

    // Compute the eigenvector
    V CS;
    i_type SGN2;

    if (DF >= ZERO )
    {
         CS         = DF + RT;
         SGN2       = 1;
    }
    else
    {
         CS         = DF - RT;
         SGN2       = -1;
    };

    V ACS           = abs( CS );

    if (ACS > AB )
    {
         V CT       = -TB / CS;
         SN1        = ONE / std::sqrt(ONE + CT*CT );
         CS1        = CT * SN1;
    }
    else
    {
        if ( AB == ZERO )
        {
            CS1     = ONE;
            SN1     = ZERO;
        }
        else
        {
            V TN    = -CS / TB;
            CS1     = ONE / std::sqrt(ONE + TN*TN );
            SN1     = TN*CS1;
        };
    };

    if (SGN1 == SGN2 )
    {
        V TN        = CS1;
        CS1         = -SN1;
        SN1         = TN;
    };
};

template BLAS_EXT_EXPORT void
lapack::laev2<d_type>(const d_type& A, const d_type& B, const d_type& C, d_type& RT1, d_type& RT2, 
                      d_type& CS1, d_type& SN1 );

template BLAS_EXT_EXPORT void
lapack::laev2<s_type>(const s_type& A, const s_type& B, const s_type& C, s_type& RT1, s_type& RT2, 
                      s_type& CS1, s_type& SN1 );


};};