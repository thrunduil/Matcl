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

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid_real<void,V>::type
lapack::lasv2(V F, V G, V H, V& SSMIN, V& SSMAX, V& SNR, V& CSR, V& SNL, V& CSL)
{
    const V ZERO    = V(0.0);
    const V ONE     = V(1.0);
    const V TWO     = V(2.0);
    const V FOUR    = V(4.0);
    const V HALF    = V(0.5);

    V EPS       = lapack::lamch<V>("E");
    V FT        = F;
    V FA        = abs(FT);
    V HT        = H;
    V HA        = abs(H);

    //  PMAX points to the maximum absolute element of matrix
    //  PMAX = 1 if F largest in absolute values
    //  PMAX = 2 if G largest in absolute values
    //  PMAX = 3 if H largest in absolute values

    i_type PMAX    = 1;
    bool SWAP       = ( HA > FA );

    if (SWAP)
    {        
        PMAX        = 3;

        std::swap(FT, HT);
        std::swap(FA, HA);

        // Now FA >= HA
    };

    V GT            = G;
    V GA            = abs( GT );

    V  CLT, CRT, SLT, SRT;

    if (GA == ZERO )
    {
        // Diagonal matrix
        SSMIN   = HA;
        SSMAX   = FA;
        CLT     = ONE;
        CRT     = ONE;
        SLT     = ZERO;
        SRT     = ZERO;
    }
    else
    {
        bool GASMAL = true;

        if (GA > FA )
        {
            PMAX            = 2;

            if ( ( FA / GA ) < EPS )
            {
                // Case of very large GA
                GASMAL      = false;
                SSMAX       = GA;

                if ( HA > ONE )
                    SSMIN   = FA / ( GA / HA );
                else
                    SSMIN   = ( FA / GA )*HA;

               CLT          = ONE;
               SLT          = HT / GT;
               SRT          = ONE;
               CRT          = FT / GT;
            };
        };

        if (GASMAL == true)
        {
            // Normal case
            V D             = FA - HA;
            V L;

            if (D == FA)
            {
                // Copes with infinite F or H
                L           = ONE;
            }
            else
            {
                L           = D / FA;
            };

            // Note that 0 <= L <= 1
            V M             = GT / FT;

            // Note that abs(M) <= 1/macheps
            V T             = TWO - L;

            // Note that T >= 1
            V MM            = M*M;
            V TT            = T*T;
            V S             = std::sqrt( TT+MM );

            // Note that 1 <= S <= 1 + 1/macheps

            V R;
            if ( L == ZERO )
               R            = abs( M );
            else
               R            = std::sqrt( L*L+MM );

            // Note that 0 <= R <= 1 + 1/macheps
            V A             = HALF * ( S+R );

            // Note that 1 <= A <= 1 + abs(M)
            SSMIN           = HA / A;
            SSMAX           = FA * A;

            if (MM == ZERO )
            {
                // Note that M is very tiny
                if ( L == ZERO )
                    T       = copysign(TWO, FT) * copysign(ONE, GT);
                else
                    T       = GT / copysign(D, FT) + M / T;
            }
            else
            {
                T           = (M / (S + T) + M / (R + L)) * (ONE + A);
            };

            L               = std::sqrt(T * T + FOUR);
            CRT             = TWO / L;
            SRT             = T / L;
            CLT             = ( CRT+SRT*M ) / A;
            SLT             = ( HT / FT )*SRT / A;
        };
    };

    if (SWAP)
    {
        CSL     = SRT;
        SNL     = CRT;
        CSR     = SLT;
        SNR     = CLT;
    }
    else
    {
        CSL     = CLT;
        SNL     = SLT;
        CSR     = CRT;
        SNR     = SRT;
    };

    V TSIGN;

    // Correct signs of SSMAX and SSMIN
    if (PMAX == 1)
        TSIGN   = copysign(ONE, CSR) * copysign(ONE, CSL) * copysign(ONE, F);
    else if (PMAX == 2)
        TSIGN   = copysign(ONE, SNR) * copysign(ONE, CSL) * copysign(ONE, G);
    else if (PMAX == 3)
        TSIGN   = copysign(ONE, SNR) * copysign(ONE, SNL) * copysign(ONE, H);

    SSMAX   = copysign(SSMAX, TSIGN );
    SSMIN   = copysign(SSMIN, TSIGN * copysign(ONE, F) * copysign(ONE, H) );
};

template void BLAS_EXT_EXPORT
lapack::lasv2<d_type>(d_type F, d_type G, d_type H, d_type& SSMIN, d_type& SSMAX, 
                          d_type& SNR, d_type& CSR, d_type& SNL, d_type& CSL );
template void BLAS_EXT_EXPORT
lapack::lasv2<s_type>(s_type F, s_type G, s_type H, s_type& SSMIN, s_type& SSMAX, s_type& SNR, s_type& CSR, 
                      s_type& SNL, s_type& CSL );

}};