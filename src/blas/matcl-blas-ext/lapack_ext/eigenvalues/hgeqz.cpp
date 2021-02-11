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
#include "blas/matcl-blas-ext/lapack_ext/utils/optim_params.h"

#include <algorithm>

namespace matcl { namespace lapack
{
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid_real<void,V>::type
hgeqz2(const char* JOB, const char* COMPQ, const char* COMPZ, i_type N, i_type ILO, i_type IHI,
    V* H, i_type LDH, V* T, i_type LDT, V* ALPHAR, V* ALPHAI, V* BETA, V* Q, i_type LDQ, V* Z, 
    i_type LDZ, V* WORK, i_type LWORK, i_type& INFO)
{
    bool    ILSCHR, ILQ, ILZ, LQUERY;
    i_type  ISCHUR, ICOMPQ, ICOMPZ;

    //     Decode JOB, COMPQ, COMPZ
    if (JOB[0] == 'e' || JOB[0] == 'E')
    {
        ILSCHR  = false;
        ISCHUR  = 1;
    }
    else if (JOB[0] == 's' || JOB[0] == 'S')
    {
        ILSCHR  = true;
        ISCHUR  = 2;
    }
    else 
    {
        ISCHUR  = 0;
    };

    if (COMPQ[0] == 'n' || COMPQ[0] == 'N')
    {
         ILQ    = false;
         ICOMPQ = 1;
    }
    else if (COMPQ[0] == 'v' || COMPQ[0] == 'V')
    {
        ILQ     = true;
        ICOMPQ  = 2;
    }
    else if (COMPQ[0] == 'i' || COMPQ[0] == 'I')
    {
         ILQ    = true;
         ICOMPQ = 3;
    }
    else
    {
         ICOMPQ = 0;
    }

    if (COMPZ[0] == 'n' || COMPZ[0] == 'N')
    {
        ILZ     = false;
        ICOMPZ  = 1;
    }
    else if (COMPZ[0] == 'v' || COMPZ[0] == 'V')
    {
         ILZ    = true;
         ICOMPZ = 2;
    }
    else if (COMPZ[0] == 'i' || COMPZ[0] == 'I')
    {
         ILZ    = true;
         ICOMPZ = 3;
    }
    else
    {
         ICOMPZ = 0;
    };

    // Check Argument Values
    INFO        = 0;

    WORK[0]     = V(std::max( 1, N ));
    LQUERY      = ( LWORK == -1 );
    if( ISCHUR == 0 )
         INFO   = -1;
    else if (ICOMPQ == 0 )
         INFO   = -2;
    else if (ICOMPZ == 0 )
         INFO   = -3;
    else if (N < 0 )
         INFO   = -4;
    else if (ILO < 1 )
         INFO   = -5;
    else if (IHI > N || IHI < ILO-1 )
         INFO   = -6;
    else if (LDH < N )
         INFO   = -8;
    else if (LDT < N )
         INFO   = -10;
    else if(LDQ < 1 || (ILQ && LDQ < N ) )
         INFO   = -15;
    else if(LDZ < 1 || (ILZ && LDZ < N ) )
         INFO   = -17;
    else if( LWORK < std::max( 1, N ) && LQUERY == false )
         INFO   = -19;

    if (INFO != 0 )
        return;
    else if (LQUERY)
        return;

    // Quick return if possible
    if ( N <= 0 )
    {
        WORK[0]     = V(1.0 );
        return;
    };

    // Initialize Q and Z
    if (ICOMPQ == 3 )
        lapack::laset("Full", N, N, V(0.0), V(1.0), Q, LDQ );
    if (ICOMPZ == 3 )
        lapack::laset("Full", N, N, V(0.0), V(1.0), Z, LDZ );

    // Machine Constants
    i_type IN   = IHI + 1 - ILO;
    V SAFMIN    = lamch<V>("S");
    V SAFMAX    = V(1.0) / SAFMIN;
    V ULP       = lamch<V>("E") * lamch<V>("B");
    V ANORM     = lapack::lanhs( "F", IN, H + (ILO-1) +(ILO-1)*LDH, LDH, WORK );
    V BNORM     = lapack::lanhs( "F", IN, T + (ILO-1) +(ILO-1)*LDT, LDT, WORK );
    V ATOL      = std::max(SAFMIN, ULP*ANORM );
    V BTOL      = std::max( SAFMIN, ULP*BNORM );
    V ASCALE    = V(1.0) / std::max( SAFMIN, ANORM );
    V BSCALE    = V(1.0) / std::max( SAFMIN, BNORM );
    V SAFETY    = 1.0E+2;
    V HALF      = 0.5;

    // Set Eigenvalues IHI+1:N
    for (i_type J = IHI + 1; J <= N; ++J)
    {
        if (T[J-1 +(J-1)*LDT] <= V(0.0) )
        {
            if (ILSCHR)
            {
                for (i_type JR = 1; JR <= J; ++JR)
                {
                    H[JR-1 +(J-1)*LDH] = -H[JR-1 +(J-1)*LDH];
                    T[JR-1 +(J-1)*LDT] = -T[JR-1 +(J-1)*LDT];
                };
            }
            else
            {
               H[J-1 +(J-1)*LDH]    = -H[J-1 +(J-1)*LDH];
               T[J-1 +(J-1)*LDT]    = -T[J-1 +(J-1)*LDT];
            };

            if (ILZ)
            {
                for (i_type JR = 1; JR <= N; ++JR)
                {
                    Z[JR-1 +(J-1)*LDZ]  = -Z[JR-1 +(J-1)*LDZ];
                };
            };
        };

        ALPHAR[J-1]     = H[J-1 +(J-1)*LDH];
        ALPHAI[J-1]     = V(0.0);
        BETA[J-1]       = T[J-1 +(J-1)*LDT];
    };

    // If IHI < ILO, skip QZ steps
    if (IHI < ILO )
        goto lab_380;

    //     MAIN QZ ITERATION LOOP
    //
    //     Initialize dynamic indices
    //
    //    Eigenvalues ILAST+1:N have been found.
    //       Column operations modify rows IFRSTM:whatever.
    //        Row operations modify columns whatever:ILASTM.
    //
    //     If only eigenvalues are being computed, then
    //        IFRSTM is the row of the last splitting row above row ILAST;
    //        this is always at least ILO.
    //    IITER counts iterations since the last eigenvalue was found,
    //       to tell when to use an extraordinary shift.
    //    MAXIT is the maximum number of QZ sweeps allowed.
    i_type ILAST        = IHI;
    i_type IFRSTM, ILASTM, IFIRST;

    if(ILSCHR)
    {
         IFRSTM         = 1;
         ILASTM         = N;
    }
    else
    {
         IFRSTM         = ILO;
         ILASTM         = IHI;
    };
    i_type IITER        = 0;
    V ESHIFT            = V(0.0);
    i_type MAXIT        = 30*(IHI - ILO + 1);

    for (i_type JITER = 1; JITER <= MAXIT; ++JITER)
    {
        // Split the matrix if possible.
        //       Two tests:
        //           1: H(j,j-1)=0  or  j=ILO
        //           2: T(j,j)=0
        if (ILAST == ILO )
        {            
            // Special case: j=ILAST
            goto lab_80;
        }
        else
        {
            if (abs( H[ILAST-1 + (ILAST-1-1)*LDH] ) <= ATOL )
            {
                H[ILAST-1 + (ILAST-1-1)*LDH]    = V(0.0);
                goto lab_80;
            };
        };

        if (abs(T[ILAST-1 + (ILAST-1)*LDT]) <= BTOL)
        {
            T[ILAST-1 + (ILAST-1)*LDT]  = V(0.0);
            goto lab_70;
        };
    
        // General case: j<ILAST
        for (i_type J = ILAST - 1; J >= ILO; --J)
        {
            i_type ILAZRO   = false;
            i_type ILAZR2   = false;

            // Test 1: for H(j,j-1)=0 or j=ILO        
            if (J == ILO )
            {
                ILAZRO  = true;
            }
            else
            {
                if (abs( H[J-1 + (J-1-1)*LDH] ) <= ATOL )
                {
                    H[J-1 + (J-1-1)*LDH]    = V(0.0);
                    ILAZRO  = true;
                }
                else
                {
                    ILAZRO  = false;
                };
            };

            // Test 2: for T(j,j)=0
            if (abs( T[J-1 + (J-1)*LDT] ) < BTOL )
            {
                T[J-1 +(J-1)*LDT]   = V(0.0);

                // Test 1a: Check for 2 consecutive small subdiagonals in A
                ILAZR2      = false;
                if (ILAZRO == false)
                {
                    V TEMP      = abs( H[J-1 + (J-1-1)*LDH] );
                    V TEMP2     = abs( H[J-1 + (J-1)*LDH] );
                    V TEMPR     = std::max(TEMP, TEMP2 );

                    if (TEMPR < V(1.0) && TEMPR != V(0.0) )
                    {
                        TEMP    = TEMP / TEMPR;
                        TEMP2   = TEMP2 / TEMPR;
                    };
                    if (TEMP * (ASCALE * abs( H[J+1-1 + (J-1)*LDH] ) ) <= TEMP2 * ( ASCALE*ATOL ) )
                        ILAZR2  = true;
                };

                // If both tests pass (1 & 2), i.e., the leading diagonal
                // element of B in the block is zero, split a 1x1 block off
                // at the top. (I.e., at the J-th row/column) The leading
                // diagonal element of the remainder can also be zero, so
                // this may have to be done repeatedly.
                if (ILAZRO || ILAZR2 )
                {                                
                    for (i_type JCH = J; JCH <= ILAST - 1; ++JCH)
                    {
                        V TEMP      = H[JCH-1 +(JCH-1)*LDH];
                        
                        V C, S;

                        lapack::lartg<V>(TEMP, H[JCH+1-1 + (JCH-1)*LDH], &C, &S, &H[JCH-1 + (JCH-1)*LDH] );
                        H[JCH+1-1 + (JCH-1)*LDH]    = V(0.0);

                        lapack::rot(ILASTM-JCH, &H[JCH-1 + (JCH+1-1)*LDH], LDH, &H[JCH+1-1 + (JCH+1-1)*LDH], LDH, C, S );
                        lapack::rot(ILASTM-JCH, &T[JCH-1 + (JCH+1-1)*LDT], LDT, &T[JCH+1-1 + (JCH+1-1)*LDT], LDT, C, S );

                        if ( ILQ )
                            lapack::rot( N, &Q[1-1 + (JCH-1)*LDQ], 1, &Q[1-1 + (JCH+1-1)*LDQ], 1, C, S );

                        if (ILAZR2)
                            H[JCH-1 + (JCH-1-1)*LDH]    = H[JCH-1 + (JCH-1-1)*LDH]*C;

                        ILAZR2  = false;

                        if ( abs( T[JCH+1-1 + (JCH+1-1)*LDT] ) >= BTOL )
                        {
                            if ( JCH+1 >= ILAST )
                            {
                                goto lab_80;
                            }
                            else
                            {
                                IFIRST  = JCH + 1;
                                goto lab_110;
                            };
                        };
                        
                        T[JCH+1-1 + (JCH+1-1)*LDT]  = V(0.0);
                    };
                        
                    goto lab_70;
                }
                else
                {
                    // Only test 2 passed -- chase the zero to T(ILAST,ILAST)
                    // Then process as in the case T(ILAST,ILAST)=0
                    for (i_type JCH = J; JCH <= ILAST - 1; ++JCH)
                    {
                        V TEMP  = T[JCH-1 + (JCH+1-1)*LDT];
                        V C, S;
                        lapack::lartg<V>(TEMP, T[JCH+1-1 + (JCH+1-1)*LDT], &C, &S, &T[JCH-1 + (JCH+1-1)*LDT]);

                        T[JCH+1-1 + (JCH+1-1)*LDT] = V(0.0);

                        if ( JCH < ILASTM-1 )
                            lapack::rot(ILASTM-JCH-1, &T[JCH-1 + (JCH+2-1)*LDT], LDT, &T[JCH+1-1+(JCH+2-1)*LDT], 
                                        LDT, C, S );

                        lapack::rot(ILASTM-JCH+2, &H[JCH-1 + (JCH-1-1)*LDH], LDH, &H[JCH+1-1 +(JCH-1-1)*LDH],
                                    LDH, C, S );

                        if (ILQ)
                            lapack::rot(N, &Q[1-1 + (JCH-1)*LDQ], 1, &Q[1-1 + (JCH+1-1)*LDQ], 1, C, S );

                        TEMP    = H[JCH+1-1 + (JCH-1)*LDH];
                        lapack::lartg<V>(TEMP, H[JCH+1-1 + (JCH-1-1)*LDH], &C, &S, &H[JCH+1-1 + (JCH-1)*LDH] );

                        H[JCH+1-1 + (JCH-1-1)*LDH]  = V(0.0);

                        lapack::rot<V>(JCH+1-IFRSTM, &H[IFRSTM-1 + (JCH-1)*LDH], 1, &H[IFRSTM-1 + (JCH-1-1)*LDH], 
                                       1, C, S );

                        lapack::rot<V>(JCH-IFRSTM, &T[IFRSTM-1 + (JCH-1)*LDT], 1, &T[IFRSTM-1 + (JCH-1-1)*LDT], 
                                       1, C, S );

                        if (ILZ)
                            lapack::rot<V>( N, &Z[1-1 + (JCH-1)*LDZ], 1, &Z[1-1 + (JCH-1-1)*LDZ], 1, C, S );
                    };
                    goto lab_70;
                };
            }
            else if( ILAZRO )
            {
                // Only test 1 passed -- work on J:ILAST
                IFIRST      = J;
                goto lab_110;
            };

            // Neither test passed -- try next J
       };
    
        // (Drop-through is "impossible")   
        INFO            = N + 1;
        goto lab_420;

        // T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
        //        1x1 block.

      lab_70:
        V TEMP  = H[ILAST-1 + (ILAST-1)*LDH];
        V C, S;
        lapack::lartg<V>(TEMP, H[ILAST-1 + (ILAST-1-1)*LDH], &C, &S, &H[ILAST-1 + (ILAST-1)*LDH] );

        H[ILAST-1  + (ILAST-1-1)*LDH]   = V(0.0);
        lapack::rot<V>(ILAST - IFRSTM, &H[IFRSTM-1 + (ILAST-1)*LDH], 1, &H[IFRSTM-1 + (ILAST-1-1)*LDH], 1, C, S );
        lapack::rot<V>(ILAST-IFRSTM, &T[IFRSTM-1 + (ILAST-1)*LDT], 1, &T[IFRSTM-1 + (ILAST-1-1)*LDT], 1, C, S );
             
        if( ILZ )
            lapack::rot<V>( N, &Z[1-1 + (ILAST-1)*LDZ], 1, &Z[1-1 + (ILAST-1)*LDZ], 1, C, S );
    
        // H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI,
        //                        and BETA
    
      lab_80:
        if (T[ILAST-1 + (ILAST-1)*LDT] < V(0.0))
        {
            if (ILSCHR)
            {
                for (i_type J = IFRSTM; J <= ILAST; ++J)
                {
                    H[J-1 + (ILAST-1)*LDH]  = -H[J-1 + (ILAST-1)*LDH];
                    T[J-1 + (ILAST-1)*LDT]  = -T[J-1 + (ILAST-1)*LDT];
                };
            }
            else
            {
                H[ILAST-1 + (ILAST-1)*LDH]  = -H[ILAST-1 + (ILAST-1)*LDH];
                T[ILAST-1 + (ILAST-1)*LDT]  = -T[ILAST-1 + (ILAST-1)*LDT];
            };

            if (ILZ)
            {
                for (i_type J = 1; J <= N; ++J)
                {
                    Z[J-1 + (ILAST-1)*LDZ]  = -Z[J-1 + (ILAST-1)*LDZ];
                };
            };
        };

        ALPHAR[ILAST-1]     = H[ILAST-1 + (ILAST-1)*LDH];
        ALPHAI[ILAST-1]     = V(0.0);
        BETA[ILAST-1]       = T[ILAST-1 + (ILAST-1)*LDT];
    
        // Go to next block -- exit if finished.
        ILAST               = ILAST - 1;
        if (ILAST < ILO )
            goto lab_380;
    
        // Reset counters
        IITER               = 0;
        ESHIFT              = V(0.0);

        if (ILSCHR == false)
        {
            ILASTM          = ILAST;
            if (IFRSTM > ILAST)
                IFRSTM      = ILO;
        };

        goto lab_350;

        //        QZ step
        //
        //        This iteration only involves rows/columns IFIRST:ILAST. We
        //        assume IFIRST < ILAST, and that the diagonal of B is non-zero.

      lab_110:
        IITER           = IITER + 1;

        if (ILSCHR == false)
            IFRSTM      = IFIRST;

        //        Compute single shifts.
        //
        //        At this point, IFIRST < ILAST, and the diagonal elements of
        //        T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
        //        magnitude)

        V S1, S2, WR, WR2, WI;

        if ( ( IITER / 10 )*10 == IITER )
        {
            // Exceptional shift.  Chosen for no particularly good reason.
            // (Single shift only.)
            if ( V(MAXIT) * SAFMIN * abs( H[ILAST-1 + (ILAST-1-1)*LDH] ) < abs( T[ILAST-1-1 + (ILAST-1-1)*LDT] ) )
            {
                ESHIFT = H[ILAST-1 + (ILAST-1-1)*LDH] / T[ILAST-1-1 + (ILAST-1-1)*LDT];
            }
            else
            {
                ESHIFT  = ESHIFT + V(1.0) / ( SAFMIN * V( MAXIT ) );
            };
            S1      = V(1.0);
            WR      = ESHIFT;
        }
        else
        {
            //  Shifts based on the generalized eigenvalues of the
            //  bottom-right 2x2 block of A and B. The first eigenvalue
            //  returned by DLAG2 is the Wilkinson shift (AEP p.512),

            lapack::lag2( &H[ILAST-1-1 + (ILAST-1-1)*LDH], LDH, &T[ILAST-1-1 + (ILAST-1-1)*LDT], LDT, 
                         SAFMIN * SAFETY, S1, S2, WR, WR2, WI );
    
            if ( abs( (WR/S1) * T[ILAST-1 + (ILAST-1)*LDT] - H[ILAST-1 + (ILAST-1)*LDH] )
                    > abs( (WR2/S2) * T[ILAST-1 + (ILAST-1)*LDT] - H[ILAST-1 + (ILAST-1)*LDH] ) )
            {
                TEMP    = WR;
                WR      = WR2;
                WR2     = TEMP;
                TEMP    = S1;
                S1      = S2;
                S2      = TEMP;
            };

            //TEMP        = std::max( S1, SAFMIN * std::max( std::max( V(1.0), abs( WR )), abs( WI ) ) );

            if (WI != V(0.0) )
                goto lab_200;
        };

        //        Fiddle with shift to avoid overflow   
        TEMP            = std::min( ASCALE, V(1.0) ) * ( HALF*SAFMAX );
        V SCALE;

        if ( S1 > TEMP )
            SCALE       = TEMP / S1;
        else
            SCALE       = V(1.0);

        TEMP            = std::min(BSCALE, V(1.0) ) * ( HALF*SAFMAX );

        if ( abs( WR ) > TEMP )
            SCALE       = std::min( SCALE, TEMP / abs( WR ) );

        S1              = SCALE*S1;
        WR              = SCALE*WR;

        i_type ISTART;

        // Now check for two consecutive small subdiagonals.
        for (i_type J = ILAST - 1; J >= IFIRST + 1; --J)
        {
            ISTART      = J;
            V TEMP      = abs( S1 * H[J-1 + (J-1-1)*LDH] );
            V TEMP2     = abs( S1 * H[J-1 +(J-1)*LDH] - WR * T[J-1 + (J-1)*LDT] );
            V TEMPR     = std::max(TEMP, TEMP2);
        
            if ( TEMPR < V(1.0) && TEMPR != V(0.0) )
            {
                TEMP    = TEMP / TEMPR;
                TEMP2   = TEMP2 / TEMPR;
            };
        
            if ( abs( ASCALE * H[J+1-1 +(J-1)*LDH] * TEMP ) <= ASCALE * ATOL * TEMP2 )
                goto lab_130;
        };

        ISTART          = IFIRST;

      lab_130:
        //  Do an implicit single-shift QZ sweep.
        //
        //  Initial Q

        V TEMPR;
        TEMP            = S1 * H[ISTART - 1 + (ISTART-1)*LDH] - WR * T[ISTART-1 + (ISTART-1)*LDT];
        V TEMP2         = S1 * H[ISTART+1-1 + (ISTART-1)*LDH];
        lapack::lartg(TEMP, TEMP2, &C, &S, &TEMPR );

        // Sweep
        for (i_type J = ISTART; J <= ILAST - 1; ++J)
        {
            if ( J > ISTART )
            {
                TEMP    = H[J-1 + (J-1-1)*LDH];
                lapack::lartg(TEMP, H[J+1-1 + (J-1-1)*LDH], &C, &S, &H[J-1 + (J-1-1)*LDH] );
                H[J+1-1 + (J-1-1)*LDH] = V(0.0);
            };

            for (i_type JC = J; JC <= ILASTM; ++JC)
            {
                TEMP                    = C*H[J-1 +(JC-1)*LDH] + S*H[J+1-1 +(JC-1)*LDH];
                H[J+1-1 + (JC-1)*LDH]   = -S*H[J-1 +(JC-1)*LDH] + C*H[J+1-1 + (JC-1)*LDH];
                H[J-1 + (JC-1)*LDH]     = TEMP;

                TEMP2                   = C*T[J-1 + (JC-1)*LDT] + S*T[J+1-1 + (JC-1)*LDT];
                T[J+1-1 + (JC-1)*LDT]   = -S*T[J-1 + (JC-1)*LDT] + C*T[J+1-1 + (JC-1)*LDT];
                T[J-1 + (JC-1)*LDT]     = TEMP2;
            };

            if (ILQ)
            {
                for (i_type JR = 1; JR <= N; ++JR)
                {
                    TEMP                    = C*Q[JR-1+(J-1)*LDQ] + S*Q[JR-1+(J+1-1)*LDQ];
                    Q[JR-1 + (J+1-1)*LDQ]   = -S*Q[JR-1 + (J-1)*LDQ] + C*Q[JR-1+(J+1-1)*LDQ];
                    Q[JR-1 + (J-1)*LDQ]     = TEMP;
                };
            };

            TEMP                        = T[J+1-1 + (J+1-1)*LDT];
            lapack::lartg(TEMP, T[J+1-1 + (J-1)*LDT], &C, &S, &T[J+1-1 + (J+1-1)*LDT] );
            T[J+1-1 + (J-1)*LDT]        = V(0.0);

            for (i_type JR = IFRSTM; JR <= std::min( J+2, ILAST ); ++JR)
            {
                TEMP                    = C*H[JR-1 + (J+1-1)*LDH] + S*H[JR-1 + (J-1)*LDH];
                H[JR-1 + (J-1)*LDH]     = -S*H[JR-1+ (J+1-1)*LDH] + C*H[JR-1 + (J-1)*LDH];
                H[JR-1 + (J+1-1)*LDH]   = TEMP;
            };

            for (i_type JR = IFRSTM; JR <= J; ++JR)
            {
                TEMP                    = C*T[JR-1 + (J+1-1)*LDT] + S*T[JR-1 + (J-1)*LDT];
                T[JR-1 + (J-1)*LDT]     = -S*T[JR-1+ (J+1-1)*LDT] + C*T[JR-1 + (J-1)*LDT];
                T[JR-1 + (J+1-1)*LDT]   = TEMP;
            };

            if (ILZ)
            {
                for (i_type JR = 1; JR <= N; ++JR)
                {
                    TEMP                = C*Z[JR-1 + (J+1-1)*LDZ] + S*Z[JR-1 + (J-1)*LDZ];
                    Z[JR-1 + (J-1)*LDZ] = -S*Z[JR-1 + (J+1-1)*LDZ] + C*Z[JR-1 + (J-1)*LDZ];
                    Z[JR-1 + (J+1-1)*LDZ] = TEMP;
                };
            };
        };

        goto lab_350;

        //        Use Francis double-shift
        //
        //        Note: the Francis double-shift should work with real shifts,
        //              but only if the block is at least 3x3.
        //              This code may break if this point is reached with
        //              a 2x2 block with real eigenvalues.

      lab_200:

        if (IFIRST+1 == ILAST )
        {
            //           Special case -- 2x2 block with complex eigenvectors
            //
            //           Step 1: Standardize, that is, rotate so that
            //
            //                       ( B11  0  )
            //                   B = (         )  with B11 non-negative.
            //                       (  0  B22 )

            V B22, B11, SR, CR, SL, CL;

            lapack::lasv2( T[ILAST-1-1 + (ILAST-1-1)*LDT], T[ILAST-1-1 +(ILAST-1)*LDT],
                           T[ILAST-1 + (ILAST-1)*LDT], B22, B11, SR, CR, SL, CL );

            if (B11 < V(0.0) )
            {
                CR      = -CR;
                SR      = -SR;
                B11     = -B11;
                B22     = -B22;
            };

            lapack::rot(ILASTM+1-IFIRST, &H[ILAST-1-1+(ILAST-1-1)*LDH], LDH,
                            &H[ILAST-1+ (ILAST-1-1)*LDH], LDH, CL, SL );

            lapack::rot(ILAST+1-IFRSTM, &H[IFRSTM-1+ (ILAST-1-1)*LDH], 1,
                            &H[IFRSTM-1+ (ILAST-1)*LDH], 1, CR, SR );

            if (ILAST < ILASTM )
                lapack::rot(ILASTM-ILAST, &T[ILAST-1-1+ (ILAST+1-1)*LDT], LDT,
                            &T[ILAST-1+ (ILAST+1-1)*LDT], LDT, CL, SL );

            if (IFRSTM < ILAST-1 )
                lapack::rot(IFIRST-IFRSTM, &T[IFRSTM-1+ (ILAST-1-1)*LDT], 1,
                            &T[IFRSTM-1+ (ILAST-1)*LDT], 1, CR, SR );

            if (ILQ)
                lapack::rot(N, &Q[1-1+ (ILAST-1-1)*LDQ], 1, &Q[1-1+ (ILAST-1)*LDQ], 1, CL, SL );

            if (ILZ)
                lapack::rot(N, &Z[1-1+ (ILAST-1-1)*LDZ], 1, &Z[1-1+ (ILAST-1)*LDZ], 1, CR, SR );

            T[ILAST-1-1+ (ILAST-1-1)*LDT]   = B11;
            T[ILAST-1-1+ (ILAST-1)*LDT]     = V(0.0);
            T[ILAST-1+ (ILAST-1-1)*LDT]     = V(0.0);
            T[ILAST-1+ (ILAST-1)*LDT]       = B22;

            // If B22 is negative, negate column ILAST
            if (B22 < V(0.0) )
            {
                for (i_type J = IFRSTM; J <= ILAST; ++J)
                {
                    H[J-1+ (ILAST-1)*LDH]   = -H[J-1+ (ILAST-1)*LDH];
                    T[J-1+ (ILAST-1)*LDT]   = -T[J-1+ (ILAST-1)*LDT];
                };

                if (ILZ)
                {
                    for (i_type J = 1; J <= N; ++J)
                    {
                        Z[J-1+ (ILAST-1)*LDZ]   = -Z[J-1+ (ILAST-1)*LDZ];
                    };
                };

                B22 = -B22;
            };

            // Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.)
            //
            // Recompute shift

            lapack::lag2(&H[ILAST-1-1+ (ILAST-1-1)*LDH], LDH, &T[ILAST-1-1+ (ILAST-1-1)*LDT], 
                         LDT, SAFMIN*SAFETY, S1, TEMP, WR, TEMP2, WI );

            // If standardization has perturbed the shift onto real line,
            // do another (real single-shift) QR step.
            if (WI == V(0.0))
                goto lab_350;

            V S1INV     = V(1.0) / S1;

            // Do EISPACK (QZVAL) computation of alpha and beta
            V A11       = H[ILAST-1-1+ (ILAST-1-1)*LDH];
            V A21       = H[ILAST-1+ (ILAST-1-1)*LDH];
            V A12       = H[ILAST-1-1+ (ILAST-1)*LDH];
            V A22       = H[ILAST-1+ (ILAST-1)*LDH];

            //           Compute complex Givens rotation on right
            //           (Assume some element of C = (sA - wB) > unfl )
            //                            __
            //           (sA - wB) ( CZ   -SZ )
            //                     ( SZ    CZ )
            V C11R      = S1*A11 - WR*B11;
            V C11I      = -WI*B11;
            V C12       = S1*A12;
            V C21       = S1*A21;
            V C22R      = S1*A22 - WR*B22;
            V C22I      = -WI*B22;

            V CZ, SZR, SZI;

            if (abs( C11R ) + abs(C11I) + abs(C12) > abs(C21) + abs(C22R) + abs(C22I) )
            {
                V T1    = lapack::lapy3( C12, C11R, C11I );
                CZ      = C12 / T1;
                SZR     = -C11R / T1;
                SZI     = -C11I / T1;
            }
            else
            {
                CZ      = lapack::lapy2( C22R, C22I );
                if (CZ <= SAFMIN )
                {
                    CZ  = V(0.0);
                    SZR = V(1.0);
                    SZI = V(0.0);
                }
                else
                {
                    TEMPR   = C22R / CZ;
                    V TEMPI = C22I / CZ;
                    V T1    = lapack::lapy2( CZ, C21 );
                    CZ      = CZ / T1;
                    SZR     = -C21*TEMPR / T1;
                    SZI     = C21*TEMPI / T1;
                };
            };

            //           Compute Givens rotation on left
            //
            //           (  CQ   SQ )
            //           (  __      )  A or B
            //           ( -SQ   CQ )

            V AN        = abs( A11 ) + abs( A12 ) + abs( A21 ) + abs( A22 );
            V BN        = abs( B11 ) + abs( B22 );
            V WABS      = abs( WR ) + abs( WI );

            V CQ, SQR, SQI;
            if ( S1 * AN > WABS * BN )
            {
                CQ      = CZ*B11;
                SQR     = SZR*B22;
                SQI     = -SZI*B22;
            }
            else
            {
                V A1R   = CZ*A11 + SZR*A12;
                V A1I   = SZI*A12;
                V A2R   = CZ*A21 + SZR*A22;
                V A2I   = SZI*A22;
                CQ      = lapack::lapy2( A1R, A1I );

                if (CQ <= SAFMIN )
                {
                    CQ  = V(0.0);
                    SQR = V(1.0);
                    SQI = V(0.0);
                }
                else
                {
                    TEMPR   = A1R / CQ;
                    V TEMPI = A1I / CQ;
                    SQR     = TEMPR*A2R + TEMPI*A2I;
                    SQI     = TEMPI*A2R - TEMPR*A2I;
                };
            };

            V T1    = lapy3( CQ, SQR, SQI );
            CQ      = CQ / T1;
            SQR     = SQR / T1;
            SQI     = SQI / T1;

            // Compute diagonal elements of QBZ
            TEMPR   = SQR*SZR - SQI*SZI;
            V TEMPI = SQR*SZI + SQI*SZR;
            V B1R   = CQ*CZ*B11 + TEMPR*B22;
            V B1I   = TEMPI*B22;
            V B1A   = lapy2( B1R, B1I );
            V B2R   = CQ*CZ*B22 + TEMPR*B11;
            V B2I   = -TEMPI*B11;
            V B2A   = lapy2( B2R, B2I );

            // Normalize so beta > 0, and Im( alpha1 ) > 0
            BETA[ILAST-1-1]     = B1A;
            BETA[ILAST-1]       = B2A;
            ALPHAR[ILAST-1-1]   = ( WR*B1A )*S1INV;
            ALPHAI[ILAST-1-1]   = ( WI*B1A )*S1INV;
            ALPHAR[ILAST-1]     = ( WR*B2A )*S1INV;
            ALPHAI[ILAST-1]     = -( WI*B2A )*S1INV;

            // Step 3: Go to next block -- exit if finished.
            ILAST               = IFIRST - 1;
            if (ILAST < ILO )
                goto lab_380;

            // Reset counters
            IITER               = 0;
            ESHIFT              = V(0.0);
            if (ILSCHR == false)
            {
                ILASTM          = ILAST;
                if (IFRSTM > ILAST)
                    IFRSTM      = ILO;
            };

            goto lab_350;
        }
        else
        {
            //           Usual case: 3x3 or larger block, using Francis implicit
            //                       double-shift
            //
            //                                    2
            //           Eigenvalue equation is  w  - c w + d = 0,
            //
            //                                         -1 2        -1
            //           so compute 1st column of  (A B  )  - c A B   + d
            //           using the formula in QZIT (from EISPACK)
            //
            //           We assume that the block is at least 3x3

            V AD11      = ( ASCALE * H[ILAST-1-1+ (ILAST-1-1)*LDH] ) / ( BSCALE * T[ILAST-1-1+ (ILAST-1-1)*LDT]);
            V AD21      = ( ASCALE * H[ILAST-1+ (ILAST-1-1)*LDH] ) / ( BSCALE * T[ILAST-1-1+ (ILAST-1-1)*LDT]);
            V AD12      = ( ASCALE * H[ILAST-1-1+ (ILAST-1)*LDH] ) / ( BSCALE * T[ILAST-1+ (ILAST-1)*LDT]);
            V AD22      = ( ASCALE * H[ILAST-1+ (ILAST-1)*LDH] ) / ( BSCALE * T[ILAST-1+ (ILAST-1)*LDT]);
            V U12       = T[ILAST-1-1+ (ILAST-1)*LDT] / T[ILAST-1+ (ILAST-1)*LDT];

            V AD11L     = ( ASCALE * H[IFIRST-1+ (IFIRST-1)*LDH] ) / ( BSCALE * T[IFIRST-1+ (IFIRST-1)*LDT] );
            V AD21L     = ( ASCALE * H[IFIRST+1-1+ (IFIRST-1)*LDH] ) / ( BSCALE * T[IFIRST-1+ (IFIRST-1)*LDT] );
            V AD12L     = ( ASCALE * H[IFIRST-1+ (IFIRST+1-1)*LDH] ) / ( BSCALE * T[IFIRST+1-1+ (IFIRST+1-1)*LDT] );
            V AD22L     = ( ASCALE * H[IFIRST+1-1+ (IFIRST+1-1)*LDH] ) / ( BSCALE * T[IFIRST+1-1+ (IFIRST+1-1)*LDT] );
            V AD32L     = ( ASCALE * H[IFIRST+2-1+ (IFIRST+1-1)*LDH] ) / ( BSCALE * T[IFIRST+1-1+ (IFIRST+1-1)*LDT] );
            V U12L      = T[IFIRST-1+ (IFIRST+1-1)*LDT] / T[IFIRST+1-1+ (IFIRST+1-1)*LDT];

            V VV[3];

            VV[1-1]     = ( AD11-AD11L )*( AD22-AD11L ) - AD12*AD21 + AD21*U12*AD11L + ( AD12L-AD11L*U12L )*AD21L;
            VV[2-1]     = ( ( AD22L-AD11L )-AD21L*U12L-( AD11-AD11L )- ( AD22-AD11L )+AD21*U12 )*AD21L;
            VV[3-1]     = AD32L*AD21L;

            ISTART      = IFIRST;
            V TAU;
            lapack::larfg<V>( 3, &VV[0], &VV[1], 1, &TAU );

            VV[0]       = V(1.0);

            //           Sweep
            for (i_type J = ISTART; J <= ILAST - 2; ++J)
            {
                // All but last elements: use 3x3 Householder transforms.
                // Zero (j-1)st column of A

                if ( J > ISTART )
                {
                    VV[0]   = H[J-1+ (J-1-1)*LDH];
                    VV[1]   = H[J+1-1+ (J-1-1)*LDH];
                    VV[2]   = H[J+2-1+ (J-1-1)*LDH];
        
                    lapack::larfg( 3, &H[J-1+(J-1-1)*LDH], &VV[1], 1, &TAU );

                    VV[0]                   = V(1.0);
                    H[J+1-1+(J-1-1)*LDH]    = V(0.0);
                    H[J+2-1+(J-1-1)*LDH]    = V(0.0);
                };

                for (i_type JC = J; JC <= ILASTM; ++JC)
                {
                    TEMP                = TAU * ( H[J-1+(JC-1)*LDH] + VV[1] * H[J+1-1+(JC-1)*LDH] 
                                                 + VV[2]* H[J+2-1+(JC-1)*LDH] );
                    H[J-1+(JC-1)*LDH]   = H[J-1+(JC-1)*LDH] - TEMP;
                    H[J+1-1+(JC-1)*LDH] = H[J+1-1+(JC-1)*LDH] - TEMP*VV[1];
                    H[J+2-1+(JC-1)*LDH] = H[J+2-1+(JC-1)*LDH] - TEMP*VV[2];

                    TEMP2               = TAU * ( T[J-1+(J-1)*LDT] + VV[1] * T[J+1-1+(JC-1)*LDT]
                                                 + VV[2] * T[J+2-1+(JC-1)*LDT] );
                    T[J-1+(JC-1)*LDT]   = T[J-1+(JC-1)*LDT] - TEMP2;
                    T[J+1-1+(JC-1)*LDT] = T[J+1-1+(JC-1)*LDT] - TEMP2 * VV[1];
                    T[J+2-1+(JC-1)*LDT] = T[J+2-1+(JC-1)*LDT] - TEMP2 * VV[2];
                };

                if (ILQ)
                {
                    for (i_type JR = 1; JR <= N; ++JR)
                    {
                        TEMP    = TAU * ( Q[JR-1+ (J-1)*LDQ] + VV[1] * Q[JR-1+ (J+1-1)*LDQ]
                                         + VV[2] * Q[JR-1+ (J+2-1)*LDQ] );
                        Q[JR-1+ (J-1)*LDQ]  = Q[JR-1+ (J-1)*LDQ] - TEMP;
                        Q[JR-1+ (J+1-1)*LDQ] = Q[JR-1+ (J+1-1)*LDQ] - TEMP*VV[1];
                        Q[JR-1+ (J+2-1)*LDQ] = Q[JR-1+ (J+2-1)*LDQ] - TEMP*VV[2];
                    };
                };

                // Zero j-th column of B (see DLAGBC for details)
                // Swap rows to pivot

                bool ILPIVT = false;
                TEMP        = std::max( abs(T[J+1-1+ (J+1-1)*LDT]), abs(T[J+1-1+(J+2-1)*LDT] ) );
                TEMP2       = std::max( abs(T[J+2-1+ (J+1-1)*LDT] ), abs(T[J+2-1+(J+2-1)*LDT] ) );

                V U1, U2, W11, W21, W12, W22;

                if (std::max(TEMP, TEMP2) < SAFMIN )
                {
                    SCALE   = V(0.0);
                    U1      = V(1.0);
                    U2      = V(1.0);

                    goto lab_250;
                }
                else if ( TEMP >= TEMP2 )
                {
                    W11     = T[J+1-1+(J+1-1)*LDT];
                    W21     = T[J+2-1+(J+1-1)*LDT];
                    W12     = T[J+1-1+(J+2-1)*LDT];
                    W22     = T[J+2-1+(J+2-1)*LDT];
                    U1      = T[J+1-1+(J-1)*LDT];
                    U2      = T[J+2-1+(J-1)*LDT];
                }
                else
                {
                    W21     = T[J+1-1+(J+1-1)*LDT];
                    W11     = T[J+2-1+(J+1-1)*LDT];
                    W22     = T[J+1-1+(J+2-1)*LDT];
                    W12     = T[J+2-1+(J+2-1)*LDT];
                    U2      = T[J+1-1+(J-1)*LDT];
                    U1      = T[J+2-1+(J-1)*LDT];
                };

                // Swap columns if nec.
                if (abs(W12) > abs( W11 ) )
                {
                    ILPIVT  = true;
                    TEMP    = W12;
                    TEMP2   = W22;
                    W12     = W11;
                    W22     = W21;
                    W11     = TEMP;
                    W21     = TEMP2;
                };

                // LU-factor
                TEMP        = W21 / W11;
                U2          = U2 - TEMP*U1;
                W22         = W22 - TEMP*W12;
                W21         = V(0.0);

                // Compute SCALE
                SCALE       = V(1.0);
                if ( abs( W22 ) < SAFMIN )
                {
                    SCALE   = V(0.0);
                    U2      = V(1.0);
                    U1      = -W12 / W11;
                    goto lab_250;
                };

                if (abs(W22) < abs(U2) )
                    SCALE   = abs(W22/U2);

                if (abs(W11) < abs(U1) )
                    SCALE   = std::min( SCALE, abs(W11/U1) );

                // Solve
                U2          = ( SCALE*U2 ) / W22;
                U1          = ( SCALE*U1-W12*U2 ) / W11;

              lab_250:
                if (ILPIVT)
                {
                    TEMP    = U2;
                    U2      = U1;
                    U1      = TEMP;
                };

                // Compute Householder Vector
                V T1        = sqrt( SCALE*SCALE + U1*U1 + U2*U2 );
                TAU         = V(1.0) + SCALE / T1;
                V VS        = -V(1.0) / (SCALE + T1);
                VV[0]       = V(1.0);
                VV[1]       = VS*U1;
                VV[2]       = VS*U2;

                // Apply transformations from the right.
                for (i_type JR = IFRSTM; JR <= std::min( J+3, ILAST ); ++JR)
                {
                    TEMP    = TAU * ( H[JR-1+ (J-1)*LDH] + VV[1] * H[JR-1+(J+1-1)*LDH] 
                                     + VV[2] * H[JR-1 + (J+2-1)*LDH] );
                    H[JR-1+(J-1)*LDH]   = H[JR-1+(J-1)*LDH] - TEMP;
                    H[JR-1+(J+1-1)*LDH] = H[JR-1+(J+1-1)*LDH] - TEMP*VV[1];
                    H[JR-1+(J+2-1)*LDH] = H[JR-1+(J+2-1)*LDH] - TEMP*VV[2];
                };

                for (i_type JR = IFRSTM; JR <= J + 2; ++JR)
                {
                    TEMP    = TAU * ( T[JR-1+(J-1)*LDT] + VV[1] * T[JR-1+(J+1-1)*LDT] 
                                     + VV[2]* T[JR-1 + (J+2-1)*LDT] );
                    T[JR-1+(J-1)*LDT]   = T[JR-1+(J-1)*LDT] - TEMP;
                    T[JR-1+(J+1-1)*LDT] = T[JR-1+(J+1-1)*LDT] - TEMP*VV[1];
                    T[JR-1+(J+2-1)*LDT] = T[JR-1+(J+2-1)*LDT] - TEMP*VV[2];
                };

                if (ILZ)
                {
                    for (i_type JR = 1; JR <= N; ++JR)
                    {
                        TEMP    = TAU * ( Z[JR-1+(J-1)*LDZ] + VV[1] * Z[JR-1 + (J+1-1)*LDZ]
                                         + VV[2] * Z[JR-1 + (J+2-1)*LDZ] );
                        Z[JR-1+(J-1)*LDZ]   = Z[JR-1+(J-1)*LDZ] - TEMP;
                        Z[JR-1+(J+1-1)*LDZ] = Z[JR-1+(J+1-1)*LDZ] - TEMP*VV[1];
                        Z[JR-1+(J+2-1)*LDZ] = Z[JR-1+(J+2-1)*LDZ] - TEMP*VV[2];
                    };
                };

                T[J+1-1+ (J-1)*LDT] = V(0.0);
                T[J+2-1+ (J-1)*LDT] = V(0.0);
            };

            // Last elements: Use Givens rotations
            // Rotations from the left
            i_type J        = ILAST - 1;
            TEMP            = H[J-1 +(J-1-1)*LDH];
            lapack::lartg<V>(TEMP, H[J+1-1 + (J-1-1)*LDH], &C, &S, &H[J-1+(J-1-1)*LDH] );
            H[J+1-1 +(J-1-1)*LDH] = V(0.0);

            for (i_type JC = J; JC <= ILASTM; ++JC)
            {
                TEMP        = C*H[J-1+(JC-1)*LDH] + S*H[J+1-1+(JC-1)*LDH];
                H[J+1-1+ (JC-1)*LDH]= -S*H[J-1+(JC-1)*LDH] + C*H[J+1-1+(JC-1)*LDH];
                H[J-1+ (JC-1)*LDH]  = TEMP;

                TEMP2       = C*T[J-1+(JC-1)*LDT] + S*T[J+1-1+(JC-1)*LDT];
                T[J+1-1 +(JC-1)*LDT]= -S*T[J-1+(JC-1)*LDT] + C*T[J+1-1+(JC-1)*LDT];
                T[J-1+(JC-1)*LDT]   = TEMP2;
            };

            if (ILQ)
            {
                for (i_type JR = 1; JR <= N; ++JR)
                {
                    TEMP    = C*Q[JR-1+(J-1)*LDQ] + S*Q[JR-1+(J+1-1)*LDQ];
                    Q[JR-1+(J+1-1)*LDQ]     = -S*Q[JR-1+(J-1)*LDQ] + C*Q[JR-1+(J+1-1)*LDQ];
                    Q[JR-1+(J-1)*LDQ]       = TEMP;
                };
            };

            // Rotations from the right.
            TEMP            = T[J+1-1+(J+1-1)*LDT];
            lapack::lartg<V>( TEMP, T[J+1-1+(J-1)*LDT], &C, &S, &T[J+1-1+(J+1-1)*LDT] );
            T[J+1-1+(J-1)*LDT]  = V(0.0);

            for (i_type JR = IFRSTM; JR <= ILAST; ++JR)
            {
                TEMP        = C*H[JR-1+(J+1-1)*LDH] + S*H[JR-1+(J-1)*LDH];
                H[JR-1+ (J-1)*LDH]  = -S*H[JR-1+(J+1-1)*LDH] + C*H[JR-1+(J-1)*LDH];
                H[JR-1+ (J+1-1)*LDH]= TEMP;
            };

            for (i_type JR = IFRSTM; JR <= ILAST - 1; ++JR)
            {
                TEMP        = C*T[JR-1+(J+1-1)*LDT] + S*T[JR-1+(J-1)*LDT];
                T[JR-1+(J-1)*LDT]   = -S*T[JR-1+(J+1-1)*LDT] + C*T[JR-1+(J-1)*LDT];
                T[JR-1+(J+1-1)*LDT] = TEMP;
            };

            if ( ILZ )
            {
                for (i_type JR = 1; JR <= N; ++JR)
                {
                    TEMP    = C*Z[JR-1+(J+1-1)*LDZ] + S*Z[JR-1+(J-1)*LDZ];
                    Z[JR-1+(J-1)*LDZ]   = -S*Z[JR-1+(J+1-1)*LDZ] + C*Z[JR-1+(J-1)*LDZ];
                    Z[JR-1+(J+1-1)*LDZ] = TEMP;
                };
            };
        };

        // End of iteration loop
      lab_350:
        ;
    };

    // Drop-through = non-convergence
    INFO            = ILAST;
    goto lab_420;

    // Successful completion of all QZ steps
  lab_380:
  
    // Set Eigenvalues 1:ILO-1
    for (i_type J = 1; J <= ILO - 1; ++J)
    {
        if (T[J-1 +(J-1)*LDT] <= V(0.0) )
        {
            if (ILSCHR)
            {
                for (i_type JR = 1; JR <= J; ++JR)
                {
                    H[JR-1 +(J-1)*LDH] = -H[JR-1 +(J-1)*LDH];
                    T[JR-1 +(J-1)*LDT] = -T[JR-1 +(J-1)*LDT];
                };
            }
            else
            {
               H[J-1 +(J-1)*LDH]    = -H[J-1 +(J-1)*LDH];
               T[J-1 +(J-1)*LDT]    = -T[J-1 +(J-1)*LDT];
            };

            if (ILZ)
            {
                for (i_type JR = 1; JR <= N; ++JR)
                {
                    Z[JR-1 +(J-1)*LDZ]  = -Z[JR-1 +(J-1)*LDZ];
                };
            };
        };

        ALPHAR[J-1]     = H[J-1 +(J-1)*LDH];
        ALPHAI[J-1]     = V(0.0);
        BETA[J-1]       = T[J-1 +(J-1)*LDT];
    };
    
    // Normal Termination
    INFO        = 0;

    // Exit (other than argument error) -- return optimal workspace size

  lab_420:
    WORK[1]     = V(N);
      
    return;
};

template BLAS_EXT_EXPORT void
hgeqz2<d_type>(const char* job, const char* compq, const char* compz, i_type n, i_type ilo, i_type ihi,
    d_type* h, i_type ldh, d_type* t, i_type ldt, d_type* alphar, d_type* alphai, d_type* beta, d_type* q, 
    i_type ldq, d_type* z, i_type ldz, d_type* work, i_type lwork, i_type& info);

template BLAS_EXT_EXPORT void
hgeqz2<s_type>(const char* job, const char* compq, const char* compz, i_type n, i_type ilo, i_type ihi,
    s_type* h, i_type ldh, s_type* t, i_type ldt, s_type* alphar, s_type* alphai, s_type* beta, s_type* q, 
    i_type ldq, s_type* z, i_type ldz, s_type* work, i_type lwork, i_type& info);


}};

