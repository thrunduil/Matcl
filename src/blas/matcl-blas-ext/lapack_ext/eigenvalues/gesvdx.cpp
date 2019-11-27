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

#if 0
TODO: not finished

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-blas-lapack/blas/details/blas_utils.h"

#include <algorithm>

namespace matcl { namespace lapack
{

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::gesvdx(const char* JOBU, const char* JOBVT, const char* RANGE, i_type M, i_type N, V* A, 
       i_type LDA, typename details::real_type<V>::type VL, typename details::real_type<V>::type VU, 
       i_type IL, i_type IU, i_type& NS, typename details::real_type<V>::type* S, V* U, i_type LDU, 
       V* VT, i_type LDVT, V* WORK, i_type LWORK, typename details::real_type<V>::type* RWORK, 
       i_type* IWORK, i_type& INFO)
{
    using VR        = typename details::real_type<V>::type;
    bool is_complex = details::is_complex<V>::value;

    V ZERO          = V(0.0);
    VR RZERO        = VR(0.0);
    V ONE           = V(1.0);
    VR RONE         = VR(1.0);

    // Test the input arguments.
    NS              = 0;
    INFO            = 0;
    VR ABSTOL       = 2*lamch<VR>("S");
    i_type LQUERY   = ( LWORK == -1 );
    i_type MINMN    = std::min( M, N );

    bool WANTU      = JOBU[0] == 'V' || JOBU[0] == 'v';
    bool WANTVT     = JOBVT[0] == 'V' || JOBVT[0] == 'v';

    const char* JOBZ;

    if (WANTU || WANTVT )
        JOBZ        = "V";
    else
        JOBZ        = "N";

    bool ALLS       = RANGE[0] == 'A' || RANGE[0] == 'a';
    bool VALS       = RANGE[0] == 'V' || RANGE[0] == 'v';
    bool INDS       = RANGE[0] == 'I' || RANGE[0] == 'i';

    if ( !(JOBU[0] == 'V' || JOBU[0] == 'v') && !(JOBU[0] == 'N' || JOBU[0] == 'n') )
        INFO        = -1;
    else if (!(JOBVT[0] == 'V' || JOBVT[0] == 'v') && !(JOBVT[0] == 'N' || JOBVT[0] == 'n') )
        INFO        = -2;
    else if ( !( ALLS || VALS || INDS ) )
        INFO        = -3;
    else if ( M < 0 )
        INFO        = -4;
    else if ( N < 0 )
        INFO        = -5;
    else if ( M > LDA )
        INFO        = -7;
    else if ( MINMN > 0 )
    {
        if ( VALS )
        {
            if ( VL < RZERO )
                INFO    = -8;
            else if ( VU <= VL )
                INFO    = -9;
        }
        else if ( INDS )
        {
            if ( IL < 1 || IL > std::max( 1, MINMN ) )
                INFO    = -10;
            else if ( IU < std::min( MINMN, IL ) || IU > MINMN )
                INFO    = -11;
        };
        
        if ( INFO == 0 )
        {
            if ( WANTU && LDU < M )
                INFO    = -15;
            else if ( WANTVT && LDVT < MINMN )
                INFO    = -16;
        };
    };

    // Compute workspace
    // (Note: Comments in the code beginning "Workspace:" describe the
    // minimal amount of workspace needed at that point in the code,
    // as well as the preferred amount for good performance.
    // NB refers to the optimal block size for the immediately
    // following subroutine, as returned by ilaenv.)

    i_type MINWRK;
    i_type MAXWRK;
    i_type MNTHR;

    const char* name_GESVD;
    const char* name_GEQRF;
    const char* name_GEBRD;
    const char* name_GELQF;

    if (std::is_same<V,d_type>::value)
    {
        name_GESVD  = "DGESVD";
        name_GEQRF  = "DGEQRF";
        name_GEBRD  = "DGEBRD";
        name_GELQF  = "DGELQF";
    }
    else if (std::is_same<V,s_type>::value)
    {
        name_GESVD  = "SGESVD";
        name_GEQRF  = "SGEQRF";
        name_GEBRD  = "SGEBRD";
        name_GELQF  = "SGELQF";
    }

    if ( INFO == 0 )
    {
        MINWRK              = 1;
        MAXWRK              = 1;

        if (MINMN > 0)
        {
            if ( M >= N )
            {
                char type[2];
                type[0]         = JOBU[0];
                type[1]         = JOBVT[0];

                MNTHR           = ilaenv( 6, name_GESVD, type, M, N, 0, 0 );

                if ( M >= MNTHR )
                {
                    // Path 1 (M much larger than N)
                    if (is_complex == false)
                    {
                        MAXWRK = N*(N*2+16) +  N*ilaenv( 1, name_GEQRF, " ", M, N, -1, -1 );
                        MAXWRK = std::max( MAXWRK, N*(N*2+21) + 2*N*ilaenv( 1, name_GEBRD, " ", N, N, -1, -1 ) );
                        MINWRK = N*(N*2+22);
                    }
                    else
                    {
                        MAXWRK = N + N*ilaenv( 1, name_GEQRF, " ", M, N, -1, -1 );
                        MAXWRK = std::max( MAXWRK, N*N + N + 2*N*ilaenv( 1, name_GEBRD, " ", N, N, -1, -1 ) );
                        MINWRK = N*(N+4);
                    };
                }
                else
                {
                    // Path 2 (M at least N, but not much larger)
                    if (is_complex == false)
                    {
                        MAXWRK = N*(N*2+20) + ( M+N )*ilaenv( 1, name_GEBRD, " ", M, N, -1, -1 );
                        MINWRK = N*(N*2+21) + M;
                    }
                    else
                    {
                        MAXWRK = 2*N + ( M+N )*ilaenv( 1, name_GEBRD, " ", M, N, -1, -1 );
                        MINWRK = 2*N + M;
                    };
                };
            }
            else
            {
                char type[2];
                type[0]         = JOBU[0];
                type[1]         = JOBVT[0];

                MNTHR           = ilaenv( 6, name_GESVD, type, M, N, 0, 0 );

                if ( N >= MNTHR )
                {
                    // Path 1t (N much larger than M)
                    if (is_complex == false)
                    {
                        MAXWRK = M*(M*2+16) + M*ilaenv( 1, name_GELQF, " ", M, N, -1, -1 );
                        MAXWRK = std::max( MAXWRK, M*(M*2+21) + 2*M*ilaenv( 1, name_GEBRD, " ", M, M, -1, -1 ) );
                        MINWRK = M*(M*2+22);
                    }
                    else
                    {
                        MAXWRK = M + M*ilaenv( 1, name_GELQF, " ", M, N, -1, -1 );
                        MAXWRK = std::max( MAXWRK, M*M + M + 2*M*ilaenv( 1, name_GEBRD, " ", M, M, -1, -1 ) );
                        MINWRK = M*(M+4);
                    };
                }
                else
                {
                    //Path 2t (N greater than M, but not much larger)
                    if (is_complex == false)
                    {
                        MAXWRK = M*(M*2+20) + ( M+N )*ilaenv( 1, name_GEBRD, " ", M, N, -1, -1 );
                        MINWRK = M*(M*2+21) + N;
                    }
                    else
                    {
                        MAXWRK = M*(M*2+19) + ( M+N )*ilaenv( 1, name_GEBRD, " ", M, N, -1, -1 );
                        MINWRK = 2*M + N;
                    };
                };
            };
        };
        
        MAXWRK      = std::max( MAXWRK, MINWRK );
        WORK[1-1]   = VR(MAXWRK);

        if ( LWORK < MINWRK && !LQUERY )
            INFO    = -19;
    };

    if ( INFO != 0 )
        return;
    else if ( LQUERY )
        return;

    // Quick return if possible
    if ( M == 0 || N == 0 )
        return;

    // Set singular values indices accord to RANGE.
    const char* RNGTGK;
    i_type ILTGK;
    i_type IUTGK;

    if ( ALLS )
    {
        RNGTGK  = "I";
        ILTGK   = 1;
        IUTGK   = std::min( M, N );
    }
    else if ( INDS )
    {
        RNGTGK  = "I";
        ILTGK   = IL;
        IUTGK   = IU;
    }
    else      
    {
        RNGTGK  = "V";
        ILTGK   = 0;
        IUTGK   = 0;
    };

    // Get machine constants
    VR EPS      = lamch<VR>("P");
    VR SMLNUM   = sqrt(lamch<VR>("S") ) / EPS;
    VR BIGNUM   = RONE / SMLNUM;

    // Scale A if max element outside range [SMLNUM,BIGNUM]
    VR ANRM     = lange<V>("M", M, N, A, LDA, nullptr);
    i_type ISCL = 0;

    if ( ANRM > RZERO && ANRM < SMLNUM )
    {
        ISCL    = 1;
        lapack::lascl<V>("G", 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
    }
    else if ( ANRM > BIGNUM )
    {
        ISCL    = 1;
        lapack::lascl<V>("G", 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
    };

    if ( M >= N )
    {
        // A has at least as many rows as columns. If A has sufficiently
        // more rows than columns, first reduce A using the QR
        // decomposition.
        if ( M >= MNTHR )
        {
            // Path 1 (M much larger than N):
            // A = Q * R = Q * ( QB * B * PB**T )
            //           = Q * ( QB * ( UB * S * VB**T ) * PB**T )
            // U = Q * QB * UB; V**T = VB**T * PB**T
            // 
            // Compute A=Q*R
            // (Workspace: need 2*N, prefer N+N*NB)

            i_type ITAU     = 1;
            i_type ITEMP    = ITAU + N;
            lapack::geqrf<V>( M, N, A, LDA, WORK + ITAU-1, WORK + ITEMP-1, LWORK-ITEMP+1, &INFO );
              
            // Copy R into WORK and bidiagonalize it:
            // (Workspace: need N*N+5*N, prefer N*N+4*N+2*N*NB)
            i_type IQRF     = ITEMP;
            i_type ID       = IQRF + N*N;
            i_type IE       = ID + N;
            i_type ITAUQ    = IE + N;
            i_type ITAUP    = ITAUQ + N;
            ITEMP           = ITAUP + N;

            lapack::lacpy<V>("U", N, N, A, LDA, WORK + IQRF-1, N );
            lapack::laset<V>("L", N-1, N-1, ZERO, ZERO, WORK + IQRF+1-1, N );
            lapack::gebrd<V>( N, N, WORK + IQRF-1, N, WORK + ID-1, WORK + IE-1, 
                WORK + ITAUQ-1, WORK + ITAUP-1, WORK + ITEMP-1, LWORK-ITEMP+1, INFO );

/*
*  
*           Copy R into WORK and bidiagonalize it:
*           (Workspace: need N*N+3*N, prefer N*N+N+2*N*NB)
*   
            ITAUQ = ITEMP + N*N
            ITAUP = ITAUQ + N
            ITEMP = ITAUP + N
            ID = 1   
            IE = ID + N
            ITGKZ = IE + N
            CALL CLACPY( 'U', N, N, A, LDA, WORK( IQRF ), N )
            CALL CLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                   WORK( IQRF+1 ), N )
            CALL CGEBRD( N, N, WORK( IQRF ), N, RWORK( ID ), 
     $                   RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ),
     $                   WORK( ITEMP ), LWORK-ITEMP+1, INFO )
            ITEMP = ITGKZ + N*(N*2+1)
*/
            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 14*N + 2*N*(N+1))            
            i_type ITGKZ    = ITEMP;
            ITEMP           = ITGKZ + 2*N*N + 2*N;
            
            lapack::bdsvdx<V>("U", JOBZ, RNGTGK, N, WORK + ID-1, WORK + IE-1, VL, VU, ILTGK, IUTGK, NS, S,
                              WORK + ITGKZ-1, N*2, WORK + ITEMP-1, IWORK, INFO);
            
            // If needed, compute left singular vectors.
            if ( WANTU )
            {                
                i_type J    = ITGKZ;
                for (i_type I = 1; I <= NS; ++I)
                {
                    copy<V>(N, WORK + J-1, 1, U + 1-1 + (I-1)*LDU, 1 );
                    J       = J + N*2;
                };
                             
                lapack::laset<V>("A", M-N, NS, ZERO, ZERO, U + N+1-1 + (1-1)*LDU, LDU );
    
                // Call DORMBR to compute QB*UB.
                // (Workspace in WORK[ITEMP-1]: need N, prefer N*NB)
                lapack::ormbr<V>("Q", "L", "N", N, NS, N, WORK + IQRF-1, N, WORK + ITAUQ-1, U, LDU, 
                                 WORK + ITEMP-1, LWORK-ITEMP+1, INFO );

                // Call DORMQR to compute Q*(QB*UB).
                // (Workspace in WORK[ITEMP-1]: need N, prefer N*NB)
                lapack::ormqr<V>("L", "N", M, NS, N, A, LDA, WORK + ITAU-1, U, LDU, WORK + ITEMP-1,
                                LWORK-ITEMP+1, INFO );
            };

            // If needed, compute right singular vectors.
            if ( WANTVT) 
            {
                i_type J        = ITGKZ + N;
                for (i_type I = 1; I <= NS; ++I)
                {
                    copy<V>(N, WORK + J-1, 1, VT + I-1 +(1-1)*LDVT, LDVT );
                    J           = J + N*2;
                };
    
                // Call DORMBR to compute VB**T * PB**T
                // (Workspace in WORK[ITEMP-1]: need N, prefer N*NB)
                lapack::ormbr<V>("P", "R", "T", NS, N, N, WORK + IQRF-1, N, WORK + ITAUP-1, VT, LDVT, 
                                 WORK + ITEMP-1, LWORK-ITEMP+1, INFO );
            };            

            /*
*
*           Solve eigenvalue problem TGK*Z=Z*S.
*           (Workspace: need 2*N*N+14*N)          
*                   
            CALL SBDSVDX( 'U', JOBZ, RNGTGK, N, RWORK( ID ),
     $                    RWORK( IE ), VL, VU, ILTGK, IUTGK, NS, S,
     $                    RWORK( ITGKZ ), N*2, RWORK( ITEMP ),
     $                    IWORK, INFO)
*
*           If needed, compute left singular vectors.
*
            IF( WANTU ) THEN
               K = ITGKZ
               DO I = 1, NS
                  DO J = 1, N
                     U( J, I ) = CMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  END DO
                  K = K + N
               END DO
               CALL CLASET( 'A', M-N, N, CZERO, CZERO, U( N+1,1 ), LDU )
*
*              Call CUNMBR to compute QB*UB.
*              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
*
               CALL CUNMBR( 'Q', 'L', 'N', N, NS, N, WORK( IQRF ), N, 
     $                      WORK( ITAUQ ), U, LDU, WORK( ITEMP ), 
     $                      LWORK-ITEMP+1, INFO )
*
*              Call CUNMQR to compute Q*(QB*UB).
*              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
*
               CALL CUNMQR( 'L', 'N', M, NS, N, A, LDA, 
     $                      WORK( ITAU ), U, LDU, WORK( ITEMP ),
     $                      LWORK-ITEMP+1, INFO )
            END IF  
*      
*           If needed, compute right singular vectors.
*
            IF( WANTVT) THEN
               K = ITGKZ + N
               DO I = 1, NS
                  DO J = 1, N
                     VT( I, J ) = CMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  END DO
                  K = K + N
               END DO
*
*              Call CUNMBR to compute VB**T * PB**T
*              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
*
               CALL CUNMBR( 'P', 'R', 'C', NS, N, N, WORK( IQRF ), N, 
     $                      WORK( ITAUP ), VT, LDVT, WORK( ITEMP ),
     $                      LWORK-ITEMP+1, INFO )
            END IF
*/
        }
        else
        {
            // Path 2 (M at least N, but not much larger)
            // Reduce A to bidiagonal form without QR decomposition
            // A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            // U = QB * UB; V**T = VB**T * PB**T
            // 
            // Bidiagonalize A
            // (Workspace: need 4*N+M, prefer 4*N+(M+N)*NB)

            i_type ID       = 1;
            i_type IE       = ID + N;
            i_type ITAUQ    = IE + N;
            i_type ITAUP    = ITAUQ + N;
            i_type ITEMP    = ITAUP + N;

            lapack::gebrd<V>( M, N, A, LDA, WORK + ID-1, WORK + IE-1, WORK + ITAUQ-1, WORK + ITAUP-1, 
                             WORK + ITEMP-1, LWORK-ITEMP+1, INFO );

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 14*N + 2*N*(N+1))           
            i_type ITGKZ    = ITEMP;
            ITEMP           = ITGKZ + 2*N*N + 2*N;
            lapack::bdsvdx<V>("U", JOBZ, RNGTGK, N, WORK + ID-1, WORK + IE-1, VL, VU, ILTGK, IUTGK, NS, 
                              S, WORK + ITGKZ-1, N*2, WORK + ITEMP-1, IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU )
            {
                i_type J        = ITGKZ;
                for (i_type I = 1; I <= NS; ++I)
                {
                    copy<V>( N, WORK + J-1, 1, U + 1-1 +(I-1)*LDU, 1 );
                    J           = J + N*2;
                };

                lapack::laset<V>("A", M-N, NS, ZERO, ZERO, U + N+1-1 + (1-1)*LDU, LDU );

                // Call DORMBR to compute QB*UB.
                // (Workspace in WORK[ITEMP-1]: need N, prefer N*NB)
                i_type IERR;
                lapack::ormbr<V>("Q", "L", "N", M, NS, N, A, LDA, WORK + ITAUQ-1, U, LDU, WORK + ITEMP-1, 
                    LWORK-ITEMP+1, IERR );
            };

            // If needed, compute right singular vectors.

            if ( WANTVT)
            {
                i_type J        = ITGKZ + N;
                for (i_type I = 1; I <= NS; ++I)
                {
                    copy<V>( N, WORK + J-1, 1, VT + I-1 + (1-1)*LDVT, LDVT);
                    J           = J + N*2;
                };
    
                // Call DORMBR to compute VB**T * PB**T
                // (Workspace in WORK[ITEMP-1]: need N, prefer N*NB)
                i_type IERR;
                lapack::ormbr<V>("P", "R", "T", NS, N, N, A, LDA, WORK + ITAUP-1, VT, LDVT, 
                                 WORK + ITEMP-1, LWORK-ITEMP+1, IERR );
            };
        };

        /*
         IF( M.GE.MNTHR ) THEN
*
         ELSE
*
*           Path 2 (M at least N, but not much larger)
*           Reduce A to bidiagonal form without QR decomposition
*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
*           U = QB * UB; V**T = VB**T * PB**T
*
*           Bidiagonalize A
*           (Workspace: need 2*N+M, prefer 2*N+(M+N)*NB)
*
            ITAUQ = 1
            ITAUP = ITAUQ + N
            ITEMP = ITAUP + N          
            ID = 1
            IE = ID + N
            ITGKZ = IE + N
            CALL CGEBRD( M, N, A, LDA, RWORK( ID ), RWORK( IE ), 
     $                   WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ),
     $                   LWORK-ITEMP+1, INFO )
            ITEMP = ITGKZ + N*(N*2+1)
*
*           Solve eigenvalue problem TGK*Z=Z*S.
*           (Workspace: need 2*N*N+14*N)          
*                     
            CALL SBDSVDX( 'U', JOBZ, RNGTGK, N, RWORK( ID ),
     $                    RWORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, 
     $                    RWORK( ITGKZ ), N*2, RWORK( ITEMP ), 
     $                    IWORK, INFO)
*
*           If needed, compute left singular vectors.
*
            IF( WANTU ) THEN
               K = ITGKZ
               DO I = 1, NS
                  DO J = 1, N              
                     U( J, I ) = CMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  END DO
                  K = K + N
               END DO
               CALL CLASET( 'A', M-N, N, CZERO, CZERO, U( N+1,1 ), LDU )
*
*              Call CUNMBR to compute QB*UB.
*              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
*   
               CALL CUNMBR( 'Q', 'L', 'N', M, NS, N, A, LDA, 
     $                      WORK( ITAUQ ), U, LDU, WORK( ITEMP ), 
     $                      LWORK-ITEMP+1, IERR )
            END IF  
*      
*           If needed, compute right singular vectors.
*
            IF( WANTVT) THEN
               K = ITGKZ + N
               DO I = 1, NS
                  DO J = 1, N
                     VT( I, J ) = CMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  END DO
                  K = K + N
               END DO
*
*              Call CUNMBR to compute VB**T * PB**T
*              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
*
               CALL CUNMBR( 'P', 'R', 'C', NS, N, N, A, LDA, 
     $                      WORK( ITAUP ), VT, LDVT, WORK( ITEMP ),
     $                      LWORK-ITEMP+1, IERR )
            END IF
         END IF                             
*/
    }
    else
    {
        // A has more columns than rows. If A has sufficiently more
        // columns than rows, first reduce A using the LQ decomposition.

        if ( N >= MNTHR )
        {
            // Path 1t (N much larger than M):
            // A = L * Q = ( QB * B * PB**T ) * Q 
            //           = ( QB * ( UB * S * VB**T ) * PB**T ) * Q
            // U = QB * UB ; V**T = VB**T * PB**T * Q
            // 
            // Compute A=L*Q
            // (Workspace: need 2*M, prefer M+M*NB)

            i_type ITAU     = 1;
            i_type ITEMP    = ITAU + M;

            lapack::gelqf(M, N, A, LDA, WORK + ITAU-1, WORK + ITEMP-1, LWORK-ITEMP+1, INFO );

            // Copy L into WORK and bidiagonalize it:
            //     (Workspace in WORK[ITEMP-1]: need M*M+5*N, prefer M*M+4*M+2*M*NB)
            i_type ILQF     = ITEMP;
            i_type ID       = ILQF + M*M;
            i_type IE       = ID + M;
            i_type ITAUQ    = IE + M;
            i_type ITAUP    = ITAUQ + M;
            ITEMP           = ITAUP + M;

            lapack::lacpy<V>("L", M, M, A, LDA, WORK + ILQF-1, M );
            lapack::laset<V>("U", M-1, M-1, ZERO, ZERO, WORK + ILQF+M-1, M );
            lapack::gebrd<V>( M, M, WORK + ILQF-1, M, WORK + ID-1, WORK + IE-1, WORK + ITAUQ-1, 
                             WORK + ITAUP-1, WORK + ITEMP-1, LWORK-ITEMP+1, INFO );

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*M*M+14*M)          
            i_type ITGKZ    = ITEMP;
            ITEMP           = ITGKZ + 2*M*M + 2*M;

            lapack::bdsvdx<V>("U", JOBZ, RNGTGK, M, WORK + ID-1, WORK + IE-1, VL, VU, ILTGK, IUTGK, NS, 
                              S, WORK + ITGKZ-1, M*2, WORK + ITEMP-1, IWORK, INFO);

            // If needed, compute left singular vectors.
            if ( WANTU )
            {
                i_type J    = ITGKZ;
                for (i_type I = 1; I <= NS; ++I)
                {
                    copy<V>( M, WORK + J-1, 1, U + 1-1 +(I-1)*LDU, 1 );
                    J       = J + M*2;
                };
    
                // Call DORMBR to compute QB*UB.
                // (Workspace in WORK[ITEMP-1]: need M, prefer M*NB)
                lapack::ormbr<V>( "Q", "L", "N", M, NS, M, WORK + ILQF-1, M, WORK + ITAUQ-1, U, LDU, 
                                 WORK + ITEMP-1, LWORK-ITEMP+1, INFO );
            };

            // If needed, compute right singular vectors.
            if ( WANTVT)
            {
                i_type J    = ITGKZ + M;
                for (i_type I = 1; I <= NS; ++I)
                {
                    copy<V>( M, WORK + J-1, 1, VT + I-1 + (1-1)*LDVT, LDVT);
                    J       = J + M*2;
                };
                lapack::laset<V>("A", M, N-M, ZERO, ZERO, VT + 1-1 + (M+1-1)*LDVT, LDVT );
    
                // Call DORMBR to compute (VB**T)*(PB**T)
                // (Workspace in WORK[ITEMP-1]: need M, prefer M*NB)
                lapack::ormbr<V>( "P", "R", "T", NS, M, M, WORK + ILQF-1, M, WORK + ITAUP-1, 
                                 VT, LDVT, WORK + ITEMP-1, LWORK-ITEMP+1, INFO );
    
                // Call DORMLQ to compute ((VB**T)*(PB**T))*Q.
                // (Workspace in WORK[ITEMP-1]: need M, prefer M*NB)
                lapack::ormlq("R", "N", NS, N, M, A, LDA, WORK + ITAU-1, VT, LDVT, WORK + ITEMP-1,
                    LWORK-ITEMP+1, INFO );
            };
        }
        else
        {
            // Path 2t (N greater than M, but not much larger)
            // Reduce to bidiagonal form without LQ decomposition
            // A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            // U = QB * UB; V**T = VB**T * PB**T           
            // 
            // Bidiagonalize A
            // (Workspace: need 4*M+N, prefer 4*M+(M+N)*NB)

            i_type ID       = 1;
            i_type IE       = ID + M;
            i_type ITAUQ    = IE + M;
            i_type ITAUP    = ITAUQ + M;
            i_type ITEMP    = ITAUP + M;

            lapack::gebrd<V>(M, N, A, LDA, WORK + ID-1, WORK + IE-1, WORK + ITAUQ-1, 
                             WORK + ITAUP-1, WORK + ITEMP-1, LWORK-ITEMP+1, INFO );

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*M*M+14*M)          

            i_type ITGKZ    = ITEMP;
            ITEMP           = ITGKZ + 2*M*M + 2*M;
            lapack::bdsvdx<V>("L", JOBZ, RNGTGK, M, WORK + ID-1, WORK + IE-1, VL, VU, ILTGK, IUTGK, 
                              NS, S, WORK + ITGKZ-1, M*2, WORK + ITEMP-1, IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU )
            {
                i_type J    = ITGKZ;
                for (i_type I = 1; I <= NS; ++I)
                {
                    copy<V>( M, WORK + J-1, 1, U + 1-1 + (I-1)*LDU, 1 );
                    J       = J + M*2;
                };
    
                // Call DORMBR to compute QB*UB.
                // (Workspace in WORK[ITEMP-1]: need M, prefer M*NB)
                lapack::ormbr<V>("Q", "L", "N", M, NS, N, A, LDA, WORK + ITAUQ-1, U, LDU, 
                                 WORK + ITEMP-1, LWORK-ITEMP+1, INFO );
            };

            // If needed, compute right singular vectors.
            if ( WANTVT)
            {
                i_type J        = ITGKZ + M;
                for (i_type I = 1; I <= NS; ++I)
                {
                    copy<V>( M, WORK + J-1, 1, VT + I-1 + (1-1)*LDVT, LDVT);
                    J           = J + M*2;
                };
            
                lapack::laset<V>("A", M, N-M, ZERO, ZERO, VT + 1-1 + (M+1-1)*LDVT, LDVT );

                // Call DORMBR to compute VB**T * PB**T
                // (Workspace in WORK[ITEMP-1]: need M, prefer M*NB)
                lapack::ormbr<V>("P", "R", "T", NS, N, M, A, LDA, WORK+ITAUP-1, VT, LDVT, 
                                 WORK + ITEMP-1, LWORK-ITEMP+1, INFO );
            };
        };

        /*
*
*        A has more columns than rows. If A has sufficiently more
*        columns than rows, first reduce A using the LQ decomposition.
*
         IF( N.GE.MNTHR ) THEN
*
*           Path 1t (N much larger than M):
*           A = L * Q = ( QB * B * PB**T ) * Q 
*                     = ( QB * ( UB * S * VB**T ) * PB**T ) * Q
*           U = QB * UB ; V**T = VB**T * PB**T * Q
*
*           Compute A=L*Q
*           (Workspace: need 2*M, prefer M+M*NB)
*
            ITAU = 1
            ITEMP = ITAU + M
            CALL CGELQF( M, N, A, LDA, WORK( ITAU ), WORK( ITEMP ),
     $                   LWORK-ITEMP+1, INFO )

*           Copy L into WORK and bidiagonalize it:
*           (Workspace in WORK( ITEMP ): need M*M+3*M, prefer M*M+M+2*M*NB)
*
            ILQF = ITEMP       
            ITAUQ = ILQF + M*M
            ITAUP = ITAUQ + M
            ITEMP = ITAUP + M
            ID = 1
            IE = ID + M
            ITGKZ = IE + M
            CALL CLACPY( 'L', M, M, A, LDA, WORK( ILQF ), M )
            CALL CLASET( 'U', M-1, M-1, CZERO, CZERO, 
     $                   WORK( ILQF+M ), M )
            CALL CGEBRD( M, M, WORK( ILQF ), M, RWORK( ID ),
     $                   RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), 
     $                   WORK( ITEMP ), LWORK-ITEMP+1, INFO )
            ITEMP = ITGKZ + M*(M*2+1)
*
*           Solve eigenvalue problem TGK*Z=Z*S.
*           (Workspace: need 2*M*M+14*M)          
*
            CALL SBDSVDX( 'U', JOBZ, RNGTGK, M, RWORK( ID ),
     $                    RWORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, 
     $                    RWORK( ITGKZ ), M*2, RWORK( ITEMP ), 
     $                    IWORK, INFO)
*
*           If needed, compute left singular vectors.
*
            IF( WANTU ) THEN
               K = ITGKZ
               DO I = 1, NS
                  DO J = 1, M
                     U( J, I ) = CMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  END DO
                  K = K + M
               END DO
*
*              Call CUNMBR to compute QB*UB.
*              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
*
               CALL CUNMBR( 'Q', 'L', 'N', M, NS, M, WORK( ILQF ), M, 
     $                      WORK( ITAUQ ), U, LDU, WORK( ITEMP ), 
     $                      LWORK-ITEMP+1, INFO )
            END IF  
*      
*           If needed, compute right singular vectors.
*
            IF( WANTVT) THEN
               K = ITGKZ + M
               DO I = 1, NS
                  DO J = 1, M
                     VT( I, J ) = CMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  END DO
                  K = K + M
               END DO
               CALL CLASET( 'A', M, N-M, CZERO, CZERO, 
     $                      VT( 1,M+1 ), LDVT )
*
*              Call CUNMBR to compute (VB**T)*(PB**T)
*              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
*
               CALL CUNMBR( 'P', 'R', 'C', NS, M, M, WORK( ILQF ), M, 
     $                      WORK( ITAUP ), VT, LDVT, WORK( ITEMP ),
     $                      LWORK-ITEMP+1, INFO )
*
*              Call CUNMLQ to compute ((VB**T)*(PB**T))*Q.
*              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
*
               CALL CUNMLQ( 'R', 'N', NS, N, M, A, LDA, 
     $                      WORK( ITAU ), VT, LDVT, WORK( ITEMP ),
     $                      LWORK-ITEMP+1, INFO )
            END IF  
         ELSE
*
*           Path 2t (N greater than M, but not much larger)
*           Reduce to bidiagonal form without LQ decomposition
*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
*           U = QB * UB; V**T = VB**T * PB**T           
*
*           Bidiagonalize A
*           (Workspace: need 2*M+N, prefer 2*M+(M+N)*NB)
*            
            ITAUQ = 1
            ITAUP = ITAUQ + M
            ITEMP = ITAUP + M
            ID = 1
            IE = ID + M
            ITGKZ = IE + M
            CALL CGEBRD( M, N, A, LDA, RWORK( ID ), RWORK( IE ), 
     $                   WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ),
     $                   LWORK-ITEMP+1, INFO )
            ITEMP = ITGKZ + M*(M*2+1)
*
*           Solve eigenvalue problem TGK*Z=Z*S.
*           (Workspace: need 2*M*M+14*M)          
*          
            CALL SBDSVDX( 'L', JOBZ, RNGTGK, M, RWORK( ID ), 
     $                    RWORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, 
     $                    RWORK( ITGKZ ), M*2, RWORK( ITEMP ), 
     $                    IWORK, INFO)
* 
*           If needed, compute left singular vectors.
*
            IF( WANTU ) THEN
               K = ITGKZ
               DO I = 1, NS
                  DO J = 1, M
                     U( J, I ) = CMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  END DO
                  K = K + M
               END DO
*
*              Call CUNMBR to compute QB*UB.
*              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
*
               CALL CUNMBR( 'Q', 'L', 'N', M, NS, N, A, LDA, 
     $                      WORK( ITAUQ ), U, LDU, WORK( ITEMP ), 
     $                      LWORK-ITEMP+1, INFO )
            END IF  
*      
*           If needed, compute right singular vectors.
*
            IF( WANTVT) THEN
               K = ITGKZ + M
               DO I = 1, NS
                  DO J = 1, M
                     VT( I, J ) = CMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  END DO
                  K = K + M
               END DO
               CALL CLASET( 'A', M, N-M, CZERO, CZERO, 
     $                      VT( 1,M+1 ), LDVT )
*
*              Call CUNMBR to compute VB**T * PB**T
*              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
*
               CALL CUNMBR( 'P', 'R', 'C', NS, N, M, A, LDA, 
     $                      WORK( ITAUP ), VT, LDVT, WORK( ITEMP ),
     $                      LWORK-ITEMP+1, INFO )
            END IF 
         END IF
*/
    };

    // Undo scaling if necessary
    if ( ISCL == 1 )
    {
        if ( ANRM > BIGNUM )
            lapack::lascl<V>( "G", 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, INFO );
        if ( ANRM < SMLNUM )
            lapack::lascl<V>( "G", 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, INFO );
    };

    // Return optimal workspace in WORK(1)
    WORK[1-1]   = VR( MAXWRK );
    return;
}

/*
template BLAS_EXT_EXPORT void
gesvdx<d_type>(const char* JOBU, const char* JOBVT, const char* RANGE, i_type M, i_type N, d_type* A, 
       i_type LDA, d_type VL, d_type VU, i_type IL, i_type IU, i_type& NS, d_type* S, d_type* U, i_type LDU, 
       d_type* VT, i_type LDVT, d_type* WORK, i_type LWORK, d_type* RWORK, i_type* IWORK, i_type& INFO);

template BLAS_EXT_EXPORT void
gesvdx<s_type>(const char* JOBU, const char* JOBVT, const char* RANGE, i_type M, i_type N, s_type* A, 
       i_type LDA, s_type VL, s_type VU, i_type IL, i_type IU, i_type& NS, s_type* S, s_type* U, i_type LDU, 
       s_type* VT, i_type LDVT, s_type* WORK, i_type LWORK, s_type* RWORK, i_type* IWORK, i_type& INFO);

template BLAS_EXT_EXPORT void
gesvdx<c_type>(const char* JOBU, const char* JOBVT, const char* RANGE, i_type M, i_type N, c_type* A, 
       i_type LDA, s_type VL, s_type VU, i_type IL, i_type IU, i_type& NS, s_type* S, c_type* U, i_type LDU, 
       c_type* VT, i_type LDVT, c_type* WORK, i_type LWORK, s_type* RWORK, i_type* IWORK, i_type& INFO);
*/

}};
#endif