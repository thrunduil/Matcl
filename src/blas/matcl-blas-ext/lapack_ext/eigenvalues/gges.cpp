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
#include "blas/matcl-blas-ext/lapack_ext/utils/utils.h"

#include <algorithm>

#pragma warning(push)
#pragma warning(disable:4100) //unreferenced formal parameter

namespace matcl { namespace lapack
{

template<class V, bool Is_compl = details::is_complex<V>::value>
struct gges2_alg
{
    using VR = typename details::real_type<V>::type;
    using VC = typename details::complex_type<V>::type;

    static void rescale_eigs(bool ILASCL, bool ILBSCL, i_type N, VC* ALPHA, V* BETA, VR SAFMAX, VR SAFMIN,
                    VR ANRMTO, VR ANRM, VR BNRMTO, VR BNRM, const V* A, i_type LDA, const V* B, i_type LDB)
    {
        // Check if unscaling would cause over/underflow, if so, rescale
        //   (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of
        //    B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I)

        VR* ALPHAR  = reinterpret_cast<VR*>(ALPHA);
        VR* ALPHAI  = ALPHAR + N;

        VR WORK[1];

        if (ILASCL)
        {        
            for (i_type I = 0; I < N; ++I)
            {
                if (ALPHAI[I] != VR(0) )
                {
                    if ( ALPHAR[I] / SAFMAX > (ANRMTO / ANRM )
                        || SAFMIN / ALPHAR[I] > ( ANRM / ANRMTO ) )
                    {
                        WORK[0]     = std::abs( A[I + I*LDA] / ALPHAR[I] );
                        BETA[I]     = BETA[I] * WORK[0];
                        ALPHAR[I]   = ALPHAR[I] * WORK[0];
                        ALPHAI[I]   = ALPHAI[I] * WORK[0];
                    }
                    else if ( ALPHAI[I] / SAFMAX > ANRMTO / ANRM
                             || SAFMIN / ALPHAI[I] > ANRM / ANRMTO )
                    {
                        WORK[0]     = std::abs( A[I + (I+1)*LDA] / ALPHAI[I] );
                        BETA[I]     = BETA[I] * WORK[0];
                        ALPHAR[I]   = ALPHAR[I] * WORK[0];
                        ALPHAI[I]   = ALPHAI[I] * WORK[0];
                    };
                };
            };
        };

        if (ILBSCL)
        {
            for (i_type I = 0; I < N; ++I)
            {
                if ( ALPHAI[I] != VR(0.0) )
                {
                    if ( BETA[I] / SAFMAX > BNRMTO / BNRM
                        || SAFMIN / BETA[I] > BNRM / BNRMTO )
                    {
                        WORK[0]     = std::abs( B[I + I*LDB] / BETA[I] );
                        BETA[I]     = BETA[I] * WORK[0];
                        ALPHAR[I]   = ALPHAR[I] * WORK[0];
                        ALPHAI[I]   = ALPHAI[I] * WORK[0];
                    };
                };
            };
        };
    };

    static void undo_scale_alpha(VR ANRMTO, VR ANRM, i_type N, VC* ALPHA)
    {
        VR* ALPHAR  = reinterpret_cast<VR*>(ALPHA);
        VR* ALPHAI  = ALPHAR + N;

        i_type IERR;
        lapack::lascl( "G", 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR );
        lapack::lascl( "G", 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR );
    }

    static void make_complex_alpha(VC* ALPHA, i_type N, V* WORK)
    {
        // make ALPHA complex

        VR* ALPHAR  = reinterpret_cast<VR*>(ALPHA);
        VR* ALPHAI  = ALPHAR + N;

        VR* WORK_R      = WORK;
        VR* WORK_I      = WORK_R + N;

        for (i_type I = 0; I < N; ++I)
            WORK_R[I]   = ALPHAR[I];

        for (i_type I = 0; I < N; ++I)
            WORK_I[I]   = ALPHAI[I];

        for (i_type I = 0; I < N; ++I)
            ALPHA[I]    = VC(WORK_R[I], WORK_I[I]);
    };

    static void hgeqz(const char* job, const char* compq, const char* compz, i_type n, i_type ilo, i_type ihi,
                V* h, i_type ldh, V* t, i_type ldt, VC* alpha, V* beta, V* q, i_type ldq, V* z, 
                i_type ldz, V* work, i_type lwork, VR* rwork, i_type& info)
    {
        (void)rwork;

        VR* ALPHAR  = reinterpret_cast<VR*>(alpha);
        VR* ALPHAI  = ALPHAR + n;

        lapack::hgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, ALPHAR, ALPHAI, beta, q, ldq, z, 
                ldz, work, lwork, info);
    };
};

template<class V>
struct gges2_alg<V,true>
{
    using VR = typename details::real_type<V>::type;

    static void rescale_eigs(bool ILASCL, bool ILBSCL, i_type N, V* ALPHA, V* BETA, VR SAFMAX, VR SAFMIN,
                    VR ANRMTO, VR ANRM, VR BNRMTO, VR BNRM, const V* A, i_type LDA, const V* B, i_type LDB)
    {};

    static void undo_scale_alpha(VR ANRMTO, VR ANRM, i_type N, V* ALPHA)
    {};

    static void make_complex_alpha(V* ALPHA, i_type N, V* WORK)
    {};

    static void hgeqz(const char* job, const char* compq, const char* compz, i_type n, i_type ilo, i_type ihi,
                V* h, i_type ldh, V* t, i_type ldt, V* alpha, V* beta, V* q, i_type ldq, V* z, 
                i_type ldz, V* work, i_type lwork, VR* rwork, i_type& info)
    {
        lapack::hgeqz<V>(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, 
                ldz, work, lwork, rwork, info);
    };
};


template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::gges2(const char *JOBVSL, const char *JOBVSR, i_type N, V* A, i_type LDA, 
    V* B, i_type LDB, typename details::complex_type<V>::type* ALPHA, V* BETA, V* VSL, i_type LDVSL, V* VSR,
    i_type LDVSR, V* WORK, i_type LWORK, i_type* IWORK, i_type LIWORK, i_type& INFO)
{
    using VR = typename details::real_type<V>::type;
    using VC = typename details::complex_type<V>::type;

    bool is_V_real  = details::is_complex<V>::value == false;

    // Decode the input arguments

    i_type  IJOBVL, IJOBVR, MINWRK, MAXWRK, IERR;
    bool    ILVSL, ILVSR;

    if ( JOBVSL[0] == 'N' || JOBVSL[0] == 'n')
    {
        IJOBVL  = 1;
        ILVSL   = false;
    }
    else if ( JOBVSL[0] == 'V' || JOBVSL[0] == 'v')
    {
        IJOBVL  = 2;
        ILVSL   = true;
    }
    else
    {
        IJOBVL  = -1;
        ILVSL   = false;
    };

    if ( JOBVSR[0] == 'N' || JOBVSR[0] == 'n' )
    {
        IJOBVR  = 1;
        ILVSR   = false;
    }
    else if ( JOBVSR[0] == 'V' || JOBVSR[0] == 'v' )
    {
        IJOBVR  = 2;
        ILVSR   = true;
    }
    else
    {
        IJOBVR  = -1;
        ILVSR   = false;
    };

    const char* ctrans_name = is_V_real? "T" : "C";
    const char* qtriu_name  = is_V_real? "H" : "U";

    // Test the input arguments
    INFO        = 0;
    bool LQUERY = ( LWORK == -1 || LIWORK == -1);

    if ( IJOBVL <= 0 )
        INFO    = -1;
    else if ( IJOBVR <= 0 )
        INFO    = -2;
    else if ( N < 0 )
        INFO    = -3;
    else if ( LDA < std::max( 1, N ) )
        INFO    = -4;
    else if ( LDB < std::max( 1, N ) )
        INFO    = -7;
    else if ( LDVSL < 1 || ( ILVSL && LDVSL < N ) )
        INFO    = -11;
    else if ( LDVSR < 1 || ( ILVSR && LDVSR < N ) )
        INFO    = -13;

    // Compute workspace
    // (Note: Comments in the code beginning "Workspace:" describe the
    //  minimal amount of workspace needed at that point in the code,
    //  as well as the preferred amount for good performance.
    //  NB refers to the optimal block size for the immediately
    //  following subroutine, as returned by ILAENV.)

    i_type lrwork;
    i_type lwork_hess;
    i_type liwork_hess;

    if ( INFO == 0 )
    {
        if ( N > 0 )
        {
            i_type NB_QR, NB_ORMQR, NB_ORGQR;
            
            if (std::is_same<V,d_type>::value)
            {
                NB_QR       = lapack::ilaenv( 1, "DGEQRF", " ", N, 1, N, 0 );
                NB_ORMQR    = lapack::ilaenv( 1, "DORMQR", " ", N, 1, N, -1 );
                NB_ORGQR    = lapack::ilaenv( 1, "DORGQR", " ", N, 1, N, -1 );
            }
            else if (std::is_same<V,s_type>::value)
            {
                NB_QR       = lapack::ilaenv( 1, "SGEQRF", " ", N, 1, N, 0 );
                NB_ORMQR    = lapack::ilaenv( 1, "SORMQR", " ", N, 1, N, -1 );
                NB_ORGQR    = lapack::ilaenv( 1, "SORGQR", " ", N, 1, N, -1 );
            }
            else if (std::is_same<V,c_type>::value)
            {
                NB_QR       = lapack::ilaenv( 1, "CGEQRF", " ", N, 1, N, 0 );
                NB_ORMQR    = lapack::ilaenv( 1, "CUNMQR", " ", N, 1, N, -1 );
                NB_ORGQR    = lapack::ilaenv( 1, "CUNGQR", " ", N, 1, N, -1 );
            }
            else
            {
                NB_QR       = lapack::ilaenv( 1, "ZGEQRF", " ", N, 1, N, 0 );
                NB_ORMQR    = lapack::ilaenv( 1, "ZUNMQR", " ", N, 1, N, -1 );
                NB_ORGQR    = lapack::ilaenv( 1, "ZUNGQR", " ", N, 1, N, -1 );
            }

            lrwork          = is_V_real? 8*N : 8*N / 2;

            i_type MIN_QR   = N;
            i_type MAX_QR   = N * NB_QR;
            i_type MIN_ORM  = N;
            i_type MAX_ORM  = N * NB_ORMQR;
            i_type MIN_ORG  = N;
            i_type MAX_ORG  = N * NB_ORGQR;
            i_type MIN_QZ   = N;
            i_type MAX_QZ   = N;

            i_type WORK_S   = N;

            MINWRK          = WORK_S;
            MAXWRK          = WORK_S;

            MINWRK          = std::max(MINWRK, WORK_S + MIN_QR );
            MINWRK          = std::max(MINWRK, WORK_S + MIN_ORM );                

            MAXWRK          = std::max( MAXWRK, WORK_S + MAX_QR );
            MAXWRK          = std::max( MAXWRK, WORK_S + MAX_ORM );                

            if ( ILVSL )
            {
                MINWRK      = std::max(MINWRK, MIN_ORG );
                MAXWRK      = std::max( MAXWRK, WORK_S + MAX_ORG );            
            };

            MINWRK          = std::max(MINWRK, MIN_QZ );
            MAXWRK          = std::max( MAXWRK, MAX_QZ );

            V work_query;
            i_type iwork_query;

            gghrd2(JOBVSL, "I", N, 1, N, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, 
                       &work_query, -1, &iwork_query, -1, IERR );

            lwork_hess      = (i_type)real(work_query);
            liwork_hess     = iwork_query;

            MINWRK          = std::max(MINWRK, lwork_hess + WORK_S);
            MAXWRK          = std::max(MAXWRK, lwork_hess + WORK_S);

            MINWRK          = MINWRK + lrwork;
            MAXWRK          = MAXWRK + lrwork;
        }
        else
        {
            MINWRK  = 1;
            MAXWRK  = 1;
        };

        WORK[0]     = VR(MAXWRK);
        IWORK[0]    = liwork_hess;

        if (LWORK < MINWRK && !LQUERY )
            INFO    = -15;
        if (LIWORK < liwork_hess && !LQUERY )
            INFO    = -17;
    };

    if ( INFO != 0 )
        return;
    else if (LQUERY )
        return;

    // Quick return if possible
    if ( N == 0 )
        return;

    // Get machine constants
    VR EPS      = lapack::lamch<V>( "P" );
    VR SAFMIN   = lapack::lamch<V>( "S" );
    VR SAFMAX   = VR(1.0) / SAFMIN;

    lapack::labad( SAFMIN, SAFMAX );
    VR SMLNUM   = std::sqrt( SAFMIN ) / EPS;
    VR BIGNUM   = VR(1.0) / SMLNUM;

    V* WORK_SAV = WORK;
    VR* RWORK   = reinterpret_cast<VR*>(WORK);
    WORK        = WORK + lrwork;

    // Scale A if max element outside range [SMLNUM,BIGNUM]
    VR ANRM     = lapack::lange( "M", N, N, A, LDA, RWORK );
    bool ILASCL = false;

    VR ANRMTO;

    if (ANRM > VR(0.0) &&  ANRM < SMLNUM )
    {
        ANRMTO  = SMLNUM;
        ILASCL  = true;
    }
    else if ( ANRM > BIGNUM )
    {
        ANRMTO  = BIGNUM;
        ILASCL  = true;
    };
    
    if ( ILASCL )
        lapack::lascl("G", 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR );

    // Scale B if max element outside range [SMLNUM,BIGNUM]
    VR BNRM     = lapack::lange("M", N, N, B, LDB, RWORK );
    VR BNRMTO;
    bool ILBSCL = false;

    if (BNRM > VR(0.0) && BNRM < SMLNUM )
    {
        BNRMTO  = SMLNUM;
        ILBSCL  = true;
    }
    else if ( BNRM > BIGNUM )
    {
        BNRMTO  = BIGNUM;
        ILBSCL  = true;
    };
    
    if ( ILBSCL )
         lapack::lascl("G", 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR );

    // Permute the matrix to make it more nearly triangular
    // (Real Workspace: need 6*N + 2*N space for storing balancing factors)
    i_type ILEFT    = 1;
    i_type IRIGHT   = N + 1;
    i_type IRWRK    = IRIGHT + N;
    i_type ILO, IHI;

    //lapack::ggba2("P", N, A, LDA, B, LDB, ILO, IHI, RWORK + ILEFT - 1, RWORK + IRIGHT-1, 
    //        RWORK + IRWRK - 1, IERR );

    const char* bal_char    = "P";
    lapack::ggbal2(bal_char, N, A, LDA, B, LDB, ILO, IHI, RWORK + ILEFT - 1, RWORK + IRIGHT-1, IERR );

    i_type DOHESS   = check_hess_type<V>(N, A, LDA, B, LDB, ILO, IHI);

    // Initialize VSL

    if ( ILVSL )
        lapack::laset( "Full", N, N, V(0.0), V(1.0), VSL, LDVSL );

    // Initialize VSR
    if (ILVSR && DOHESS == 0)
        lapack::laset("Full", N, N, V(0.0), V(1.0), VSR, LDVSR );

    i_type ITAU     = 1;
    i_type IWRK     = ITAU;

    if (DOHESS > 0)
    {
        if (DOHESS == 2)
        {
            // Reduce B to triangular form (QR decomposition of B)
            //   (Workspace: need N, prefer N*NB)

            i_type IROWS    = IHI + 1 - ILO;
            i_type ICOLS    = N + 1 - ILO;
            IWRK            = ITAU + IROWS;

            lapack::geqrf<V>(IROWS, ICOLS, B + (ILO-1) + (ILO-1)*LDB, LDB, WORK + ITAU-1, WORK + IWRK - 1, 
                    LWORK+1-IWRK, &IERR );

            // Apply the orthogonal transformation to matrix A
            //     (Workspace: need N, prefer N*NB)
            lapack::ormqr("L", ctrans_name, IROWS, ICOLS, IROWS, B+ (ILO-1) + (ILO-1)*LDB, LDB, WORK + ITAU - 1, 
                    A + (ILO-1) + (ILO-1)*LDA, LDA, WORK + IWRK - 1, LWORK+1-IWRK, IERR);

            // Initialize VSL
            //   (Workspace: need N, prefer N*NB)

            if ( ILVSL )
            {
                if (IROWS > 1 )
                {
                    lapack::lacpy("L", IROWS-1, IROWS-1, B + (ILO+1-1) + (ILO-1)*LDB, LDB,
                            VSL + (ILO+1-1) + (ILO-1)*LDVSL, LDVSL );
                };

                lapack::orgqr<V>(IROWS, IROWS, IROWS, VSL + (ILO-1) + (ILO-1)*LDVSL, LDVSL, WORK + ITAU-1, 
                        WORK + IWRK - 1, LWORK+1-IWRK, &IERR );
            };
        };

        // Reduce to generalized Hessenberg form
        //   (Workspace: need lwork_hess)
        //   (integer Workspace: need liwork_hess)

        lapack::gghrd2(JOBVSL, "I", N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, 
                       WORK + IWRK - 1, lwork_hess, IWORK, liwork_hess, IERR );
    };

    // Perform QZ algorithm, computing Schur vectors if desired
    // (Workspace: need N)
    // (Real Workspace: need N if complex)

    IWRK        = ITAU;
    VR* ALPHAR  = reinterpret_cast<VR*>(ALPHA);
    VR* ALPHAI  = ALPHAR + N;

    (void)ALPHAI;

    gges2_alg<V>::hgeqz("S", JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VSL, 
                        LDVSL, VSR, LDVSR, WORK + IWRK-1, LWORK+1-IWRK, RWORK + IRWRK - 1, IERR );

    if (IERR != 0 )
    {
        if (IERR > 0 && IERR < N )
            INFO    = IERR;
        else if (IERR > N && IERR <= 2*N )
            INFO    = IERR - N;
        else
            INFO    = N + 1;
        
        goto lab_50;
    };

    // Apply back-permutation to VSL and VSR
    //   (Workspace: none needed)

    if (ILVSL)
        lapack::ggbak( bal_char, "L", N, ILO, IHI, RWORK + ILEFT - 1, RWORK + IRIGHT - 1, N, VSL, LDVSL, IERR );

    if (ILVSR)
        lapack::ggbak( bal_char, "R", N, ILO, IHI, RWORK + ILEFT - 1, RWORK + IRIGHT - 1, N, VSR, LDVSR, IERR );

    gges2_alg<V>::rescale_eigs(ILASCL, ILBSCL, N, ALPHA, BETA, SAFMAX, SAFMIN, ANRMTO, ANRM, BNRMTO, BNRM,
            A, LDA, B, LDB);

    // Undo scaling
    if ( ILASCL )
    {
        lapack::lascl( qtriu_name, 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR );
        gges2_alg<V>::undo_scale_alpha(ANRMTO, ANRM, N, ALPHA);
    };

    if (ILBSCL)
    {
        lapack::lascl( "U", 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR );
        lapack::lascl( "G", 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );
    };

    gges2_alg<V>::make_complex_alpha(ALPHA, N, WORK_SAV);

  lab_50:
    WORK_SAV[0] = VR(MAXWRK);
    IWORK[0]    = liwork_hess;
};

template void BLAS_EXT_EXPORT
gges2<d_type>(const char *JOBVSL, const char *JOBVSR, i_type N, d_type *a, i_type LDA, 
        d_type *b, i_type LDB, z_type *alpha, d_type * beta, d_type * vsl, i_type LDVSL, d_type *vsr,
        i_type LDVSR, d_type *work, i_type LWORK, i_type* IWORK, i_type LIWORK, i_type& INFO);

template void BLAS_EXT_EXPORT
gges2<s_type>(const char *JOBVSL, const char *JOBVSR, i_type N, s_type *a, i_type LDA, 
        s_type *b, i_type LDB, c_type *alpha, s_type * beta, s_type * vsl, i_type LDVSL, s_type *vsr,
        i_type LDVSR, s_type *work, i_type LWORK, i_type* IWORK, i_type LIWORK, i_type& INFO);

template void BLAS_EXT_EXPORT
gges2<c_type>(const char *JOBVSL, const char *JOBVSR, i_type N, c_type *a, i_type LDA, 
        c_type *b, i_type LDB, c_type *alpha, c_type * beta, c_type * vsl, i_type LDVSL, c_type *vsr,
        i_type LDVSR, c_type *work, i_type LWORK, i_type* IWORK, i_type LIWORK, i_type& INFO);

template void BLAS_EXT_EXPORT
gges2<z_type>(const char *JOBVSL, const char *JOBVSR, i_type N, z_type *a, i_type LDA, 
        z_type *b, i_type LDB, z_type *alpha, z_type * beta, z_type * vsl, i_type LDVSL, z_type *vsr,
        i_type LDVSR, z_type *work, i_type LWORK, i_type* IWORK, i_type LIWORK, i_type& INFO);

}};

#pragma warning(pop)