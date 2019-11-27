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

#include <algorithm>

namespace matcl { namespace lapack
{

template<class V>
typename details::real_type<V>::type a2sum(i_type N, const V* A, i_type inc)
{
    using VR = typename details::real_type<V>::type;

    VR ret  = VR(0.0);

    for (i_type i = 0; i < N; ++i)
    {
        ret += abs2(A[0]);
        A   += inc;
    };

    return ret;
};

template<class V, bool Is_compl = details::is_complex<V>::value>
struct ldexp_helper
{
    static V eval(const V& val, i_type e)
    {
        return ldexp(val,e);
    };
};

template<class V>
struct ldexp_helper<V,true>
{
    static V eval(const V& val, i_type e)
    {
        using VR = typename details::real_type<V>::type;

        VR re   = ldexp(real(val),e);
        VR im   = ldexp(imag(val),e);
        return V(re,im);
    };
};

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
ggbal2(const char *JOB, i_type N, V* A, i_type LDA, V* B, i_type LDB, i_type& ILO, i_type& IHI, 
        typename details::real_type<V>::type* LSCALE, typename details::real_type<V>::type* RSCALE, i_type& INFO)
{
    using VR    = typename details::real_type<V>::type;

    //limit number of iterations in scaling
    static const i_type max_iter    = 5;

    // Test the input parameters
    INFO            = 0;

    if (JOB[0] != 'n' && JOB[0] != 'N' && JOB[0] != 'p' && JOB[0] != 'P' && JOB[0] != 's'
        && JOB[0] != 'S' && JOB[0] != 'b' && JOB[0] != 'B')
    {
        INFO        = -1;
    }
    else if (N < 0 )
        INFO        = -2;
    else if (LDA < std::max( 1, N ) )
        INFO        = -4;
    else if (LDB < std::max( 1, N ) )
         INFO       = -6;

    if ( INFO != 0 )
        return;

    // Quick return if possible
    if ( N == 0 )
    {
         ILO        = 1;
         IHI        = N;
         return;
    };

    if (N == 1 )
    {
        ILO         = 1;
        IHI         = N;
        LSCALE[0]   = VR(1.0);
        RSCALE[0]   = VR(1.0);
        return;
    };

    if (JOB[0] == 'n' || JOB[0] == 'N')
    {
        ILO         = 1;
        IHI         = N;

        for (i_type I = 0; I < N; ++I)
        {
            LSCALE[I]   = VR(1.0);
            RSCALE[I]   = VR(1.0);
        };

        return;
    };

    i_type K    = 1;
    i_type L    = N;

    i_type      LM1 = 0, I = 0, J = 0, JP1 = 0, M = 0, IFLOW = 0, IP1 = 0;

    V ZERO      = V(0.0);

    if ( JOB[0] == 's' || JOB[0] == 'S') 
        goto lab_190;
    else
        goto lab_30;

    // Permute the matrices A and B to isolate the eigenvalues.

    //     Find row with one nonzero in columns 1 through L
  lab_20:
    L           = LM1;

    if ( L != 1 )
        goto lab_30;

    RSCALE[0]   = VR(1.0);
    LSCALE[0]   = VR(1.0);
    goto lab_190;

  lab_30:
    LM1         = L - 1;
    for (I = L; I >= 1; --I)
    {
        for (J = 1; J <= LM1; ++J)
        {
            JP1     = J + 1;

            if ( A[I-1+(J-1)*LDA] != ZERO || B[I-1+(J-1)*LDB] != ZERO )
                goto lab_50;
        };
        
        J       = L;
        goto lab_70;

      lab_50:
        for (J = JP1; J <= L; ++J)
        {
            if ( A[I-1+(J-1)*LDA] != ZERO || B[I-1+(J-1)*LDB] != ZERO )
                goto lab_80;
        };

        J       = JP1 - 1;

      lab_70:
        M       = L;
        IFLOW   = 1;

        goto lab_160;

      lab_80:
        ;
    };

    goto lab_100;

    // Find column with one nonzero in rows K through N
  lab_90:
    K           = K + 1;

  lab_100:
    for( J = K; J <= L; ++J)
    {
        for (I = K; I <= LM1; ++I)
        {
            IP1     = I + 1;

            if( A[I-1+(J-1)*LDA] != ZERO || B[I-1+(J-1)*LDB] != ZERO )
                goto lab_120;
        };
         
        I           = L;
        goto lab_140;

      lab_120:
        for (I = IP1; I <= L; ++I)
        {
            if ( A[I-1+(J-1)*LDA] != ZERO || B[I-1+(J-1)*LDB] != ZERO )
                goto lab_150;
        };

        I           = IP1 - 1;

      lab_140:
         M          = K;
         IFLOW      = 2;
         goto lab_160;

       lab_150:
         ;
    };
    
    goto lab_190;

    // Permute rows M and I
  lab_160:

    LSCALE[M-1]     = VR(I);
    if ( I != M )
    {
        lapack::swap( N-K+1, A + (I-1) + (K-1)*LDA, LDA, A + M-1 +(K-1)*LDA, LDA );
        lapack::swap( N-K+1, B + (I-1) + (K-1)*LDB, LDB, B + M-1 + (K-1)* LDB, LDB );
    };

    // Permute columns M and J
    RSCALE[M-1]     = VR(J);
    if ( J != M )
    {
        lapack::swap( L, A+(J-1)*LDA, 1, A+(M-1)*LDA, 1 );
        lapack::swap( L, B+(J-1)*LDB, 1, B+(M-1)*LDB, 1 );
    };

    if (IFLOW == 1)
        goto lab_20;
    else 
        goto lab_90;

  lab_190:
    ;

    ILO     = K;
    IHI     = L;

    if (JOB[0] == 'P' || JOB[0] == 'p')
    {
        for (I = ILO; I <= IHI; ++I)
        {
            LSCALE[I-1] = VR(1.0);
            RSCALE[I-1] = VR(1.0);
        };
        return;
    };

    if (ILO == IHI )
        return;

    i_type KB       = IHI - ILO + 1;

    for (I = ILO; I <= IHI; ++I)
    {
        LSCALE[I-1] = VR(1.0);
        RSCALE[I-1] = VR(1.0);
    };

    for (i_type iter = 1; iter <= max_iter; ++iter)
    {
        i_type emax     = 0;
        i_type emin     = 0;

        for (i_type I = ILO - 1; I < IHI; ++I)
        {
            // scale the rows of abs2(A) + abs2(M) to have approximate row sum 1
            VR d        = a2sum<V>(KB, A + I + (ILO-1)*LDA, LDA)
                        + a2sum<V>(KB, B + I + (ILO-1)*LDB, LDB);

            if (d == VR(0.0))
                continue;

            i_type e    = -lround(log2(d)/2);
            emax        = std::max(emax,e);
            emin        = std::min(emin,e);

            LSCALE[I]   = ldexp_helper<VR>::eval(LSCALE[I], e);

            if (e != VR(0.0) )
            {
                //scale rows of A and B
                V* A_loc    = A + (ILO-1)*LDA;
                V* B_loc    = B + (ILO-1)*LDB;
                for (i_type k = 0; k < N - ILO + 1; ++k)
                {
                    A_loc[I]    = ldexp_helper<V>::eval(A_loc[I],e);
                    B_loc[I]    = ldexp_helper<V>::eval(B_loc[I],e);

                    A_loc   += LDA;
                    B_loc   += LDB;
                };
            };
        }        

        for (i_type I = ILO - 1; I < IHI; ++I)
        {
            // scale the columns of abs2(A) + abs2(M) to have approximate column sum 1

            VR d        = lapack::a2sum<V>(KB, A + (ILO-1) + I * LDA, 1)
                        + lapack::a2sum<V>(KB, B + (ILO-1) + I * LDB, 1);

            if (d == VR(0.0))
                continue;

            i_type e    = -lround(log2(d)/2);
            emax        = std::max(emax,e);
            emin        = std::min(emin,e);

            RSCALE[I]   = ldexp_helper<VR>::eval(RSCALE[I], e);

            if (e != VR(0.0) )
            {
                //scale columns of A and B
                V* A_loc    = A + I * LDA;
                V* B_loc    = B + I * LDB;
                
                for (i_type k = 0; k < IHI; ++k)
                {
                    A_loc[k]    = ldexp_helper<V>::eval(A_loc[k],e);
                    B_loc[k]    = ldexp_helper<V>::eval(B_loc[k],e);
                };
            };
        };

        // Stop if norms are all between 1/2 and 2
        if (emax <= emin+2)
            break;
    };
};

template BLAS_EXT_EXPORT void
ggbal2<d_type>(const char *job, i_type n, d_type* a, i_type lda, d_type* b, i_type ldb, i_type& ilo,
        i_type& ihi, d_type* lscale, d_type* rscale, i_type& info);
template BLAS_EXT_EXPORT void
ggbal2<s_type>(const char *job, i_type n, s_type* a, i_type lda, s_type* b, i_type ldb, i_type& ilo,
        i_type& ihi, s_type* lscale, s_type* rscale, i_type& info);
template BLAS_EXT_EXPORT void
ggbal2<c_type>(const char *job, i_type n, c_type* a, i_type lda, c_type* b, i_type ldb, i_type& ilo,
        i_type& ihi, s_type* lscale, s_type* rscale, i_type& info);
template BLAS_EXT_EXPORT void
ggbal2<z_type>(const char *job, i_type n, z_type* a, i_type lda, z_type* b, i_type ldb, i_type& ilo,
        i_type& ihi, d_type* lscale, d_type* rscale, i_type& info);

}};