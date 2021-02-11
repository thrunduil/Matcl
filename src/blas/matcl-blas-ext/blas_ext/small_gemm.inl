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

#include "matcl-blas-lapack/blas/details/blas_simd_complex.h"
#include "blas/matcl-blas-ext/blas_ext/blas_simd.h"
#include "matcl-simd/simd.h"
#include <algorithm>

namespace matcl { namespace lapack
{

// evaluate A = A * S, S is small matrix
template<class V, bool Use_simd = true>
struct gemm_M_S4
{
    static void eval(i_type M, V* A, i_type LDA, const V* B, i_type LDB)
    {
        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        if (M == 0)
            return;

        if (M < VS)
            return gemm_M_S4<V,false>::eval(M, A, LDA, B, LDB);

        V* A0               = A + 0 * LDA;
        V* A1               = A + 1 * LDA;
        V* A2               = A + 2 * LDA;
        V* A3               = A + 3 * LDA;

        simd_type b00       = simd_type::broadcast(B + 0);
        simd_type b10       = simd_type::broadcast(B + 1);
        simd_type b20       = simd_type::broadcast(B + 2);
        simd_type b30       = simd_type::broadcast(B + 3);

        simd_type b01       = simd_type::broadcast(B + 0 + 1*LDB);
        simd_type b11       = simd_type::broadcast(B + 1 + 1*LDB);
        simd_type b21       = simd_type::broadcast(B + 2 + 1*LDB);
        simd_type b31       = simd_type::broadcast(B + 3 + 1*LDB);

        simd_type b02       = simd_type::broadcast(B + 0 + 2*LDB);
        simd_type b12       = simd_type::broadcast(B + 1 + 2*LDB);
        simd_type b22       = simd_type::broadcast(B + 2 + 2*LDB);
        simd_type b32       = simd_type::broadcast(B + 3 + 2*LDB);

        simd_type b03       = simd_type::broadcast(B + 0 + 3*LDB);
        simd_type b13       = simd_type::broadcast(B + 1 + 3*LDB);
        simd_type b23       = simd_type::broadcast(B + 2 + 3*LDB);
        simd_type b33       = simd_type::broadcast(B + 3 + 3*LDB);

        i_type size_1       = M / VS;

        for (i_type i = 0; i < size_1; ++i)
        {
            simd_type a0    = simd_type::load(A0 + i * VS, std::false_type());
            simd_type a1    = simd_type::load(A1 + i * VS, std::false_type());
            simd_type a2    = simd_type::load(A2 + i * VS, std::false_type());
            simd_type a3    = simd_type::load(A3 + i * VS, std::false_type());

            simd_type c0    = (a0 * b00 + a1 * b10) + (a2*b20 + a3*b30);
            simd_type c1    = (a0 * b01 + a1 * b11) + (a2*b21 + a3*b31);
            simd_type c2    = (a0 * b02 + a1 * b12) + (a2*b22 + a3*b32);
            simd_type c3    = (a0 * b03 + a1 * b13) + (a2*b23 + a3*b33);

            c0.store(A0 + i * VS, std::false_type());
            c1.store(A1 + i * VS, std::false_type());
            c2.store(A2 + i * VS, std::false_type());
            c3.store(A3 + i * VS, std::false_type());
        };

        i_type pos          = size_1 * VS;

        if (pos != M)
            return gemm_M_S4<V,false>::eval(M - pos, A + pos, LDA, B, LDB);
    };
};

// evaluate A = A * S, S is small matrix
template<class V>
struct gemm_M_S4<V,false>
{
    static void eval(i_type M, V* A, i_type LDA, const V* B, i_type LDB)
    {
        //gemm<V>("No","No", M, N, K, V(1.0), A, LDA, B, LDB, V(0.0), C, LDC);

        V* A0               = A + 0 * LDA;
        V* A1               = A + 1 * LDA;
        V* A2               = A + 2 * LDA;
        V* A3               = A + 3 * LDA;

        V b00               = B[0 + 0*LDB];
        V b10               = B[1 + 0*LDB];
        V b20               = B[2 + 0*LDB];
        V b30               = B[3 + 0*LDB];

        V b01               = B[0 + 1*LDB];
        V b11               = B[1 + 1*LDB];
        V b21               = B[2 + 1*LDB];
        V b31               = B[3 + 1*LDB];

        V b02               = B[0 + 2*LDB];
        V b12               = B[1 + 2*LDB];
        V b22               = B[2 + 2*LDB];
        V b32               = B[3 + 2*LDB];

        V b03               = B[0 + 3*LDB];
        V b13               = B[1 + 3*LDB];
        V b23               = B[2 + 3*LDB];
        V b33               = B[3 + 3*LDB];

        for (i_type i = 0; i < M; ++i)
        {
            V a0            = A0[i];
            V a1            = A1[i];
            V a2            = A2[i];
            V a3            = A3[i];

            V c0            = a0 * b00 + a1 * b10 + a2 * b20 + a3 * b30;
            V c1            = a0 * b01 + a1 * b11 + a2 * b21 + a3 * b31;
            V c2            = a0 * b02 + a1 * b12 + a2 * b22 + a3 * b32;
            V c3            = a0 * b03 + a1 * b13 + a2 * b23 + a3 * b33;

            A0[i]           = c0;
            A1[i]           = c1;
            A2[i]           = c2;
            A3[i]           = c3;
        };
    };
};

// evaluate A = A * S, S is small matrix
template<class V, bool Use_simd = true>
struct gemm_M_S3
{
    static void eval(i_type M, V* A, i_type LDA, const V* B, i_type LDB)
    {
        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        if (M == 0)
            return;

        if (M < VS)
            return gemm_M_S3<V,false>::eval(M, A, LDA, B, LDB);

        V* A0               = A + 0 * LDA;
        V* A1               = A + 1 * LDA;
        V* A2               = A + 2 * LDA;

        simd_type b00       = simd_type::broadcast(B + 0);
        simd_type b10       = simd_type::broadcast(B + 1);
        simd_type b20       = simd_type::broadcast(B + 2);

        simd_type b01       = simd_type::broadcast(B + 0 + 1*LDB);
        simd_type b11       = simd_type::broadcast(B + 1 + 1*LDB);
        simd_type b21       = simd_type::broadcast(B + 2 + 1*LDB);

        simd_type b02       = simd_type::broadcast(B + 0 + 2*LDB);
        simd_type b12       = simd_type::broadcast(B + 1 + 2*LDB);
        simd_type b22       = simd_type::broadcast(B + 2 + 2*LDB);

        i_type size_1       = M / VS;

        for (i_type i = 0; i < size_1; ++i)
        {
            simd_type a0    = simd_type::load(A0 + i * VS, std::false_type());
            simd_type a1    = simd_type::load(A1 + i * VS, std::false_type());
            simd_type a2    = simd_type::load(A2 + i * VS, std::false_type());

            simd_type c0    = (a0 * b00 + a1 * b10) + a2*b20;
            simd_type c1    = (a0 * b01 + a1 * b11) + a2*b21;
            simd_type c2    = (a0 * b02 + a1 * b12) + a2*b22;

            c0.store(A0 + i * VS, std::false_type());
            c1.store(A1 + i * VS, std::false_type());
            c2.store(A2 + i * VS, std::false_type());
        };

        i_type pos          = size_1 * VS;

        if (pos != M)
            return gemm_M_S3<V,false>::eval(M - pos, A + pos, LDA, B, LDB);
    };
};

// evaluate A = A * S, S is small matrix
template<class V>
struct gemm_M_S3<V,false>
{
    static void eval(i_type M, V* A, i_type LDA, const V* B, i_type LDB)
    {
        //gemm<V>("No","No", M, N, K, V(1.0), A, LDA, B, LDB, V(0.0), C, LDC);

        V* A0               = A + 0 * LDA;
        V* A1               = A + 1 * LDA;
        V* A2               = A + 2 * LDA;

        V b00               = B[0 + 0*LDB];
        V b10               = B[1 + 0*LDB];
        V b20               = B[2 + 0*LDB];

        V b01               = B[0 + 1*LDB];
        V b11               = B[1 + 1*LDB];
        V b21               = B[2 + 1*LDB];

        V b02               = B[0 + 2*LDB];
        V b12               = B[1 + 2*LDB];
        V b22               = B[2 + 2*LDB];

        for (i_type i = 0; i < M; ++i)
        {
            V a0            = A0[i];
            V a1            = A1[i];
            V a2            = A2[i];

            V c0            = a0 * b00 + a1 * b10 + a2 * b20;
            V c1            = a0 * b01 + a1 * b11 + a2 * b21;
            V c2            = a0 * b02 + a1 * b12 + a2 * b22;

            A0[i]           = c0;
            A1[i]           = c1;
            A2[i]           = c2;
        };
    };
};

// evaluate A = S * A, S is small matrix
template<class V, bool Use_simd = simd::has_simd_size<V,4>::value>
struct gemm_S4_M
{
    static void eval(i_type N, const V* S, i_type LDS, V* A, i_type LDA)
    {
        using simd_type = typename simd::default_simd_type_elems<V,4>::type;

        static const i_type VS  = simd_type::vector_size;

        simd_type S0        = simd_type::load(S + 0 * LDS, std::false_type());
        simd_type S1        = simd_type::load(S + 1 * LDS, std::false_type());
        simd_type S2        = simd_type::load(S + 2 * LDS, std::false_type());
        simd_type S3        = simd_type::load(S + 3 * LDS, std::false_type());

        for (i_type i = 0; i < N; ++i)
        {
            simd_type a0    = simd_type::broadcast(A+0);
            simd_type a1    = simd_type::broadcast(A+1);
            simd_type a2    = simd_type::broadcast(A+2);
            simd_type a3    = simd_type::broadcast(A+3);

            simd_type C     = (S0 * a0 + S1 * a1) + (S2 * a2 + S3 * a3);

            C.store(A,std::false_type());

            A               += LDA;
        }
    }
};

template<class V>
struct gemm_S4_M<V,false>
{
    static void eval(i_type N, const V* B, i_type LDB, V* A, i_type LDA)
    {
        V b00               = B[0 + 0*LDB];
        V b10               = B[1 + 0*LDB];
        V b20               = B[2 + 0*LDB];
        V b30               = B[3 + 0*LDB];

        V b01               = B[0 + 1*LDB];
        V b11               = B[1 + 1*LDB];
        V b21               = B[2 + 1*LDB];
        V b31               = B[3 + 1*LDB];

        V b02               = B[0 + 2*LDB];
        V b12               = B[1 + 2*LDB];
        V b22               = B[2 + 2*LDB];
        V b32               = B[3 + 2*LDB];

        V b03               = B[0 + 3*LDB];
        V b13               = B[1 + 3*LDB];
        V b23               = B[2 + 3*LDB];
        V b33               = B[3 + 3*LDB];

        for (i_type i = 0; i < N; ++i)
        {
            V a0            = A[0];
            V a1            = A[1];
            V a2            = A[2];
            V a3            = A[3];

            V c0            = b00 * a0 + b01 * a1 + b02 * a2 + b03 * a3;
            V c1            = b10 * a0 + b11 * a1 + b12 * a2 + b13 * a3;
            V c2            = b20 * a0 + b21 * a1 + b22 * a2 + b23 * a3;
            V c3            = b30 * a0 + b31 * a1 + b32 * a2 + b33 * a3;

            A[0]            = c0;
            A[1]            = c1;
            A[2]            = c2;
            A[3]            = c3;

            A               += LDA;
        }
    }
};

template<class V>
struct gemm_S3_M
{
    static void eval(i_type N, const V* B, i_type LDB, V* A, i_type LDA)
    {
        V b00               = B[0 + 0*LDB];
        V b10               = B[1 + 0*LDB];
        V b20               = B[2 + 0*LDB];

        V b01               = B[0 + 1*LDB];
        V b11               = B[1 + 1*LDB];
        V b21               = B[2 + 1*LDB];

        V b02               = B[0 + 2*LDB];
        V b12               = B[1 + 2*LDB];
        V b22               = B[2 + 2*LDB];

        for (i_type i = 0; i < N; ++i)
        {
            V a0            = A[0];
            V a1            = A[1];
            V a2            = A[2];

            V c0            = b00 * a0 + b01 * a1 + b02 * a2;
            V c1            = b10 * a0 + b11 * a1 + b12 * a2;
            V c2            = b20 * a0 + b21 * a1 + b22 * a2;

            A[0]            = c0;
            A[1]            = c1;
            A[2]            = c2;

            A               += LDA;
        }
    }
};

};};