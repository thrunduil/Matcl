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

#pragma once

#include "matcl-blas-lapack/blas/blas.h"
#include "blas/matcl-blas-ext/timer.h"
#include <iostream>

#include "matcl-blas-lapack/blas/details/blas_simd_complex.h"
#include "matcl-simd/simd.h"

namespace matcl { namespace lapack
{

//-----------------------------------------------------------------------
//          gemm no_trans x no_trans 2x2 kernel
//-----------------------------------------------------------------------
template<class V>
struct gemm2_NN_M2_22_B0
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B
        // A = Mx2, B = 2*2, C = Mx2

        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        i_type size_1       = M / VS;

        const V* A0         = A;
        const V* A1         = A + LDA;

        V* C0               = C;
        V* C1               = C + LDC;

        V b00_v             = alpha*B[0 + 0 * LDB];
        V b10_v             = alpha*B[1 + 0 * LDB];
        V b01_v             = alpha*B[0 + 1 * LDB];
        V b11_v             = alpha*B[1 + 1 * LDB];

        if (M < VS)
            goto lab_tail;

        simd_type b00       = simd_type::broadcast(b00_v);
        simd_type b10       = simd_type::broadcast(b10_v);
        simd_type b01       = simd_type::broadcast(b01_v);
        simd_type b11       = simd_type::broadcast(b11_v);

        for (i_type k = 0; k < size_1; ++k) 
        {            
            simd_type a0    = simd_type::load(A0, std::false_type());
            simd_type a1    = simd_type::load(A1, std::false_type());
            
            simd_type c     = fma_f(a0, b00, a1*b10);
            c.store(C0, std::false_type());

            c               = fma_f(a0, b01, a1*b11);
            c.store(C1, std::false_type());

            C0              = C0 + VS;
            C1              = C1 + VS;
            A0              = A0 + VS;
            A1              = A1 + VS;
        };

        i_type pos          = size_1 * VS;

        if (pos == M)
            return;

        M                   = M - pos;

    lab_tail:

        for (i_type k = 0; k < M; ++k) 
        {
            V a0        = A0[k];
            V a1        = A1[k];

            V c0        = a0 * b00_v + a1 * b10_v;
            V c1        = a0 * b01_v + a1 * b11_v;

            C0[k]       = c0;
            C1[k]       = c1;
        };
    };
};

template<class V>
struct gemm2_NN_M2_22_B1
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B + C
        // A = Mx2, B = 2*2, C = Mx2

        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        i_type size_1       = M / VS;

        const V* A0         = A;
        const V* A1         = A + LDA;

        V* C0               = C;
        V* C1               = C + LDC;

        V b00_v             = alpha*B[0 + 0 * LDB];
        V b10_v             = alpha*B[1 + 0 * LDB];
        V b01_v             = alpha*B[0 + 1 * LDB];
        V b11_v             = alpha*B[1 + 1 * LDB];

        if (M < VS)
            goto lab_tail;

        simd_type b00       = simd_type::broadcast(b00_v);
        simd_type b10       = simd_type::broadcast(b10_v);
        simd_type b01       = simd_type::broadcast(b01_v);
        simd_type b11       = simd_type::broadcast(b11_v);

        for (i_type k = 0; k < size_1; ++k) 
        {            
            simd_type a0    = simd_type::load(A0, std::false_type());
            simd_type a1    = simd_type::load(A1, std::false_type());
            
            simd_type c     = simd_type::load(C0, std::false_type());
            c               = fma_f(a0, b00, fma_f(a1, b10, c));
            c.store(C0, std::false_type());

            c               = simd_type::load(C1, std::false_type());
            c               = fma_f(a0, b01, fma_f(a1, b11, c));
            c.store(C1, std::false_type());

            C0              = C0 + VS;
            C1              = C1 + VS;
            A0              = A0 + VS;
            A1              = A1 + VS;
        };

        i_type pos          = size_1 * VS;

        if (pos == M)
            return;

        M                   = M - pos;

    lab_tail:

        for (i_type k = 0; k < M; ++k) 
        {
            V a0        = A0[k];
            V a1        = A1[k];

            V c0        = a0 * b00_v + a1 * b10_v;
            V c1        = a0 * b01_v + a1 * b11_v;

            C0[k]       = C0[k] + c0;
            C1[k]       = C1[k] + c1;
        };
    };
};

template<class V>
struct gemm2_NN_M2_22_B
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B + beta * C
        // A = Mx2, B = 2*2, C = Mx2

        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        i_type size_1       = M / VS;

        V* C0               = C;
        V* C1               = C + LDC;

        if (M < VS)
            goto lab_tail;

        simd_type b         = simd_type::broadcast(beta);

        for (i_type k = 0; k < size_1; ++k) 
        {            
            simd_type c     = simd_type::load(C0, std::false_type());
            c               = b * c;
            c.store(C0, std::false_type());

            c               = simd_type::load(C1, std::false_type());
            c               = b * c;
            c.store(C1, std::false_type());

            C0              = C0 + VS;
            C1              = C1 + VS;
        };

        i_type pos          = size_1 * VS;

        goto lab_mult;

        M                   = M - pos;

    lab_tail:

        for (i_type k = 0; k < M; ++k) 
        {
            C0[k]       = beta * C0[k];
            C1[k]       = beta * C1[k];
        };

    lab_mult:
        gemm2_NN_M2_22_B1<V>::eval(M, alpha, A, LDA, B, LDB, C, LDC);

    };
};

template<class V>
struct gemm2_NN_M2_22
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        // C = alpha*A*B + beta*C
        // A = Mx2, B = 2*2, C = Mx2

        if (beta == V(0.0))
            return gemm2_NN_M2_22_B0<V>::eval(M, alpha, A, LDA, B, LDB, C, LDC);
        else if (beta == V(1.0))
            return gemm2_NN_M2_22_B1<V>::eval(M, alpha, A, LDA, B, LDB, C, LDC);
        else
            return gemm2_NN_M2_22_B<V>::eval(M, alpha, A, LDA, B, LDB, beta, C, LDC);
    };
};

//-----------------------------------------------------------------------
//          gemm no_trans x no_trans 2x1 kernel
//-----------------------------------------------------------------------
template<class V>
struct gemm2_NN_M2_21_B1
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B + C
        // A = Mx2, B = 2*1, C = Mx1

        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        i_type size_1       = M / VS;

        const V* A0         = A;
        const V* A1         = A + LDA;

        V* C0               = C;
        V* C1               = C + LDC;

        V b00_v             = alpha*B[0 + 0 * LDB];
        V b10_v             = alpha*B[1 + 0 * LDB];

        if (M < VS)
            goto lab_tail;

        simd_type b00       = simd_type::broadcast(b00_v);
        simd_type b10       = simd_type::broadcast(b10_v);

        for (i_type k = 0; k < size_1; ++k) 
        {            
            simd_type a0    = simd_type::load(A0, std::false_type());
            simd_type a1    = simd_type::load(A1, std::false_type());
            
            simd_type c     = simd_type::load(C0, std::false_type());
            c               = fma_f(a0, b00, fma_f(a1, b10, c));
            c.store(C0, std::false_type());

            C0              = C0 + VS;
            A0              = A0 + VS;
            A1              = A1 + VS;
        };

        i_type pos          = size_1 * VS;

        if (pos == M)
            return;

        M                   = M - pos;

    lab_tail:

        for (i_type k = 0; k < M; ++k) 
        {
            V a0        = A0[k];
            V a1        = A1[k];

            V c0        = a0 * b00_v + a1 * b10_v;

            C0[k]       = C0[k] + c0;
        };
    };
};
template<class V>
struct gemm2_NN_M2_21_B0
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B
        // A = Mx2, B = 2*1, C = Mx1

        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        i_type size_1       = M / VS;

        const V* A0         = A;
        const V* A1         = A + LDA;

        V* C0               = C;
        V* C1               = C + LDC;

        V b00_v             = alpha*B[0 + 0 * LDB];
        V b10_v             = alpha*B[1 + 0 * LDB];

        if (M < VS)
            goto lab_tail;

        simd_type b00       = simd_type::broadcast(b00_v);
        simd_type b10       = simd_type::broadcast(b10_v);

        for (i_type k = 0; k < size_1; ++k) 
        {            
            simd_type a0    = simd_type::load(A0, std::false_type());
            simd_type a1    = simd_type::load(A1, std::false_type());
            
            simd_type c     = fma_f(a0, b00, a1 *b10);
            c.store(C0, std::false_type());

            C0              = C0 + VS;
            A0              = A0 + VS;
            A1              = A1 + VS;
        };

        i_type pos          = size_1 * VS;

        if (pos == M)
            return;

        M                   = M - pos;

    lab_tail:

        for (i_type k = 0; k < M; ++k) 
        {
            V a0        = A0[k];
            V a1        = A1[k];

            V c0        = a0 * b00_v + a1 * b10_v;
            C0[k]       = c0;
        };
    };
};
template<class V>
struct gemm2_NN_M2_21_B
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B + C
        // A = Mx2, B = 2*1, C = Mx1

        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        i_type size_1       = M / VS;

        const V* A0         = A;
        const V* A1         = A + LDA;

        V* C0               = C;
        V* C1               = C + LDC;

        V b00_v             = alpha*B[0 + 0 * LDB];
        V b10_v             = alpha*B[1 + 0 * LDB];

        if (M < VS)
            goto lab_tail;

        simd_type b00       = simd_type::broadcast(b00_v);
        simd_type b10       = simd_type::broadcast(b10_v);
        simd_type b         = simd_type::broadcast(beta);

        for (i_type k = 0; k < size_1; ++k) 
        {            
            simd_type a0    = simd_type::load(A0, std::false_type());
            simd_type a1    = simd_type::load(A1, std::false_type());
            
            simd_type c     = simd_type::load(C0, std::false_type());
            c               = fma_f(a0, b00, fma_f(a1, b10, b*c));
            c.store(C0, std::false_type());

            C0              = C0 + VS;
            A0              = A0 + VS;
            A1              = A1 + VS;
        };

        i_type pos          = size_1 * VS;

        if (pos == M)
            return;

        M                   = M - pos;

    lab_tail:

        for (i_type k = 0; k < M; ++k) 
        {
            V a0        = A0[k];
            V a1        = A1[k];

            V c0        = a0 * b00_v + a1 * b10_v;

            C0[k]       = C0[k] + beta*c0;
        };
    };
};
template<class V>
struct gemm2_NN_M2_21
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B + beta*C
        // A = Mx2, B = 2*2, C = Mx2

        if (beta == V(0.0))
            return gemm2_NN_M2_21_B0<V>::eval(M, alpha, A, LDA, B, LDB, C, LDC);
        else if (beta == V(1.0))
            return gemm2_NN_M2_21_B1<V>::eval(M, alpha, A, LDA, B, LDB, C, LDC);
        else
            return gemm2_NN_M2_21_B<V>::eval(M, alpha, A, LDA, B, LDB, beta, C, LDC);
    };
};

//-----------------------------------------------------------------------
//          gemm no_trans x no_trans 1x2 kernel
//-----------------------------------------------------------------------
template<class V>
struct gemm2_NN_M1_12_B1
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B + C
        // A = Mx1, B = 1*2, C = Mx2

        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        i_type size_1       = M / VS;

        const V* A0         = A;

        V* C0               = C;
        V* C1               = C + LDC;

        V b00_v             = alpha*B[0 + 0 * LDB];
        V b01_v             = alpha*B[0 + 1 * LDB];

        if (M < VS)
            goto lab_tail;

        simd_type b00       = simd_type::broadcast(b00_v);
        simd_type b01       = simd_type::broadcast(b01_v);

        for (i_type k = 0; k < size_1; ++k) 
        {            
            simd_type a0    = simd_type::load(A0, std::false_type());
            
            simd_type c     = simd_type::load(C0, std::false_type());
            c               = fma_f(a0, b00, c);
            c.store(C0, std::false_type());

            c               = simd_type::load(C1, std::false_type());
            c               = fma_f(a0, b01, c);
            c.store(C1, std::false_type());

            C0              = C0 + VS;
            C1              = C1 + VS;
            A0              = A0 + VS;
        };

        i_type pos          = size_1 * VS;

        if (pos == M)
            return;

        M                   = M - pos;

    lab_tail:

        for (i_type k = 0; k < M; ++k) 
        {
            V a0        = A0[k];

            V c0        = a0 * b00_v;
            V c1        = a0 * b01_v;

            C0[k]       = C0[k] + c0;
            C1[k]       = C1[k] + c1;
        };
    };
};

template<class V>
struct gemm2_NN_M1_12_B0
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B
        // A = Mx1, B = 1*2, C = Mx2

        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        i_type size_1       = M / VS;

        const V* A0         = A;

        V* C0               = C;
        V* C1               = C + LDC;

        V b00_v             = alpha*B[0 + 0 * LDB];
        V b01_v             = alpha*B[0 + 1 * LDB];

        if (M < VS)
            goto lab_tail;

        simd_type b00       = simd_type::broadcast(b00_v);
        simd_type b01       = simd_type::broadcast(b01_v);

        for (i_type k = 0; k < size_1; ++k) 
        {            
            simd_type a0    = simd_type::load(A0, std::false_type());
            
            simd_type c     = a0 * b00;
            c.store(C0, std::false_type());

            c               = a0 * b01;
            c.store(C1, std::false_type());

            C0              = C0 + VS;
            C1              = C1 + VS;
            A0              = A0 + VS;
        };

        i_type pos          = size_1 * VS;

        if (pos == M)
            return;

        M                   = M - pos;

    lab_tail:

        for (i_type k = 0; k < M; ++k) 
        {
            V a0        = A0[k];

            V c0        = a0 * b00_v;
            V c1        = a0 * b01_v;

            C0[k]       = c0;
            C1[k]       = c1;
        };
    };
};

template<class V>
struct gemm2_NN_M1_12_B
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B + beta * C
        // A = Mx1, B = 1*2, C = Mx2

        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        i_type size_1       = M / VS;

        const V* A0         = A;

        V* C0               = C;
        V* C1               = C + LDC;

        V b00_v             = alpha*B[0 + 0 * LDB];
        V b01_v             = alpha*B[0 + 1 * LDB];

        if (M < VS)
            goto lab_tail;

        simd_type b00       = simd_type::broadcast(b00_v);
        simd_type b01       = simd_type::broadcast(b01_v);
        simd_type b         = simd_type::broadcast(beta);

        for (i_type k = 0; k < size_1; ++k) 
        {            
            simd_type a0    = simd_type::load(A0, std::false_type());
            
            simd_type c     = simd_type::load(C0, std::false_type());
            c               = fma_f(a0, b00, c*b);
            c.store(C0, std::false_type());

            c               = simd_type::load(C1, std::false_type());
            c               = fma_f(a0, b01, c*b);
            c.store(C1, std::false_type());

            C0              = C0 + VS;
            C1              = C1 + VS;
            A0              = A0 + VS;
        };

        i_type pos          = size_1 * VS;

        if (pos == M)
            return;

        M                   = M - pos;

    lab_tail:

        for (i_type k = 0; k < M; ++k) 
        {
            V a0        = A0[k];

            V c0        = a0 * b00_v;
            V c1        = a0 * b01_v;

            C0[k]       = beta*C0[k] + c0;
            C1[k]       = beta*C1[k] + c1;
        };
    };
};

template<class V>
struct gemm2_NN_M1_12
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B + beta*C
        // A = Mx2, B = 2*2, C = Mx2

        if (beta == V(0.0))
            return gemm2_NN_M1_12_B0<V>::eval(M, alpha, A, LDA, B, LDB, C, LDC);
        else if (beta == V(1.0))
            return gemm2_NN_M1_12_B1<V>::eval(M, alpha, A, LDA, B, LDB, C, LDC);
        else
            return gemm2_NN_M1_12_B<V>::eval(M, alpha, A, LDA, B, LDB, beta, C, LDC);
    };
};

//-----------------------------------------------------------------------
//          gemm no_trans x no_trans 1x1 kernel
//-----------------------------------------------------------------------
template<class V>
struct gemm2_NN_M1_11_B1
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B + C
        // A = Mx1, B = 1*1, C = Mx1

        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        i_type size_1       = M / VS;

        const V* A0         = A;
        V* C0               = C;

        V b00_v             = alpha*B[0 + 0 * LDB];

        if (M < VS)
            goto lab_tail;

        simd_type b00       = simd_type::broadcast(b00_v);

        for (i_type k = 0; k < size_1; ++k) 
        {            
            simd_type a0    = simd_type::load(A0, std::false_type());
            
            simd_type c     = simd_type::load(C0, std::false_type());
            c               = fma_f(a0, b00, c);
            c.store(C0, std::false_type());

            C0              = C0 + VS;
            A0              = A0 + VS;
        };

        i_type pos          = size_1 * VS;

        if (pos == M)
            return;

        M                   = M - pos;

    lab_tail:

        for (i_type k = 0; k < M; ++k) 
        {
            V a0        = A0[k];
            V c0        = a0 * b00_v;
            C0[k]       = C0[k] + c0;
        };
    };
};
template<class V>
struct gemm2_NN_M1_11_B0
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B
        // A = Mx1, B = 1*1, C = Mx1

        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        i_type size_1       = M / VS;

        const V* A0         = A;
        V* C0               = C;

        V b00_v             = alpha*B[0 + 0 * LDB];

        if (M < VS)
            goto lab_tail;

        simd_type b00       = simd_type::broadcast(b00_v);

        for (i_type k = 0; k < size_1; ++k) 
        {            
            simd_type a0    = simd_type::load(A0, std::false_type());
            
            simd_type c     = fma_f(a0, b00, c);
            c.store(C0, std::false_type());

            C0              = C0 + VS;
            A0              = A0 + VS;
        };

        i_type pos          = size_1 * VS;

        if (pos == M)
            return;

        M                   = M - pos;

    lab_tail:

        for (i_type k = 0; k < M; ++k) 
        {
            V a0        = A0[k];
            V c0        = a0 * b00_v;
            C0[k]       = c0;
        };
    };
};
template<class V>
struct gemm2_NN_M1_11_B
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B + C
        // A = Mx1, B = 1*1, C = Mx1

        using simd_type = typename simd::default_simd_type<V>::type;

        static const i_type VS  = simd_type::vector_size;

        i_type size_1       = M / VS;

        const V* A0         = A;
        V* C0               = C;

        V b00_v             = alpha*B[0 + 0 * LDB];

        if (M < VS)
            goto lab_tail;

        simd_type b00       = simd_type::broadcast(b00_v);
        simd_type b         = simd_type::broadcast(beta);

        for (i_type k = 0; k < size_1; ++k) 
        {            
            simd_type a0    = simd_type::load(A0, std::false_type());
            
            simd_type c     = simd_type::load(C0, std::false_type());
            c               = fma_f(a0, b00, c*b);
            c.store(C0, std::false_type());

            C0              = C0 + VS;
            A0              = A0 + VS;
        };

        i_type pos          = size_1 * VS;

        if (pos == M)
            return;

        M                   = M - pos;

    lab_tail:

        for (i_type k = 0; k < M; ++k) 
        {
            V a0        = A0[k];
            V c0        = a0 * b00_v;
            C0[k]       = beta*C0[k] + c0;
        };
    };
};
template<class V>
struct gemm2_NN_M1_11
{
    static void eval(i_type M, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        // C = alpha*A*B + beta*C
        // A = Mx1, B = 1*1, C = Mx1

        if (beta == V(0.0))
            return gemm2_NN_M1_11_B0<V>::eval(M, alpha, A, LDA, B, LDB, C, LDC);
        else if (beta == V(1.0))
            return gemm2_NN_M1_11_B1<V>::eval(M, alpha, A, LDA, B, LDB, C, LDC);
        else
            return gemm2_NN_M1_11_B<V>::eval(M, alpha, A, LDA, B, LDB, beta, C, LDC);
    };
};

//-----------------------------------------------------------------------
//          gemm no_trans x no_trans
//-----------------------------------------------------------------------
template<class V>
struct gemm2_NN
{
    static void eval(i_type M, i_type N, i_type K, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, 
                     V beta, V* C, i_type LDC)
    {
        // C = alpha*A*B + beta*C

        i_type K2   = K / 2;
        i_type N2   = N / 2;

        if (K2 == 0)
        {
            for (i_type j = 0; j < N2; ++j)
            {
                gemm2_NN_M1_12<V>::eval(M, alpha, A, LDA, B, LDB, beta, C, LDC);

                B       += LDB * 2;
                C       += LDC * 2;
            };

            if (N2*2 != N)
                gemm2_NN_M1_11<V>::eval(M, alpha, A, LDA, B, LDB, beta, C, LDC);

            return;
        };

        for (i_type j = 0; j < N2; ++j)
        {
            gemm2_NN_M2_22<V>::eval(M, alpha, A, LDA, B, LDB, beta, C, LDC);
        
            for (i_type k = 1; k < K2; ++k)
                gemm2_NN_M2_22_B1<V>::eval(M, alpha, A + 2*k*LDA, LDA, B + 2*k, LDB, C, LDC);

            if (K2*2 != K)
                gemm2_NN_M1_12_B1<V>::eval(M, alpha, A + 2*K2*LDA, LDA, B + 2*K2, LDB, C, LDC);

            B       += LDB * 2;
            C       += LDC * 2;
        };

        if (N2*2 != N)
        {
            gemm2_NN_M2_21<V>::eval(M, alpha, A, LDA, B, LDB, beta, C, LDC);
        
            for (i_type k = 1; k < K2; ++k)
                gemm2_NN_M2_21_B1<V>::eval(M, alpha, A + 2*k*LDA, LDA, B + 2*k, LDB, C, LDC);

            if (K2*2 != K)
                gemm2_NN_M1_11_B1<V>::eval(M, alpha, A + 2*K2*LDA, LDA, B + 2*K2, LDB, C, LDC);
        };
    }
};

//-----------------------------------------------------------------------
//          make_trans
//-----------------------------------------------------------------------
template<class V>
void make_trans_22(V* BT, const V* B, i_type LDB)
{
    BT[0 + 0*2] = B[0 + 0*LDB];
    BT[1 + 0*2] = B[0 + 1*LDB];
    BT[0 + 1*2] = B[1 + 0*LDB];
    BT[1 + 1*2] = B[1 + 1*LDB];
};
template<class V>
void make_trans_12(V* BT, const V* B, i_type LDB)
{
    BT[0 + 0*2] = B[0 + 0*LDB];
    BT[0 + 1*2] = B[1 + 0*LDB];
};
template<class V>
void make_trans_21(V* BT, const V* B, i_type LDB)
{
    BT[0 + 0*2] = B[0 + 0*LDB];
    BT[1 + 0*2] = B[0 + 1*LDB];
};
template<class V>
void make_trans_11(V* BT, const V* B, i_type LDB)
{
    BT[0 + 0*2] = B[0 + 0*LDB];
};

//-----------------------------------------------------------------------
//          gemm no_trans x trans
//-----------------------------------------------------------------------
template<class V>
struct gemm2_NT
{
    static void eval(i_type M, i_type N, i_type K, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, 
                     V beta, V* C, i_type LDC)
    {
        // C = alpha*A*B' + beta*C

        i_type K2   = K / 2;
        i_type N2   = N / 2;

        V BT[2*2];

        if (K2 == 0)
        {
            for (i_type j = 0; j < N2; ++j)
            {
                make_trans_12(BT, B, LDB);
                gemm2_NN_M1_12<V>::eval(M, alpha, A, LDA, BT, 2, beta, C, LDC);

                B       += 2;
                C       += LDC * 2;
            };

            if (N2*2 != N)
            {
                make_trans_11(BT, B, LDB);
                gemm2_NN_M1_11<V>::eval(M, alpha, A, LDA, BT, 2, beta, C, LDC);
            }

            return;
        };

        for (i_type j = 0; j < N2; ++j)
        {
            make_trans_22(BT, B, LDB);
            gemm2_NN_M2_22<V>::eval(M, alpha, A, LDA, BT, 2, beta, C, LDC);
        
            for (i_type k = 1; k < K2; ++k)
            {
                make_trans_22(BT, B + 2*k*LDB, LDB);
                gemm2_NN_M2_22_B1<V>::eval(M, alpha, A + 2*k*LDA, LDA, BT, 2, C, LDC);
            }

            if (K2*2 != K)
            {
                make_trans_12(BT, B + 2*K2*LDB, LDB);
                gemm2_NN_M1_12_B1<V>::eval(M, alpha, A + 2*K2*LDA, LDA, BT, 2, C, LDC);
            }

            B       += 2;
            C       += LDC * 2;
        };

        if (N2*2 != N)
        {
            make_trans_21(BT, B, LDB);
            gemm2_NN_M2_21<V>::eval(M, alpha, A, LDA, BT, 2, beta, C, LDC);
        
            for (i_type k = 1; k < K2; ++k)
            {
                make_trans_21(BT, B + 2*k*LDB, LDB);
                gemm2_NN_M2_21_B1<V>::eval(M, alpha, A + 2*k*LDA, LDA, BT, 2, C, LDC);
            }

            if (K2*2 != K)
            {
                make_trans_11(BT, B + 2*K2*LDB, LDB);
                gemm2_NN_M1_11_B1<V>::eval(M, alpha, A + 2*K2*LDA, LDA, BT, 2, C, LDC);
            }
        };
    }
};

//-----------------------------------------------------------------------
//          gemm no_trans x no_trans band
//-----------------------------------------------------------------------
template<class V>
struct gemm_NN_band
{
    static void eval(i_type LD_A, i_type UD_A, i_type LD_B, i_type UD_B, 
                    i_type M, i_type N, i_type K, V alpha, const V* A, 
                    i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        // C = alpha*A*B + beta*C, where A has LD_A subdiagonals and UD_A superdiagonals
        // and B has LD_B subdiagonals and UD_B superdiagonals

        i_type K2   = K / 2;
        i_type N2   = N / 2;

        if (K2 == 0)
        {
            for (i_type j = 0; j < N2; ++j)
            {
                gemm2_NN_M1_12<V>::eval(M, alpha, A, LDA, B, LDB, beta, C, LDC);

                B       += LDB * 2;
                C       += LDC * 2;
            };

            if (N2*2 != N)
                gemm2_NN_M1_11<V>::eval(M, alpha, A, LDA, B, LDB, beta, C, LDC);

            return;
        };

        for (i_type j = 0; j < N2; ++j)
        {
            gemm2_NN_M2_22<V>::eval(M, alpha, A, LDA, B, LDB, beta, C, LDC);

            i_type k1   = std::max(j - (UD_B + 1)/2, 1);
            i_type k2   = std::min(j + (LD_B + 1)/2 + 1, K2);

            for (i_type k = k1; k < k2; ++k)
            {
                i_type F    = std::max(2*k + 1 - UD_A, 1);
                i_type L    = std::min(2*k + 1 + LD_A + 1, M);
                i_type M2   = L - F + 1;
                gemm2_NN_M2_22_B1<V>::eval(M2, alpha, A + F - 1 + 2*k*LDA, LDA, B + 2*k, LDB, C + F - 1, LDC);
            }

            if (K2*2 != K && k2 == K2)
            {
                i_type F    = std::max(2*K2 + 1 - UD_A, 1);
                i_type L    = std::min(2*K2 + 1 + LD_A + 1, M);
                i_type M2   = L - F + 1;
                gemm2_NN_M1_12_B1<V>::eval(M2, alpha, A + F - 1 + 2*K2*LDA, LDA, B + 2*K2, LDB, C + F - 1, LDC);
            };

            B       += LDB * 2;
            C       += LDC * 2;
        };

        if (N2*2 != N)
        {
            gemm2_NN_M2_21<V>::eval(M, alpha, A, LDA, B, LDB, beta, C, LDC);

            i_type k1   = std::max(N2 - (UD_B + 1)/2, 1);
            i_type k2   = std::min(N2 + (LD_B + 1)/2 + 1, K2);

            for (i_type k = k1; k < k2; ++k)
            {
                i_type F    = std::max(2*k + 1 - UD_A, 1);
                i_type L    = std::min(2*k + 1 + LD_A + 1, M);
                i_type M2   = L - F + 1;
                gemm2_NN_M2_21_B1<V>::eval(M2, alpha, A + F - 1 + 2*k*LDA, LDA, B + 2*k, LDB, C + F - 1, LDC);
            };

            if (K2*2 != K && k2 == K2)
            {
                i_type F    = std::max(2*K2 + 1 - UD_A, 1);
                i_type L    = std::min(2*K2 + 1 + LD_A + 1, M);
                i_type M2   = L - F + 1;
                gemm2_NN_M1_11_B1<V>::eval(M2, alpha, A + F - 1 + 2*K2*LDA, LDA, B + 2*K2, LDB, C + F - 1, LDC);
            };
        };
    }
};

//-----------------------------------------------------------------------
//          gemm no_trans x trans band
//-----------------------------------------------------------------------
template<class V>
struct gemm_NT_band
{
    static void eval(i_type LD_B, i_type UD_B, i_type M, i_type N, i_type K, V alpha, const V* A, 
                     i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        // C = alpha*A*B' + beta*C

        i_type K2   = K / 2;
        i_type N2   = N / 2;

        V BT[2*2];

        if (K2 == 0)
        {
            for (i_type j = 0; j < N2; ++j)
            {
                make_trans_12(BT, B, LDB);
                gemm2_NN_M1_12<V>::eval(M, alpha, A, LDA, BT, 2, beta, C, LDC);

                B       += 2;
                C       += LDC * 2;
            };

            if (N2*2 != N)
            {
                make_trans_11(BT, B, LDB);
                gemm2_NN_M1_11<V>::eval(M, alpha, A, LDA, BT, 2, beta, C, LDC);
            }

            return;
        };

        for (i_type j = 0; j < N2; ++j)
        {
            make_trans_22(BT, B, LDB);
            gemm2_NN_M2_22<V>::eval(M, alpha, A, LDA, BT, 2, beta, C, LDC);

            i_type k1   = std::max(j - (UD_B + 1)/2, 1);
            i_type k2   = std::min(j + (LD_B + 1)/2 + 1, K2);

            for (i_type k = k1; k < k2; ++k)
            {
                make_trans_22(BT, B + 2*k*LDB, LDB);
                gemm2_NN_M2_22_B1<V>::eval(M, alpha, A + 2*k*LDA, LDA, BT, 2, C, LDC);
            }

            if (K2*2 != K && k2 == K2)
            {
                make_trans_12(BT, B + 2*K2*LDB, LDB);
                gemm2_NN_M1_12_B1<V>::eval(M, alpha, A + 2*K2*LDA, LDA, BT, 2, C, LDC);
            }

            B       += 2;
            C       += LDC * 2;
        };

        if (N2*2 != N)
        {
            make_trans_21(BT, B, LDB);
            gemm2_NN_M2_21<V>::eval(M, alpha, A, LDA, BT, 2, beta, C, LDC);
        
            i_type k1   = std::max(N2 - (UD_B + 1)/2, 1);
            i_type k2   = std::min(N2 + (LD_B + 1)/2 + 1, K2);

            for (i_type k = k1; k < k2; ++k)
            {
                make_trans_21(BT, B + 2*k*LDB, LDB);
                gemm2_NN_M2_21_B1<V>::eval(M, alpha, A + 2*k*LDA, LDA, BT, 2, C, LDC);
            }

            if (K2*2 != K && k2 == K2)
            {
                make_trans_11(BT, B + 2*K2*LDB, LDB);
                gemm2_NN_M1_11_B1<V>::eval(M, alpha, A + 2*K2*LDA, LDA, BT, 2, C, LDC);
            }
        };
    }
};

};}
