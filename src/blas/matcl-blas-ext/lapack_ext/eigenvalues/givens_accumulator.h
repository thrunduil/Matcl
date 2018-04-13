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
#include "blas/matcl-blas-ext/blas_ext/gemm2_NN.inl"
#include "matcl-blas-ext/blas_concurrency.h"

namespace matcl { namespace lapack
{

template<class V, bool Is_compl = details::is_complex<V>::value>
struct gemm2
{
    static void eval(const char* tr_A, const char* tr_B, i_type M, i_type N, i_type K, 
                     V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        if (M == 0 || N == 0 || K == 0)
            return;

        //return lapack::gemm(tr_A, tr_B, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);

        if ((tr_A[0] == 'N' || tr_A[0] == 'n') && (tr_B[0] == 'N' || tr_B[0] == 'n'))
        {
            gemm2_NN<V>::eval(M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
        }
        else if ((tr_A[0] == 'N' || tr_A[0] == 'n') && !(tr_B[0] == 'N' || tr_B[0] == 'n'))
        {
            gemm2_NT<V>::eval(M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
        }
        else
        {
            return lapack::gemm(tr_A, tr_B, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
        };
    };
};

template<class V>
struct gemm2<V,true>
{
    static void eval(const char* tr_A, const char* tr_B, i_type M, i_type N, i_type K, 
                     V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        return lapack::gemm(tr_A, tr_B, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    };
};

template<class V, bool Is_compl = details::is_complex<V>::value>
struct gemm_band_l
{
    static void eval(i_type LD, i_type UD, i_type M, i_type N, i_type K, 
                     V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        //return lapack::gemm("No", "No", M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);

        if (M == 0 || N == 0 || K == 0)
            return;

        gemm_NN_band<V>::eval(LD, UD, K, N, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    };
};
template<class V>
struct gemm_band_l<V, true>
{
    static void eval(i_type LD, i_type UD, i_type M, i_type N, i_type K, 
                     V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V beta, V* C, i_type LDC)
    {
        return lapack::gemm("No", "No", M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    };
};

template<class V, bool Is_compl = details::is_complex<V>::value>
struct gemm_band_r
{
    static void eval(i_type LD_A, i_type UD_A, i_type LD_B, i_type UD_B, const char* tr_A, const char* tr_B, 
                     i_type M, i_type N, i_type K, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, V
                     beta, V* C, i_type LDC)
    {
        static const bool is_compl  = details::is_complex<V>::value;

        if (M == 0 || N == 0 || K == 0)
            return;

        //return lapack::gemm(tr_A, tr_B, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);

        if ((tr_A[0] == 'N' || tr_A[0] == 'n') && (tr_B[0] == 'N' || tr_B[0] == 'n'))
        {
            gemm_NN_band<V>::eval(LD_A, UD_A, LD_B, UD_B, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
        }
        else if ((tr_A[0] == 'N' || tr_A[0] == 'n') && !(tr_B[0] == 'N' || tr_B[0] == 'n'))
        {
            gemm_NT_band<V>::eval(LD_B, UD_B, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
        }
        else
        {
            return lapack::gemm(tr_A, tr_B, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
        };
    };
};

template<class V>
struct gemm_band_r<V,true>
{
    static void eval(i_type LD_A, i_type UD_A, i_type LD_B, i_type UD_B, const char* tr_A, const char* tr_B, 
                     i_type M, i_type N, i_type K, V alpha, const V* A, i_type LDA, const V* B, i_type LDB, 
                     V beta, V* C, i_type LDC)
    {
        return lapack::gemm(tr_A, tr_B, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    };
};

template<class V, bool Is_compl = details::is_complex<V>::value>
struct get_trans_char
{
    static const char* eval(const char* trans)
    {
        if (trans[0] == 'N' || trans[0] == 'n')
            return "no";
        else
            return "trans";
    };
};
template<class V>
struct get_trans_char<V, true>
{
    static const char* eval(const char* trans)
    {
        return trans;
    };
};

template<class V>
class givens_accumulator
{
    private:
        using VR    = typename details::real_type<V>::type;

    private:
        i_type      m_N_max;
        i_type      m_N_current;
        i_type      m_NB;
        i_type      m_NBM_left;
        i_type      n_workspaces_max;
        i_type      n_workspaces;
        i_type      m_NB_last;
        i_type      n_blocks;
        V*          m_WORK;
        V*          m_WORK_mult;
        bool        m_left;

    public:
        givens_accumulator(i_type N, i_type NB, i_type NBM, i_type n_threads);

        i_type      work_size() const;
        void        set_work(V* WORK);

        void        accumulate(bool left, i_type N, i_type NJ, const VR* SEQ_C, const V* SEQ_S, 
                               const i_type* SEQ_I1, const i_type* SEQ_I2, const i_type* SEQ_J);

        void        apply_left_left(i_type M, V* A, i_type LDA);
        void        apply_left_right(i_type M, i_type& UD, V* A, i_type LDA);
        void        apply_right_right(i_type M, i_type& UD, V* A, i_type LDA);

    private:
        V*          get_U(i_type num) const;
        V*          get_U_last() const;
        i_type      get_offset(i_type bl) const;
        i_type      get_size(i_type bl) const;
        void        initialize(V* U, i_type K);
        void        get_block_position(i_type i1, i_type i2, i_type j, i_type& block_num, i_type& offset) const;
        void        add_rotation(bool left, V* U, i_type K, i_type i1, i_type i2, i_type j, VR c, V s);
        void        set_block_params(i_type N);
        void        apply_right(const char* trans, i_type M, i_type& UD, V* A, i_type LDA);
        void        gemm_struct_left(i_type size, i_type NB_mult, const V* U, i_type LDU, const V* A, i_type LDA, 
                        V* C, i_type LDC) const;
        void        gemm_struct_right(i_type UD, const char* trand, i_type LEN, i_type size, const V* A, i_type LDA,
                        const V* U, i_type LDU, V* C, i_type LDC) const;
        void        apply_left_left_block(i_type i0, i_type i1, V* A, i_type LDA, matcl::workspace_pool& wp) const;
        void        apply_right_block(i_type i0, i_type i1, const char* tr, i_type UD, V* A, i_type LDA,
                        workspace_pool& wp) const;
        i_type      get_mult_workspace_size() const;
        i_type      get_mult_workspace_block() const;
        V*          get_workspace_left(i_type i) const;
        V*          get_workspace_right(i_type i) const;
        i_type      max_rows_per_task() const;
};

template<class V>
givens_accumulator<V>::givens_accumulator(i_type N, i_type NB, i_type NBM, i_type n_threads)
    :m_NB(NB), m_WORK(nullptr), m_WORK_mult(nullptr), m_N_max(N), m_NBM_left(NBM), m_left(false)
    ,n_workspaces_max(n_threads), m_N_current(N)
{    
    set_block_params(N);
};

template<class V>
void givens_accumulator<V>::set_block_params(i_type N)
{
    n_blocks        = std::max((N-1)/m_NB - 1, 0);
    m_NB_last       = N - 1 - n_blocks*m_NB;
    m_N_current     = N;
    n_workspaces    = std::max(1,std::min(n_blocks, n_workspaces_max));
};

template<class V>
i_type givens_accumulator<V>::work_size() const
{
    i_type lmwork   = get_mult_workspace_size();        
    i_type lwork    = 4 * m_NB * m_NB * n_blocks + m_NB_last * m_NB_last;
    return lwork + lmwork;
};

template<class V>
i_type givens_accumulator<V>::get_mult_workspace_size() const
{
    return std::max(m_NBM_left, m_N_max - 1 + n_workspaces)
                * get_mult_workspace_block();
};

template<class V>
i_type givens_accumulator<V>::max_rows_per_task() const
{
    return (m_N_max - 1)/n_workspaces + 1;
};

template<class V>
i_type givens_accumulator<V>::get_mult_workspace_block() const
{
    return std::max(m_NB_last, 2*m_NB);
};

template<class V>
V* givens_accumulator<V>::get_workspace_left(i_type i) const
{
    i_type work_size    = get_mult_workspace_block() * m_NBM_left;
    i_type off          = i * work_size;

    //i_type max_work_size    = m_WORK - m_WORK_mult;
    //if ((i+1) * work_size > max_work_size)
    //    throw;

    return m_WORK_mult + off;
};

template<class V>
V* givens_accumulator<V>::get_workspace_right(i_type i) const
{
    i_type work_size    = get_mult_workspace_block() * max_rows_per_task();
    i_type off          = i * work_size;

    //i_type max_work_size    = m_WORK - m_WORK_mult;
    //if ((i+1) * work_size > max_work_size)
    //    throw;

    return m_WORK_mult + off;
};

template<class V>
void givens_accumulator<V>::set_work(V* WORK)
{
    i_type lmwork   = get_mult_workspace_size();
    m_WORK_mult     = WORK;
    m_WORK          = WORK + lmwork;
};

template<class V>
V* givens_accumulator<V>::get_U(i_type num) const
{
    if (num < 0)
        return get_U_last();
    else
        return m_WORK + num * m_NB * m_NB * 4;
};

template<class V>
V* givens_accumulator<V>::get_U_last() const
{
    return m_WORK + n_blocks * m_NB * m_NB * 4;
};

template<class V>
void givens_accumulator<V>::initialize(V* U, i_type K)
{
    lapack::laset("F", K, K, V(0.0), V(1.0), U, K);
};

template<class V>
void givens_accumulator<V>::accumulate(bool left, i_type N, i_type NJ, const VR* SEQ_C, const V* SEQ_S, 
                        const i_type* SEQ_I1, const i_type* SEQ_I2, const i_type* SEQ_J)
{
    m_left = left;

    // set block parameters
    set_block_params(N);

    // initialize blocks to unit matrices
    for (i_type i = 0; i < n_blocks; ++i)
    {
        V* U    = get_U(i);
        initialize(U, 2 * m_NB);
    };

    V* Ul       = get_U_last();
    initialize(Ul, m_NB_last);

    //accumulate rotations
    for (i_type j = 0; j < NJ; ++j)
    {
        i_type k1       = SEQ_J[j];
        i_type k2       = SEQ_J[j + 1];
        for (i_type i = k1; i < k2; ++i)
        {
            i_type i1   = SEQ_I1[i];
            i_type i2   = SEQ_I2[i];
            VR c        = SEQ_C[i];
            V s         = SEQ_S[i];

            i_type num, off;
            get_block_position(i1, i2, j + 1, num, off);

            V* U        = get_U(num);
            i_type K    = num < 0 ? m_NB_last : 2 * m_NB;
            add_rotation(left, U, K, i1 - off, i2 - off, i+1, c, s);
        };
    };
};

template<class V>
i_type givens_accumulator<V>::get_size(i_type bl) const
{
    if (bl < 0)
        return m_NB_last;
    else
        return m_NB * 2;
};

template<class V>
i_type givens_accumulator<V>::get_offset(i_type bl) const
{
    i_type off_last = m_N_current - m_NB_last;

    if (bl < 0)
        return off_last;
    else
        return off_last - (bl + 1) * m_NB;
};

template<class V>
void givens_accumulator<V>::get_block_position(i_type i1, i_type i2, i_type j, i_type& block_num, i_type& offset) const
{
    i_type min_last = m_N_current - m_NB_last + 1;

    i_type ind_min  = std::min(i1, i2) - (j-1);
    i_type ind_max  = std::max(i1, i2);

    if (ind_min >= min_last)
    {
        block_num   = -1;
        offset      = min_last - 1;
        return;
    };

    i_type ind_min2 = min_last - ind_min - 1;
    block_num       = ind_min2/m_NB;
    offset          = min_last - 1 - (block_num+1) * m_NB;
};

template<class V>
void givens_accumulator<V>::add_rotation(bool left, V* U, i_type K, i_type i1, i_type i2, i_type j, VR c, V s)
{
    V val_j         = V(VR(j));

    V s1            = s;
    V s2            = - conj(s);

    if (left == true)
    {
        for (i_type i = 0; i < K; ++i)
        {
            V x         = U[i1-1];
            V y         = U[i2-1];

            V x2        = c  * x + s1 * y;
            V y2        = s2 * x + c  * y;

            U[i1-1]     = x2;
            U[i2-1]     = y2;

            U           += K;
        };
    }
    else
    {
        V* U1           = U + (i1 - 1) * K;
        V* U2           = U + (i2 - 1) * K;

        for (i_type i = 0; i < K; ++i)
        {
            V x         = U1[i];
            V y         = U2[i];

            V x2        = c  * x + s1 * y;
            V y2        = s2 * x + c  * y;

            U1[i]       = x2;
            U2[i]       = y2;
        };
    };
};

template<class V>
void givens_accumulator<V>::apply_left_right(i_type M, i_type& UD, V* A, i_type LDA)
{
    const char* tr      = get_trans_char<V>::eval("conj tr");
    return apply_right(tr, M, UD, A, LDA);
};

template<class V>
void givens_accumulator<V>::apply_right_right(i_type M, i_type& UD, V* A, i_type LDA)
{
    return apply_right("no", M, UD, A, LDA);
};

template<class V>
void givens_accumulator<V>::apply_left_left_block(i_type i0, i_type i1, V* A, i_type LDA, workspace_pool& wp) const
{
    auto work_ptr           = wp.get();
    V* work                 = reinterpret_cast<V*>(work_ptr.get());

    for (i_type i = i0; i < i1;)
    {
        V* U                = get_U(-1);
        i_type off          = get_offset(-1);
        i_type size         = get_size(-1);
        i_type NB_mult      = std::min(m_NBM_left, i1 - i);

        lapack::lacpy("F", size, NB_mult, A + off, LDA, work, size);
        lapack::gemm2<V>::eval("No", "No", size, NB_mult, size, V(1.0), U, size, work, size, V(0.0), A + off, LDA);

        for (i_type bl = 0; bl < n_blocks; ++bl)
        {
            U               = get_U(bl);
            off             = get_offset(bl);
            size            = get_size(bl);

            lapack::lacpy("F", size, NB_mult, A + off, LDA, work, size);
            gemm_struct_left(size, NB_mult, U, size, work, size, A + off, LDA);
        };

        i               += m_NBM_left;
        A               += m_NBM_left * LDA;
    };
};

template<class V>
void givens_accumulator<V>::apply_right_block(i_type i0, i_type i1, const char* tr, i_type UD, V* A, i_type LDA,
                                              workspace_pool& wp) const
{
    auto work_ptr       = wp.get();
    V* work             = reinterpret_cast<V*>(work_ptr.get());

    V* U                = get_U(-1);
    i_type off          = get_offset(-1);
    i_type size         = get_size(-1);
    i_type NB_mult      = i1 - i0;

    //first and last row positions in column off + 1 and off + size
    i_type F            = off + 1 - UD;
    i_type L            = i1;

    //first and last row in given band
    i_type FP           = std::max(F, i0 + 1);
    i_type LP           = std::min(L, i1);

    i_type LEN          = LP - FP + 1;

    if (LEN > 0)
    {
        lapack::lacpy("F", LEN, size, A + FP - 1 + off * LDA, LDA, work, NB_mult);
        lapack::gemm2<V>::eval("No", tr, LEN, size, size, V(1.0), work, NB_mult, U, size, V(0.0), 
                        A + FP - 1 + off * LDA, LDA);
    };

    for (i_type bl = 0; bl < n_blocks; ++bl)
    {
        U               = get_U(bl);
        off             = get_offset(bl);
        size            = get_size(bl);

        //first and last row positions in column off + 1 and off + size
        i_type F        = off + 1 - UD;
        i_type L        = i1;

        //first and last row in given band
        i_type FP       = std::max(F, i0 + 1);
        i_type LP       = std::min(L, i1);

        LEN             = LP - FP + 1;

        if (LEN > 0)
        {
            lapack::lacpy("F", LEN, size, A + FP - 1 + off * LDA, LDA, work, NB_mult);
            gemm_struct_right(UD + m_NB + i0, tr, LEN, size, work, NB_mult, U, size, A + FP - 1 + off * LDA, LDA);
        };
    };
};

template<class V>
void givens_accumulator<V>::apply_left_left(i_type M, V* A, i_type LDA)
{
    auto comp_f         = [this](i_type i0, i_type i1, V* A, i_type LDA, workspace_pool& wp) -> void
                        { 
                            return this->apply_left_left_block(i0, i1, A, LDA, wp);
                        };

    i_type n_threads    = n_workspaces;

    workspace_pool wp;

    for (i_type i = 0; i < n_threads; ++i)
        wp.add(get_workspace_left(i));

    if (n_threads == 1)
    {
        comp_f(0, M, A, LDA, wp);
        return;
    };

    i_type n_bl_per_task    = (n_blocks-1)/n_threads + 1;

    tast_group tg;

    for (i_type i = 0; i < M;)
    {
        i_type i0       = i;
        i_type i1       = std::min(i + n_bl_per_task * m_NBM_left, M);
        auto comp_0     = [i0, i1, A, LDA, comp_f, &wp]()-> void {return comp_f(i0, i1, A, LDA, wp); };

        tg.add(comp_0);

        i               += m_NBM_left * n_bl_per_task;
        A               += m_NBM_left * n_bl_per_task * LDA;
    };

    tg.wait();
};

template<class V>
void givens_accumulator<V>::apply_right(const char* tr, i_type M, i_type& UD, V* A, i_type LDA)
{    
    auto comp_f         = [this](i_type i0, i_type i1, const char* tr, i_type UD, V* A, 
                                    i_type LDA, workspace_pool& wp) -> void
                        { 
                            return this->apply_right_block(i0, i1, tr, UD, A, LDA, wp);
                        };

    i_type n_threads    = n_workspaces;

    workspace_pool wp;

    for (i_type i = 0; i < n_threads; ++i)
        wp.add(get_workspace_right(i));

    if (n_threads == 1)
    {
        comp_f(0, M, tr, UD, A, LDA, wp);
        UD                  += m_NB;
        return;
    };

    i_type n_bl_per_task    = (M-1)/n_threads + 1;

    tast_group tg;

    i_type pos          = 0;
    i_type UD2          = UD;
    for (i_type i = 0; i < M; ++ pos)
    {
        i_type i0       = i;
        i_type i1       = std::min(i + n_bl_per_task, M);
        auto comp_0     = [i0, i1, tr, UD2, A, LDA, comp_f, &wp]()-> void 
                        { return comp_f(i0, i1, tr, UD2, A, LDA, wp); };

        //comp_0();
        tg.add(comp_0);

        i               += n_bl_per_task;
    };

    tg.wait();

    UD                  += m_NB;
};

template<class V>
void givens_accumulator<V>::gemm_struct_left(i_type size, i_type NB_mult, const V* U, i_type LDU, const V* A, i_type LDA, 
                V* C, i_type LDC) const
{
    lapack::gemm_band_l<V>::eval(size/2, size/2, size, NB_mult, size, V(1.0), U, LDU, A, LDA, V(0.0), C, LDC);
    return;
};

template<class V>
void givens_accumulator<V>::gemm_struct_right(i_type UD, const char* trans, i_type LEN, i_type size, const V* A, 
        i_type LDA, const V* U, i_type LDU, V* C, i_type LDC) const
{
    lapack::gemm_band_r<V>::eval(LEN, UD, size/2, size/2, "No", trans, LEN, size, size, V(1.0), A, LDA, U, LDU,
                                 V(0.0), C, LDC);
    return;
};

};}
