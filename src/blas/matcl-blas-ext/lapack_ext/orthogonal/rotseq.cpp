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

#include "matcl-blas-lapack/blas/details/blas_simd_complex.h"
#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "blas/matcl-blas-ext/blas_ext/small_gemm.inl"
#include "matcl-blas-lapack/blas_ext/rot.inl"
#include "blas/matcl-blas-ext/lapack_ext/utils/optim_params.h"
#include "matcl-blas-ext/blas_concurrency.h"
#include "blas/matcl-blas-ext/timer.h"

#include <algorithm>
#include <atomic>
#include <iostream>

namespace matcl { namespace lapack
{

//maximum size of accumulated rotations; cannot be changed without code modifications
static const i_type ROT_BLOCK_SIZE          = 4;

static const i_type min_length_accumulate   = optim_params::rotseq_min_length_to_acc;
static const i_type panel_col_size_max      = optim_params::rotseq_panel_col_size_max;
static const i_type panel_col_size_min      = optim_params::rotseq_panel_col_size_min;

enum class op_s_type
{
    type_s, type_min_conj, type_min, type_conj
};

static op_s_type get_op_s(const char* TRANS, bool min_conj)
{
    if (TRANS[0] == 'n' || TRANS[0] == 'N')
        return min_conj? op_s_type::type_min_conj : op_s_type::type_s;
    else if (TRANS[0] == 'z' || TRANS[0] == 'Z')
        return min_conj? op_s_type::type_min : op_s_type::type_conj;
    else if (TRANS[0] == 't' || TRANS[0] == 'T')
        return min_conj? op_s_type::type_s : op_s_type::type_min_conj;
    else
        return min_conj? op_s_type::type_conj : op_s_type::type_min;
};

template<class V>
V convert_s(const V& s, op_s_type op_s)
{
    switch (op_s)
    {
        case op_s_type::type_s:
            return s;
        case op_s_type::type_min_conj:
            return -conj(s);
        case op_s_type::type_min:
            return -s;
        case op_s_type::type_conj:
            return conj(s);
    };
    return s;
};

static void calc_interval_right(i_type& F, i_type& L, i_type M, i_type ind1, i_type ind2, i_type& LDIAGS, 
                                i_type& UDIAGS, i_type& LDIAGS_0, i_type& UDIAGS_0, i_type& MIN_MOD, i_type& MAX_MOD,
                                bool& need_init)
{
    i_type IMIN         = std::min(ind1, ind2);
    i_type IMAX         = std::max(ind1, ind2);

    i_type F0_1         = std::max(IMIN - UDIAGS, 1);
    i_type F0_2         = std::max(std::min(IMIN - UDIAGS_0, MIN_MOD), 1);

    i_type L0_1         = std::min(IMAX + LDIAGS, M);
    i_type L0_2         = std::min(std::max(IMAX + LDIAGS_0, MAX_MOD), M);

    //check if estimation based on current LDIAGS and UDIAGS is better
    //than based on modified interval
    if ((F0_1 > F0_2 && L0_1 <= L0_2) || (F0_1 == F0_2 && L0_1 < L0_2))
    {
        //reset interval and update diags estimator
        LDIAGS_0        = LDIAGS;
        UDIAGS_0        = UDIAGS;

        MIN_MOD         = M;
        MAX_MOD         = 0;
        need_init       = true;
    };

    F                   = std::max(F0_1, F0_2);
    L                   = std::min(L0_1, L0_2);

    if (F > L)
        return;

    // update estimation of ldiags and udiags
    UDIAGS      = std::max(UDIAGS, IMAX - F);
    LDIAGS      = std::max(LDIAGS, L - IMIN);

    // update bounds on modified rows
    if (need_init == true)
    {
        MIN_MOD         = F;
        MAX_MOD         = L;
        need_init       = false;
    }
    else
    {
        MIN_MOD         = std::min(MIN_MOD,F);
        MAX_MOD         = std::max(MAX_MOD,L);
    };
};

static void calc_interval_right(i_type& LEN, i_type& F, i_type M, i_type len, i_type j, i_type j_step, 
            const i_type* IND1, const i_type* IND2, i_type& LDIAGS, i_type& UDIAGS, i_type& LDIAGS_0, 
            i_type& UDIAGS_0, i_type& MIN_MOD, i_type& MAX_MOD, bool& need_init)
{
    F                   = M + 1;
    i_type L            = 0;
    for (i_type i = 0; i < len; ++i, j += j_step)
    {
        i_type Fl, Ll;

        i_type ind1     = IND1[j];
        i_type ind2     = IND2[j];
        calc_interval_right(Fl, Ll, M, ind1, ind2, LDIAGS, UDIAGS, LDIAGS_0, UDIAGS_0, MIN_MOD, MAX_MOD, need_init);

        F               = std::min(F, Fl);
        L               = std::max(L, Ll);
    };

    if (F > L)
        LEN             = 0;
    else
        LEN             = L - F + 1;
};

static void calc_interval_left(i_type& F, i_type& L, i_type N, i_type ind1, i_type ind2, i_type& LDIAGS, 
                                i_type& UDIAGS, i_type& LDIAGS_0, i_type& UDIAGS_0, i_type& MIN_MOD, i_type& MAX_MOD,
                                bool& need_init)
{
    // calc interval to modify and its length
    i_type IMIN         = std::min(ind1, ind2);
    i_type IMAX         = std::max(ind1, ind2);
   
    i_type F0_1         = std::max(IMIN - LDIAGS, 1);
    i_type F0_2         = std::max(std::min(IMIN - LDIAGS_0, MIN_MOD), 1);

    i_type L0_1         = std::min(IMAX + UDIAGS, N);
    i_type L0_2         = std::min(std::max(IMAX + UDIAGS_0, MAX_MOD), N);

    //check if estimation based on current LDIAGS and UDIAGS is better
    //than based on modified interval
    if ((F0_1 > F0_2 && L0_1 <= L0_2) || (F0_1 == F0_2 && L0_1 < L0_2))
    {
        //reset interval and update diags estimator
        LDIAGS_0        = LDIAGS;
        UDIAGS_0        = UDIAGS;

        MIN_MOD         = N;
        MAX_MOD         = 0;
        need_init       = true;
    };

    F                   = std::max(F0_1, F0_2);
    L                   = std::min(L0_1, L0_2);    

    if (F > L)
        return;

    // update estimation of ldiags and udiags
    LDIAGS              = std::max(LDIAGS, IMAX - F);
    UDIAGS              = std::max(UDIAGS, L - IMIN);

    // update bounds on modified columns
    if (need_init == true)
    {
        MIN_MOD         = F;
        MAX_MOD         = L;
        need_init       = false;
    }
    else
    {
        MIN_MOD         = std::min(MIN_MOD,F);
        MAX_MOD         = std::max(MAX_MOD,L);
    };
};

static void calc_interval_left(i_type& LEN, i_type& F, i_type N, i_type len, i_type j, i_type j_step, 
            const i_type* IND1, const i_type* IND2, i_type& LDIAGS, i_type& UDIAGS, i_type& LDIAGS_0, 
            i_type& UDIAGS_0, i_type& MIN_MOD, i_type& MAX_MOD, bool& need_init)
{
    F                   = N + 1;
    i_type L            = 0;
    for (i_type i = 0; i < len; ++i, j += j_step)
    {
        i_type Fl, Ll;

        i_type ind1     = IND1[j];
        i_type ind2     = IND2[j];
        calc_interval_left(Fl, Ll, N, ind1, ind2, LDIAGS, UDIAGS, LDIAGS_0, UDIAGS_0, MIN_MOD, MAX_MOD, need_init);

        F               = std::min(F, Fl);
        L               = std::max(L, Ll);
    };

    if (F > L)
        LEN             = 0;
    else
        LEN             = L - F + 1;
};

static void check_rotation_index(i_type ind1, i_type ind2, i_type N, i_type& INFO)
{
    if (ind1 <= 0 || ind1 > N)
    {
        INFO        = -7;
    };
    if (ind2 <= 0 || ind2 > N)
    {
        INFO        = -8;
    };
};

static void accumulate_rotations_select_l(i_type& ind_1, i_type& ind_2, i_type& len, i_type j, 
                i_type j_step, i_type max_len, i_type N, i_type N_vec, const i_type* IND1, 
                const i_type* IND2, i_type& INFO)
{
    if (N_vec < min_length_accumulate)
    {
        len             = 1;
        ind_1           = IND1[j];
        ind_2           = IND2[j];

        check_rotation_index(ind_1, ind_2, N, INFO);
        return;
    };

    len                 = 0;
    i_type min_ind      = N;
    i_type max_ind      = 0;
    i_type j0           = j;

    for (i_type i = 0; i < max_len; ++i, j += j_step)
    {
        // get rotation info
        i_type ind1     = IND1[j];
        i_type ind2     = IND2[j];

        // check rotation
        check_rotation_index(ind1, ind2, N, INFO);
        if (INFO != 0)
            return;

        i_type min_sav  = min_ind;
        i_type max_sav  = max_ind;
        min_ind         = std::min(std::min(ind1, ind2), min_ind);
        max_ind         = std::max(std::max(ind1, ind2), max_ind);

        i_type dif      = max_ind - min_ind + 1;

        if (dif > ROT_BLOCK_SIZE)
        {
            if (len <= 1)
            {
                len     = 1;
                ind_1   = IND1[j0];
                ind_2   = IND2[j0];
                return;
            }
            else
            {
                ind_1   = min_sav;
                ind_2   = max_sav;
                return;
            };
        };

        ++len;
    };

    if (len <= 1)
    {
        len             = 1;
        ind_1           = IND1[j0];
        ind_2           = IND2[j0];
    }
    else
    {
        ind_1           = min_ind;
        ind_2           = max_ind;
    };
    return;
};
static void accumulate_rotations_select_r(i_type& ind_1, i_type& ind_2, i_type& ind_3, i_type& len, 
                bool& rev1, bool& rev2, i_type j, i_type j_step, i_type max_len, i_type N, 
                const i_type* IND1, const i_type* IND2, i_type& INFO)
{
    len                 = 0;
    rev1                = false;
    rev2                = false;
    ind_3               = -1;
    i_type j0           = j;

    i_type ind11        = IND1[j0];
    i_type ind21        = IND2[j0];

    ind_1               = ind11;
    ind_2               = ind21;

    check_rotation_index(ind11, ind21, N, INFO);

    if (max_len < 2)
    {
        len             = 1;
        return;
    };

    i_type ind12        = IND1[j0 + j_step];
    i_type ind22        = IND2[j0 + j_step];
    
    check_rotation_index(ind12, ind22, N, INFO);

    i_type ind1_min     = std::min(ind11, ind21);
    i_type ind1_max     = std::max(ind11, ind21);
    i_type ind2_min     = std::min(ind12, ind22);
    i_type ind2_max     = std::max(ind12, ind22);
   
    if (ind1_max == ind2_min)
    {
        len             = 2;
        ind_1           = ind1_min;
        ind_2           = ind1_max;
        ind_3           = ind2_max;
    }
    else if (ind1_min == ind2_max)
    {
        len             = 2;
        ind_1           = ind1_max;
        ind_2           = ind1_min;
        ind_3           = ind2_min;
    }
    else
    {
        len             = 1;
        return;
    };

    if (ind_1 != ind11)
        rev1            = true;

    if (ind_2 != ind12)
        rev2            = true;
};

template<class VM, class V, class VR>
static void accumulate_rotations_collect_left(VM* mat, i_type& SIZE, i_type ind1, i_type ind2, i_type len, i_type j, 
                i_type j_step, const i_type* IND1, const i_type* IND2, const VR* C, const V* S, op_s_type op_s)
{
    SIZE                = (ind2 - ind1) + 1;

    for (i_type i = 0; i < SIZE*SIZE; ++i)
        mat[i]          = V(0.0);

    for (i_type i = 0; i < SIZE*SIZE; i += SIZE + 1)
        mat[i]          = V(1.0);

    j                   = j + (len-1)*j_step;

    for (i_type i = len - 1; i >= 0; --i, j -= j_step)
    {
        VR c            = C[j];
        V s             = convert_s(S[j],op_s);
        i_type i1       = IND1[j] - ind1;
        i_type i2       = IND2[j] - ind1;

        lapack::rot<VM>(SIZE, mat + i1 * SIZE, 1, mat + i2 * SIZE, 1, c, -conj(s));
    };
};

template<class V>
void apply_rotations_left(i_type LEN, i_type SIZE, V* R, i_type LDR, const V* M, i_type LDM)
{    
    if (SIZE == 3)
        gemm_S3_M<V>::eval(LEN, M, SIZE, R, LDR);
    else
        gemm_S4_M<V>::eval(LEN, M, SIZE, R, LDR);
};

template<class V1, class V2>
static void apply_seq_1(i_type i_start, i_type i_step, i_type K, i_type M, const typename details::real_type<V1>::type* C, 
                        const V1* S, const i_type* IND1, const i_type* IND2, i_type& LDIAGS, i_type& UDIAGS, 
                        V2* R, i_type LDR, op_s_type op_s, i_type& INFO)
{
    using VR    = typename details::real_type<V1>::type;
    using VR2   = typename details::real_type<V2>::type;

    for (i_type i = 0, j = i_start; i < K; ++i, j += i_step)
    {
        i_type ind1, ind2;

        //get working interval
        ind1            = IND1[j];
        ind2            = IND2[j];

        if (LDIAGS < M - 1)
        {
            i_type ind_min  = std::min(ind1, ind2);
            i_type ind_max  = std::max(ind1, ind2);

            if (ind_min - 1 > LDIAGS)
                continue;

            LDIAGS      = std::max(LDIAGS, ind_max - 1);
        };

        VR2 c           = VR2(C[j]);
        V2 s            = V2(convert_s(S[j],op_s));

        // apply rotations
        V2 x            = R[ind1-1];
        V2 y            = R[ind2-1];
        V2 x2           = c * x + s * y;
        V2 y2           = c * y - conj(s) * x;

        R[ind1-1]       = x2;
        R[ind2-1]       = y2;
    };
};

static i_type sum_arith(i_type start, i_type steps0)
{
    if (start <= 0 || steps0 <= 0)
        return 0;

    i_type steps        = std::min(start, steps0);
    i_type end          = start - steps;
    return ((start + end) * steps)/2;
};
static i_type calc_nz(i_type M, i_type N, i_type LD, i_type UD)
{
    i_type nz   = M * N 
                - sum_arith(M - 1 - LD, std::min(M - 1, N - 1))
                - sum_arith(N - 1 - UD, std::min(M - 1, N - 1));
    return nz;
};

static void get_thread_block_params(i_type& n_blocks, i_type& NB, i_type& n_threads, i_type M, i_type N0, i_type LD, 
                                    i_type UD, i_type NB_max, i_type NB_min)
{
    i_type nz               = calc_nz(M, N0, LD, UD);
    i_type N_eff            = nz / M;
    i_type max_threads      = get_num_threads(domain::matcl);
    i_type exp_threads_max  = (N_eff - 1)/NB_max + 1;
    i_type exp_threads_min  = std::max(N_eff/NB_min, 1);
    
    if (max_threads == 1 || exp_threads_min == 1)
        n_threads           = 1;
    else
        n_threads           = std::min(max_threads, (exp_threads_max + exp_threads_min)/2);

    i_type n_blocks_NB      = (N_eff - 1)/n_threads + 1;
    n_blocks                = (n_blocks_NB - 1)/NB_max + 1;
    NB                      = (n_blocks_NB - 1)/n_blocks + 1;
};

template<class V1, class V2>
static void eval_size_r(bool type_r, bool tr_no, bool tr_conj, i_type K, const char* TRANS, bool tr_tr, bool tr_ctr,
                        i_type& LDIAGS, i_type& UDIAGS, i_type M, i_type N, const typename details::real_type<V1>::type* C, 
                        const V1* S, const i_type* IND1, const i_type* IND2, V2* R, i_type LDR, i_type& INFO)
{
    using VR    = typename details::real_type<V1>::type;
    using VR2   = typename details::real_type<V2>::type;

    i_type i_start, i_step;
    op_s_type op_s;

    if (type_r == true)
    {
        if (tr_no == true || tr_conj == true)
        {
            i_start     = 0;
            i_step      = 1;
        }
        else
        {
            i_start     = K - 1;
            i_step      = -1;
        };

        op_s            = get_op_s(TRANS, false);
    }
    else
    {
        if (tr_tr == true || tr_ctr == true)
        {
            i_start     = 0;
            i_step      = 1;
        }
        else
        {
            i_start     = K - 1;
            i_step      = -1;
        };

        op_s            = get_op_s(TRANS, true);
    };

    i_type LDIAGS_0     = LDIAGS;
    i_type UDIAGS_0     = UDIAGS;    
    i_type LDIAGS_s     = LDIAGS;
    i_type UDIAGS_s     = UDIAGS;       

    // uninitialized bound on modified rows
    i_type MIN_MOD      = M;
    i_type MAX_MOD      = 0;
    bool need_init      = true;

    for (i_type i = 0, j = i_start; i < K;)
    {
        i_type ind1, ind2, ind3, len;
        bool rev1, rev2;

        //get working interval
        accumulate_rotations_select_r(ind1, ind2, ind3, len, rev1, rev2, j, i_step, K - i, N, 
                                        IND1, IND2, INFO);

        if (INFO != 0 || len == 0)
            return;

        // calc interval to modify and its length
        i_type F;
        i_type LEN;
        calc_interval_right(LEN, F, M, len, j, i_step, IND1, IND2, LDIAGS, UDIAGS, LDIAGS_0, 
                            UDIAGS_0, MIN_MOD, MAX_MOD, need_init);

        if (LEN > 0)
        {
            if (len == 1)
            {
                VR c        = C[j];
                V1 s        = convert_s(S[j],op_s);

                // apply rotations
                lapack::rot(LEN, R + F-1 + (ind1-1)*LDR, 1, R + F-1 + (ind2-1)*LDR, 1, VR2(c), V2(s));
            }
            else
            {
                //len = 2;
                VR c1       = C[j];
                V1 s1       = convert_s(S[j],op_s);
                VR c2       = C[j + i_step];
                V1 s2       = convert_s(S[j + i_step],op_s);

                if (rev1)
                    s1      = -conj(s1);
                if (rev2)
                    s2      = -conj(s2);

                // apply rotations
                lapack::rot_2<V2>::eval(LEN, R + F-1 + (ind1-1)*LDR, R + F-1 + (ind2-1)*LDR, 
                        R + F-1 + (ind3-1)*LDR, VR2(c1), V2(s1), VR2(c2), V2(s2));
            };
        };

        i                   += len;
        j                   += i_step*len;
    };

    return;
};

template<class V1, class V2>
static void eval_size_l(bool type_r, bool tr_no, bool tr_conj, i_type K, const char* TRANS, bool tr_tr, bool tr_ctr,
                        i_type& LDIAGS, i_type& UDIAGS, i_type M, i_type N, const typename details::real_type<V1>::type* C, 
                        const V1* S, const i_type* IND1, const i_type* IND2, V2* R, i_type LDR, i_type& INFO)
{
    using VR    = typename details::real_type<V1>::type;
    using VR2   = typename details::real_type<V2>::type;

    i_type i_start, i_step;
    op_s_type op_s;

    if (type_r == false)
    {
        if (tr_no == true || tr_conj == true)
        {
            i_start     = 0;
            i_step      = 1;
        }
        else
        {
            i_start     = K - 1;
            i_step      = -1;
        };

        op_s            = get_op_s(TRANS, false);
    }
    else
    {
        if (tr_tr == true || tr_ctr == true)
        {
            i_start     = 0;
            i_step      = 1;
        }
        else
        {
            i_start     = K - 1;
            i_step      = -1;
        };

        op_s            = get_op_s(TRANS, true);
    };

    if (N == 1)
        return apply_seq_1(i_start, i_step, K, M, C, S, IND1, IND2, LDIAGS, UDIAGS, R, LDR, op_s, INFO);

    i_type LDIAGS_s     = LDIAGS;
    i_type UDIAGS_s     = UDIAGS;

    auto task           = [N, i_start, i_step, op_s, K, M, C, S, IND1, IND2, LDIAGS_s, UDIAGS_s]
                        (i_type l1, i_type l2, i_type NB, V2* R, i_type LDR, i_type& LDIAGS_r, i_type& UDIAGS_r)->void
    {
        LDIAGS_r            = LDIAGS_s;
        UDIAGS_r            = UDIAGS_s;
        i_type INFO         = 0;
        R                   += l1 * LDR;

        V2 mat[ROT_BLOCK_SIZE * ROT_BLOCK_SIZE];

        for (i_type l = l1; l < l2;)
        {
            // uninitialized bound on modified columns
            i_type MIN_MOD  = N;
            i_type MAX_MOD  = 0;
            i_type NBJ      = (l2 - l)/NB >= 2 ? NB : l2 - l;

            i_type LDIAGS   = LDIAGS_s + l;
            i_type UDIAGS   = UDIAGS_s - l;
            i_type LDIAGS_0 = LDIAGS;
            i_type UDIAGS_0 = UDIAGS;
            bool need_init  = true;                

            for (i_type i = 0, j = i_start; i < K;)
            {
                i_type ind1, ind2, len;

                //get working interval
                accumulate_rotations_select_l(ind1, ind2, len, j, i_step, K - i, M, NBJ, IND1, IND2, INFO);

                if (INFO != 0 || len == 0)
                    return;

                // calc interval to modify and its length
                i_type LEN, F;
                calc_interval_left(LEN, F, NBJ, len, j, i_step, IND1, IND2, LDIAGS, UDIAGS, LDIAGS_0, UDIAGS_0,
                                    MIN_MOD, MAX_MOD, need_init);

                if (LEN > 0)
                {
                    if (len == 1)
                    {
                        VR c        = C[j];
                        V1 s        = convert_s(S[j],op_s);

                        // apply rotations
                        lapack::rot(LEN, R + ind1-1 + (F-1)*LDR, LDR, R + ind2-1 + (F-1)*LDR, LDR, VR2(c), V2(s));
                    }
                    else
                    {
                        i_type SIZE;
                        accumulate_rotations_collect_left(mat, SIZE, ind1, ind2, len, j, i_step, IND1, IND2, C, S, op_s);
                        apply_rotations_left(LEN, SIZE, R + ind1-1 + (F-1)*LDR, LDR, mat, SIZE);
                    };
                };

                i               += len;
                j               += i_step*len;                    
            };

            LDIAGS_r            = std::max(LDIAGS_r, LDIAGS - l);
            UDIAGS_r            = std::max(UDIAGS_r, UDIAGS + l);

            R                   += NBJ * LDR;
            l                   += NBJ;
        };
    };
        
    i_type NB_max               = panel_col_size_max;
    i_type NB_min               = panel_col_size_min;
        
    i_type n_blocks, NB, n_threads;

    get_thread_block_params(n_blocks, NB, n_threads, M, N, LDIAGS, UDIAGS, NB_max, NB_min);

    if (n_threads == 1)
    {
        task(0, N, NB, R, LDR, LDIAGS, UDIAGS);
        return;
    };

    using atomic_int        = std::atomic<int>;

    atomic_int LDIAGS_r     = LDIAGS;
    atomic_int UDIAGS_r     = UDIAGS;
    i_type exp_nz           = n_blocks*NB*M;

    tast_group tg;

    for (i_type i = 0, j = 0; i < N; ++j)
    {
        i_type i1           = i;
        i_type i2           = std::min(i + n_blocks * NB, N);
        i_type block_nz     = calc_nz(M, i2 - i1 + 1, LDIAGS + i, UDIAGS - i);

        while(double(block_nz)/double(exp_nz) < 0.75 && i2 < N)
        {
            i2              = std::min(i2 + NB, N);
            block_nz        = calc_nz(M, i2 - i1 + 1, LDIAGS + i, UDIAGS - i);
        }
            
        auto comp_0         = [i1, i2, task, NB, R, LDR, &LDIAGS_r, &UDIAGS_r]()-> void 
        { 
            i_type LDIAGS_l;
            i_type UDIAGS_l;
            task(i1, i2, NB, R, LDR, LDIAGS_l, UDIAGS_l); 

            i_type LDIAGS_r0    = LDIAGS_r;
            i_type UDIAGS_r0    = UDIAGS_r;

            LDIAGS_r0           = std::max(LDIAGS_r0, LDIAGS_l);
            UDIAGS_r0           = std::max(UDIAGS_r0, UDIAGS_l);

            LDIAGS_r            = LDIAGS_r0;
            UDIAGS_r            = UDIAGS_r0;
        };
            
        tg.add(comp_0);

        i                       = i2;
    };

    tg.wait();

    LDIAGS                  = LDIAGS_r;
    UDIAGS                  = UDIAGS_r;
};

template<class V1, class V2> BLAS_EXT_EXPORT
typename details::enable_if_valid2<void,V1, V2>::type
lapack::rotseq(const char* TYPE, const char* SIDE, const char* TRANS, i_type K, 
    const typename details::real_type<V1>::type* C, const V1* S, const i_type* IND1, const i_type* IND2, 
    i_type M, i_type N, i_type& LDIAGS, i_type& UDIAGS, V2* R, i_type LDR, i_type& INFO)
{
    using VR    = typename details::real_type<V1>::type;
    using VR2   = typename details::real_type<V2>::type;

    bool type_l     = (TYPE[0] == 'l' || TYPE[0] == 'L');
    bool type_r     = (TYPE[0] == 'r' || TYPE[0] == 'R');
    bool side_l     = (SIDE[0] == 'l' || SIDE[0] == 'L');
    bool side_r     = (SIDE[0] == 'r' || SIDE[0] == 'R');
    bool tr_no      = (TRANS[0] == 'n' || TRANS[0] == 'N');
    bool tr_conj    = (TRANS[0] == 'z' || TRANS[0] == 'Z');
    bool tr_tr      = (TRANS[0] == 't' || TRANS[0] == 'T');
    bool tr_ctr     = (TRANS[0] == 'c' || TRANS[0] == 'C');

    // test arguments    
    INFO        = 0;

    if (type_l == false && type_r == false)
        INFO    = -1;
    else if (side_r == false && side_l == false)
        INFO    = -2;
    else if (tr_no == false && tr_conj == false && tr_tr == false && tr_ctr == false)
        INFO    = -3;
    else if (K < 0)
        INFO    = -4;
    else if (M < 0)
        INFO    = -9;
    else if (N < 0)
        INFO    = -10;
    else if (LDIAGS < 0)
        INFO    = -11;
    else if (UDIAGS < 0)
        INFO    = -12;
    else if (LDR < std::max(M, 1))
        INFO    = -14;

    if (INFO != 0)
        return;    

    //fast exit
    if (M == 0 || N == 0 || K == 0)
        return;

    if (side_r == true)
    {
        return eval_size_r(type_r, tr_no, tr_conj, K, TRANS, tr_tr, tr_ctr, LDIAGS, UDIAGS, M, N, 
                           C, S, IND1, IND2, R, LDR, INFO);
    }
    else if (side_r == false)
    {
        return eval_size_l(type_r, tr_no, tr_conj, K, TRANS, tr_tr, tr_ctr, LDIAGS, UDIAGS, M, N, 
                           C, S, IND1, IND2, R, LDR, INFO);
    };
};

template void BLAS_EXT_EXPORT
rotseq<d_type, d_type>(const char* TYPE, const char* SIDE, const char* TRANS, i_type K, 
       const d_type* C, const d_type* S, const i_type* IND1, const i_type* IND2,
       i_type M, i_type N, i_type& LDIAGS, i_type& UDIAGS, d_type* X, i_type LDX, i_type& INFO);

template void BLAS_EXT_EXPORT
rotseq<s_type, s_type>(const char* TYPE, const char* SIDE, const char* TRANS, i_type K, 
       const s_type* C, const s_type* S, const i_type* IND1, const i_type* IND2,
       i_type M, i_type N, i_type& LDIAGS, i_type& UDIAGS, s_type* X, i_type LDX, i_type& INFO);

template void BLAS_EXT_EXPORT
rotseq<s_type, d_type>(const char* TYPE, const char* SIDE, const char* TRANS, i_type K, 
       const s_type* C, const s_type* S, const i_type* IND1, const i_type* IND2,
       i_type M, i_type N, i_type& LDIAGS, i_type& UDIAGS, d_type* X, i_type LDX, i_type& INFO);

template void BLAS_EXT_EXPORT
rotseq<c_type, c_type>(const char* TYPE, const char* SIDE, const char* TRANS, i_type K, 
       const s_type* C, const c_type* S, const i_type* IND1, const i_type* IND2,
       i_type M, i_type N, i_type& LDIAGS, i_type& UDIAGS, c_type* X, i_type LDX, i_type& INFO);

template void BLAS_EXT_EXPORT
rotseq<z_type, z_type>(const char* TYPE, const char* SIDE, const char* TRANS, i_type K, 
       const d_type* C, const z_type* S, const i_type* IND1, const i_type* IND2,
       i_type M, i_type N, i_type& LDIAGS, i_type& UDIAGS, z_type* X, i_type LDX, i_type& INFO);

template void BLAS_EXT_EXPORT
rotseq<c_type, z_type>(const char* TYPE, const char* SIDE, const char* TRANS, i_type K, 
       const s_type* C, const c_type* S, const i_type* IND1, const i_type* IND2,
       i_type M, i_type N, i_type& LDIAGS, i_type& UDIAGS, z_type* X, i_type LDX, i_type& INFO);

template void BLAS_EXT_EXPORT
rotseq<d_type, z_type>(const char* TYPE, const char* SIDE, const char* TRANS, i_type K, 
       const d_type* C, const d_type* S, const i_type* IND1, const i_type* IND2,
       i_type M, i_type N, i_type& LDIAGS, i_type& UDIAGS, z_type* X, i_type LDX, i_type& INFO);

template void BLAS_EXT_EXPORT
rotseq<s_type, z_type>(const char* TYPE, const char* SIDE, const char* TRANS, i_type K, 
       const s_type* C, const s_type* S, const i_type* IND1, const i_type* IND2,
       i_type M, i_type N, i_type& LDIAGS, i_type& UDIAGS, z_type* X, i_type LDX, i_type& INFO);

template void BLAS_EXT_EXPORT
rotseq<s_type, c_type>(const char* TYPE, const char* SIDE, const char* TRANS, i_type K, 
       const s_type* C, const s_type* S, const i_type* IND1, const i_type* IND2,
       i_type M, i_type N, i_type& LDIAGS, i_type& UDIAGS, c_type* X, i_type LDX, i_type& INFO);

};};