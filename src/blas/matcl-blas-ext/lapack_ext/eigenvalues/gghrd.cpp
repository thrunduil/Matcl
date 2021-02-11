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
#include "blas/matcl-blas-ext/lapack_ext/eigenvalues/givens_accumulator.h"
#include "matcl-blas-ext//blas_concurrency.h"
#include "blas/matcl-blas-ext/lapack_ext/utils/utils.h"

#include <algorithm>

namespace matcl { namespace lapack
{

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
gghrd2_s(const char* COMPQ, const char *COMPZ, i_type N, i_type ILO, i_type IHI, V *A, i_type LDA, 
           V* B, i_type LDB, V* Q, i_type LDQ, V* Z, i_type LDZ, i_type &info)
{
    using VR    = details::real_type<V>::type;

    const V ZERO   = 0.;
    const V ONE    = 1.;

    bool ILQ    = false;
    bool ILZ    = false;
    int ICOMPQ  = 0;
    int ICOMPZ  = 0;

    // Decode C
    if( COMPQ[0] == 'N' || COMPQ[0] == 'n')
    {
        ILQ     = false;
        ICOMPQ  = 1;
    }
    else if (COMPQ[0] == 'V' || COMPQ[0] == 'v' )
    {
        ILQ     = true;
        ICOMPQ  = 2;
    }
    else if (COMPQ[0] == 'I' || COMPQ[0] == 'i')
    {
        ILQ     = true;
        ICOMPQ  = 3;
    }
    else
    {
        ICOMPQ  = 0;
    };

    // Decode COMPZ
    if (COMPZ[0] == 'N' || COMPZ[0] == 'n')
    {
        ILZ     = false;
        ICOMPZ  = 1;
    }
    else if (COMPZ[0] == 'V' || COMPZ[0] == 'v')
    {
        ILZ     = true;
        ICOMPZ  = 2;
    }
    else if (COMPZ[0] == 'I' || COMPZ[0] == 'i')
    {
        ILZ     = true;
        ICOMPZ  = 3;
    }
    else
    {
        ICOMPZ  = 0;
    };

    // Test the input parameters.
    int INFO = 0;

    if (ICOMPQ <= 0 )
        INFO    = -1;
    else if (ICOMPZ <= 0)
        INFO    = -2;
    else if (N < 0)
        INFO    = -3;
    else if (ILO < 1)
        INFO    = -4;
    else if (IHI > N || IHI < ILO-1 )
        INFO    = -5;
    else if (LDA < lapack::maximum( 1, N ) )
        INFO    = -7;
    else if (LDB < lapack::maximum( 1, N ) )
        INFO    = -9;
    else if ( ( ILQ && LDQ < N ) || LDQ < 1 )
        INFO    = -11;
    else if ( ( ILZ && LDZ < N ) || LDZ < 1 )
        INFO    = -13;

    if ( INFO != 0 )
    {
        info   = INFO;
        return;
    };

    // Initialize Q and Z if desired.
    if (ICOMPQ == 3)
        lapack::laset("Full", N, N, ZERO, ONE, Q, LDQ );

    if (ICOMPZ == 3)
        lapack::laset("Full", N, N, ZERO, ONE, Z, LDZ );

    // Quick return if possible
    if (N <= 1)
        return;

    // Zero out lower triangle of B
    V* tmp_B = B;

    for (i_type JCOL = 0; JCOL < N - 1; ++JCOL)
    {
         for (i_type JROW = JCOL + 1; JROW < N; ++JROW)
            tmp_B[JROW] = ZERO;

         tmp_B  = tmp_B + LDB;
    };

    // Reduce A and B

    VR C;
    V  S;

    for (i_type JCOL = ILO - 1; JCOL <= IHI - 1 - 2; ++JCOL)
    {        
        for (i_type JROW = IHI - 1; JROW >= JCOL + 2; --JROW)
        {

            // Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)
            V TEMP          = A[JROW - 1 + JCOL * LDA];

            lapack::lartg<V>(TEMP, A[JROW + JCOL * LDA], &C, &S, A + JROW - 1 + JCOL * LDA );

            if (C != VR(1.0))
            {
                A[JROW  + JCOL * LDA] = ZERO;

                lapack::rot<V>( N - JCOL - 1, A + JROW - 1 + (JCOL + 1) * LDA, LDA, 
                                A + JROW + (JCOL + 1) * LDA, LDA, C, S);
                lapack::rot<V>( N + 2 - JROW - 1, B + JROW - 1 + (JROW - 1) * LDB, LDB, 
                                B + JROW + (JROW - 1) * LDB, LDB, C, S );
                      
                if (ILQ)
                    lapack::rot<V>(N, Q + (JROW - 1) * LDQ, 1, Q + JROW * LDQ, 1, C, conj(S));

                // Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)
                TEMP            = B[JROW + JROW * LDB];
            
                lapack::lartg<V>(TEMP, B[JROW + (JROW - 1) * LDB], &C, &S, B + JROW + JROW * LDB );

                B[JROW + (JROW - 1) * LDB] = ZERO;

                lapack::rot<V>( IHI, A + JROW * LDA, 1, A + (JROW - 1) * LDA, 1, C, S );
                lapack::rot<V>(JROW, B + JROW * LDB, 1, B + (JROW - 1) * LDB, 1, C, S );

                if (ILZ)
                    lapack::rot<V>(N, Z + JROW * LDZ, 1, Z + (JROW - 1) * LDZ, 1, C, S );
            };
        };
    };

    return;
};

class profiler
{
    private:
        //blas::timer m_t;
        blas::null_timer m_t;
        double      m_t_acc;
        double      m_t_left;
        double      m_t_right;
        double      m_t_q;
        double      m_t_l2_pA;
        double      m_t_l2_pB;
        double      m_t_l2_A;
        double      m_t_l2_B;

    public:
        profiler()                  : m_t_acc(0.0), m_t_left(0.0), m_t_right(0.0), m_t_q(0.0), m_t_l2_pA(0.0)
                                    , m_t_l2_pB(0.0), m_t_l2_A(0.0), m_t_l2_B(0.0){};
        void start_acc()            { m_t.tic(); };
        void end_acc()              { m_t_acc += m_t.toc(); };
        void start_left()           { m_t.tic(); };
        void end_left()             { m_t_left += m_t.toc(); };
        void start_q()              { m_t.tic(); };
        void end_q()                { m_t_q += m_t.toc(); };
        void start_right()          { m_t.tic(); };
        void end_right()            { m_t_right += m_t.toc(); };

        void start_l2_prepare_A()   { m_t.tic(); };
        void end_l2_prepare_A()     { m_t_l2_pA += m_t.toc(); };
        void start_l2_B()           { m_t.tic(); };
        void end_l2_B()             { m_t_l2_B += m_t.toc(); };
        void start_l2_B_anih()      { m_t.tic(); };
        void end_l2_B_anih()        { m_t_l2_pB += m_t.toc(); };
        void start_l2_A()           { m_t.tic(); };
        void end_l2_A()             { m_t_l2_A += m_t.toc(); };

        void print()
        {
            double l2   = m_t_l2_pA + m_t_l2_pB + m_t_l2_A + m_t_l2_B;
            std::cout   << "lev 2: " << l2 << ", left: " << m_t_left << ", right: " << m_t_right 
                        << ", q: " << m_t_q << ", acc: " << m_t_acc << "\n";
            std::cout   << "lev 2 prep A: " << m_t_l2_pA << ", prep B: " << m_t_l2_pB
                        << ", A: " << m_t_l2_A << ", B: " << m_t_l2_B << "\n";
        };
};

// minimal SLEN_A: (IHI - ILO)*NB - NB * (NB + 1) / 2
// minimal SLEN_B: (IHI - ILO)*NB - NB * (NB + 1) / 2
template<class V, class VR = typename details::real_type<V>::type>
void gghrd_lev2(i_type N, i_type ILO, i_type IHI, i_type& NB, V *A, i_type LDA, V* B, i_type LDB, 
        VR* SEQ_AC, V* SEQ_AS, i_type* SEQ_AI1, i_type* SEQ_AI2, i_type* SEQ_AJ, i_type& SLEN_A,
        VR* SEQ_BC, V* SEQ_BS, i_type* SEQ_BI1, i_type* SEQ_BI2, i_type* SEQ_BJ, i_type& SLEN_B,
        i_type &INFO, profiler& pr)
{
    INFO                = 0;

    // Quick return if possible
    if (IHI - ILO <= 1)
        return;

    i_type first_col    = ILO - 1;
    i_type last_col     = std::min(first_col + NB, IHI - 2);
    i_type col_index    = 0;
    NB                  = last_col - first_col;
    SLEN_A              = 0;
    SLEN_B              = 0;

    V* A_col            = A + (ILO - 1) * LDA;

    VR c;
    V s, r;

    i_type LD, UD;
    
    SEQ_AJ[0]           = 0;
    SEQ_BJ[0]           = 0;

    for (i_type col = first_col; col < last_col; ++col, ++col_index)
    {
        i_type K                = IHI - 1 - col;

        //prepare columns of A
        pr.start_l2_prepare_A();
        if (col_index > 0)
        {
            LD                  = N;
            UD                  = 1;
            lapack::rotseq("left", "left", "no", SLEN_A, SEQ_AC, SEQ_AS, SEQ_AI1, SEQ_AI2, N - first_col, 
                           1, LD, UD, &A_col[first_col], LDA, INFO);

            if (INFO != 0)
                return;
        };                

        // anihilate current column of A
        i_type pos_Ait          = SLEN_A;

        for (i_type row = IHI - 1; row >= col + 2; --row)
        {
            lartg<V>(A_col[row - 1], A_col[row], &c, &s,  &r);
        
            SEQ_AC[SLEN_A]      = c;
            SEQ_AS[SLEN_A]      = s;
            SEQ_AI1[SLEN_A]     = row - first_col;
            SEQ_AI2[SLEN_A]     = row + 1 - first_col;
            ++SLEN_A;

            A_col[row - 1]      = r;
            A_col[row]          = V(0.0);
        };        
        
        pr.end_l2_prepare_A();

        // apply rotations to B
        LD                      = 0;
        UD                      = N;
        i_type SLEN_Ait         = SLEN_A - pos_Ait;

        pr.start_l2_B();
        lapack::rotseq("left", "left", "no", SLEN_Ait, SEQ_AC + pos_Ait, SEQ_AS + pos_Ait, SEQ_AI1 + pos_Ait, 
                       SEQ_AI2 + pos_Ait, IHI - first_col, IHI - first_col, LD, UD, &B[first_col + first_col * LDB], 
                       LDB, INFO);
        pr.end_l2_B();

        if (INFO != 0)
            return;
        
        // anihilate subdiagonals of B        
        i_type pos_Bit          = SLEN_B;
        i_type SLEN_BH          = 0;

        pr.start_l2_B_anih();
        lapack::huundr(K, 1, col_index + 1, &B[first_col + (col + 1) * LDB], LDB, SEQ_BC + SLEN_B, SEQ_BS + SLEN_B, 
                        SEQ_BI1 + SLEN_B, SEQ_BI2 + SLEN_B, SLEN_BH, INFO);
        pr.end_l2_B_anih();

        if (INFO != 0)
            return;

        SLEN_B                  += SLEN_BH;

        //apply rotations to A
        pr.start_l2_A();

        LD                      = N - first_col;
        UD                      = K;
        lapack::rotseq("right", "right", "no", SLEN_BH, SEQ_BC + pos_Bit, SEQ_BS + pos_Bit, SEQ_BI1 + pos_Bit, 
                       SEQ_BI2 + pos_Bit, IHI - first_col, K, LD, UD, &A[first_col + (col + 1) * LDA], LDA, INFO);
        pr.end_l2_A();

        if (INFO != 0)
            return;

        //update sequence indices
        for (i_type i = pos_Bit; i < SLEN_B; ++i)
        {
            SEQ_BI1[i]      += col + 1 - first_col;
            SEQ_BI2[i]      += col + 1 - first_col;
        };

        //move to next column
        A_col                   += LDA;
        SEQ_AJ[col_index+1]     = SLEN_A;
        SEQ_BJ[col_index+1]     = SLEN_B;
    };

    return;
};

template<class V>
void clear_lt(i_type M, i_type N, V* B, i_type LDB)
{
    for (i_type JCOL = 0; JCOL < N; ++JCOL)
    {
         for (i_type JROW = JCOL + 1; JROW < M; ++JROW)
            B[JROW]     = V(0.0);

         B              = B + LDB;
    };
};

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::gghrd2(const char* COMPQ, const char *COMPZ, i_type N, i_type ILO, i_type IHI, V *A, i_type LDA, 
           V* B, i_type LDB, V* Q, i_type LDQ, V* Z2, i_type LDZ2, V* WORK, i_type LWORK, 
           i_type* IWORK, i_type LIWORK, i_type &info)
{
    using VR            = typename details::real_type<V>::type;
    bool is_compl       = details::is_complex<V>::value;

    // block size
    i_type NB           = optim_params::gghrd_block();
    i_type n_thread     = get_num_threads(domain::matcl);
    i_type NBM          = optim_params::gghrd_mult_block();

    bool accumulate_l   = true;
    bool accumulate_r   = true;

    bool ILQ            = false;
    bool ILZ            = false;
    int ICOMPQ          = 0;
    int ICOMPZ          = 0;

    // Decode C
    if( COMPQ[0] == 'N' || COMPQ[0] == 'n')
    {
        ILQ             = false;
        ICOMPQ          = 1;
    }
    else if (COMPQ[0] == 'V' || COMPQ[0] == 'v' )
    {
        ILQ             = true;
        ICOMPQ          = 2;
    }
    else if (COMPQ[0] == 'I' || COMPQ[0] == 'i')
    {
        ILQ             = true;
        ICOMPQ          = 3;
    }
    else
    {
        ICOMPQ          = 0;
    };

    // Decode COMPZ
    if (COMPZ[0] == 'N' || COMPZ[0] == 'n')
    {
        ILZ             = false;
        ICOMPZ          = 1;
    }
    else if (COMPZ[0] == 'V' || COMPZ[0] == 'v' )
    {
        ILZ             = true;
        ICOMPZ          = 2;
    }
    else if (COMPZ[0] == 'I' || COMPZ[0] == 'i')
    {
        ILZ             = true;
        ICOMPZ          = 3;
    }
    else
    {
        ICOMPZ          = 0;
    };

    // Test the input parameters.
    int INFO = 0;

    if (ICOMPQ <= 0 )
        INFO            = -1;
    else if (ICOMPZ <= 0)
        INFO            = -2;
    else if (N < 0)
        INFO            = -3;
    else if (ILO < 1)
        INFO            = -4;
    else if (IHI > N || IHI < ILO-1 )
        INFO            = -5;
    else if (LDA < lapack::maximum( 1, N ) )
        INFO            = -7;
    else if (LDB < lapack::maximum( 1, N ) )
        INFO            = -9;
    else if ( ( ILQ && LDQ < N ) || LDQ < 1 )
        INFO            = -11;
    else if ( ( ILZ && LDZ2 < N ) || LDZ2 < 1 )
        INFO            = -13;

    if ( INFO != 0 )
    {
        info            = INFO;
        return;
    };    

    if (IHI - ILO - 1 < NB)
    {
        if (LWORK == -1 || LIWORK == -1)
        {
            WORK[0]     = V(1.0);
            IWORK[0]    = 1;
            return;
        }
        else
        {
            return gghrd2_s(COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, LDQ, Z2, LDZ2, info);
        };
    };

    // initialize givens accumulator
    givens_accumulator<V> gac(N, NB, NBM, n_thread);

    i_type max_it       = (IHI - ILO - 2) / NB + 1;    
    i_type SLEN_A       = (IHI - ILO) * NB - NB * (NB + 1) / 2;
    i_type SLEN_B       = SLEN_A;
    i_type SLEN_AR      = is_compl ? SLEN_A/2 + 1 : SLEN_A;
    i_type SLEN_BR      = SLEN_AR;
    i_type LWORK_AC     = (accumulate_l || accumulate_r)? gac.work_size() : 0;

    // calculation of minimum work size    
    i_type LWORK_min;
    i_type LIWORK_min;

    if (IHI - ILO - 1 <= 0)
    {
        LWORK_min       = 1;
        LIWORK_min      = 1;
    }
    else
    {
        LWORK_min       = SLEN_A + SLEN_AR + SLEN_B + SLEN_BR + LWORK_AC;
        LIWORK_min      = 2*SLEN_A + 2*SLEN_B + 2 * NB;
    };

    WORK[0]             = V(VR(LWORK_min));
    IWORK[0]            = LIWORK_min;

    if (LWORK == -1 || LIWORK == -1)
        return;

    if (LWORK < LWORK_min)
    {
        INFO            = -15;
        return;
    };
    if (LIWORK < LIWORK_min)
    {
        INFO            = -17;
        return;
    };

    i_type LD_Q = N, LD_Z = N, UD_Z = N, UD_Q = N;

    // Initialize Q and Z if desired.
    if (ICOMPQ == 3)
    {
        LD_Q            = 0;
        UD_Q            = 0;
        lapack::laset("Full", N, N, V(0.0), V(1.0), Q, LDQ );
    }

    if (ICOMPZ == 3)
    {
        LD_Z            = 0;
        UD_Z            = 0;
        lapack::laset("Full", N, N, V(0.0), V(1.0), Z2, LDZ2 );
    };

    // Quick return if possible
    if (N <= 1)
        return;

    // Zero out lower triangle of B
    clear_lt(N, N, B, LDB);

    // Quick return if possible
    if (IHI - ILO + 1 <= 2)
        return;

    //prepare work arrays
    V* WORK_s           = WORK;
    i_type* IWORK_s     = IWORK;

    V* SEQ_AS           = WORK;
    WORK                = WORK + SLEN_A;

    VR* SEQ_AC          = reinterpret_cast<VR*>(WORK);
    WORK                = WORK + SLEN_AR;    

    V* SEQ_BS           = WORK;
    WORK                = WORK + SLEN_B;

    VR* SEQ_BC          = reinterpret_cast<VR*>(WORK);
    WORK                = WORK + SLEN_BR;    

    V* WORK_AC          = WORK;
    WORK                = WORK + LWORK_AC;

    i_type* SEQ_AI1     = IWORK;
    IWORK               = IWORK + SLEN_A;

    i_type* SEQ_AI2     = IWORK;
    IWORK               = IWORK + SLEN_A;

    i_type* SEQ_AJ      = IWORK;
    IWORK               = IWORK + NB + 1;

    i_type* SEQ_BI1     = IWORK;
    IWORK               = IWORK + SLEN_B;

    i_type* SEQ_BI2     = IWORK;
    IWORK               = IWORK + SLEN_B;

    i_type* SEQ_BJ      = IWORK;
    IWORK               = IWORK + NB + 1;

    i_type LD_N         = N;
    i_type UD_N         = N;
    i_type LD, UD;

    gac.set_work(WORK_AC);
    
    profiler pr;

    LD_Z                += ILO - 1;
    UD_Z                -= ILO - 1;
    LD_Q                += ILO - 1;
    UD_Q                -= ILO - 1;

    UD_Z                = std::max(UD_Z, 0);
    UD_Q                = std::max(UD_Q, 0);

    for (i_type i = 0; i < max_it; ++i)
    {
        i_type NBJ      = NB;

        //process NB columns        
        gghrd_lev2<V>(N, ILO, IHI, NBJ, A, LDA, B, LDB, SEQ_AC, SEQ_AS, SEQ_AI1, SEQ_AI2, SEQ_AJ, SLEN_A,
                   SEQ_BC, SEQ_BS, SEQ_BI1, SEQ_BI2, SEQ_BJ, SLEN_B, info, pr);

        pr.start_acc();
        // accumulate A rotations
        if (accumulate_l)
            gac.accumulate(true, IHI - ILO + 1, NBJ, SEQ_AC, SEQ_AS, SEQ_AI1, SEQ_AI2, SEQ_AJ);
        pr.end_acc();

        pr.start_left();

        // apply A rotations to remaining part of A from the left
        if (accumulate_l)
        {
            gac.apply_left_left(N - ILO + 1 - NBJ, &A[ILO-1 + (ILO-1+NBJ)*LDA], LDA);
        }
        else
        {
            LD              = N;
            UD              = N;
            lapack::rotseq("left", "left", "no", SLEN_A, SEQ_AC, SEQ_AS, SEQ_AI1, SEQ_AI2, 
                           IHI - ILO + 1, N - ILO + 1 - NBJ, LD, UD, &A[ILO-1 + (ILO-1+NBJ)*LDA], LDA, INFO);
        };
        
        // apply A rotations to remaining part of B from the left
        if (IHI < N)
        {
            if (accumulate_l)
            {
                gac.apply_left_left(N - IHI, &B[ILO-1 + IHI*LDB], LDB);
            }
            else
            {
                LD          = N;
                UD          = N;
                lapack::rotseq("left", "left", "no", SLEN_A, SEQ_AC, SEQ_AS, SEQ_AI1, SEQ_AI2, 
                        IHI - ILO + 1, N - IHI, LD, UD, &B[ILO-1 + IHI*LDB], LDB, INFO);
            };
        };

        pr.end_left();

        pr.start_q();
        // apply A rotations to Q from the right
        if (ILQ)
        {
            if (accumulate_l)
            {
                gac.apply_left_right(N, UD_Q, &Q[(ILO-1) * LDQ], LDQ);
                UD_Q        = std::max(UD_Q - NBJ, 0);
            }
            else
            {
                lapack::rotseq("left", "right", "conj", SLEN_A, SEQ_AC, SEQ_AS, SEQ_AI1, SEQ_AI2, 
                    N, IHI - ILO + 1, LD_Q, UD_Q, &Q[(ILO-1) * LDQ], LDQ, INFO);

                //update bandwidth estimation for next block
                LD_Q        = std::min(LD_Q + NBJ, N);
                UD_Q        = std::max(UD_Q - NBJ, 0);
            };
        };

        pr.end_q();

        pr.start_acc();
        // accumulate B rotations
        if (accumulate_r)
            gac.accumulate(false, IHI - ILO + 1, NBJ, SEQ_BC, SEQ_BS, SEQ_BI1, SEQ_BI2, SEQ_BJ);    
        pr.end_acc();

        pr.start_right();

        // apply B rotations to Z 
        if (ILZ)
        {
            if (accumulate_r)
            {
                gac.apply_right_right(N, UD_Z, &Z2[(ILO-1)*LDZ2], LDZ2);
                UD_Z        = std::max(UD_Z - NBJ, 0);
            }
            else
            {
                lapack::rotseq("right", "right", "no", SLEN_B, SEQ_BC, SEQ_BS, SEQ_BI1, SEQ_BI2, 
                       N, IHI - ILO + 1, LD_Z, UD_Z, &Z2[(ILO-1)*LDZ2], LDZ2, INFO);

                //update bandwidth estimation for next block
                LD_Z        = std::min(LD_Z + NBJ, N);
                UD_Z        = std::max(UD_Z - NBJ, 0);
            };
        };
        
        if (ILO > 1)
        {
            // apply B rotations from the right to upper part of A
            if (accumulate_r)
            {
                gac.apply_right_right(ILO - 1, UD_N, &A[(ILO-1)*LDA], LDA);
            }
            else
            {
                LD          = N;
                UD          = N;
                lapack::rotseq("right", "right", "no", SLEN_B, SEQ_BC, SEQ_BS, SEQ_BI1, SEQ_BI2, 
                        ILO - 1, IHI - ILO + 1, LD, UD, &A[(ILO-1)*LDA], LDA, INFO);
            };

            // apply B rotations from the right to upper part of B
            if (accumulate_r)
            {
                gac.apply_right_right(ILO - 1, UD_N, &B[(ILO-1)*LDB], LDB);
            }
            else
            {
                LD          = N;
                UD          = N;
                lapack::rotseq("right", "right", "no", SLEN_B, SEQ_BC, SEQ_BS, SEQ_BI1, SEQ_BI2, 
                               ILO - 1, IHI - ILO + 1, LD, UD, &B[(ILO-1)*LDB], LDB, INFO);
            };
        };

        pr.end_right();
        ILO             += NBJ;
    };

    WORK_s[0]       = V(VR(LWORK_min));
    IWORK_s[0]      = LIWORK_min;

    //pr.print();    
};

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::gghrdf(const char *COMPQ, const char *COMPZ, i_type N, V* A, i_type LDA, 
    V* B, i_type LDB, V* Q, i_type LDQ, V* Z, i_type LDZ, V* WORK, i_type LWORK, i_type* IWORK, 
    i_type LIWORK, i_type& INFO)
{
    using VR = typename details::real_type<V>::type;
    using VC = typename details::complex_type<V>::type;

    bool is_V_real  = details::is_complex<V>::value == false;

    bool ILQ        = false;
    bool ILZ        = false;
    i_type ICOMPQT  = 0;
    i_type ICOMPZT  = 0;

    if( COMPQ[0] == 'N' || COMPQ[0] == 'n')
    {
        ILQ             = false;
    }
    else if (COMPQ[0] == 'V' || COMPQ[0] == 'v' )
    {
        ILQ             = true;
    }
    else
    {
        ICOMPQT         = 1;
    };

    // Decode COMPZ
    if (COMPZ[0] == 'N' || COMPZ[0] == 'n')
    {
        ILZ             = false;
    }
    else if (COMPZ[0] == 'V' || COMPZ[0] == 'v' )
    {
        ILZ             = true;
    }
    else
    {
        ICOMPZT         = 1;
    };

    i_type  MINWRK, MAXWRK, IERR;

    const char* ctrans_name = is_V_real? "T" : "C";
    const char* qtriu_name  = is_V_real? "H" : "U";

    // Test the input arguments
    INFO        = 0;
    bool LQUERY = ( LWORK == -1 || LIWORK == -1);

    // Test the input parameters.
    if (ICOMPQT < 0 )
        INFO            = -1;
    else if (ICOMPZT < 0)
        INFO            = -2;
    else if (N < 0)
        INFO            = -3;
    else if (LDA < lapack::maximum( 1, N ) )
        INFO            = -5;
    else if (LDB < lapack::maximum( 1, N ) )
        INFO            = -7;
    else if ( ( ILQ && LDQ < N ) || LDQ < 1 )
        INFO            = -9;
    else if ( ( ILZ && LDZ < N ) || LDZ < 1 )
        INFO            = -11;

    if ( INFO != 0 )
    {
        return;
    };    

    // Compute workspace
    // (Note: Comments in the code beginning "Workspace:" describe the
    //  minimal amount of workspace needed at that point in the code,
    //  as well as the preferred amount for good performance.
    //  NB refers to the optimal block size for the immediately
    //  following subroutine, as returned by ILAENV.)

    i_type lrwork;
    i_type lwork_hess;
    i_type liwork_hess;

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

        if ( ILQ )
        {
            MINWRK      = std::max(MINWRK, MIN_ORG );
            MAXWRK      = std::max( MAXWRK, WORK_S + MAX_ORG );            
        };

        MINWRK          = std::max(MINWRK, MIN_QZ );
        MAXWRK          = std::max( MAXWRK, MAX_QZ );

        V work_query;
        i_type iwork_query;

        gghrd2(COMPQ, COMPZ, N, 1, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, 
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
        INFO    = -13;
    if (LIWORK < liwork_hess && !LQUERY )
        INFO    = -15;

    if ( INFO != 0 )
        return;
    else if (LQUERY )
        return;

    // Quick return if possible
    if ( N <= 1 )
        return;

    V* WORK_SAV = WORK;
    VR* RWORK   = reinterpret_cast<VR*>(WORK);
    WORK        = WORK + lrwork;
        
    // Permute the matrix to make it more nearly triangular
    // (Real Workspace: need 6*N + 2*N space for storing balancing factors)
    i_type ILEFT    = 1;
    i_type IRIGHT   = N + 1;
    i_type IRWRK    = IRIGHT + N;
    i_type ILO, IHI;

    const char* bal_char    = "P";
    lapack::ggbal2(bal_char, N, A, LDA, B, LDB, ILO, IHI, RWORK + ILEFT - 1, RWORK + IRIGHT-1, IERR );

    i_type DOHESS       = check_hess_type<V>(N, A, LDA, B, LDB, ILO, IHI);

    const char* JOBQ    = ILQ ? ((DOHESS == 2) ? "V" : "I") : "N";
    const char* JOBZ    = ILZ ? "I" : "N";

    if ( ILZ && DOHESS == 0)
        lapack::laset( "Full", N, N, V(0.0), V(1.0), Z, LDZ );

    if (ILQ && (DOHESS == 0 || DOHESS == 2))
        lapack::laset("Full", N, N, V(0.0), V(1.0), Q, LDQ );

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

            // Initialize Q
            //   (Workspace: need N, prefer N*NB)

            if ( ILQ )
            {
                if (IROWS > 1 )
                {
                    lapack::lacpy("L", IROWS-1, IROWS-1, B + (ILO+1-1) + (ILO-1)*LDB, LDB,
                            Q + (ILO+1-1) + (ILO-1)*LDQ, LDQ );
                };

                lapack::orgqr<V>(IROWS, IROWS, IROWS, Q + (ILO-1) + (ILO-1)*LDQ, LDQ, WORK + ITAU-1, 
                        WORK + IWRK - 1, LWORK+1-IWRK, &IERR );
            };
        };

        // Reduce to generalized Hessenberg form
        //   (Workspace: need lwork_hess)
        //   (integer Workspace: need liwork_hess)

        lapack::gghrd2(JOBQ, JOBZ, N, ILO, IHI, A, LDA, B, LDB, Q, LDQ, Z, LDZ, 
                       WORK + IWRK - 1, lwork_hess, IWORK, liwork_hess, IERR );
    };

    // Apply back-permutation to VSL and VSR
    //   (Workspace: none needed)

    if (ILQ)
        lapack::ggbak( bal_char, "L", N, ILO, IHI, RWORK + ILEFT - 1, RWORK + IRIGHT - 1, N, Q, LDQ, IERR );

    if (ILZ)
        lapack::ggbak( bal_char, "R", N, ILO, IHI, RWORK + ILEFT - 1, RWORK + IRIGHT - 1, N, Z, LDZ, IERR );

    WORK_SAV[0] = VR(MAXWRK);
    IWORK[0]    = liwork_hess;
};

template void BLAS_EXT_EXPORT
gghrd2<s_type>(const char* COMPQ, const char *COMPZ, i_type N, i_type ILO, i_type IHI, s_type *A, i_type LDA, 
           s_type* B, i_type LDB, s_type* Q, i_type LDQ, s_type* Z, i_type LDZ, s_type* WORK, i_type LWORK, 
           i_type* IWORK, i_type LIWORK, i_type& info);

template void BLAS_EXT_EXPORT
gghrd2<d_type>(const char* COMPQ, const char *COMPZ, i_type N, i_type ILO, i_type IHI, d_type *A, i_type LDA, 
           d_type* B, i_type LDB, d_type* Q, i_type LDQ, d_type* Z, i_type LDZ, d_type* WORK, i_type LWORK, 
           i_type* IWORK, i_type LIWORK, i_type& info);

template void BLAS_EXT_EXPORT
gghrd2<c_type>(const char* COMPQ, const char *COMPZ, i_type N, i_type ILO, i_type IHI, c_type *A, i_type LDA, 
           c_type* B, i_type LDB, c_type* Q, i_type LDQ, c_type* Z, i_type LDZ, c_type* WORK, i_type LWORK, 
           i_type* IWORK, i_type LIWORK, i_type& info);

template void BLAS_EXT_EXPORT
gghrd2<z_type>(const char* COMPQ, const char *COMPZ, i_type N, i_type ILO, i_type IHI, z_type *A, i_type LDA, 
           z_type* B, i_type LDB, z_type* Q, i_type LDQ, z_type* Z, i_type LDZ, z_type* WORK, i_type LWORK, 
           i_type* IWORK, i_type LIWORK, i_type& info);

template void BLAS_EXT_EXPORT
gghrdf<s_type>(const char* COMPQ, const char *COMPZ, i_type N, s_type *A, i_type LDA, 
           s_type* B, i_type LDB, s_type* Q, i_type LDQ, s_type* Z, i_type LDZ, s_type* WORK, i_type LWORK, 
           i_type* IWORK, i_type LIWORK, i_type& info);
template void BLAS_EXT_EXPORT
gghrdf<d_type>(const char* COMPQ, const char *COMPZ, i_type N, d_type *A, i_type LDA, 
           d_type* B, i_type LDB, d_type* Q, i_type LDQ, d_type* Z, i_type LDZ, d_type* WORK, i_type LWORK, 
           i_type* IWORK, i_type LIWORK, i_type& info);

template void BLAS_EXT_EXPORT
gghrdf<c_type>(const char* COMPQ, const char *COMPZ, i_type N, c_type *A, i_type LDA, 
           c_type* B, i_type LDB, c_type* Q, i_type LDQ, c_type* Z, i_type LDZ, c_type* WORK, i_type LWORK, 
           i_type* IWORK, i_type LIWORK, i_type& info);

template void BLAS_EXT_EXPORT
gghrdf<z_type>(const char* COMPQ, const char *COMPZ, i_type N, z_type *A, i_type LDA, 
           z_type* B, i_type LDB, z_type* Q, i_type LDQ, z_type* Z, i_type LDZ, z_type* WORK, i_type LWORK, 
           i_type* IWORK, i_type LIWORK, i_type& info);
};};
