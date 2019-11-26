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
#include "tgsen_helpers.h"
#include "blas/matcl-blas-ext/lapack_ext/utils/optim_params.h"
#include "matcl-blas-lapack/blas/details/config_blas_lib.h"

#include <vector>
#include <cassert>

namespace matcl { namespace lapack
{

//-----------------------------------------------------------------
//                  MAIN ROUTINES
//-----------------------------------------------------------------

#define assert2 assert

#define MAGIC_NUMBER 12345.

template<class V, bool Is_compl = details::is_complex<V>::value>
struct compute_eigenvalues_tgsen3
{
    static void eval(i_type n, i_type k, V* a, i_type lda, V* b, i_type ldb, 
                    i_type wantq, i_type wantz, V* Q, i_type LDQ, V* Z, i_type LDZ, 
                    typename details::complex_type<V>::type* alpha, V* beta)
    {
        using VR    = typename details::real_type<V>::type;

        // Compute generalized eigenvalues of reordered pair (A, B) and
        // normalize the generalized Schur form.

        // Get machine constants
        VR SMLNUM       = lamch<VR>( "S" );

        V work[8];

        bool PAIR = false;
        for (i_type K = 1; K <= n; ++K)
        {
            if ( PAIR )
                PAIR = false;
            else
            {
                if ( K < n )
                {
                    if ( a[K + (K-1)*lda] != 0. )
                    {
                        PAIR = true;
                    };
                };

                if ( PAIR )
                {
                    // Compute the eigenvalue(s) at position K.

                    work[0] = a[K-1 + (K-1)*lda];
                    work[1] = a[K + (K-1)*lda];
                    work[2] = a[K-1 + K*lda];
                    work[3] = a[K + K*lda];
                    work[4] = b[K-1 + (K-1)*ldb];
                    work[5] = b[K + (K-1)*ldb];
                    work[6] = b[K-1 + K*ldb];
                    work[7] = b[K + K*ldb];

                    i_type ld_tmp = 2;

                    V a1, a2, a3;
                    lag2<V>( work, ld_tmp, work + 4, ld_tmp, SMLNUM, beta[K-1],
                                         beta[K], a1, a2, a3 );
                    alpha[K-1].real(a1);
                    alpha[K].real(a2);
                    alpha[K-1].imag(a3);
                    alpha[K].imag(-a3);
                }
                else
                {
                    if ( b[K-1 + (K-1)*ldb] < 0. )
                    {
                        // If B(K,K) is negative, make it positive
                        for (i_type I = 1; I <= n; ++I)
                        {
                            a[K-1 + (I-1)*lda] = -a[K-1 + (I-1)*lda];
                            b[K-1 + (I-1)*ldb] = -b[K-1 + (I-1)*ldb];
                        };
                        if (wantq)
                        {
                            for (i_type I = 1; I <= k; ++I)
                            {
                                Q[I-1 + (K-1)*LDQ] = -Q[I-1 + (K-1)*LDQ];
                            };
                        };
                    };

                    alpha[K-1].real(a[K-1 + (K-1)*lda]);
                    alpha[K-1].imag(0.);

                    beta[K-1]   = b[K-1 + (K-1)*ldb];
                };
            };
        };
    };
};

template<class V>
struct compute_eigenvalues_tgsen3<V,true>
{
    static void eval(i_type n, i_type k, V* a, i_type lda, V* b, i_type ldb, 
                    i_type wantq, i_type wantz, V* Q, i_type LDQ, V* Z, i_type LDZ, 
                    typename details::complex_type<V>::type* alpha, V* beta)
    {
        using VR    = typename details::real_type<V>::type;

        VR SMLNUM       = lamch<VR>( "S" );

        // Compute generalized eigenvalues of reordered pair (A, B) and
        // normalize the generalized Schur form.

        for (i_type K = 1; K <= n; ++K)
        {
            //     If B(K,K) is complex, make it real and positive (normalization
            //     of the generalized Schur form) and Store the generalized
            //     eigenvalues of reordered pair (A, B)
            V bkk       = b[K-1 + (K-1)*ldb];
            VR DSCALE   = abs( bkk );

            if ( DSCALE > SMLNUM )
            {
                V TEMP1     = conj( bkk / DSCALE );
                V TEMP2     = bkk / DSCALE;

                b[K-1 + (K-1)*ldb]  = DSCALE;

                lapack::scal(n-K, TEMP1, b + K-1 + (K+1-1)*ldb, ldb);
                lapack::scal(n-K+1, TEMP1, a + K-1 + (K-1)*lda, lda);

                if (wantq)
                {
                    for (i_type I = 1; I <= k; ++I)
                    {
                        Q[I-1 + (K-1)*LDQ] = Q[I-1 + (K-1)*LDQ] * TEMP2;
                    };
                };
            }
            else
            {
                b[K-1 + (K-1)*ldb]  = V(0.0);
            };

            alpha[K-1]  = a[K-1 + (K-1)*lda];
            beta[K-1]   = b[K-1 + (K-1)*ldb];
        };
    };
};

template<class V>
static bool compute_block(i_type& BS, i_type BB, i_type n, i_type K, i_type max_eig_in_window, i_type max_window_size,
              i_type* sel_vec, V* a, i_type lda, V* b, i_type ldb, V* Q, i_type LDQ, 
              V* Z, i_type LDZ, i_type wantq, i_type wantz, matrix_lock* a_lock, matrix_lock* b_lock,
              V* work, i_type& lwork)
{    
    using VC = typename details::complex_type<V>::type;

    i_type LWMIN_q      = (max_window_size+1) * (max_window_size + 1);
    i_type LWMIN_z      = LWMIN_q;
    i_type LWMIN_A      = LWMIN_q;
    i_type LWMIN_B      = LWMIN_q;
    i_type NBQ          = lapack::minimum(K,2 * (max_window_size+1));
    i_type LWMIN_m      = (max_window_size+1) * NBQ;
    i_type LWMIN_AB     = max_window_size+1;
    i_type LWMIN_l      = 4*(max_window_size+1)+16;
    i_type COMPL_SIZE   = sizeof(VC)/sizeof(V);
    i_type LWMIN        = LWMIN_q + LWMIN_z + lapack::maximum(LWMIN_m, 2*LWMIN_A + COMPL_SIZE*LWMIN_AB 
                                                              + LWMIN_AB + LWMIN_l);

    if (lwork == -1)
    {
        lwork           = LWMIN;
        lwork           = lwork + 1;
        return true;
    };

    if (BS > 1)
    {
        //complex eigenvalues cannot be split
        assert2(is_equal(a[BS - 1 + (BS - 2)*lda], V(0)) == true);
    };
    i_type BS_save      = BS;

    V* work_test        = work + LWMIN;

    V* work_q           = work; work = work + LWMIN_q;
    V* work_z           = work; work = work + LWMIN_z;
    V* work_mat         = work;

    V* work_A           = work; work = work + LWMIN_A;
    V* work_B           = work; work = work + LWMIN_B;
    V* alpha_re         = work; work = work + LWMIN_AB * COMPL_SIZE;
    V* beta             = work; work = work + LWMIN_AB;
    V* work_LAP         = work; work = work + LWMIN_l;   
   
    work_test[0]        = MAGIC_NUMBER;

    VC* alpha           = reinterpret_cast<VC*>(alpha_re);

    assert2(sel_vec[BS-1] == 0);

    bool success        = true;

    while (BS <= BB)
    {
        // find working region
        i_type n_eigs   = 0;
        i_type BE       = find_block_end(BS,BB, max_eig_in_window, sel_vec, n_eigs);

        if (n_eigs == 0)
        {
            // nothing to do
            break;
        };

        //reorder eigenvalues in windows
        i_type WE       = BE;
        i_type WS       = BE;

        while (WS != BS)
        {
            WS          = find_window_start(BS,WE,max_window_size,sel_vec,a,lda);

            // reorder current window
            i_type wsize            = WE - WS + 1;

            i_type sel_1 = sel_vec[WS - 1];
            i_type sel_2 = sel_vec[WE - 1];

            assert2(wsize <= max_window_size+1);
            assert2(sel_1 == 0 && sel_2 != 0);

            // prepare q, z matrices
            make_id(work_q, wsize);
            make_id(work_z, wsize);

            // use unblocked version
            i_type n_reorder        = 0;
            i_type tmp_info         = 0;
            i_type tmp_liwork       = 0;

            // operate on local workpsace to avoid cache misses
            lacpy("N", wsize, wsize, a + WS - 1 + (WS - 1) * lda, lda, work_A, wsize);
            lacpy("N", wsize, wsize, b + WS - 1 + (WS - 1) * ldb, ldb, work_B, wsize);

            tgsen<V>(0,1,1,sel_vec + WS - 1, wsize, work_A, wsize, work_B, wsize, 
                alpha, beta, work_q, wsize, work_z, wsize, &n_reorder, 
                nullptr, nullptr, nullptr, work_LAP, LWMIN_l, &tmp_liwork, 1, &tmp_info);

            lacpy("N", wsize, wsize, work_A, wsize, a + WS - 1 + (WS - 1) * lda, lda);
            lacpy("N", wsize, wsize, work_B, wsize, b + WS - 1 + (WS - 1) * ldb, ldb);

            if ( tmp_info > 0 )
            {
                // Swap is rejected: exit.
                success             = false;
                goto lab_exit;
            };

            bool calc               = false;
            bool calc_A1            = calc;
            bool calc_A2            = calc;
            bool calc_B1            = calc;
            bool calc_B2            = calc;
            bool calc_Q             = (wantq == false);
            bool calc_Z             = (wantz == false);

            while (!calc_A1 || !calc_A2 || !calc_B1 || !calc_B2 || !calc_Q || !calc_Z)
            {
                // apply delayed rotations
                // matrix TA
                if (n - WS + 1 - wsize <= 0)
                {
                    calc_A1 = true;
                }
                else if (calc_A1 == false)
                {
                    access al = a_lock->get_row();

                    if (al.is_granted() == true)
                    {
                        lacpy("N", wsize, n - WS + 1 - wsize, a + WS - 1 + (WS - 1 + wsize) * lda, 
                            lda, work_mat,wsize);
                        gemm<V>("T", "N", wsize, n - WS + 1 - wsize, wsize, 1.0, work_q, wsize, 
                                work_mat, wsize, 0.0, a + WS - 1 + (WS - 1 + wsize) * lda, lda);

                        calc_A1 = true;
                    };
                };
                if (WS - 1 <= 0)
                {
                    calc_A2 = true;
                }
                else if (calc_A2 == false)
                {
                    access al = a_lock->get_col();

                    if (al.is_granted() == true)
                    {
                        lacpy("N", WS - 1, wsize, a + (WS - 1) * lda, lda, work_mat, WS - 1);
                        gemm<V>("N", "N", WS - 1 , wsize, wsize, 1.0, work_mat, WS - 1, 
                                work_z, wsize, 0.0, a + (WS - 1) * lda, lda);

                        calc_A2 = true;
                    };
                };

                // matrix TB
                if (n - WS + 1 - wsize <= 0)
                {
                    calc_B1 = true;
                }
                else if (calc_B1 == false)
                {
                    access bl = b_lock->get_row();

                    if (bl.is_granted() == true)
                    {
                        lacpy("N", wsize, n - WS + 1 - wsize, b + WS - 1 + (WS - 1 + wsize) * ldb, 
                            ldb, work_mat,wsize);
                        gemm<V>("T", "N", wsize, n - WS + 1 - wsize, wsize, 1.0, work_q, wsize, 
                                work_mat, wsize, 0.0, b + WS - 1 + (WS - 1 + wsize) * ldb, ldb);

                        calc_B1 = true;
                    };
                };
                if (WS - 1 <= 0)
                {
                    calc_B2 = true;
                }
                else if (calc_B2 == false)
                {
                    access bl = b_lock->get_col();

                    if (bl.is_granted() == true)
                    {
                        lacpy("N", WS - 1, wsize, b + (WS - 1) * ldb, ldb, work_mat, WS - 1);
                        gemm<V>("N", "N", WS - 1 , wsize, wsize, 1.0, work_mat, WS - 1, 
                                work_z, wsize, 0.0, b + (WS - 1) * ldb, ldb);

                        calc_B2 = true;
                    };
                };

                if (calc_Q == false)
                {
                    // matrix Q
                    for (i_type j = 0; j < K; j += NBQ)
                    {
                        i_type JNBQ = lapack::minimum(NBQ, K - j);
                        lacpy("N", JNBQ, wsize, Q + j + (WS - 1)*LDQ, LDQ, work_mat, JNBQ);
                        gemm<V>("N", "N", JNBQ, wsize, wsize, 1.0, work_mat, JNBQ, work_q, wsize, 
                            0.0, Q + j + (WS - 1)*LDQ, LDQ);
                    };

                    calc_Q = true;
                    continue;
                };

                if (calc_Z == false)
                {
                    // matrix Z
                    for (i_type j = 0; j < K; j += NBQ)
                    {
                        i_type JNBQ = lapack::minimum(NBQ, K - j);
                        lacpy("N", JNBQ, wsize, Z + j + (WS - 1)*LDZ, LDZ, work_mat, JNBQ);
                        gemm<V>("N", "N", JNBQ, wsize, wsize, 1.0, work_mat, JNBQ, work_q, wsize, 
                            0.0, Z + j + (WS - 1)*LDZ, LDZ);
                    };

                    calc_Z = true;
                    continue;
                };
            };

            // move to next window
            WE              = WS + n_reorder - 1;

            analyse_select2<V>::eval(n_reorder, sel_vec + WS - 1, a + WS - 1 + (WS - 1)*lda, lda);

            assert2(WS >= BS_save);
            assert2(WS + n_reorder - 1 <= BB);
        };

        // finish processing current block of eigenvalues
        for (i_type K = BS + n_eigs; K <= BE; ++K)
        {
            sel_vec[K-1]    = 0;
        };
        BS                  = BS + n_eigs;

        if (BS > 1)
        {
            //complex eigenvalues cannot be split
            assert2(is_equal(a[BS - 1 + (BS - 2)*lda], V(0)) == true);
        };

        assert2(BS + n_eigs >= BS_save);
        assert2(BE <= BB);
        assert2(sel_vec[BS-1] == 0);
    };

  lab_exit:

    assert2(is_equal(work_test[0], V(MAGIC_NUMBER)) == true);

    return success;
};

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
matcl::lapack::tgsen3(i_type wantq,i_type wantz,const i_type *select,i_type n, i_type K,
        V *a,i_type lda,V *b,i_type ldb, typename details::complex_type<V>::type *alpha, V *beta, 
        V *Q, i_type LDQ, V* Z, i_type LDZ,i_type *m, V *work,i_type lwork, i_type *info)
{
    using VR    = typename details::real_type<V>::type;

    //parameters of the algorithm
    const i_type max_eig_in_window  = optim_params::dtgsen3_max_eig_in_window;
    const i_type max_window_size    = optim_params::dtgsen3_max_window_size;
    const i_type max_threads        = get_num_threads(domain::matcl);
    //const i_type max_threads        = 1;

    i_type INFO = 0;

    if ( n < 0 )
        INFO    = -4;
    else if ( K < 0 )
        INFO    = -5;
    else if (lda < lapack::maximum( 1, n ) )
        INFO    = -7;
    else if (ldb < lapack::maximum( 1, n ) )
        INFO    = -9;
    else if ( LDQ < 1 || ( wantq && LDQ < K ) )
        INFO    = -14;
    else if (LDZ < 1 || ( wantz && LDZ < K ) )
        INFO    = -16;
    
    if (INFO != 0 )
    {
        *info = INFO;
        return;
    };

    if (n <= max_window_size && K == n)
    {
        i_type tmp_iwork = 0;
        tgsen<V>(0,wantq,wantz,select, n, a, lda, b, ldb, alpha, beta, Q, LDQ, Z, LDZ, m, 
            nullptr,nullptr,nullptr, work, lwork, &tmp_iwork, 1, info);

        return;
    };
    
    bool LQUERY = ( lwork == -1 );   

    // query workspace
    i_type LWMIN    = -1;    
    i_type BS_dum   = 1;

    compute_block<V>(BS_dum, 0, n, K, max_eig_in_window, max_window_size,nullptr, a, lda, b, ldb, Q, LDQ, 
                  Z, LDZ, wantq, wantz, nullptr, nullptr, nullptr, LWMIN);

    // Set M to the dimension of the specified pair of deflating subspaces.
    //TODO: change this to workspace request
    std::vector<i_type> sel_vec(n);        

    i_type M            = analyse_select<V>::eval(n, select, a, lda, sel_vec.data());

    // Collect the selected blocks at the top-left corner of (A, B).
    i_type start_pos    = ignore_sorted_eigenvalues(sel_vec.data(), 1, n);
    i_type tot_eigs     = 0;

    find_block_end(start_pos, n, n, sel_vec.data(), tot_eigs);

    i_type eigs_thread  = lapack::maximum(tot_eigs / max_threads / 2, max_eig_in_window * 2);
    i_type n_threads    = lapack::maximum(tot_eigs / eigs_thread, 1);
    n_threads           = lapack::minimum(n_threads, max_threads);

    i_type LIWMIN   = 1;
    work[0]         = V(VR(LWMIN * n_threads));

    if( lwork < LWMIN && LQUERY == false )
        INFO        = -18;

    if ( INFO != 0 )
    {
        *info       = INFO;
        return;
    }
    else if ( LQUERY )
        return;

    // Quick return if possible.
    if ( M == n || M == 0 )
    {
        compute_eigenvalues_tgsen3<V>::eval(n, K, a, lda, b, ldb, wantq, wantz,
                                            Q, LDQ, Z, LDZ, alpha, beta);
        return;
    };

    i_type IERR         = 0;

    matrix_lock         a_lock;
    matrix_lock         b_lock;
    workspace_manager<V>   m_workspace(work, n_threads,LWMIN);    

    auto comp_f         = [n,K,max_eig_in_window,max_window_size, &sel_vec,a,lda,b,ldb,Q,LDQ,Z,LDZ,
                            wantq, wantz, &a_lock, &b_lock, &lwork
                          ]
                        (i_type& BS, i_type BB, V* work) -> bool
                        { 
                            return compute_block<V>(BS, BB, n, K,max_eig_in_window, max_window_size, 
                                    sel_vec.data(), a, lda, b, ldb, Q, LDQ, Z, LDZ, wantq, wantz, 
                                    &a_lock, &b_lock, work, lwork);
                        };


    bool res = false;

    if (n_threads == 1)
    {
        //use non-threaded version
        i_type BS       = start_pos;
        res             = comp_f(BS, n, work);
        *m              = BS - 1;
    }
    else
    {
        task_manager<V>    m_tasks(&m_workspace);    

        // initialize tasks
        i_type BS_0     = start_pos;
        while (BS_0 < n)
        {
            i_type dum  = 0;        
            i_type BS_1 = BS_0;
            i_type BB_1 = find_block_end(BS_1, n, eigs_thread, sel_vec.data(), dum);
            BB_1        = ignore_sorted_eigenvalues(sel_vec.data(), BB_1, n) - 1;
            BS_0        = BB_1 + 1;

            if (BS_1 < BB_1)
                m_tasks.initial_task(BS_1,BB_1);
        };

        res             = m_tasks.eval(comp_f);
        *m              = m_tasks.number_selected();
    };

    if (res == false)
    {
        // Swap is rejected: exit.
        *info = 1;
    };

    compute_eigenvalues_tgsen3<V>::eval(n, K, a, lda, b, ldb, wantq, wantz, Q, 
                                        LDQ, Z, LDZ, alpha, beta);

    return;
};

template BLAS_EXT_EXPORT void
matcl::lapack::tgsen3<d_type>(i_type wantq,i_type wantz,const i_type *select,i_type n, i_type K,
        d_type *a,i_type lda,d_type *b,i_type ldb, z_type *alpha, d_type *beta, 
        d_type *Q, i_type LDQ, d_type* Z, i_type LDZ,i_type *m, d_type *work,i_type lwork, i_type *info);

template BLAS_EXT_EXPORT void
matcl::lapack::tgsen3<s_type>(i_type wantq,i_type wantz,const i_type *select,i_type n, i_type K,
        s_type *a,i_type lda,s_type *b,i_type ldb, c_type *alpha, s_type *beta, 
        s_type *Q, i_type LDQ, s_type* Z, i_type LDZ,i_type *m, s_type *work,i_type lwork, i_type *info);

template BLAS_EXT_EXPORT void
matcl::lapack::tgsen3<z_type>(i_type wantq,i_type wantz,const i_type *select,i_type n, i_type K,
        z_type *a,i_type lda,z_type *b,i_type ldb, z_type *alpha, z_type *beta, 
        z_type *Q, i_type LDQ, z_type* Z, i_type LDZ,i_type *m, z_type *work,i_type lwork, i_type *info);

template BLAS_EXT_EXPORT void
matcl::lapack::tgsen3<c_type>(i_type wantq,i_type wantz,const i_type *select,i_type n, i_type K,
        c_type *a,i_type lda,c_type *b,i_type ldb, c_type *alpha, c_type *beta, 
        c_type *Q, i_type LDQ, c_type* Z, i_type LDZ,i_type *m, c_type *work,i_type lwork, i_type *info);

};};