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
#include "matcl-blas-lapack/blas/config_blas.h"
#include "tgsen_helpers.h"
#include "blas/matcl-blas-ext/lapack_ext/utils/optim_params.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-core/error/exception.h"

#include <vector>
#include <cassert>

namespace matcl { namespace lapack
{

#define assert2 matcl_assert

#define MAGIC_NUMBER 12345.

template<class V, bool Is_compl = details::is_complex<V>::value>
struct compute_eigenvalues_trsen3
{
    static void eval(i_type n, V* t, i_type ldt, typename details::complex_type<V>::type* w)
    {
        using VC = typename details::complex_type<V>::type;

        for(i_type k = 0; k < n; ++k)
            w[k]    = VC(t[k + k*ldt], 0.);

        for (i_type k = 0; k < n-1; ++k)
        {
            if( t[k+1 + k*ldt] != 0.)
            {
                V t1        = t[k + (k+1) * ldt];
                V t2        = t[k+1 + k * ldt];
                V t3        = sqrt(abs(t1)) * sqrt(abs(t2));
                w[k].imag(t3);
                w[k+1].imag(-t3);
            };
        };
    };
};

template<class V>
struct compute_eigenvalues_trsen3<V,true>
{
    static void eval(i_type n, V* t, i_type ldt, typename details::complex_type<V>::type* w)
    {
        using VC = typename details::complex_type<V>::type;

        for(i_type k = 0; k < n; ++k)
            w[k]    = t[k + k*ldt];
    };
};

template<class V>
static bool compute_block(i_type& BS, i_type BB, i_type n, i_type K, i_type max_eig_in_window, 
            i_type max_window_size, i_type* sel_vec, V* t, i_type ldt, V* Q, i_type LDQ, 
            i_type wantq, matrix_lock* t_lock, V* work, i_type& lwork)
{    
    using VC = typename details::complex_type<V>::type;

    i_type LWMIN_q      = (max_window_size+1) * (max_window_size + 1);
    i_type LWMIN_t      = LWMIN_q;
    i_type NBQ          = lapack::minimum(K,2 * (max_window_size+1));
    i_type LWMIN_m      = (max_window_size+1) * n;
    i_type LWMIN_AB     = max_window_size+1;
    i_type LWMIN_l      = max_window_size+1;
    i_type COMPL_SIZE   = sizeof(VC)/sizeof(V);
    i_type LWMIN        = LWMIN_q + lapack::maximum(LWMIN_m, LWMIN_t + COMPL_SIZE*LWMIN_AB + LWMIN_l);

    if (lwork == -1)
    {
        lwork           = LWMIN;
        lwork           = lwork + 1;
        return true;
    };

    bool is_V_real          = details::is_complex<V>::value == false;
    const char* ctrans_name = is_V_real? "T" : "C";

    if (BS > 1)
    {
        //complex eigenvalues cannot be split
        assert2(is_equal(t[BS - 1 + (BS - 2)*ldt], V(0)) == true, "");
    };

    i_type BS_save      = BS;

    V* work_test        = work + LWMIN;
    V* work_q           = work; work = work + LWMIN_q;

    V* work_mat_T1      = work; //wsize x (n - WS + 1 - wsize) <= N x wsize
    V* work_mat_T2      = work; // (WS - 1) x wsize <= N x wsize
    V* work_mat_Q       = work; // <= NBQ x wsize <= N x wsize

    V* work_T           = work; work = work + LWMIN_t;
    V* wre              = work; work = work + LWMIN_AB * COMPL_SIZE;
    V* work_LAP         = work;

    VC* w               = reinterpret_cast<VC*>(wre);

    work_test[0]        = MAGIC_NUMBER;

    assert2(sel_vec[BS-1] == 0, "");

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
            WS          = find_window_start(BS,WE,max_window_size,sel_vec,t,ldt);

            // reorder current window
            i_type wsize            = WE - WS + 1;

            i_type sel_1 = sel_vec[WS - 1];
            i_type sel_2 = sel_vec[WE - 1];

            assert2(sel_1 == 0 && sel_2 != 0, "");

            // prepare q matrix
            make_id(work_q, wsize);

            // use unblocked version
            i_type n_reorder        = 0;
            i_type tmp_info         = 0;
            i_type tmp_liwork       = 0;

            // operate on local workpsace to avoid cache misses
            lacpy("N", wsize, wsize, t + WS - 1 + (WS - 1) * ldt, ldt, work_T, wsize);

            trsen<V>("N", "V", sel_vec + WS - 1, wsize, work_T, wsize, work_q, wsize,
                    w, &n_reorder, nullptr, nullptr, work_LAP, LWMIN_l, &tmp_liwork, 1, &tmp_info);

            lacpy("N", wsize, wsize, work_T, wsize, t + WS - 1 + (WS - 1) * ldt, ldt);

            if ( tmp_info > 0 )
            {
                // Swap is rejected: exit.
                return false;
            };

            bool calc               = false;
            bool calc_T1            = calc;
            bool calc_T2            = calc;
            bool calc_Q             = (wantq == false);

            while (!calc_T1 || !calc_T2 || !calc_Q)
            {
                // apply delayed rotations
                // matrix TA
                if (n - WS + 1 - wsize <= 0)
                {
                    calc_T1 = true;
                }
                else if (calc_T1 == false)
                {
                    access tl = t_lock->get_row();

                    if (tl.is_granted() == true)
                    {
                        lacpy("N", wsize, n - WS + 1 - wsize, t + WS - 1 + (WS - 1 + wsize) * ldt, 
                            ldt, work_mat_T1, wsize);
                        gemm<V>(ctrans_name, "N", wsize, n - WS + 1 - wsize, wsize, 1.0, work_q, wsize, 
                                work_mat_T1, wsize, 0.0, t + WS - 1 + (WS - 1 + wsize) * ldt, ldt);

                        calc_T1 = true;
                    };
                };
                if (WS - 1 <= 0)
                {
                    calc_T2 = true;
                }
                else if (calc_T2 == false)
                {
                    access tl = t_lock->get_col();

                    if (tl.is_granted() == true)
                    {
                        lacpy("N", WS - 1, wsize, t + (WS - 1) * ldt, ldt, work_mat_T2, WS - 1);
                        gemm<V>("N", "N", WS - 1 , wsize, wsize, 1.0, work_mat_T2, WS - 1, 
                                work_q, wsize, 0.0, t + (WS - 1) * ldt, ldt);

                        calc_T2 = true;
                    };
                };

                if (calc_Q == false)
                {
                    // matrix Q
                    for (i_type j = 0; j < K; j += NBQ)
                    {
                        i_type JNBQ = lapack::minimum(NBQ, K - j);
                        
                        lacpy("N", JNBQ, wsize, Q + j + (WS - 1)*LDQ, LDQ, work_mat_Q, JNBQ);
                        gemm<V>("N", "N", JNBQ, wsize, wsize, 1.0, work_mat_Q, JNBQ, work_q, wsize, 
                            0.0, Q + j + (WS - 1)*LDQ, LDQ);
                    };

                    calc_Q = true;
                    continue;
                };
            };

            // move to next window
            WE              = WS + n_reorder - 1;

            analyse_select2<V>::eval(n_reorder, sel_vec + WS - 1, t + WS - 1 + (WS - 1)*ldt, ldt);

            assert2(WS >= BS_save, "");
            assert2(WS + n_reorder - 1 <= BB, "");
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
            assert2(is_equal(t[BS - 1 + (BS - 2)*ldt], V(0)) == true, "");
        };

        assert2(BS + n_eigs >= BS_save, "");
        assert2(BE <= BB, "");
        assert2(sel_vec[BS-1] == 0, "");
    };

    assert2(is_equal(work_test[0], V(MAGIC_NUMBER)) == true, "");

    return true;
};

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::trsen3(const char* compq, const i_type *select, i_type n, i_type K, V* t, i_type ldt, V* Q, i_type LDQ, 
        typename details::complex_type<V>::type *w, i_type* m, V* work, i_type lwork, i_type* info)
{
    using VR    = typename details::real_type<V>::type;

    //parameters of the algorithm
    const i_type max_eig_in_window  = optim_params::dtrsen3_max_eig_in_window;
    const i_type max_window_size    = optim_params::dtrsen3_max_window_size;
    const i_type max_threads        = get_num_threads(domain::matcl);
    //const i_type max_threads        = 1;

    bool wantq 	    = (compq[0] == 'V' || compq[0] == 'v');

    i_type INFO = 0;

    if (compq[0] != 'N' && compq[0] != 'n' && compq[0] != 'V' && compq[0] != 'v')
    {
        INFO    = -1;
    }
    else if ( n < 0 )
        INFO    = -3;
    else if ( K < 0 )
        INFO    = -4;
    else if (ldt < lapack::maximum( 1, n ) )
        INFO    = -6;
    else if ( LDQ < 1 || ( wantq && LDQ < K ) )
        INFO    = -8;
    
    if (INFO != 0 )
    {
        *info = INFO;
        return;
    };

    if (n <= max_window_size && K == n)
    {
        i_type tmp_iwork = 0;

        lapack::trsen<V>("N", compq, select, n, t, ldt, Q, LDQ, w, m, nullptr, nullptr, 
            work, lwork, &tmp_iwork, 1, info);

        return;
    };

    bool LQUERY     = ( lwork == -1 );

    // query workspace
    i_type LWMIN    = -1;    
    i_type BS_dum   = 1;

    compute_block<V>(BS_dum, 0, n, K, max_eig_in_window, max_window_size,nullptr, t, ldt, Q, LDQ, 
                  wantq, nullptr, nullptr, LWMIN);

    // Set M to the dimension of the specified pair of deflating subspaces.
    using workspace     = matcl::pod_workspace<i_type>;

    workspace sel_vec   = workspace(n + 1);

    i_type M            = analyse_select<V>::eval(n, select, t, ldt, sel_vec.ptr());

    // Collect the selected blocks at the top-left corner of (A, B).
    i_type start_pos    = ignore_sorted_eigenvalues(sel_vec.ptr(), 1, n);
    i_type tot_eigs     = 0;

    find_block_end(start_pos, n, n, sel_vec.ptr(), tot_eigs);

    i_type eigs_thread  = lapack::maximum(tot_eigs / max_threads / 2, max_eig_in_window * 2);
    i_type n_threads    = lapack::maximum(tot_eigs / eigs_thread, 1);
    n_threads           = lapack::minimum(n_threads, max_threads);

    i_type LIWMIN   = 1;
    work[0]         = V(VR(LWMIN * n_threads));

    if( lwork < LWMIN && LQUERY == false )
        INFO        = -12;

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
        compute_eigenvalues_trsen3<V>::eval(n, t, ldt, w);
        return;
    };
	
    i_type IERR         = 0;

    using work_manager  = workspace_manager<V>;

    matrix_lock         t_lock;
    work_manager        m_workspace(work, n_threads,LWMIN);    

    auto comp_f         = [n,K,max_eig_in_window,max_window_size, &sel_vec,t,ldt,Q,LDQ,
                            wantq, &t_lock, &lwork
                          ]
                        (i_type& BS, i_type BB, V* work) -> bool
                        { 
                            return compute_block<V>(BS, BB, n, K, max_eig_in_window, max_window_size, 
                                    sel_vec.ptr(), t, ldt, Q, LDQ, wantq, &t_lock, work, lwork);
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
        task_manager<V> m_tasks(&m_workspace);

        // initialize tasks
        i_type BS_0     = start_pos;
        while (BS_0 < n)
        {
            i_type dum  = 0;        
            i_type BS_1 = BS_0;
            i_type BB_1 = find_block_end(BS_1, n, eigs_thread, sel_vec.ptr(), dum);
            BB_1        = ignore_sorted_eigenvalues(sel_vec.ptr(), BB_1, n) - 1;
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
    
    compute_eigenvalues_trsen3<V>::eval(n, t, ldt, w);

    return;
};

template BLAS_EXT_EXPORT void
lapack::trsen3<d_type>(const char* compq, const i_type *select, i_type n, i_type k, d_type* t, i_type ldt, 
        d_type* q, i_type ldq, z_type* w, i_type* m, d_type* work, i_type lwork, i_type* info);

template BLAS_EXT_EXPORT void
lapack::trsen3<s_type>(const char* compq, const i_type *select, i_type n, i_type k, s_type* t, i_type ldt, 
        s_type* q, i_type ldq, c_type* w, i_type* m, s_type* work, i_type lwork, i_type* info);

template BLAS_EXT_EXPORT void
lapack::trsen3<c_type>(const char* compq, const i_type *select, i_type n, i_type k, c_type* t, i_type ldt, 
        c_type* q, i_type ldq, c_type* w, i_type* m, c_type* work, i_type lwork, i_type* info);

template BLAS_EXT_EXPORT void
lapack::trsen3<z_type>(const char* compq, const i_type *select, i_type n, i_type k, z_type* t, i_type ldt, 
        z_type* q, i_type ldq, z_type* w, i_type* m, z_type* work, i_type lwork, i_type* info);

};};