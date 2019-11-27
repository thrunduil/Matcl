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

#include "tgsen_helpers.h"
#include <iostream>

#include <cassert>

namespace matcl { namespace lapack
{

#define assert2 assert

//-----------------------------------------------------------------
//                  analyse_select
//-----------------------------------------------------------------
template<class V, bool Is_compl = details::is_complex<V>::value>
struct analyse_select_impl
{
    static i_type eval(const i_type n, const i_type* select, const V* a, const i_type lda, i_type* sel_vec)
    {
        bool PAIR   = false;
        i_type M    = 0;

        for (i_type K = 0; K < n; ++K)
        {        
            if (PAIR == true)
            {
                PAIR        = false;
                sel_vec[K]  = 1;
            }
            else
            {
                if ( K < n - 1 )
                {
                    if( a[K+1 + K * lda] == 0. )
                    {
                        if (select[K])
                        {
                            M           = M + 1;
                            sel_vec[K]  = 1;
                        }
                        else
                        {
                            sel_vec[K]  = 0;
                        };
                    }
                    else
                    {                    
                        if (select[K] || select[K+1] )
                        {
                            PAIR        = true;
                            M           = M + 2;
                            sel_vec[K]  = 2;
                        }
                        else
                        {
                            sel_vec[K]  = 0;
                        };
                    };
                }
                else
                {
                    if (select[n-1])
                    {
                        sel_vec[K]  = 1;
                        M           = M + 1;
                    }
                    else
                    {
                        sel_vec[K]  = 0;
                    };
                };
            };
        };

        return M;
    };
};
template<class V>
struct analyse_select_impl<V, true>
{
    static i_type eval(const i_type n, const i_type* select, const V* a, const i_type lda, i_type* sel_vec)
    {
        i_type M    = 0;

        for (i_type K = 0; K < n; ++K)
        {        
            if ( K < n - 1 )
            {
                if (select[K])
                {
                    M           = M + 1;
                    sel_vec[K]  = 1;
                }
                else
                {
                    sel_vec[K]  = 0;
                };
            }
            else
            {
                if (select[n-1])
                {
                    sel_vec[K]  = 1;
                    M           = M + 1;
                }
                else
                {
                    sel_vec[K]  = 0;
                };
            };
        };

        return M;
    };
};

template<class V>
i_type analyse_select<V>
    ::eval(const i_type n, const i_type* select, const V* a, const i_type lda, i_type* sel_vec)
{
    return analyse_select_impl<V>::eval(n,select,a,lda,sel_vec);
};

template<class V, bool Is_compl = details::is_complex<V>::value>
struct analyse_select2_impl
{
    static i_type eval(const i_type n, i_type* sel_vec, const V* a, const i_type lda)
    {
        bool PAIR   = false;
        i_type M    = 0;

        for (i_type K = 0; K < n; ++K)
        {        
            if (PAIR == true)
            {
                PAIR        = false;
                sel_vec[K]  = 1;
            }
            else
            {
                if ( K < n - 1 )
                {
                    if( a[K+1 + K * lda] == 0. )
                    {
                        M               = M + 1;
                        sel_vec[K]      = 1;
                    }
                    else
                    {                    
                        PAIR            = true;
                        M               = M + 2;
                        sel_vec[K]      = 2;
                    };
                }
                else
                {
                    sel_vec[K]          = 1;
                    M                   = M + 1;
                };
            };
        };

        return M;
    };
};
template<class V>
struct analyse_select2_impl<V,true>
{
    static i_type eval(const i_type n, i_type* sel_vec, const V* a, const i_type lda)
    {
        i_type M    = 0;

        for (i_type K = 0; K < n; ++K)
        {        
            if ( K < n - 1 )
            {
                M               = M + 1;
                sel_vec[K]      = 1;
            }
            else
            {
                sel_vec[K]      = 1;
                M               = M + 1;
            };
        };

        return M;
    };
};

template<class V>
i_type analyse_select2<V>::eval(const i_type n, i_type* sel_vec, const V* a, const i_type lda)
{
    return analyse_select2_impl<V>::eval(n,sel_vec,a,lda);
};

template struct analyse_select<d_type>;
template struct analyse_select<s_type>;
template struct analyse_select<c_type>;
template struct analyse_select<z_type>;

template struct analyse_select2<d_type>;
template struct analyse_select2<s_type>;
template struct analyse_select2<c_type>;
template struct analyse_select2<z_type>;

//-----------------------------------------------------------------
//                  find_block_end
//-----------------------------------------------------------------
i_type find_block_end(i_type BS, i_type n, int max_eig_in_window, const i_type* sel_vec, 
                      i_type& eig_in_block)
{
    i_type n_eigs  = 0;
    i_type BE      = n;

    for (i_type K = BS; K <= n; ++K)
    {
        i_type SEL      = sel_vec[K-1];
        if (SEL != 0)
        {
            n_eigs      = n_eigs + 1;
            BE          = K;
            
            if (n_eigs >= max_eig_in_window)
            {
                if (SEL == 2)
                {
                    n_eigs  = n_eigs + 1;
                    BE      = K+1;
                };

                break;
            };
        };
    };

    eig_in_block = n_eigs;

    return BE;
};

//-----------------------------------------------------------------
//                  make_id
//-----------------------------------------------------------------
template<class V>
void lapack::make_id(V* mat, i_type n)
{
    memset(mat, 0, n * n * sizeof(V));

    for (i_type k = 0; k < n; ++k)
    {
        mat[k]  = V(1.0);
        mat     = mat + n;
    };
};

template void lapack::make_id<d_type>(d_type* mat, i_type n);
template void lapack::make_id<s_type>(s_type* mat, i_type n);
template void lapack::make_id<c_type>(c_type* mat, i_type n);
template void lapack::make_id<z_type>(z_type* mat, i_type n);

//-----------------------------------------------------------------
//                  ignore_sorted_eigenvalues
//-----------------------------------------------------------------
i_type ignore_sorted_eigenvalues(const i_type* sel_vec, i_type s, i_type n)
{
    i_type start_pos = s;

    for (i_type K = s; K <= n; ++K)
    {
        if (sel_vec[K-1] != 0)
        {
            start_pos   = start_pos + 1;
        }
        else
            break;
    };
    return start_pos;
};

//-----------------------------------------------------------------
//                  find_window_start
//-----------------------------------------------------------------
template<class V, bool Is_compl = details::is_complex<V>::value>
struct find_window_start_impl
{
    static i_type eval(i_type BS, i_type WE, i_type max_window_size, const i_type* sel_vec,
                         const V* a, const i_type lda)
    {
        i_type WS   = lapack::maximum(BS, WE - max_window_size + 1);

        while(sel_vec[WS-1] != 0)
        {
            // do not reorder eigenvalues in top-left part of the window
            WS      = WS + 1;
        };

        if (WS > 1)
        {
            if( a[WS-1 + (WS-2) * lda] != 0. )
            {
                // do not split complex eigenvalues
                WS  = WS - 1;
            };
        };

        return WS;
    };
};
template<class V>
struct find_window_start_impl<V,true>
{
    static i_type eval(i_type BS, i_type WE, i_type max_window_size, const i_type* sel_vec,
                         const V* a, const i_type lda)
    {
        i_type WS   = lapack::maximum(BS, WE - max_window_size + 1);

        while(sel_vec[WS-1] != 0)
        {
            // do not reorder eigenvalues in top-left part of the window
            WS      = WS + 1;
        };

        return WS;
    };
};

template<class V>
i_type lapack::find_window_start(i_type BS, i_type WE, i_type max_window_size, const i_type* sel_vec,
                         const V* a, const i_type lda)
{
    return find_window_start_impl<V>::eval(BS, WE, max_window_size, sel_vec, a, lda);
};

template
i_type lapack::find_window_start<d_type>(i_type BS, i_type WE, i_type max_window_size, const i_type* sel_vec,
                         const d_type* a, const i_type lda);
template
i_type lapack::find_window_start<s_type>(i_type BS, i_type WE, i_type max_window_size, const i_type* sel_vec,
                         const s_type* a, const i_type lda);
template
i_type lapack::find_window_start<c_type>(i_type BS, i_type WE, i_type max_window_size, const i_type* sel_vec,
                         const c_type* a, const i_type lda);
template
i_type lapack::find_window_start<z_type>(i_type BS, i_type WE, i_type max_window_size, const i_type* sel_vec,
                         const z_type* a, const i_type lda);

//-----------------------------------------------------------------
//                  get_density
//-----------------------------------------------------------------
template<class V>
d_type lapack::get_density(const V* mat, const i_type n)
{
    i_type n2   = n * n;
    i_type nz   = 0;

    //calculate nonzeros
    for (i_type k = 0; k < n2; ++k)
    {
        if (mat[k] != V(0.))
            ++nz;
    };

    d_type density = d_type(nz) / d_type(n2);

    return density;
};

template d_type lapack::get_density<d_type>(const d_type* mat, const i_type n);
template d_type lapack::get_density<s_type>(const s_type* mat, const i_type n);
template d_type lapack::get_density<c_type>(const c_type* mat, const i_type n);
template d_type lapack::get_density<z_type>(const z_type* mat, const i_type n);

//-----------------------------------------------------------------
//                  access
//-----------------------------------------------------------------
access::access()
    :m_owner(nullptr), m_rows_lock(false)
{};
access::access(matrix_lock* lock, bool rows)
    :m_owner(lock), m_rows_lock(rows)
{};
access::~access()
{
    if (m_owner)
        m_owner->release(m_rows_lock);
};

access::access(access&& other)    
{
    m_owner         = other.m_owner;
    m_rows_lock     = other.m_rows_lock;
    other.m_owner   = nullptr;
};

bool access::is_granted() const
{
    return m_owner != nullptr;
};

//-----------------------------------------------------------------
//                  matrix_lock
//-----------------------------------------------------------------
matrix_lock::matrix_lock()
{
    m_locked_rows   = 0;
    m_locked_cols   = 0;
    m_mutex.v_.clear();
};

access matrix_lock::get_row()
{
    scoped_lock lock(m_mutex);

    if (m_locked_cols.load() != 0)
    {        
        return access();
    }
    else
    {
        m_locked_rows += 1;
        return access(this,true);
    };
};
access matrix_lock::get_col()
{
    scoped_lock lock(m_mutex);

    if (m_locked_rows.load() != 0)
    {        
        return access();
    }
    else
    {
        m_locked_cols += 1;
        return access(this,false);
    };
};
void matrix_lock::release(bool row)
{
    scoped_lock lock(m_mutex);

    if (row)
        m_locked_rows -= 1;
    else
        m_locked_cols -= 1;
};

//-----------------------------------------------------------------
//                  workspace_manager
//-----------------------------------------------------------------
template<class V>
workspace_manager<V>::workspace_manager(V* work, i_type n_threads, i_type work_size)
{
    for (int i = 0; i < n_threads; ++i)
    {
        m_work.push_back(work);
        work = work + work_size;
    };
    m_mutex_vec.v_.clear();
};
template<class V>
V* workspace_manager<V>::get()
{
    std::unique_lock<std::mutex> lock(m_mutex);
    V* work = nullptr;

    while(work == nullptr)
    {
        scoped_lock lock(m_mutex_vec);

        if (this->m_work.size() > 0)
        {
            work = m_work.back();
            m_work.pop_back();
        };
    };

    return work;
};
template<class V>
void workspace_manager<V>::release(V* work)
{
    scoped_lock lock(m_mutex_vec);
    m_work.push_back(work);
};

template class workspace_manager<d_type>;
template class workspace_manager<s_type>;
template class workspace_manager<c_type>;
template class workspace_manager<z_type>;

//-----------------------------------------------------------------
//                  task_manager
//-----------------------------------------------------------------
template<class V>
task_manager<V>::task_manager(workspace_manager<V>* wm)
    :m_workspace(wm)
{
    m_error         = false;
    m_pending       = 0;
};

template<class V>
task_manager<V>::~task_manager()
{};

template<class V>
void  task_manager<V>::initial_task(i_type BS, i_type BB)
{
    task<V> t;
    t.BB        = BB;
    t.BS        = BS;
    t.code      = (int)m_tasks.size();
    t.done      = false;

    m_tasks.push_back(t);
};

template<class V>
void task_manager<V>::enqueue_task(i_type code)
{
    ++m_pending;

    auto task_f = [this, code]
                () -> void
                {
                    bool res        = false;
                    V* work         = this->m_workspace->get();

                    assert2(work != nullptr);
                    assert2(code < (int)this->m_tasks.size() && code >= 0);
                    
                    task<V> t       = this->m_tasks[code];

                    try
                    {                        
                        res         = this->m_function(t.BS, t.BB, work);
                    }
                    catch(std::exception& ex)
                    {
                        std::cerr << "reordering failed, reason: " << ex.what() << "\n";
                    }
                    catch(...)
                    {
                        std::cerr << "reordering failed, reason: unknown error" << "\n";
                    }

                    this->m_workspace->release(work);                    

                    this->set_result(t,res);
                    return;
                };

    concurrent_eval(task_f);
};

template<class V>
void task_manager<V>::set_result(const task<V>& t, bool is_ok)
{
    task_lock lock(m_task_mutex);

    --m_pending;

    i_type code             = t.code;

    assert2(code < (int)this->m_tasks.size() && code >= 0);

    task<V>& tc             = m_tasks[code];
    tc.done                 = true;
    tc.BS                   = t.BS;

    if (is_ok == false)
    {
        m_error             = true;
    };

    if (m_error == false)
    {
        if (code > 0)
        {
            task<V>& tp     = m_tasks[code-1];
            tp.BB           = t.BS - 1;
        
            if (tp.done == true)
            {
                tp.done     = false;
                enqueue_task(code-1);
            };
        };

        if (tc.BB != t.BB)
        {
            tc.done         = false;
            enqueue_task(code);
        };
    };

    //all tasks must finish 
    if (m_pending == 0)
    {
        m_condition_done.notify_all();
    };
};

template<class V>
bool task_manager<V>::eval(function_t f)
{    
    task_lock lock_task(m_task_mutex);

    m_function = f;

    for (size_t i = 0; i < m_tasks.size(); ++i)
    {
        //cheap operation
        enqueue_task(int(i));
    };

    //avoid spurious wake-ups
    while(m_pending > 0)
    {
        //this unlocks m_task_mutex
        m_condition_done.wait(lock_task);
    };

    return (m_error == false);
};

template<class V>
i_type task_manager<V>::number_selected() const
{
    if (m_tasks.size() == 0)
        return 0;

    return m_tasks[0].BS - 1;
};

template class task_manager<d_type>;
template class task_manager<s_type>;
template class task_manager<c_type>;
template class task_manager<z_type>;

};};