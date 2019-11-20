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
#include "matcl-blas-ext/blas_concurrency.h"
#include <boost/smart_ptr/detail/spinlock.hpp>
#include <atomic>
#include <vector>
#include <mutex>
#include <condition_variable>

#ifndef _MSC_VER
    #include <assert.h>
    #include <string.h> // memset
#endif
namespace matcl { namespace lapack
{

template<class V>
struct analyse_select
{
    static i_type eval(const i_type n, const i_type* select, const V* a, const i_type lda, 
                      i_type* sel_vec);
};

template<class V>
struct analyse_select2
{
    static i_type eval(const i_type n, i_type* sel_vec, const V* a, const i_type lda);
};

i_type  find_block_end(i_type BS, i_type n, int max_eig_in_window, const i_type* sel_vec, 
                      i_type& eig_in_block);

template<class V>
void    make_id(V* mat, i_type n);

i_type  ignore_sorted_eigenvalues(const i_type* sel_vec, i_type s, i_type n);

template<class V>
i_type  find_window_start(i_type BS, i_type WE, i_type max_window_size, const i_type* sel_vec,
                         const V* a, const i_type lda);

template<class V>
d_type  get_density(const V* mat, const i_type n);

class matrix_lock;

template<class V>
bool is_equal(const V& a, const V& b)
{
    return a == b;
};

//not shared between threads
class access
{
    private:
        matrix_lock*    m_owner;
        bool            m_rows_lock;

    public:
        ~access();
        access(access&&);

        bool    is_granted() const;

    private:
        access();
        access(matrix_lock* lock, bool rows);

        access(const access&);
        access& operator=(const access&);

        friend class matrix_lock;
};

//shared between threads
class matrix_lock
{
    private:
        using scoped_lock   = boost::detail::spinlock::scoped_lock;
        using mutex         = boost::detail::spinlock;
        using atomic_long   = std::atomic<long>;

    private:
        atomic_long m_locked_rows;
        atomic_long m_locked_cols;
        mutex       m_mutex;

    public:
        matrix_lock();

        access      get_row();
        access      get_col();

    private:
        matrix_lock(const matrix_lock&);
        matrix_lock& operator=(const matrix_lock&);

        void        release(bool row);

        friend class access;
};

//shared between threads
template<class V>
class workspace_manager
{
    private:
        using work_vec      = std::vector<V*>;
        using scoped_lock   = boost::detail::spinlock::scoped_lock;
        using mutex_vec     = boost::detail::spinlock;

    private:
        std::mutex  m_mutex;
        work_vec    m_work;
        mutex_vec   m_mutex_vec;

    public:
        workspace_manager(V* work, i_type n_threads, i_type work_size);

        V*          get();
        void        release(V* work);
};

template<class V>
class task_manager;

template<class V>
class task
{
    private:
        i_type      code;
        bool        done;

    public:
        i_type      BS;
        i_type      BB;

        friend task_manager<V>;
};

//shared between threads
template<class V>
class task_manager
{
    private:
        using task_lock     = std::unique_lock<std::mutex>;
        using condition     = std::condition_variable;
        using function_t    = std::function<bool (i_type& BS, i_type BB, V* work)>;
        using atomic_int    = std::atomic<int>;

    private:
        std::vector<task<V>>    m_tasks;

        workspace_manager<V>*   m_workspace;
        function_t              m_function;        

        std::mutex              m_task_mutex;
        condition               m_condition_done;

        atomic_int              m_pending;
        atomic_int              m_error;

        concurrent_scope        m_scope;

    public:
        task_manager(workspace_manager<V>* wm);
        ~task_manager();

        void    set_result(const task<V>& t, bool is_ok);

        void    initial_task(i_type BS, i_type BB);
        bool    eval(function_t f);
        i_type  number_selected() const;

    private:
        task_manager(const task_manager&);
        task_manager& operator=(const task_manager&);

        void    enqueue_task(i_type code);
};

};};