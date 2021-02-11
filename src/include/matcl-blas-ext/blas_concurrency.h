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

#pragma once

#include "matcl-blas-ext/config_blas_ext.h"
#include "matcl-blas-lapack/blas/blas.h"
#include <functional>
#include <atomic>
#include <condition_variable>
#include <vector>
#include <boost/smart_ptr/detail/spinlock.hpp>

#pragma warning(push)
#pragma warning(disable: 4251) // needs to have dll-interface to be used by clients of class

namespace matcl
{

enum class domain
{
    blas_kernel, matcl
};

// return number of threads, that can run together in given domain
BLAS_EXT_EXPORT int     get_num_threads(domain dom);

// return default value of number of threads, that can run together in given domain
BLAS_EXT_EXPORT int     get_default_threads(domain dom);

// set number of threads, that can run together in given domain
BLAS_EXT_EXPORT void   set_num_threads(int n, domain dom);

// return if the user of blas kernel can create threads; some blas kernels
// breaks down if are called in concurrent environment
BLAS_EXT_EXPORT bool   are_user_threads_allowed();

// evaluate given function by thread pool; function f must be thread safe and
// may not throw exceptions
BLAS_EXT_EXPORT void   concurrent_eval(const std::function<void()>& f);

// change number of threads that can be executed together; at exit
// reset previous value
class BLAS_EXT_EXPORT concurrent_scope
{
    private:
        int         m_old_kernel_threads;

    public:
        concurrent_scope();
        ~concurrent_scope();

    private:
        concurrent_scope(const concurrent_scope&) = delete;
        concurrent_scope& operator=(const concurrent_scope&) = delete;
};

// execute a set of functions concurently
class BLAS_EXT_EXPORT tast_group
{
    private:
        using atomic_int    = std::atomic<int>;
        using condition     = std::condition_variable;
        using mutex_type    = std::mutex;
        using lock_type     = std::unique_lock<mutex_type>;

    private:
        concurrent_scope    m_old_scope;
        atomic_int          m_pending_tasks;
        condition           m_condition_done;
        std::mutex          m_task_mutex;

    public:
        tast_group();
        ~tast_group();

        // add a new function to the task group
        void            add(const std::function<void()>& f);

        // wait until all functions in the task group are evaluated
        void            wait();

    private:
        tast_group(const tast_group&) = delete;
        tast_group& operator=(const tast_group&) = delete;

        void            task_finished();
};

//pool of workspaces available for use by concurrently executed tasks
class BLAS_EXT_EXPORT workspace_pool
{
    private:
        using work_vec      = std::vector<void*>;
        using mutex_type    = boost::detail::spinlock;
        using lock_type     = std::unique_lock<mutex_type>;

    private:
        work_vec            m_work;
        mutex_type          m_mutex_vec;

    public:
        // class responsible for releasing workspace at the end of its lifetine;
        // one can access stored workspace only from nontemporary instance of this
        // class; when destructor is called, stored workspace is released and cannot
        // be used any longer; workspace_ptr cannot live longer than workspace_pool
        // and must be destroyed before concurrently executed task that use this workspace
        // if completed
        class workspace_ptr
        {
            private:
                void*               m_ptr;
                workspace_pool*     m_owner;

            private:
                workspace_ptr(void* ptr, workspace_pool* owner);

                workspace_ptr(const workspace_ptr&) = delete;
                workspace_ptr& operator=(const workspace_ptr&) = delete;

                friend workspace_pool;

            public:
                // initialize workspace to nullptr
                workspace_ptr();

                // release workspace 
                ~workspace_ptr();

                // workspace pointer is movable and move assignable
                workspace_ptr(workspace_ptr&& other);
                workspace_ptr& operator=(workspace_ptr&& other);

                // get stored workspace
                void*   get() &;

                // release stored workspace; workspace cannot be used any longer
                void    release();
        };

    public:
        // create empty worskpace pool
        workspace_pool();
        ~workspace_pool();

        // add new workspace, workspace pool must be populated before use,
        // this function can be called only in the main thread; there should be 
        // as many workspaces as number of working threads; class workspace pool
        // is not responsible for freeing memory
        void            add(void* work);

        // get workspace (if available) to use by a concurrent task
        workspace_ptr   get();

    private:
        workspace_pool(const workspace_pool&) = delete;
        workspace_pool& operator=(const workspace_pool&) = delete;

        void            release(void* work);
};

};

#pragma warning(pop)