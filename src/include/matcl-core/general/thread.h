/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-core/config.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/matrix/scalar_types.h"
#include <mutex>
#include <atomic>

namespace matcl
{

// default spinlock mutex type for multi-threaded applications based
// on boost::spinlock
class MATCL_CORE_EXPORT spinlock_mutex
{
    private:
        long        m_mutex;

    public:
        spinlock_mutex();

        spinlock_mutex(const spinlock_mutex&) = delete;
        spinlock_mutex& operator=(const spinlock_mutex&) = delete;

        bool        try_lock();
        void        lock();
        void        unlock();
};

// default mutex type for single threaded applications
class nonblocking_mutex
{
    public:
        nonblocking_mutex(){};

        nonblocking_mutex(const nonblocking_mutex&) = delete;
        nonblocking_mutex& operator=(const nonblocking_mutex&) = delete;

        bool        try_lock()  {return true;}
        void        lock()      {};
        void        unlock()    {};
};

// implements interface of atomics in single threaded version of matcl
template<class V>
struct atomic_nothread
{
    public:
        atomic_nothread();
        atomic_nothread(V v);

        V                   load(std::memory_order = std::memory_order_seq_cst) const;
        void                store(V val, std::memory_order = std::memory_order_seq_cst );
                            operator V() const;

        atomic_nothread&    operator++();
        atomic_nothread&    operator--();

        atomic_nothread&    operator=(V val);

        template<class T>
        V                   fetch_add(T inc, std::memory_order = std::memory_order_seq_cst);

        template<class T>
        V                   fetch_sub(T inc, std::memory_order = std::memory_order_seq_cst);

    private:
        atomic_nothread(const atomic_nothread&)             = delete;
        atomic_nothread& operator=(const atomic_nothread&)  = delete;

    private:
        V               m_data;
};

// default blocking mutex
using mutex = std::mutex;

#if MATCL_THREADING_MODE == MATCL_MULTITHREAD_MODE
    using default_mutex             = mutex;
    using default_spinlock_mutex    = spinlock_mutex;

    template<class T>
    using atomic                    = std::atomic<T>;

    #define MATCL_THREAD_LOCAL thread_local

#else

    using default_mutex             = nonblocking_mutex;
    using default_spinlock_mutex    = nonblocking_mutex;

    template<class T>
    using atomic                    = atomic_nothread<T>;

    #define MATCL_THREAD_LOCAL

#endif

};

#include "matcl-core/details/thread.inl"