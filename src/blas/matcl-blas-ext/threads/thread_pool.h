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

#include <atomic>
#include <functional>
#include "simple_thread_pool.h"

namespace matcl { namespace details
{
  
class static_thread_pool
{
    private:
        using atomic_int        = std::atomic<int>;
        using thread_pool       = simple_thread_pool;
        using thread_pool_ptr   = std::atomic<thread_pool*>;
        using mutex_type        = boost::detail::spinlock;
        using lock_type         = boost::detail::spinlock::scoped_lock;

    private:
        atomic_int      m_number_threads;
        thread_pool_ptr m_pool;
        mutex_type      m_pool_mutex;        

    public:
        static_thread_pool();
        ~static_thread_pool();

        void            concurent_eval(const std::function<void()>& f);

        int             get_num_threads();
        int             get_default_threads();
        void            set_num_threads(int n);

    private:
        static_thread_pool(const static_thread_pool&) = delete;
        static_thread_pool& operator=(const static_thread_pool&) = delete;
};

};};