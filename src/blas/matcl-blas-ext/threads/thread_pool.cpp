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

#include "thread_pool.h"
#include "matcl-blas-ext/blas_concurrency.h"
#include <thread>

namespace matcl { namespace details
{

static_thread_pool::static_thread_pool()
    :m_pool(nullptr), m_number_threads(get_default_threads())
{
    m_pool_mutex.v_.clear();
};
static_thread_pool::~static_thread_pool()
{
    delete m_pool.load();
};

int static_thread_pool::get_num_threads()
{
    if (m_pool.load() != nullptr)
        return m_pool.load()->max_running_threads();
    else
        return m_number_threads;
};
int static_thread_pool::get_default_threads()
{
    if (are_user_threads_allowed() == true)
        return std::thread::hardware_concurrency();
    else
        return 1;
};
void static_thread_pool::set_num_threads(int n)
{
    m_number_threads = std::max<int>(n,1);

    if (m_pool.load() != nullptr)
    {
        m_pool.load()->resize(n);
    };
};
void static_thread_pool::concurent_eval(const std::function<void()>& f)
{
    {        
        lock_type m_lock(m_pool_mutex);

        if (m_pool.load() == nullptr)
            m_pool = new thread_pool(this->get_num_threads());
    };

    m_pool.load()->schedule(f);
}

}};
