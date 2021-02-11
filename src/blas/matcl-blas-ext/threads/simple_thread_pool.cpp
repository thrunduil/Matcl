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

#include "simple_thread_pool.h"
#include "matcl-blas-ext/blas_concurrency.h"
#include <iostream>

namespace matcl { namespace details
{

worker::worker(simple_thread_pool &s)
    : m_pool(s) 
{ };

void worker::operator()()
{
    using lock_type = std::unique_lock<std::mutex>;
    using task_type = std::function<void()>;

    while(true)
    {
        task_type task;

        {
            lock_type   lock(m_pool.m_queue_mutex);
             
            // look for a work item
            while(m_pool.quiting() == false
                  && (m_pool.m_running_threads > m_pool.m_max_running_threads || m_pool.m_tasks.empty()))
            { 
                // if there are none wait for notification
                m_pool.m_condition.wait(lock);
            }
 
            if(m_pool.quiting() == true) 
            {
                m_pool.remove_worker(*this);
                return;
            };
        }
        {
            //there may be many threads at this point
            lock_type   lock(m_pool.m_queue_mutex);

            if (m_pool.m_tasks.size() == 0)
                continue;

            m_pool.new_task_executed();

            // get the task from the queue
            task = m_pool.m_tasks.front();
            m_pool.m_tasks.pop_front();
        }
 
        // execute the task
        //worker cannot throw exceptions
        try
        {
            task();
        }
        catch(std::exception& ex)
        {
            std::cerr<< "task failed: " << ex.what() << "\n";
        }
        catch(...)
        {
            std::cerr<< "task failed: unknown error";
        }

        m_pool.task_finished();
    }
}

void simple_thread_pool::resize(int n)
{
    using lock_type     = std::unique_lock<std::mutex>;
    using value_type    = thread_map::value_type;

    lock_type lock(m_thread_map_mutex);

    int old_max         = m_max_threads.load();
    int new_max         = std::max(1, std::min(n,get_default_threads(domain::matcl)));

    for (int i = old_max; i < new_max; ++i)
    {        
        thread_ptr th(new std::thread(worker(*this)));
        m_workers.insert(value_type(th->get_id(), th));
    };

    m_max_threads       = std::max(new_max, old_max);
    m_max_running_threads= new_max;
};

void simple_thread_pool::remove_worker(worker& w)
{
    --m_max_threads;
};

int simple_thread_pool::max_running_threads() const
{
    return m_max_running_threads;
};
// the constructor just launches some amount of workers
simple_thread_pool::simple_thread_pool(int threads)
    :   m_max_threads(std::max(1,std::min(threads,get_default_threads(domain::matcl))) )
    , m_max_running_threads(0), m_running_threads(0), m_quiting(0)
{
    using value_type    = thread_map::value_type;

    for(int i = 0; i < m_max_threads; ++i)
    {
        thread_ptr th(new std::thread(worker(*this)));
        m_workers.insert(value_type(th->get_id(), th));
    };

    m_max_running_threads = m_max_threads.load();
}
   
bool simple_thread_pool::quiting()
{
    return m_quiting == 1;
};
void simple_thread_pool::new_task_executed()
{
    ++m_running_threads;
};
void simple_thread_pool::task_finished()
{
    --m_running_threads;
};

// the destructor joins all threads
simple_thread_pool::~simple_thread_pool()
{
    m_quiting = 1;
    {
        std::unique_lock<std::mutex> lock(m_queue_mutex);

        while(m_max_threads > 0)
            m_condition.notify_all();
    };

    for (auto pos = m_workers.begin(); pos != m_workers.end(); ++pos)
    {
        pos->second->join();
    };

    m_workers.clear();
}

void simple_thread_pool::schedule(const function_type& f)
{
    using lock_type = std::unique_lock<std::mutex>;

    lock_type   lock(m_queue_mutex);
         
    m_tasks.push_back(f);

    // wake up one thread
    m_condition.notify_one();
}

}};
