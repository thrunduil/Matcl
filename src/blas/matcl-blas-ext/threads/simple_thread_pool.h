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

#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <map>
#include <atomic>
#include <boost/shared_ptr.hpp>

namespace matcl { namespace details
{

class simple_thread_pool;
  
class worker 
{
    public:
        worker(simple_thread_pool& s);
        void operator()();

    private:
        simple_thread_pool&    m_pool;
};
  
class simple_thread_pool 
{
    private:
        using id_type       = std::thread::id;
        using thread_ptr    = boost::shared_ptr<std::thread>;
        using thread_map    = std::map<id_type, thread_ptr>;
        using function_type = std::function<void()>;
        using task_queue    = std::deque<function_type>;
        using condition_var = std::condition_variable;
        using atomic_int    = std::atomic<int>;

    public:
        simple_thread_pool(int n_threads);
        ~simple_thread_pool();

        void            schedule(const std::function<void()>& f);        
        void            resize(int n);
        int             max_running_threads() const;

    private:        
        void            remove_worker(worker& w);
        bool            quiting();
        void            new_task_executed();
        void            task_finished();

    private:
        friend class worker;
 
        thread_map      m_workers;
        task_queue      m_tasks;
        std::mutex      m_queue_mutex;
        std::mutex      m_thread_map_mutex;
        condition_var   m_condition;

        atomic_int      m_max_running_threads;
        atomic_int      m_running_threads;
        atomic_int      m_max_threads;
        atomic_int      m_quiting;
};

};};