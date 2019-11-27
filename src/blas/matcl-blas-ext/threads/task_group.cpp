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

#include "matcl-blas-ext/blas_concurrency.h"

namespace matcl
{

tast_group::tast_group()
    :m_pending_tasks(0)
{};
tast_group::~tast_group()
{};

// add a new function to the task group
void tast_group::add(const std::function<void()>& f)
{
    auto task   = [f,this]()
                {
                    f();
                    this->task_finished();
                };

    ++m_pending_tasks;
    concurrent_eval(task);
};

void tast_group::task_finished()
{
    lock_type lock(m_task_mutex);

    --m_pending_tasks;

    //all tasks must finish 
    if (m_pending_tasks == 0)
    {
        m_condition_done.notify_all();
    };
};
void tast_group::wait()
{
    lock_type lock_task(m_task_mutex);

    //avoid spurious wake-ups
    while(m_pending_tasks > 0)
    {
        //this unlocks m_task_mutex
        m_condition_done.wait(lock_task);
    };
};

};