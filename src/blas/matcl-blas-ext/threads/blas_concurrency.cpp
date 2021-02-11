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

#include "matcl-blas-ext\blas_concurrency.h"
#include "matcl-blas-lapack/blas_loader/blas_loader.h"
#include "thread_pool.h"

namespace matcl
{

//destructor of static_thread_pool cannot be called during dll 
//unloading, static pointer solves this problem with a cost of
//small memory leak
static details::static_thread_pool* m_pool = new details::static_thread_pool();

BLAS_EXT_EXPORT int get_num_threads(domain dom)
{
    switch (dom)
    {
        case domain::blas_kernel:
            return raw_blas_lapack::get_num_threads_blas_kernel();
        case domain::matcl:
            return m_pool->get_default_threads();
        default:
        {
            //unknown case
            throw;
        }
    };
};
BLAS_EXT_EXPORT int get_default_threads(domain dom)
{
    switch (dom)
    {
        case domain::blas_kernel:
            return raw_blas_lapack::get_default_threads_blas_kernel();
        case domain::matcl:
            return m_pool->get_default_threads();            
        default:
        {
            //unknown case
            throw;
        }
    };
};
BLAS_EXT_EXPORT void set_num_threads(int n, domain dom)
{
    switch (dom)
    {
        case domain::blas_kernel:
            return raw_blas_lapack::set_num_threads_blas_kernel(n);
        case domain::matcl:
            return m_pool->set_num_threads(n);
        default:
            return;
    };
};
BLAS_EXT_EXPORT bool are_user_threads_allowed()
{
    return raw_blas_lapack::are_user_threads_allowed();
};

BLAS_EXT_EXPORT void concurrent_eval(const std::function<void()>& f)
{
    m_pool->concurent_eval(f);
};

concurrent_scope::concurrent_scope()
{
    m_old_kernel_threads = get_num_threads(domain::blas_kernel);
};
concurrent_scope::~concurrent_scope()
{
    set_num_threads(m_old_kernel_threads, domain::blas_kernel);
};

};