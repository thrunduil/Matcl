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

#include "matcl-internals/algs/scatter.h"
#include "matcl-core/memory/memory.h"

namespace matcl { namespace algorithm 
{

//-----------------------------------------------------------------------------
//                      scatter_provider
//-----------------------------------------------------------------------------

//helper class to register scatter_provider as cache
struct scatter_cache : cache_registerer
{
    virtual void        release_cache() override;
};

//reuse already created scatters; the largest scatter
//is saved for later use
class scatter_provider
{
    private:
        using reg_ptr   = std::shared_ptr<scatter_cache>;

    private:
        scatter     m_scatter;
        scatter     m_zero;
        reg_ptr     m_registerer;

    public:
        scatter_provider();
        ~scatter_provider();

        //get scater of given size prepared for calling scatter::next_mark
        //at most n_increase times
        scatter     get(Integer size, Integer n_increase);

        //set content to zero
        void        clear_scatter(scatter& sc);

        //given scatter can be destroyed
        void        release(scatter& sc);

        //release unused scatter
        void        release_cache();

    private:
        scatter     new_scatter(Integer K);
};

scatter_provider::scatter_provider()
{
    using atomic_int   = matcl::atomic<int>;

    static atomic_int registered = 0;

    if (registered.fetch_add(1) == 0)
    {
        //avoid registering this cache many times
        m_registerer = reg_ptr(new scatter_cache());
        register_cache(m_registerer);
    };
};

scatter_provider::~scatter_provider()
{};

scatter scatter_provider::get(Integer size, Integer n_increase)
{
    if (m_scatter.m_size < size)
    {
        return new_scatter(size);
    }
    else if (m_scatter.m_mark + n_increase < 0)
    {
        //overflow protection
        m_scatter = m_zero;
        return new_scatter(size);
    };

    scatter ret = m_scatter;
    m_scatter = m_zero;
    return ret;
};

void scatter_provider::clear_scatter(scatter& sc)
{
    if (sc.m_ptr == nullptr)
        return;

    ::memset(sc.m_ptr,0,sc.m_size * sizeof(Integer));
    sc.m_mark       = 0;
    sc.m_tmp_mark   = 0;
};

scatter scatter_provider::new_scatter(Integer K)
{
    scatter sc;

    sc.m_ptr    = new Integer[K];
    sc.m_size   = K;
    clear_scatter(sc);

    return sc;
};

void scatter_provider::release_cache()
{
    m_scatter = m_zero;
}

void scatter_provider::release(scatter& sc)
{
    //hold the scatter for later use or destroy
    //only one scatter is saved

    if (m_scatter.m_size < sc.m_size)
    {
        m_scatter = sc;
    }
    else
    {
        delete[] sc.m_ptr;

        sc.m_size   = 0;
        sc.m_ptr    = nullptr;
        sc.m_refcount->destroy();

        return;
    };
};

static MATCL_THREAD_LOCAL
scatter_provider* m_scatter_provider = nullptr;

static scatter_provider* get_scatter_provider()
{
    if (!m_scatter_provider)
        m_scatter_provider = new scatter_provider();

    return m_scatter_provider;
};

void scatter_cache::release_cache()
{
    //release only unused memory in this thread
    if (m_scatter_provider)
        m_scatter_provider->release_cache();
};

//-----------------------------------------------------------------------------
//                      scatter
//-----------------------------------------------------------------------------
scatter scatter::get(Integer size, Integer n_increase)
{
    return get_scatter_provider()->get(size, n_increase);
};

scatter::scatter(const scatter& other)
    :m_mark(other.m_mark), m_tmp_mark(other.m_tmp_mark), m_ptr(other.m_ptr), m_size(other.m_size)
{
    m_refcount  = other.m_refcount->increase();
};

scatter& scatter::operator=(const scatter& other)
{
    other.m_refcount->increase();

    if (m_refcount->decrease())
        get_scatter_provider()->release(*this);

    m_mark      = other.m_mark;
    m_tmp_mark  = other.m_tmp_mark;
    m_ptr       = other.m_ptr;
    m_size      = other.m_size;
    m_refcount  = other.m_refcount;

    return *this;
};

scatter::~scatter()
{
    m_mark      = std::max(m_mark, m_tmp_mark);

    if (m_refcount->decrease())
        get_scatter_provider()->release(*this);
};

scatter::scatter()
{
    m_mark      = 0;
    m_tmp_mark  = 0;
    m_ptr       = nullptr;
    m_size      = 0;
    m_refcount  = refcount_str::create(1);
};

};};