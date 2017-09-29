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

#include "allocator.h"
#include "utils/impl_types.h"
#include "utils/utils.h"
#include "matcl-core/general/thread.h"
#include <vector>
#include <map>
#include <unordered_map>
#include <iostream>

namespace matcl { namespace mp { namespace details
{

namespace mmd = matcl::mp::details;

template<class T, int Size>
class bucket_vec
{
    private:
        T       m_data[Size];
        int     m_pos;

    public:
        bucket_vec()
            :m_pos(0) 
        {};
        
        int size() const
        { 
            return m_pos; 
        };

        void free_memory()
        {
            for (int i = 0; i < m_pos; ++i)
                allocator_impl::hard_free(m_data[i]);
        };        

        void make_alloc(mp_float& f, precision p)
        {
            allocator_impl::make_alloc(m_data[--m_pos], f, p);
        };

        void push_back(const T& elem)
        {
            memcpy(m_data + m_pos, elem, sizeof(elem));
            ++m_pos;
        };

        bool is_full() const
        {
            return m_pos == Size;
        }
};

struct hasher
{
    size_t operator()(size_t v) const   { return v; };
};

class allocator_impl
{
    private:
        // cache at most max_size elements in given precision group
        static const size_t max_size   = 50;

        using memory_chunk  = mmd::impl_float;
        using bucket        = bucket_vec<memory_chunk, max_size>;
        using bucket_map    = std::unordered_map<size_t, bucket, hasher>;
        using map_iter      = bucket_map::iterator;
        using map_citer     = bucket_map::const_iterator;

    private:
        bucket_map  m_buckets;

    public:
        allocator_impl();
        ~allocator_impl();

        void        malloc(mp_float& f, precision p);
        void        free(mp_float& f, precision p);

    private:
        map_iter    find_bucket_it(precision p);
        map_iter    add_bucket(precision p);
        void        remove_bucket(map_iter pos, precision p);
        void        free_memory(bucket& b);
        
        void        hard_free(mp_float& f);
        void        hard_malloc(mp_float& f, precision p);

    public:
        static void hard_free(memory_chunk data);
        static void make_alloc(memory_chunk& mem, mp_float& f, precision p);
};

allocator_impl::allocator_impl()
{};

allocator_impl::~allocator_impl()
{
    for (auto& pos : m_buckets)
        free_memory(pos.second);
};

void allocator_impl::free_memory(bucket& b)
{
    b.free_memory();
}

allocator_impl::map_iter allocator_impl::find_bucket_it(precision p)
{
    return m_buckets.find(p);
};

allocator_impl::map_iter allocator_impl::add_bucket(precision p)
{
    bucket vec;

    auto pos = m_buckets.insert(bucket_map::value_type(p, std::move(vec)));
    return pos.first;
}
void allocator_impl::remove_bucket(map_iter pos, precision p)
{
    (void)p;
    m_buckets.erase(pos);
}

void allocator_impl::malloc(mp_float& f, precision p)
{
    map_iter pos = find_bucket_it(p);

    if (pos == m_buckets.end())
    {
        add_bucket(p);
        hard_malloc(f, p);
        return;
    };

    size_t size = pos->second.size();

    if (size == 0)
    {
        hard_malloc(f, p);
        return;
    };

    pos->second.make_alloc(f, p);
};

void allocator_impl::free(mp_float& f, precision p)
{
    map_iter pos = find_bucket_it(p);

    //pos != m_buckets.end(); bucket was created when at least
    //one element of precision p was created

    if (pos->second.is_full())
        hard_free(f);
    else
        pos->second.push_back(mmd::impl_value(f));
};

void allocator_impl::hard_free(mp_float& f)
{
    mpfr_clear(mmd::impl_value(f));
};
void allocator_impl::hard_free(memory_chunk data)
{
    mpfr_clear(data);
}
void allocator_impl::hard_malloc(mp_float& f, precision p)
{
    mpfr_init2(mmd::impl_value(f), (long)p);
};
void allocator_impl::make_alloc(memory_chunk& mem, mp_float& f, precision p)
{
    //mpfr does not clear sign flags correctly and do not change 
    //sign when f is initialized later to NaN; just manually correct the sign
    mem[0]._mpfr_sign = 0;

    memcpy(&mmd::impl_value(f), mem, sizeof(memory_chunk));

    (void)p;
    //mpfr_set_prec(mmd::impl_value(f), p);
};

MATCL_THREAD_LOCAL
static allocator_impl* inst = new allocator_impl();

allocator_impl* get_alloc()
{
    return inst;
};

// initial implementation of caching in order to reduce number of malloc/free
// calls; even if all malloc/free are cached and code do mainly allocations/deallocations,
// speedup is of order 2; if some nontrivial calculations are done, speedup is at most
// 20-30%.
#define MATCL_MP_USE_CACHE 0

void mp_float_alloc::init(mp_float& f, precision prec)
{
    #if MATCL_MP_USE_CACHE == 1
        get_alloc()->malloc(f, prec);
    #else
        mpfr_init2(mmd::impl_value(f), (long)prec);
    #endif
};

void mp_float_alloc::free(mp_float& f, precision p)
{
    #if MATCL_MP_USE_CACHE == 1
        get_alloc()->free(f, p);
    #else
        (void)p;
        mpfr_clear(mmd::impl_value(f));
    #endif
}

}}}
