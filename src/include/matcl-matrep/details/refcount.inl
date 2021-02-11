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

namespace matcl { namespace details
{

template<class refcount_pool>
inline long refcount_str<refcount_pool>::get_count() const       
{ 
    long count = m_refcount.load( std::memory_order_acquire );
    return count;
};

template<class refcount_pool>
inline refcount_str<refcount_pool>* refcount_str<refcount_pool>::increase()              
{
    m_refcount.fetch_add( 1, std::memory_order_acq_rel );
    return this; 
};

template<class refcount_pool>
inline bool refcount_str<refcount_pool>::is_unique() const
{
    return get_count() == 1;
};

template<class refcount_pool>
inline bool refcount_str<refcount_pool>::decrease()
{
    long count  = m_refcount.fetch_sub( 1, std::memory_order_acq_rel ) - 1;
    return count == 0;
};

template<class refcount_pool>
inline refcount_str<refcount_pool>::refcount_str(long val)
:m_refcount(val)
{};

template<class refcount_pool>
class MATCL_MATREP_EXPORT refstruct_pool
{
    public:
        using refcount_type = refcount_str<refcount_pool>;

    private:
        refstruct_pool();
        refstruct_pool(const refstruct_pool&)               = delete;
        refstruct_pool& operator=(const refstruct_pool&)    = delete;

    public:
        static refcount_type*   create(long val);
        static void             destroy(refcount_type* ptr);
};

template<class refcount_pool>
inline refcount_str<refcount_pool>* refcount_str<refcount_pool>::create(long val)
{
    return refstruct_pool<refcount_pool>::create(val);
};

template<class refcount_pool>
inline void refcount_str<refcount_pool>::destroy()
{
    refstruct_pool<refcount_pool>::destroy(this);
};

};};
