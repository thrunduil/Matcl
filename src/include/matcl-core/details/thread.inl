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

#include "matcl-core/general/thread.h"

namespace matcl
{

// implements interface of atomics in single threaded version of matcl
template<class V>
inline atomic_nothread<V>::atomic_nothread()
    :m_data(V(0))
{};

template<class V>
inline atomic_nothread<V>::atomic_nothread(V v)
    :m_data(v)
{};

template<class V>
inline V atomic_nothread<V>::load(std::memory_order) const
{ 
    return m_data; 
};

template<class V>
inline void atomic_nothread<V>::store(V val, std::memory_order)
{ 
    m_data = val; 
};

template<class V>
inline atomic_nothread<V>::operator V() const
{ 
    return m_data; 
};

template<class V>
inline atomic_nothread<V>& atomic_nothread<V>::operator++()
{ 
    ++m_data; 
    return *this; 
};

template<class V>
inline atomic_nothread<V>& atomic_nothread<V>::operator--()
{ 
    --m_data; 
    return *this; 
};

template<class V>
inline atomic_nothread<V>& atomic_nothread<V>::operator=(V val)
{
    m_data = val; 
    return *this;
}

template<class V>
template<class T>
inline V atomic_nothread<V>::fetch_add(T inc, std::memory_order)
{ 
    V old   = m_data; 
    m_data  += inc; 
    return old; 
};

template<class V>
template<class T>
inline V atomic_nothread<V>::fetch_sub(T inc, std::memory_order)
{ 
    V old   = m_data; 
    m_data  -= inc; 
    return old; 
};

};
