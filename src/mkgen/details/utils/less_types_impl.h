/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019 - 2021
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

#include "mkgen/mkgen_fwd.h"

namespace matcl { namespace mkgen { namespace details { namespace less_impl
{

//----------------------------------------------------------------------------------
//                              less_impl
//----------------------------------------------------------------------------------

// length of null-terminated string
constexpr size_t cstrlen(const char* p) 
{
    size_t len = 0;
    while (*p) 
    {
        ++len;
        ++p;
    }
    
    return len;
}

template<class T>
struct type_name
{
    static constexpr char* name() 
    { 
        #if defined (_MSC_VER)
            return __FUNCSIG__; 
        #else
            return __PRETTY_FUNCTION__;
        #endif
    };
};

template<class T1, class T2>
constexpr bool less_impl()
{
    const char* A   = type_name<T1>::name();
    const char* B   = type_name<T2>::name();

    size_t a_len    = cstrlen(A);
    size_t b_len    = cstrlen(B);

    size_t ab_len   = (a_len < b_len) ? a_len : b_len;

    for (size_t i = 0; i < ab_len; ++i) 
    {
        if (A[i] != B[i]) 
            return A[i] < B[i];
    }

    if (a_len < b_len)
        return true;

    if (a_len > b_len)
        return false;

    return false;
}

}}}}

