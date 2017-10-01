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

#include "matcl-core/config.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/lib_functions/constants.h"

namespace matcl 
{

using matcl_int64   = int64_t;           
using matcl_uint64  = uint64_t;            
using matcl_int32   = int32_t;           
using matcl_uint32  = uint32_t;

inline Integer cast_int64(matcl_int64 n)		    { return static_cast<Integer>(n); };
inline Integer cast_int64(matcl_int32 n)		    { return n; };

inline Integer cast_int64(matcl_uint32 n)	        { return n; };
inline Integer cast_int64(matcl_uint64 n)	        { return static_cast<Integer>(n); };


inline Integer imult(Integer a, Integer b)          { return a*b;         };
inline Integer icast(Real a)                        { return (Integer)a;  };

//safe mult
inline Integer imult_s(Integer r, Integer c);

//test if mult can be performed
inline bool imult_t(Integer r, Integer c, Integer& ret);
inline bool imult_size_t_test(size_t r, size_t c, size_t& ret);

//throw on overflow
inline Integer imult_c(Integer a, Integer b);

//throw on overflow
inline size_t imult_size_t(size_t a, size_t b);

//throw on overflow
Integer icast_c(Real a);

//test if cast can be performed
bool icast_t(Real a);

inline bool same_sign(Integer x, Integer y)
{
    return (x^y) >= 0;
};

inline void integer_div(Integer a, Integer b, Integer& d, Integer& r)
{
    d = a/b;
    r = a%b;
};

//test if mult can be performed
inline bool imult_t(Integer r, Integer c, Integer& ret)
{
    ret = r * c;

    if (r != 0 && ret/r != c)
        return false;
    else
        return true;
};
inline bool imult_size_t_test(size_t r, size_t c, size_t& ret)
{
    ret = r * c;

    if (r != 0 && ret/r != c)
        return false;
    else
        return true;
};

//safe mult
inline Integer imult_s(Integer r, Integer c)
{
    Integer ret;

    if (imult_t(r,c,ret) == false)
    {
        if (r > 0 && c > 0 || r < 0 && c < 0)
            return constants::max_int();
        else 
            return constants::min_int();
    }
    else
    {
        return ret;
    };
};

MATCL_CORE_EXPORT Integer throw_overflow_int_mult(Integer a, Integer b);
MATCL_CORE_EXPORT size_t  throw_overflow_sizet_mult(size_t a, size_t b);
MATCL_CORE_EXPORT Integer throw_overflow_int_cast(Real a);

inline Integer imult_c(Integer a, Integer b)
{
    Integer ret;
    if (imult_t(a,b,ret))
        return ret;
    else
        return throw_overflow_int_mult(a, b);
};

inline size_t imult_size_t(size_t a, size_t b)
{
    size_t ret;
    if (imult_size_t_test(a,b,ret))
        return ret;
    else
        return throw_overflow_sizet_mult(a, b);
}

//test if cast can be performed
inline bool icast_t(Real a)
{
    if (a > Real(constants::max_int()) || a < Real(constants::min_int()))
        return false;
    else
        return true;
};

inline Integer icast_c(Real a)
{
    if (icast_t(a) == true) 
        return (Integer)a;
    else
        return throw_overflow_int_cast(a);
};

};
