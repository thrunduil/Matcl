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

#include "matcl-simd/config.h"

#include <stdint.h>

namespace matcl { namespace simd { namespace details
{

// integer type of the same sizeof
template<class T> struct integer_type{};
template<>        struct integer_type<float>            { using type = int32_t;};
template<>        struct integer_type<double>           { using type = int64_t;};
template<>        struct integer_type<int32_t>          { using type = int32_t;};
template<>        struct integer_type<int64_t>          { using type = int64_t;};

// unsigned integer type of the same sizeof
template<class T> struct unsigned_integer_type{};
template<>        struct unsigned_integer_type<float>   { using type = uint32_t;};
template<>        struct unsigned_integer_type<double>  { using type = uint64_t;};
template<>        struct unsigned_integer_type<int32_t> { using type = uint32_t;};
template<>        struct unsigned_integer_type<int64_t> { using type = uint64_t;};

}}}
