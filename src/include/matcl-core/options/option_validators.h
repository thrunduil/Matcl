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

#include "matcl-core/config.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/matrix/complex_type.h"
#include "matcl-core/details/options/matcl_options_details.h"

#include <functional>

namespace matcl { namespace opt
{

//--------------------------------------------------------------------------
//           PREDEFINED OPTION VALIDATORS
//--------------------------------------------------------------------------

// test whether given option has value > 0
template<class T>
auto validator_positive()               -> std::function<optional<T> (optional<T>)>;

// test whether given option has value > 0 or given value val
template<class T>
auto validator_positive_or_val(T val)   -> std::function<optional<T> (optional<T>)>;

// test whether given option has value >= 0
template<class T>
auto validator_nonnegative()            -> std::function<optional<T> (optional<T>)>;

// test whether given option has value >= 0 or given value val
template<class T>
auto validator_nonnegative_or_val(T val)-> std::function<optional<T> (optional<T>)>;

// test whether given option is in range [min, max], (min,max], [min, max),
// or (min, max)
template<class T>
auto validator_range(T min, T max, bool open_left = false, bool open_right = false)
                                        -> std::function<optional<T> (optional<T>)>;

// test whether given option is in range [min, max], (min,max], [min, max),
// or (min, max), or given value
template<class T>
auto validator_range_or_val(T val, T min, T max, bool open_left = false, bool open_right = false)
                                        -> std::function<optional<T> (optional<T>)>;

// test whether given option value is finite
template<class T>
auto validator_finite()             -> std::function<optional<T> (optional<T>)>;

// test whether given enum value represented as integer
// is in range [0, size)
auto validator_enum(Integer size)   -> std::function<optional<Integer> (optional<Integer>)>;

};}

#include "matcl-core/details/options/option_validators.inl"