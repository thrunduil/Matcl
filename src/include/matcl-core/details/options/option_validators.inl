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

#include "matcl-core/options/option_validators.h"
#include "matcl-core/error/exception_classes.h"

namespace matcl { namespace opt
{

//--------------------------------------------------------------------------
//           PREDEFINED OPTION VALIDATORS IMPLEMENTATION
//--------------------------------------------------------------------------

// test whether given option has value > 0
template<class T>
inline auto opt::validator_positive() -> std::function<optional<T> (optional<T>)>
{
    using opt_type = optional<T>;
    return [](opt_type x)->opt_type
            {
                if (x && x.value() <= T(0))
                {
                    std::ostringstream os;
                    os << "value must be positive";
                    throw error::option_validator_error(os.str());
                }
                return x;
            };
};

// test whether given option has value > 0
template<class T>
inline auto opt::validator_positive_or_val(T val) -> std::function<optional<T> (optional<T>)>
{
    using opt_type = optional<T>;
    return [val](opt_type x)->opt_type
            {
                if (x && x.value() != val && x.value() <= T(0))
                {
                    std::ostringstream os;
                    os << "value must be positive";
                    throw error::option_validator_error(os.str());
                }
                return x;
            };
};

// test whether given option has value >= 0
template<class T>
inline auto opt::validator_nonnegative() -> std::function<optional<T> (optional<T>)>
{
    using opt_type = optional<T>;
    return [](opt_type x)->opt_type
            {
                if (x && x.value() < T(0))
                {
                    std::ostringstream os;
                    os << "value must be positive or zero";
                    throw error::option_validator_error(os.str());
                }
                return x;
            };
};

// test whether given option has value >= 0
template<class T>
inline auto opt::validator_nonnegative_or_val(T val) -> std::function<optional<T> (optional<T>)>
{
    using opt_type = optional<T>;
    return [val](opt_type x)->opt_type
            {
                if (x && x.value() != val && x.value() < T(0))
                {
                    std::ostringstream os;
                    os << "value must be positive or zero";
                    throw error::option_validator_error(os.str());
                }
                return x;
            };
};

// test whether given option is in range [min, max], (min,max], [min, max),
// or (min, max)
template<class T>
inline auto opt::validator_range(T min, T max, bool open_left, bool open_right)
                                    -> std::function<optional<T> (optional<T>)>
{
    using opt_type = optional<T>;
    return [open_left, open_right, min, max](opt_type x)->opt_type
            {
                if (!x)
                    return x;

                if (open_left == false && (x.value() >= min) == false)
                {
                    std::ostringstream os;
                    os << "value out of range";
                    throw error::option_validator_error(os.str());
                }
                if (open_left == true && (x.value() > min) == false)
                {
                    std::ostringstream os;
                    os << "value out of range";
                    throw error::option_validator_error(os.str());
                }
                if (open_right == false && (x.value() <= max) == false)
                {
                    std::ostringstream os;
                    os << "value out of range";
                    throw error::option_validator_error(os.str());
                }
                if (open_right == true && (x.value() < max) == false)
                {
                    std::ostringstream os;
                    os << "value out of range";
                    throw error::option_validator_error(os.str());
                }

                return x;
            };
}

// test whether given option is in range [min, max], (min,max], [min, max),
// or (min, max)
template<class T>
inline auto opt::validator_range_or_val(T val, T min, T max, bool open_left, bool open_right)
                                    -> std::function<optional<T> (optional<T>)>
{
    using opt_type = optional<T>;
    return [val, open_left, open_right, min, max](opt_type x)->opt_type
            {
                if (!x)
                    return x;

                if (x.value() == val)
                    return x;

                if (open_left == false && (x.value() >= min) == false)
                {
                    std::ostringstream os;
                    os << "value out of range";
                    throw error::option_validator_error(os.str());
                }
                if (open_left == true && (x.value() > min) == false)
                {
                    std::ostringstream os;
                    os << "value out of range";
                    throw error::option_validator_error(os.str());
                }
                if (open_right == false && (x.value() <= max) == false)
                {
                    std::ostringstream os;
                    os << "value out of range";
                    throw error::option_validator_error(os.str());
                }
                if (open_right == true && (x.value() < max) == false)
                {
                    std::ostringstream os;
                    os << "value out of range";
                    throw error::option_validator_error(os.str());
                }

                return x;
            };
}

// test whether given option value is finite
template<class T>
inline auto opt::validator_finite() -> std::function<optional<T> (optional<T>)>
{
    using opt_type = optional<T>;
    return [](opt_type x)->opt_type
            {
                if (x && is_finite(x.value()) == false)
                {
                    std::ostringstream os;
                    os << "value must be finite";
                    throw error::option_validator_error(os.str());
                }
                return x;
            };
}

// test whether given enum value represented as integer
// is in range [0, size)
inline auto opt::validator_enum(Integer size) -> std::function<optional<Integer> (optional<Integer>)>
{
    using opt_type = optional<Integer>;
    return [size](opt_type x)->opt_type
            {
                if (x && (x.value() < 0 || x.value() >= size))
                {
                    std::ostringstream os;
                    os << "enumeration value is out out range";
                    throw error::option_validator_error(os.str());
                }
                return x;
            };
};

};}
