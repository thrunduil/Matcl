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
#include "matcl-core/details/mpl.h"

#include <iostream>

namespace matcl
{

//-----------------------------------------------------------------------
//                          GENERAL INFO
//-----------------------------------------------------------------------

// Twofold is a representation of a double precition floating point number
// z as z = value + error, where |error| <= ulp(|value|)/2 and fl(value + error)
//  = value. 
// Functions defined for twofold type delivers nearly two times higher 
// accuracy, then corresponding functions defined for double type. 
//    
// Twofold value is not rounded according to quadruple precision, for example
// a twofold number 1 + eps^3 is not rounded to 1.
//
// Twofold does not handle INF/NAN values. Results are meaningful only when 
// they are finite. When any twofold's function is called with nonfinite 
// arguments, then always nonfinite result is returned, even when correct result
// is finite (e.g. twofold_div(1, INF) does not return 0).
//
// If a function defined for twofold numbers return a twofold number 
// x = value + error with true result is res, then if value + error = res 
// (evaluated exactly), then the result is called exact. For all functions the
// following error bound is satified:
//      |res - value - error|/|res| <= k * u^2
// for some small k (at most 16), where u is the unit roundoff (u = 2^-53 for 
// double precision)

// representation of a floating point number as value + error
template<class Float_type>
class twofold
{    
    public:
        // type of underlying floating point type
        using float_type    = Float_type;

        struct uninitialized{};

    public:
        Float_type  value;
        Float_type  error;

    public:
        // set value and error to 0.0
        twofold();

        // create uninitialized twofold value
        twofold(uninitialized);

        // set value to val and error to 0.0
        explicit twofold(const Float_type& val);

        // set value to val and error to err; val and err must be normalized
        // i.e. |err| <= ulp(|val|)/2 and fl(val + err) = val.
        twofold(const Float_type& val, const Float_type &err);

    public:
        // construct normalized twofold number representing val + err; it is 
        // required, that |err| <= |val|
        static twofold      normalize_fast(const Float_type& val, const Float_type& err);

        // construct normalized twofold number representing val + err;
        // val and err can be any numbers; this is the slowes normalization
        // function
        static twofold      normalize(const Float_type& val, const Float_type& err);

        // return value + error (with is equal to value due to normalization
        // requirement)
        const Float_type&   sum() const;

        // return true if both value and error are finite (e.g. not INF and not NAN)
        // only finite values are meaningful (see GENERAL INFO for details)
        bool                is_finite() const;
};

//-----------------------------------------------------------------------
//                      ARITHMETIC FUNCTIONS
//-----------------------------------------------------------------------

// evaluate a + b for arbitrary a and b; result is exact
template<class Float_type>
twofold<Float_type> twofold_sum(const Float_type& a, const Float_type& b);

// evaluate a + b, provided |a| >= |b|; this function is faster
// than twofold_sum; result is exact
template<class Float_type>
twofold<Float_type> twofold_sum_sorted(const Float_type& a, const Float_type& b);

// evaluate a + b, provided |a| >= |b|; this function is faster operator+
// relative forward error does not exceed 2 * u^2
template<class Float_type>
twofold<Float_type> twofold_sum_sorted(const twofold<Float_type>& a, const Float_type& b);

// evaluate a + b, provided |a| >= |b|; this function is faster operator+
// relative forward error does not exceed 2 * u^2
template<class Float_type>
twofold<Float_type> twofold_sum_sorted(const Float_type& a, const twofold<Float_type>& b);

// evaluate a - b for arbitrary a and b; result is exact
template<class Float_type>
twofold<Float_type> twofold_minus(const Float_type& a, const Float_type& b);

// evaluate a - b, provided |a| >= |b|; this function is faster
// than twofold_minus; result is exact
template<class Float_type>
twofold<Float_type> twofold_minus_sorted(const Float_type& a, const Float_type& b);

// evaluate a - b, provided |a| >= |b|; this function is faster operator+
// relative forward error does not exceed 2 * u^2
template<class Float_type>
twofold<Float_type> twofold_minus_sorted(const twofold<Float_type>& a, const Float_type& b);

// evaluate a - b, provided |a| >= |b|; this function is faster operator+
// relative forward error does not exceed 2 * u^2
template<class Float_type>
twofold<Float_type> twofold_minus_sorted(const Float_type& a, const twofold<Float_type>& b);

// evaluate a * b; result is exact
template<class Float_type>
twofold<Float_type> twofold_mult(const Float_type& a, const Float_type& b);

// evaluate a / b; relative forward error does not exceed 1 * u^2
template<class Float_type>
twofold<Float_type> twofold_div(const Float_type& a, const Float_type& b);

// square root of a; relative forward error does not exceed 1 * u^2
template<class Float_type>
twofold<Float_type> twofold_sqrt(const Float_type& a);

// evaluate a + b; relative forward error does not exceed 3 * u^2
template<class Float_type>
twofold<Float_type> operator+(const twofold<Float_type>& a, const twofold<Float_type>& b);

// evaluate a + b; relative forward error does not exceed 2 * u^2
template<class Float_type>
twofold<Float_type> operator+(const twofold<Float_type>& a, const Float_type& b);

// evaluate a + b; relative forward error does not exceed 2 * u^2
template<class Float_type>
twofold<Float_type> operator+(const Float_type& a, const twofold<Float_type>& b);

// evaluate a - b; relative forward error does not exceed 3 * u^2
template<class Float_type>
twofold<Float_type> operator-(const twofold<Float_type>& a, const twofold<Float_type>& b);

// evaluate a - b; relative forward error does not exceed 2 * u^2
template<class Float_type>
twofold<Float_type> operator-(const twofold<Float_type>& a, const Float_type& b);

// evaluate a - b; relative forward error does not exceed 2 * u^2
template<class Float_type>
twofold<Float_type> operator-(const Float_type& a, const twofold<Float_type>& b);

// evaluate -a; result is exact
template<class Float_type>
twofold<Float_type> operator-(const twofold<Float_type>& a);

// evaluate a * b; relative forward error does not exceed 7 * u^2
// (6 u^2 if FMA instruction is available)
template<class Float_type>
twofold<Float_type> operator*(const twofold<Float_type>& a, const twofold<Float_type>& b);

// evaluate a * b;  relative forward error does not exceed 3 * u^2
// (2 u^2 if FMA instruction is available)
template<class Float_type>
twofold<Float_type> operator*(const twofold<Float_type>& a, const Float_type& b);

// evaluate a * b; relative forward error does not exceed 3 * u^2
// (2 u^2 if FMA instruction is available)
template<class Float_type>
twofold<Float_type> operator*(const Float_type& a, const twofold<Float_type>& b);

// evaluate a / b; relative forward error does not exceed 15 * u^2
template<class Float_type>
twofold<Float_type> operator/(const twofold<Float_type>& a, const twofold<Float_type>& b);

// evaluate a / b; relative forward error does not exceed 4 * u^2
template<class Float_type>
twofold<Float_type> operator/(const twofold<Float_type>& a, const Float_type& b);

// evaluate a / b; relative forward error does not exceed 15 * u^2
template<class Float_type>
twofold<Float_type> operator/(const Float_type& a, const twofold<Float_type>& b);

// square root of a; relative forward error does not exceed 4 * u^2
template<class Float_type>
twofold<Float_type> sqrt(const twofold<Float_type>& a);

// absolute value of a; result is exact
template<class Float_type>
twofold<Float_type> abs(const twofold<Float_type>& a);

//-----------------------------------------------------------------------
//                      ERROR RELATED FUNCTIONS
//-----------------------------------------------------------------------
// return number of distinct properly rounded representations between x and
// y; if one of arguments is NaN, then NaN is returned;
// the result is always a signed integer value stored in double type
template<class Float_type>
Float_type  float_distance(const twofold<Float_type>& x, const twofold<Float_type>& y);

// return epsilon value eps, i.e positive distance from abs(xr) to the next
// larger in magnitude properly rounded twofold number, where xr is a properly
// rounded x; rounding is (implicitly) performed according to quadruple precision
// (i.e. for double, precision is 106 bits and eps(1) = 2^(1-106) )
template<class Float_type>
Float_type  eps(const twofold<Float_type>& x);

//-----------------------------------------------------------------------
//                      IO FUNCTIONS
//-----------------------------------------------------------------------
// save to stream
template<class Float_type>
std::ostream& operator<<(std::ostream& os, const twofold<Float_type>& x);

// load from stream
template<class Float_type>
std::istream& operator>>(std::istream& os, twofold<Float_type>& x);

}

#include "matcl-core/details/float/twofold.inl"
