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
// (evaluated exactly), then the result is called exact. If error = fl(res - value),
// thenthe result is called strict. If the result x is neither strict nor exact,
// then |res - value - error| <= k * ulp(|res|) for some small k (at most 4),
// where ulp(|res|) is the unit in the last place assuming precision two times
// higher than double precision (i.e. 105 bits, when double precision is 53 
// bits)

// representation of a floating point number as value + error
class twofold
{    
    public:
        double  value;
        double  error;

    public:
        // set value and error to 0.0
        twofold();

        // set value to val and error to 0.0
        explicit twofold(const double& val);

        // set value to val and error to err; val and err must be normalized
        // i.e. |err| <= ulp(|val|)/2 and fl(val + err) = val.
        twofold(const double& val, const double &err);

        // construct normalized twofold number representing val + err; it is 
        // required, that |err| <= |val|
        static twofold  normalize_fast(const double& val, const double& err);

        // construct normalized twofold number representing val + err;
        // val and err can be any numbers; this is the slowes normalization
        // function
        static twofold  normalize(const double& val, const double& err);

        // return value + error (with is equal to value due to normalization
        // requirement)
        const double&   sum() const;

        // return true if both value and error are finite (e.g. not INF and not NAN)
        // only finite values are meaningful (see GENERAL INFO for details)
        bool            is_finite() const;
};

//-----------------------------------------------------------------------
//                      ARITHMETIC FUNCTIONS
//-----------------------------------------------------------------------

// evaluate a + b for arbitrary a and b; result is exact
twofold twofold_plus(double a, double b);

// evaluate a + b, provided |a| >= |b|; this function is faster
// than twofold_plus; result is exact
twofold twofold_plus_sorted(double a, double b);

// evaluate a - b for arbitrary a and b; result is exact
twofold twofold_minus(double a, double b);

// evaluate a - b, provided |a| >= |b|; this function is faster
// than twofold_minus; result is exact
twofold twofold_minus_sorted(double a, double b);

// evaluate a * b; result is exact
twofold twofold_mult(double a, double b);

// evaluate a / b; result is strict
twofold twofold_div(double a, double b);

// square root of a; result is strict
twofold twofold_sqrt(double a);

// evaluate a + b; result need not be strict
twofold operator+(const twofold& a, const twofold& b);

// evaluate a + b; result is strict
twofold operator+(const twofold& a, double b);

// evaluate a + b; result is strict
twofold operator+(double a, const twofold& b);

// evaluate a - b; result need not be strict
twofold operator-(const twofold& a, const twofold& b);

// evaluate a - b; result is strict
twofold operator-(const twofold& a, double b);

// evaluate a - b; result is strict
twofold operator-(double a, const twofold& b);

// evaluate -a; result is exact
twofold operator-(const twofold& a);

// evaluate a * b; result need not be strict
twofold operator*(const twofold& a, const twofold& b);

// evaluate a * b; result is strict
twofold operator*(const twofold& a, double b);

// evaluate a * b; result is strict
twofold operator*(double a, const twofold& b);

// evaluate a / b; result need not be strict
twofold operator/(const twofold& a, const twofold& b);

// evaluate a / b
twofold operator/(const twofold& a, double b);

// evaluate a / b; result need not be strict, but maximum relative
// forward error does not exceed 4 ulp (ulp calculated as if precision
// was doubled)
twofold operator/(double a, const twofold& b);

// square root of a; result need not be strict
twofold sqrt(const twofold& a);

// absolute value of a; result is exact
twofold abs(const twofold& a);

//-----------------------------------------------------------------------
//                      ERROR RELATED FUNCTIONS
//-----------------------------------------------------------------------
// return number of distinct properly rounded representations between x and
// y; if one of arguments is NaN, then NaN is returned;
// the result is always a signed integer value stored in double type
MATCL_CORE_EXPORT
double  float_distance(const twofold& x, const twofold& y);

// return epsilon value eps, i.e positive distance from abs(xr) to the next
// larger in magnitude properly rounded twofold number, where xr is a properly
// rounded x; rounding is (implicitly) performed according to quadruple precision
MATCL_CORE_EXPORT
double  eps(const twofold& x);

//-----------------------------------------------------------------------
//                      IO FUNCTIONS
//-----------------------------------------------------------------------
// save to stream
MATCL_CORE_EXPORT
std::ostream& operator<<(std::ostream& os, const twofold& x);

// load from stream
MATCL_CORE_EXPORT
std::istream& operator>>(std::istream& os, twofold& x);

//-----------------------------------------------------------------------
//                   COMPATIBILITY FUNCTIONS
//-----------------------------------------------------------------------
// evaluate a * b using the Veltkamp/Dekker algorithm; result is exact; 
// this function is equivalent to twofold_mult and is called by twofold_mult
// and other functions, when FMA instruction is not available
twofold twofold_mult_dekker(double a, double b);

// evaluate x * y + z with one rounding using the Veltkamp/Dekker algorithm
// for double precision arguments and double arithmetics for single precision
// arguments; 
// this function is equivalent to FMA instruction and is called by twofold's
// function, when FMA instruction is not available
double  fma_dekker(double x, double y, double z);
float   fma_dekker(float x, float y, float z);

}

#include "matcl-core/details/float/twofold.inl"
