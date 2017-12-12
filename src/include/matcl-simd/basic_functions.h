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

#include "matcl-simd/simd_general.h"

namespace matcl { namespace simd
{

// value representing true (all bits set to 1)
template<class T>
struct true_value
{
    static T get();
};

// value representing false (all bits set to 0)
template<class T>
struct false_value
{
    static T get();
};

//-----------------------------------------------------------------------
//                   ARITMETHIC FUNCTIONS
//-----------------------------------------------------------------------

// vector multiply x * y
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
operator*(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector division x / y
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
operator/(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector add x + y
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
operator+(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector subtract x - y
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
operator-(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector unary minus x
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
operator-(const simd<Val, Bits, Simd_tag>& x);

// vector multiply x * y; equivalent to operator*
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
mult(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector division x / y; equivalent to operator/
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
div(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector add x + y; equivalent to operator+
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
plus(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector subtract x - y; equivalent to operator-
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
minus(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector unary minus x; equivalent to unary operator-
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
uminus(const simd<Val, Bits, Simd_tag>& x);

// alternatively subtract and add elements in x and y
// i.e. form [x[0] - y[0], x[1] + y[1], x[2] - y[2], x[3] + y[3], ...]
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
sub_add(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

//res = x*y + z; use FMA instruction if available, otherwise
// evaluate according to definition (with two roundings)
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
fma_f(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z);

//res = x*y - z; use FMA instruction if available, otherwise
// evaluate according to definition (with two roundings)
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
fms_f(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z);

//res = -x*y + z; use FMA instruction if available, otherwise
// evaluate according to definition (with two roundings)
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
fnma_f(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z);

//res = -x*y - z; use FMA instruction if available, otherwise
// evaluate according to definition (with two roundings)
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
fnms_f(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z);

//res = x*y + z; use FMA instruction if available, otherwise
// this function is evaluated as using (slow) Dekker's algorithm
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
fma_a(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z);

//res = x*y - z; use FMA instruction if available, otherwise
// this function is evaluated as using (slow) Dekker's algorithm
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
fms_a(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z);

//res = -x*y + z; use FMA instruction if available, otherwise
// this function is evaluated as using (slow) Dekker's algorithm
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
fnma_a(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z);

//res = -x*y - z; use FMA instruction if available, otherwise
// this function is evaluated as using (slow) Dekker's algorithm
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
fnms_a(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z);

// sum of all elements stored in the vector x
template<class Val, int Bits, class Simd_tag>
Val sum_all(const simd<Val, Bits, Simd_tag>& x);

// absolute value
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
abs(const simd<Val, Bits, Simd_tag>& x);

// square root
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
sqrt(const simd<Val, Bits, Simd_tag>& x);

//-----------------------------------------------------------------------
//                   COMPARISON FUNCTIONS
//-----------------------------------------------------------------------

// element by element maximum
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
max(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// element by element minumum
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
min(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// test x == y; return a vector of floating point numbers containing
// true_value or false_value
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
eeq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// test x != y; return a vector of floating point numbers containing
// true_value or false_value
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
neq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// test x < y; return a vector of floating point numbers containing
// true_value or false_value
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
lt(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// test x == y; return a vector of floating point numbers containing
// true_value or false_value
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
gt(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// test x <= y; return a vector of floating point numbers containing
// true_value or false_value
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
leq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// test x >= y; return a vector of floating point numbers containing
// true_value or false_value
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
geq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// test x for NaN values; return a vector of floating point numbers containing
// true_value or false_value
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
is_nan(const simd<Val, Bits, Simd_tag>& x);

//-----------------------------------------------------------------------
//                   ROUNDING FUNCTIONS
//-----------------------------------------------------------------------
// round toward nearest integer with ties to even
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
round(const simd<Val, Bits, Simd_tag>& x);

// round toward -INF
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
floor(const simd<Val, Bits, Simd_tag>& x);

// round toward +INF
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
ceil(const simd<Val, Bits, Simd_tag>& x);

// round toward zero
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
trunc(const simd<Val, Bits, Simd_tag>& x);

//-----------------------------------------------------------------------
//                   REDUCE FUNCTIONS
//-----------------------------------------------------------------------
// return true if at least element in the vector x is NAN
template<class Val, int Bits, class Simd_tag>
bool any_nan(const simd<Val, Bits, Simd_tag>& x);

// return true if all elements in the vector x are equal to true_value,
// where x contains only true_value or false value
template<class Val, int Bits, class Simd_tag>
bool all(const simd<Val, Bits, Simd_tag>& x);

// return true if at least one element in the vector x are equal to true_value,
// where x contains only true_value or false value
template<class Val, int Bits, class Simd_tag>
bool any(const simd<Val, Bits, Simd_tag>& x);

//-----------------------------------------------------------------------
//                   BITWISE MANIPULATION
//-----------------------------------------------------------------------
// perform bitwise or operation for all elements in vectors x and y
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
bitwise_or(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// perform bitwise and operation for all elements in vectors x and y
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
bitwise_and(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// perform bitwise xor operation for all elements in vectors x and y
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
bitwise_xor(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// perform not(x[i]) and y[i] for all elements in vectors x and y
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
bitwise_andnot(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// perform not(x[i]) for all elements in a vector x
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
bitwise_not(const simd<Val, Bits, Simd_tag>& x);

// shift left all elements; floating point type is implicitely casted to unsigned
// integer type of the same tyme, then shift is performed and result is casted to
// floating point type
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
shift_left(const simd<Val, Bits, Simd_tag>& x, unsigned int count);

// (logical) shift right all elements; floating point type is implicitely casted
// to unsigned integer type of the same tyme, then shift is performed and result
// is casted to floating point type
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
shift_right(const simd<Val, Bits, Simd_tag>& x, unsigned int count);

// arithmetic shift right all elements; floating point type is implicitely casted
// to signed integer type of the same tyme, then shift is performed and result
// is casted to floating point type
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
shift_right_arithmetic(const simd<Val, Bits, Simd_tag>& x, unsigned int count);

// test whether the sign of x is negative; return a vector of floating 
// point numbers containing -0.0 (when sign bit is set) or 0.0 (otherwise)
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
signbit_base(const simd<Val, Bits, Simd_tag>& x);

//-----------------------------------------------------------------------
//                   CONDITIONAL FUNCTIONS
//-----------------------------------------------------------------------
// evaluate test ? NaN : val_false;
// test must be a vector of floating point numbers containing true_value or false_value
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
if_nan_else(const simd<Val, Bits, Simd_tag>& test, const simd<Val, Bits, Simd_tag>& val_false);

// evaluate test ? 0.0 : val_false;
// test must be a vector of floating point numbers containing true_value or false_value
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
if_zero_else(const simd<Val, Bits, Simd_tag>& test, const simd<Val, Bits, Simd_tag>& val_false);

// evaluate test ? val_true : val_false;
// test must be a vector of floating point numbers containing true_value or false_value
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
if_then_else(const simd<Val, Bits, Simd_tag>& test, const simd<Val, Bits, Simd_tag>& val_true,
             const simd<Val, Bits, Simd_tag>& val_false);

//-----------------------------------------------------------------------
//                   MISCELLANEOUS FUNCTIONS
//-----------------------------------------------------------------------

// vector of elements in reverse order
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
reverse(const simd<Val, Bits, Simd_tag>& x);

//-----------------------------------------------------------------------
//                   IO FUNCTIONS
//-----------------------------------------------------------------------

// save to stream
template<class Val, int Bits, class Simd_tag>
std::ostream& operator<<(std::ostream& os, const simd<Val, Bits, Simd_tag>& x);

// load from stream
template<class Val, int Bits, class Simd_tag>
std::istream& operator>>(std::istream& is, simd<Val, Bits, Simd_tag>& x);

}}
