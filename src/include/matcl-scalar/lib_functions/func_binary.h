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

#include "matcl-scalar/config.h"
#include "matcl-scalar/object.h"
#include "matcl-scalar/details/matfunc_helpers.h"

//TODO: simd version of min, max, cmp

namespace matcl
{

namespace md    = matcl::details;
namespace mrd   = matcl::raw::details;

//---------------------------------------------------------------------------
//                      ARITHMETIC FUNCTIONS
//---------------------------------------------------------------------------

// addition x + y
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types_promote<S1,S2>::type
                            plus(const S1& x, const S2& y);

// addition x + y; different name of the plus function
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types_promote<S1,S2>::type
                            operator+(const S1& x, const S2& y)    { return plus(x,y); };

// subtraction x - y
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types_promote<S1,S2>::type
                            minus(const S1& x, const S2& y);

// subtraction x - y; different name of the minus function
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types_promote<S1,S2>::type
                            operator-(const S1& x, const S2& y)    { return minus(x,y); };

// multiplication x * y
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types_promote<S1,S2>::type
                            mul(const S1& x, const S2& y);

// division x / y; integer division gives Real result; if x and y are
// integer values and y is zero, then result is zero
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            div(const S1& x, const S2& y);

// division x / y; integer division gives Real result; 0/0 division gives 0
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            div_0(const S1& x, const S2& y);

// division x / y; integer division gives Real result; 0/0 division gives 1
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            div_1(const S1& x, const S2& y);

// division x / y; for build-in types default divition operator is called
// (and integer division gives Integer result), otherwise equivalent to 
// div function
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            operator/(const S1& x, const S2& y)    { return div(x,y); };

// division x / y; if one of argument has floating point type then result
// is equivalent to div function; however if both arguments are integers,
// then integer division is performed; if x and y are integer values and
// y is zero, then result is zero
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types_promote<S1,S2>::type
                            idiv(const S1& x, const S2& y);

// return is x - n * y where n = floor(x / y) if y != 0 and 0 otherwise;
// is mod(x,y) is nonzero, then has the same sign as y; if signs of x and
// y are the same, then mod(x,y) = rem(x,y); not available for complex 
// arguments
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::real_unify_types_promote<S1,S2>::type
                            mod(const S1& x, const S2& y);

// return is x - n * y where n = trunc(x / y) if y != 0 and 0 otherwise
// is rem(x,y) is nonzero, then has the same sign as x; if signs of x and
// y are the same, then mod(x,y) = rem(x,y); not available for complex 
// arguments
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::real_unify_types_promote<S1,S2>::type
                            rem(const S1& x, const S2& y);

// four quadrant arctangent of x and y; -pi <= atan2(x, y) <= pi;
// not available for complex arguments
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::real_unify_types2_promote<S1,S2,Float>::type
                            atan2(const S1& x,const S2& y);

// calculate sqrt(|x|^2 + |y|^2) avoiding owerflow and underflow at
// intermediate stages of computation
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::real_unify_types2_promote<S1,S2,Float>::type
                            hypot(const S1& x,const S2& y);

// calculate power x^y
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename mrd::pow_return<S1,S2>::type
                            pow(const S1& x, const S2& y);

// calculate power x^y; both arguments are converted to complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename mrd::pow_c_return<S1,S2>::type
                            pow_c(const S1& x, const S2& y);

 // return x^y - 1; this function gives more accurate result for x and y
 // close to zero, than the pow function; not defined for complex arguments
template<class S1, class S2, 
    class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            powm1(const S1& x, const S2& y);

//---------------------------------------------------------------------------
//                      LOGICAL FUNCTIONS
//---------------------------------------------------------------------------
// logical and, i.e. x && y; both x and y are converted to boolean value
// first
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
bool                        op_and(const S1& x, const S2& y);

// logical or, i.e. x || y; both x and y are converted to boolean value
// first
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
bool                        op_or(const S1& x, const S2& y);

// logical symmetric difference; both x and y are converted to boolean
// value first
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
bool                        op_xor(const S1& x, const S2& y);

// element-wise and, i.e. x & y; for scalar types equivalent to x && y;
// for objects elem_and function is called
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            elem_and(const S1& x, const S2& y);

// element-wise and, i.e. x & y; for scalar types equivalent to x && y;
// for objects elem_and function is called
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            operator&(const S1& x, const S2& y)         { return elem_and(x, y); };

// element-wise or, i.e. x | y; for scalar types equivalent to x || y;
// for objects elem_or function is called
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            elem_or(const S1& x, const S2& y);

// element-wise or, i.e. x | y; for scalar types equivalent to x || y;
// for objects elem_or function is called
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            operator|(const S1& x, const S2& y)         { return elem_or(x, y); };

// element-wise xor, i.e. x ^ y; for scalar types equivalent to 
// (bool)x ^ (bool)y; for objects elem_xor function is called
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            elem_xor(const S1& x, const S2& y);

//---------------------------------------------------------------------------
//                      COMPARISON FUNCTIONS
//---------------------------------------------------------------------------

// for complex numbers the lexicographic ordering is defined

// equality comparison; x == y
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            eeq(const S1& x, const S2& y);

// equality comparison; x == y; different name for the eeq function
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            operator==(const S1& x, const S2& y)        { return eeq(x,y); };

// inequality comparison; x != y
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            neq(const S1& x, const S2& y);

// inequality comparison; x != y; differet name for the neq function
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            operator!=(const S1& x, const S2& y)        { return neq(x,y); };

// greater or equal comparison; x >= y
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            geq(const S1& x, const S2& y);

// greater or equal comparison; x >= y; different name for the geq function
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            operator>=(const S1& x, const S2& y)        { return geq(x,y); };

// greater than comparison; x > y
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            gt(const S1& x, const S2& y);

// greater than comparison; x > y; different name for the gt function
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            operator>(const S1& x, const S2& y)         { return gt(x,y); };

// less or equal comparison; x <= y
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            leq(const S1& x, const S2& y);

// less or equal comparison; x <= y; different name for the leq function
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            operator<=(const S1& x, const S2& y)        { return leq(x,y); };

// less than comparison; x < y
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            lt(const S1& x, const S2& y);

template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            operator<(const S1& x, const S2& y)         { return lt(x,y); };

// equality comparison; x == y; nan values are considered as equal
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            eeq_nan(const S1& x, const S2& y);

// inequality comparison; x != y; nan values are considered as equal
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::bool_or_object2<S1,S2>::type
                            neq_nan(const S1& x, const S2& y);

// return maximum value of x and y
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types_promote<S1,S2>::type
                            max(const S1& x, const S2& y);

// return minimum value of x and y
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types_promote<S1,S2>::type
                            min(const S1& x, const S2& y);

//---------------------------------------------------------------------------
//                      MISCELLANEOUS FUNCTIONS
//---------------------------------------------------------------------------

// returns a value with the magnitude of x and the sign of y;
// not available for complex arguments
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::real_unify_types_promote<S1,Float>::type
                            copysign(const S1& x, const S2& y);

// return max(x - y, 0) or NaN if x or y is NaN; not available for complex
// numbers
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2, Float>::type
                            fdim(const S1& x, const S2& y);

// return the next representable value after x in the direction of y;
// not defined for complex arguments
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::real_unify_types_promote<S1,Float>::type
                            nextafter(const S1& x, const S2& y);

// return the next representable value after x in the direction +INF;
// not available for complex arguments
template<class S1, class Enable = typename md::enable_if_scalar_ntobj<S1,void>::type>
typename md::real_unify_types_promote<S1,Float>::type
                            nextabove(const S1& x);

// return the next representable value after x in the direction -INF;
// not available for complex arguments
template<class S1, class Enable = typename md::enable_if_scalar_ntobj<S1,void>::type>
typename md::real_unify_types_promote<S1,Float>::type
                            nextbelow(const S1& x);

// return a * b + c; only one rounding at the end of computation is 
// performed
inline Real                 fma(Real a, Real b, Real c);
inline Float                fma(Float a, Float b, Float c);
inline Object               fma(const Object& a, const Object& b, const Object& c);

// return a * b - c; only one rounding at the end of computation is 
// performed
inline Real                 fms(Real a, Real b, Real c);
inline Float                fms(Float a, Float b, Float c);
inline Object               fms(const Object& a, const Object& b, const Object& c);

// compute a * b + c * d with high accuracy using Kahan's algorithm
inline Real                 dot2_ac(Real a, Real b, Real c, Real d);
inline Float                dot2_ac(Float a, Float b, Float c, Float d);
inline Object               dot2_ac(const Object& a, const Object& b, const Object& c, 
                                    const Object& d);

};

#include "matcl-scalar/details/func_binary.inl"