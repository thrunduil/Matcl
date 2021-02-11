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

#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-matrep/matrix/matrix.h"

namespace matcl
{

namespace md = matcl::details;

// complex conversion function:
// for a function fun defined for real numbers in a set D and defined
// for complex numbers, the complex conversion function (denoted by func_c)
// operating on a matrix A is defined as follows:
//     if each element of A is in the domain of D, then function func
//         is called on every element of A
//     otherwise A is converted to a complex matrix and the function func_c
//         is caled on every element of the converted matrix

//---------------------------------------------------------------------------
//                      VALUE CLASSIFICATION
//---------------------------------------------------------------------------

// return type of floating point number; not defined for complex values
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
fp_type                     fpclassify(const S& x);

// check if value is finite for each element in the matrix m of size M x N;
// return an M x N matrix with integer values
MATCL_MATREP_EXPORT Matrix  is_finite(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  is_finite(Matrix&& m);

// check if value is finite
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::bool_or_object<S>::type
                            is_finite(const S& x);

// check if value is nan for each element in the matrix m of size M x N;
// return an M x N matrix with integer values
MATCL_MATREP_EXPORT Matrix  is_nan(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  is_nan(Matrix&& m);

// check if value is infinite
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::bool_or_object<S>::type
                            is_nan(const S& x);

// check if value is infinite for each element in the matrix m of size 
// M x N; return an M x N matrix with integer values
MATCL_MATREP_EXPORT Matrix  is_inf(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  is_inf(Matrix&& m);

// check if value is infinite
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::bool_or_object<S>::type
                            is_inf(const S& x);

// for each element in the matrix m of size M x N check if value is regular
// (i.e. neither NaN, nor an infinity nor zero); for a complex value false 
// is retured if real or imaginary part is NaN or infinity or if both real 
// and imaginary part is zero; return an M x N matrix with integer values
MATCL_MATREP_EXPORT Matrix  is_regular(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  is_regular(Matrix&& m);

// check if value is regular (i.e., neither NaN, nor an infinity nor zero)
// for a complex value false is retured if real or imaginary part is
// NaN or infinity or if both real and imaginary part is zero
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::bool_or_object<S>::type
                            is_regular(const S& x);

// for each element in the matrix m of size M x N check if value is regular
// (i.e., neither NaN, nor an infinity nor zero nor subnormal);  return 
// an M x N matrix with integer values
MATCL_MATREP_EXPORT Matrix  is_normal(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  is_normal(Matrix&& m);

// check if value is regular (i.e., neither NaN, nor an infinity nor zero
// nor subnormal)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::bool_or_object<S>::type
                            is_normal(const S& x);

// for each element in the matrix m of size M x N check if value is an 
// integer; return an M x N matrix with integer values
MATCL_MATREP_EXPORT Matrix  is_int(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  is_int(Matrix&& m);

// check if value is an integer
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::bool_or_object<S>::type
                            is_int(const S& x);

// for each element in the matrix m of size M x N check if value is real;
// for a complex value return true if the imaginary part is zero; return
// an M x N matrix with integer values
MATCL_MATREP_EXPORT Matrix  is_real(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  is_real(Matrix&& m);

// check if value is real; for a complex value return true if the imaginary
// part is zero
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::bool_or_object<S>::type
                            is_real(const S& x);

// check if value is zero; not defined for a matrix
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
bool                        is_zero(const S& x);

// check if value is one; not defined for a matrix
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
bool                        is_one(const S& x);

//---------------------------------------------------------------------------
//                      LOGICAL OPERATIONS
//---------------------------------------------------------------------------

// cast to bool; equivalent to (bool)x; x must be a scalar or 1x1 matrix,
// otherwise an exception is thrown
inline bool                 cast_bool(const Matrix& x)      { return (bool)x; };
inline bool                 cast_bool(Matrix&& x)           { Matrix tmp(std::move(x)); return (bool)tmp; };

// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
bool                        cast_bool(const S& x);

// return logical negation of x, i.e !x; x must be a scalar or 1x1 matrix,
// otherwise an exception is thrown
inline bool                 op_not(const Matrix& x)         { return !x; };
inline bool                 op_not(Matrix&& x)              { Matrix tmp(std::move(x)); return !tmp; };

// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
bool                        op_not(const S& x);

// return logical negation of x; this is a different name for the op_not
// function; 
// logical negation operator and cast to bool operator for matrices are defined
// in the class Matrix; this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
bool                        operator!(const S& x);

// for each element in the matrix m of size M x N call the neg function; 
// return an M x N matrix
MATCL_MATREP_EXPORT Matrix  neg(const Matrix& A);
MATCL_MATREP_EXPORT Matrix  neg(Matrix&& A);

// call the element-wise negation operator (i.e. ~x); for scalar types this
// function is equivalent to !x (even for Integer types)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::bool_or_object<S>::type
                            neg(const S& x);

// for each element in the matrix m of size M x N call the neg function;  
// return an M x N matrix; this is a different name for the neg function
inline Matrix               operator~(const Matrix& A)      { return neg(A); };
inline Matrix               operator~(Matrix&& A)           { return neg(std::move(A)); };

// for each element in the matrix m of size M x N call the neg function; 
// return an M x N matrix; this is a different name for the neg function
inline Matrix               is_false(const Matrix& A)       { return neg(A);};
inline Matrix               is_false(Matrix&& A)            { return neg(std::move(A));};

// call the element-wise negation operator (i.e. ~x); this is a different 
// name for the neg function
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::bool_or_object<S>::type
                            is_false(const S& x);

// for each element in the matrix m of size M x N call the is_true function; 
// return an M x N matrix
MATCL_MATREP_EXPORT Matrix  is_true(const Matrix& A);
MATCL_MATREP_EXPORT Matrix  is_true(Matrix&& A);

// call the element-wise cast bool operator; for scalar types this function 
// is equivalent to (bool)x
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::bool_or_object<S>::type
                            is_true(const S& x);

//---------------------------------------------------------------------------
//                      ARITHMETIC FUNCTIONS
//---------------------------------------------------------------------------

// unary minus, i.e. -A; different name of the operator- function
MATCL_MATREP_EXPORT Matrix  uminus(const Matrix& A);
MATCL_MATREP_EXPORT Matrix  uminus(Matrix&& A);

// unary minus, i.e. -A; different name of the operator- function
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::promote_scalar<S>::type
                            uminus(const S& x);

// for each element in the matrix m of size M x N call the uminus function; 
// return an M x N matrix
inline Matrix               operator-(const Matrix& A)      { return uminus(A); };
inline Matrix               operator-(Matrix&& A)           { return uminus(std::move(A)); };

// return unary minus, i.e -x;
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable> inline
typename md::promote_scalar<S>::type
                            operator-(const S& x);

// unary plus operator; return the matrix unchanged
inline Matrix               operator+(const Matrix& A)      { return A; };
inline Matrix               operator+(Matrix&& A)           { Matrix tmp(std::move(A)); return tmp; };

// unary plus operator; always return the input value x
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable> inline
typename md::promote_scalar<S>::type
                            operator+(const S& x);

// for each element in the matrix m of size M x N call the invs function; 
// return an M x N matrix; (inv is reserved for the matrix inverse function)
MATCL_MATREP_EXPORT Matrix  invs(const Matrix& A);
MATCL_MATREP_EXPORT Matrix  invs(Matrix&& A);

// inverse of x, i.e. 1/x; if no errors occure, then mul(x, invs(x)) = 1
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S, Float>::type
                            invs(const S& x);

// inverse of x, i.e. 1/x; if no errors occur, then mmul(x, inv(x)) = 1;
// for all scalar types and most of objects this function is equivalent to
// invs (and mul is equivalent to mmul); however for OMatrix this is a 
// matrix inversion
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S, Float>::type
                            inv(const S& x);

// for each element in the matrix m of size M x N call the imaginary part
// of a value; return an M x N matrix with real values
MATCL_MATREP_EXPORT Matrix  imag(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  imag(Matrix&& m);

// return imaginary part of a value x
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::real_type_promote<S>::type
                            imag(const S& x);

// for each element in the matrix m of size M x N call the real part of a
// value; return an M x N matrix with real values
MATCL_MATREP_EXPORT Matrix  real(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  real(Matrix&& m);

// return real part of a value x
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::real_type_promote<S>::type
                            real(const S& x);

// for each element in the matrix m of size M x N call the absolute value;
// return an M x N matrix with real values
MATCL_MATREP_EXPORT Matrix  abs(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  abs(Matrix&& m);

// return the absolute value of x
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::real_type_promote<S>::type
                            abs(const S& x);

// for each element in the matrix m of size M x N call the absolute value
// squared; return an M x N matrix with real values
MATCL_MATREP_EXPORT Matrix  abs2(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  abs2(Matrix&& m);

// return the absolute value squared of x
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::real_type_promote<S>::type
                            abs2(const S& x);

// for each element in the matrix m of size M x N call the argument (polar
// angle) of a complex number; return an M x N matrix with real values
MATCL_MATREP_EXPORT Matrix  arg(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  arg(Matrix&& m);

// return argument (polar angle) of a complex number
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::real_type_promote_unify<S,Float>::type
                            arg(const S& x);

// different name of the arg function
MATCL_MATREP_EXPORT Matrix  angle(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  angle(Matrix&& m);

// different name of the arg function
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::real_type_promote_unify<S,Float>::type
                            angle(const S& x);

// for each element in the matrix m of size M x N call the conjugate 
// transposition; return an M x N matrix
MATCL_MATREP_EXPORT Matrix  conj(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  conj(Matrix&& m);

// return the conjugate transposition of x
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::promote_scalar<S>::type
                            conj(const S& x);

//---------------------------------------------------------------------------
//               EXPONENTIAL, LOGARITHMIC AND POWER FUNCTIONS
//---------------------------------------------------------------------------

// for each element in the matrix m of size M x N call the square root;
// return an M x N matrix
MATCL_MATREP_EXPORT Matrix  sqrt(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  sqrt(Matrix&& m);

// return the square root of x
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            sqrt(const S& x);

// return a complex conversion of the sqrt function
MATCL_MATREP_EXPORT Matrix  sqrt_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  sqrt_c(Matrix&& m);

// return the square root of x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            sqrt_c(const S& x);

// for each element in the matrix m of size M x N get value of sqrt(1+x)-1;
// more accurate than sqrt for x close to 0; return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  sqrt1pm1(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  sqrt1pm1(Matrix&& m);

// return sqrt(1+x) - 1; more accurate than sqrt for x close to 0
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            sqrt1pm1(const S& x);

// return a complex conversion of the sqrt1pm1 function
MATCL_MATREP_EXPORT Matrix  sqrt1pm1_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  sqrt1pm1_c(Matrix&& m);

// return the value of sqrt1pm1 function for x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            sqrt1pm1_c(const S& x);

// for each element in the matrix m of size M x N call the cubic root, for
// negative values result is also negative; not available for complex 
// values; return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  cbrt(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  cbrt(Matrix&& m);

// return the cubic root of x, for negative values result
// is also negative; not available for complex values
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            cbrt(const S& x);

// for each element in the matrix m of size M x N call the exponential 
// function; return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  exp(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  exp(Matrix&& m);

// return the exponential function, i.e. the e (Euler's number, 2.7182818)
// raised to the given power x; 
// for complex z = x + i*y, exp(z) = exp(x) * (cos(y) + i*sin(Y))
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            exp(const S& x);

// for each element in the matrix m of size M x N call the exp2 function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  exp2(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  exp2(Matrix&& m);

// return 2^x
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            exp2(const S& x);

// different name of the exp2 function
inline Matrix               pow2(const Matrix& m)   { return exp2(m); };
inline Matrix               pow2(Matrix&& m)        { return exp2(std::move(m)); };

// different name of the exp2 function
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            pow2(const S& x);

// for each element in the matrix m of size M x N call the exp10 function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  exp10(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  exp10(Matrix&& m);

// return 10^x
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            exp10(const S& x);

// different name of the exp10 function
inline Matrix               pow10(const Matrix& m)  { return exp10(m); };
inline Matrix               pow10(Matrix&& m)       { return exp10(std::move(m)); };

// different name of the exp10 function
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            pow10(const S& x);

// for each element in the matrix m of size M x N call the expm1 function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  expm1(const Matrix& x);
MATCL_MATREP_EXPORT Matrix  expm1(Matrix&& x);

// return exp(x) - 1; this function is more accurate for x near zero
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            expm1(const S& x);

// for each element in the matrix m of size M x N call the expi function; 
// return an M x N matrix with complex values
MATCL_MATREP_EXPORT Matrix  expi(const Matrix& x);
MATCL_MATREP_EXPORT Matrix  expi(Matrix&& x);

// return exp(x * 1i), where 1i is an imaginary unit
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            expi(const S& x);

// scale x by 2 raised to the power of n, returning ldexp(x,n) = x * 2^n
// this is redeclaration of function defined in matcl-scalar
template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
                            ldexp(const S1& x, Integer exp);

// function equivalent to scalbn
// this is redeclaration of function defined in matcl-scalar
template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
                            scalbn(const S1& x, Integer n);

// for each element in the matrix m of size M x N call the log function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  log(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  log(Matrix&& m);

// computes the natural (base e, Euler's number, 2.7182818) logarithm 
// of x; for complex argument z = x + i*y compute logarithm of a complex 
// value z with a branch cut along the negative real axis; if no errors
// occur, the complex natural logarithm of z is returned, in the range of
// a strip in the interval [-i pi, +i pi] along the imaginary axis and 
// mathematically unbounded along the real axis;
// Log(z) = log(|z|) + i Arg(z) = log( hypot(x,y) ) + i atan2(y,x)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            log(const S& x);

// return a complex conversion of the log function
MATCL_MATREP_EXPORT Matrix  log_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  log_c(Matrix&& m);

// return the value of log function for x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            log_c(const S& x);

// for each element in the matrix m of size M x N call the log1p function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  log1p(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  log1p(Matrix&& m);

// compute log(1+x); more accurate for m close to 0 
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            log1p(const S& x);

// return a complex conversion of the log function
MATCL_MATREP_EXPORT Matrix  log1p_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  log1p_c(Matrix&& m);

// return the value of log1p function for x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            log1p_c(const S& x);

// for each element in the matrix m of size M x N call the log2 function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  log2(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  log2(Matrix&& m);

// computes the base 2 logarithm of x
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            log2(const S& x);

// return a complex conversion of the log2 function
MATCL_MATREP_EXPORT Matrix  log2_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  log2_c(Matrix&& m);

// return the value of log2 function for x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            log2_c(const S& x);

// for each element in the matrix m of size M x N call the log10 function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  log10(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  log10(Matrix&& m);

// computes the base 10 logarithm of x
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            log10(const S& x);

// return a complex conversion of the log10 function
MATCL_MATREP_EXPORT Matrix  log10_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  log10_c(Matrix&& m);

// return the value of log10 function for x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            log10_c(const S& x);

// for each element in the matrix m of size M x N call the logb function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  logb(const Matrix& x);
MATCL_MATREP_EXPORT Matrix  logb(Matrix&& x);

// extract the value of the unbiased exponent from the floating-point argument
// x, and returns it as a floating-point value; equivalent to log2(|x|) rounded
// to integer number toward -INF; additionally for regular values 
// logb(x) = exp - 1, where exp is the exponent returned by frexp, 
// logb(0) = -inf, logb(+- inf) = inf, logb(nan) = nan; for complex argument 
// logb(|x|) is returned
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::real_unify_types_promote<S,Float>::type
                            logb(const S& x);

// for each element in the matrix m of size M x N call the ilogb function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  ilogb(const Matrix& x);
MATCL_MATREP_EXPORT Matrix  ilogb(Matrix&& x);

// cast result of logb to Integer
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::integer_or_object<S>::type
                            ilogb(const S& x);

//---------------------------------------------------------------------------
//                      TRIGONOMETRIC FUNCTIONS
//---------------------------------------------------------------------------

// for each element in the matrix m of size M x N call the sin function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  sin(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  sin(Matrix&& m);

// return the sine function of x in radians;
// for complex z = x + i*y, sin(z) = sin(x) * cosh(y) + i cos(x) * sinh(y)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            sin(const S& x);

// for each element in the matrix m of size M x N call the cos function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  cos(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  cos(Matrix&& m);

// return the cosine function of x in radians;
// for complex z = x + i*y, cos(z) = cos(x) * cosh(y) - i sin(x) * sinh(y)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            cos(const S& x);

// for each element in the matrix m of size M x N call the tan function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  tan(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  tan(Matrix&& m);

// return the tangent function of x in radians;
// for complex z = x + i*y, tan(z) = -i * tanh(i*z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            tan(const S& x);

// for each element in the matrix m of size M x N call the cot function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  cot(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  cot(Matrix&& m);

// return the cotangent function of x in radians;
// for complex z = x + i*y, cot(z) = i * coth(i*z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            cot(const S& x);

// for each element in the matrix m of size M x N call the sec function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  sec(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  sec(Matrix&& m);

// return the secant function of x in radians;
// for complex z = x + i*y, sec(z) = 1 / cos(z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            sec(const S& x);

// for each element in the matrix m of size M x N call the csc function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  csc(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  csc(Matrix&& m);

// return the cosecant function of x in radians;
// for complex z = x + i*y, csc(z) = 1 / sin(z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            csc(const S& x);

// for each element in the matrix m of size M x N call the sinh function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  sinh(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  sinh(Matrix&& m);

// return the hyperbolic sine function of x in radians;
// for complex z = x + i*y, sinh(z) = sinh(x) * cos(y) + i cosh(x) * sin(y)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            sinh(const S& x);

// for each element in the matrix m of size M x N call the cosh function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  cosh(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  cosh(Matrix&& m);

// return the hyperbolic cosine function of x in radians;
// for complex z = x + i*y, cosh(z) = cosh(x) * cos(y) + i sinh(x) * sin(y)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            cosh(const S& x);

// for each element in the matrix m of size M x N call the tanh function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  tanh(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  tanh(Matrix&& m);

// return the hyperbolic tangent function of x in radians;
// for complex z = x + i*y, tanh(z) = sinh(z) / cosh(z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            tanh(const S& x);

// for each element in the matrix m of size M x N call the coth function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  coth(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  coth(Matrix&& m);

// return the hyperbolic cotangent function of x in radians;
// for complex z = x + i*y, coth(z) = cosh(z) / sinh(z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            coth(const S& x);

// for each element in the matrix m of size M x N call the sech function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  sech(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  sech(Matrix&& m);

// return the hyperbolic secant function of x in radians;
// for complex z = x + i*y, sech(z) = 1 / cosh(z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            sech(const S& x);

// for each element in the matrix m of size M x N call the csch function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  csch(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  csch(Matrix&& m);

// return the hyperbolic cosecant function of x in radians;
// for complex z = x + i*y, csch(z) = 1 / sinh(z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            csch(const S& x);

// for each element in the matrix m of size M x N call the asin function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  asin(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  asin(Matrix&& m);

// return the inverse sine function of x, result in radians;
// for complex z:
// asin z = 1/i log(iz + sqrt(|1 - z^2|) * exp(i/2 arg(1-z^2)))
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            asin(const S& x);

// return a complex conversion of the asin function
MATCL_MATREP_EXPORT Matrix  asin_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  asin_c(Matrix&& m);

// return the value of asin function for x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            asin_c(const S& x);

// for each element in the matrix m of size M x N call the acos function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  acos(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  acos(Matrix&& m);

// return the inverse cosine function of x, result in radians
// for complex z:
// acos z = 1/i log(z + i * sqrt(|1 - z^2|) * exp(i/2 arg(1-z^2)))
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            acos(const S& x);

// return a complex conversion of the acos function
MATCL_MATREP_EXPORT Matrix  acos_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  acos_c(Matrix&& m);

// return the value of acos function for x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            acos_c(const S& x);

// for each element in the matrix m of size M x N call the atan function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  atan(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  atan(Matrix&& m);

// return the inverse tangent function of x, result in radians;
// for a complex z, atan(z) = -i atanh(iz)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            atan(const S& x);

// for each element in the matrix m of size M x N call the acot function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  acot(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  acot(Matrix&& m);

// return the inverse cotangent function of x, result in radians;
// for a complex z, acot(z) = atan(1/z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            acot(const S& x);

// for each element in the matrix m of size M x N call the asec function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  asec(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  asec(Matrix&& m);

// return the inverse secant function of x, result in radians;
// for a complex z, asec(z) = acos(1/z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            asec(const S& x);

// return a complex conversion of the asec function
MATCL_MATREP_EXPORT Matrix  asec_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  asec_c(Matrix&& m);

// return the value of asec function for x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            asec_c(const S& x);

// for each element in the matrix m of size M x N call the acsc function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  acsc(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  acsc(Matrix&& m);

// return the inverse cosecant function of x, result in radians;
// for a complex z, acsc(z) = asin(1/z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            acsc(const S& x);

// return a complex conversion of the acsc function
MATCL_MATREP_EXPORT Matrix  acsc_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  acsc_c(Matrix&& m);

// return the value of acsc function for x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            acsc_c(const S& x);

// for each element in the matrix m of size M x N call the asinh function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  asinh(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  asinh(Matrix&& m);

// return the inverse hyperbolic sine function of x, result in radians
// for a complex z, asinh(z) = i asin(-i z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            asinh(const S& x);

// for each element in the matrix m of size M x N call the acosh function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  acosh(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  acosh(Matrix&& m);

// return the inverse hyperbolic cosine function of x, result in radians
// for a complex z, acosh(z) = +-i acos(z) with posive sign if imag(acos(z)) < 0,
// thus real( acosh(z) ) >= 0
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            acosh(const S& x);

// return a complex conversion of the acosh function
MATCL_MATREP_EXPORT Matrix  acosh_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  acosh_c(Matrix&& m);

// return the value of acosh function for x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            acosh_c(const S& x);

// for each element in the matrix m of size M x N call the atanh function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  atanh(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  atanh(Matrix&& m);

// return the inverse hyperbolic tangent function of x, result in radians
// for complex z = x + i*y, atanh(z) = re + im *i, where
// re = log1p(4 * x / ( (x - 1) * (x - 1) + y^2))
// im = atan(2y, (1 - x) * (1 + x) - y^2 ) / 2
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            atanh(const S& x);

// return a complex conversion of the atanh function
MATCL_MATREP_EXPORT Matrix  atanh_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  atanh_c(Matrix&& m);

// return the value of atanh function for x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            atanh_c(const S& x);

// for each element in the matrix m of size M x N call the acoth function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  acoth(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  acoth(Matrix&& m);

// return the inverse hyperbolic cotangent function of x, result in radians;
// for a complex z, acoth(z) = atanh(1/z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            acoth(const S& x);

// return a complex conversion of the acoth function
MATCL_MATREP_EXPORT Matrix  acoth_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  acoth_c(Matrix&& m);

// return the value of acoth function for x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            acoth_c(const S& x);

// for each element in the matrix m of size M x N call the asech function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  asech(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  asech(Matrix&& m);

// return the inverse hyperbolic secant function of x, result in radians;
// for a complex z, asech(z) = acosh(1/z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            asech(const S& x);

// return a complex conversion of the asech function
MATCL_MATREP_EXPORT Matrix  asech_c(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  asech_c(Matrix&& m);

// return the value of asech function for x converted to a complex value
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
                            asech_c(const S& x);

// for each element in the matrix m of size M x N call the acsch function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  acsch(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  acsch(Matrix&& m);

// return the inverse hyperbolic cosecant function of x, result in radians;
// for a complex z, acsch(z) = asinh(1/z)
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
                            acsch(const S& x);

//---------------------------------------------------------------------------
//                      ROUNDING FUNCTIONS
//---------------------------------------------------------------------------

// for each element in the matrix m of size M x N calculate epsilon value
// eps, i.e positive distance from abs(x) to the next larger in magnitude
// floating point number of the same precision as x; return an M x N matrix
// with real values
MATCL_MATREP_EXPORT Matrix  eps(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  eps(Matrix&& m);

// return epsilon value eps, i.e positive distance from abs(x) to the next
// larger in magnitude floating point number of the same precision as x
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::real_unify_types_promote<S,Float>::type
                            eps(const S& x);

// for each element in the matrix m of size M x N call the floor function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  floor(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  floor(Matrix&& m);

// round to nearest integer towards -INF; for complex numbers
// real and imaginary part are rounded separately
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::promote_scalar<S>::type
                            floor(const S& x);

// for each element in the matrix m of size M x N call the ceil function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  ceil(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  ceil(Matrix&& m);

// round to nearest integer towards +INF; for complex numbers
// real and imaginary part are rounded separately
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::promote_scalar<S>::type
                            ceil(const S& x);

// for each element in the matrix m of size M x N call the round function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  round(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  round(Matrix&& m);

// round to nearest integer, rounding halfway cases away from zero; 
// for complex numbers real and imaginary part are rounded separately
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::promote_scalar<S>::type
                            round(const S& x);

// round to nearest integer towards zero; for complex numbers
// real and imaginary part are rounded separately
MATCL_MATREP_EXPORT Matrix  trunc(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  trunc(Matrix&& m);

// for each element in the matrix m of size M x N call the trunc function; 
// return an M x N matrix 
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::promote_scalar<S>::type
                            trunc(const S& x);

// different name of the trunc function
inline Matrix               fix(const Matrix& m)    { return trunc(m); };
inline Matrix               fix(Matrix&& m)         { return trunc(std::move(m)); };

// different name of the trunc function
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::promote_scalar<S>::type
                            fix(const S& x);

// for each element in the matrix m of size M x N call the ifloor function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  ifloor(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  ifloor(Matrix&& m);

// round to nearest integer towards -INF and cast to integer;
// not available for complex numbers
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::integer_or_object<S>::type
                            ifloor(const S& x);

// for each element in the matrix m of size M x N call the iceil function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  iceil(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  iceil(Matrix&& m);

// round to nearest integer towards +INF and cast to integer;
// not available for complex numbers
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::integer_or_object<S>::type
                            iceil(const S& x);

// for each element in the matrix m of size M x N call the iround function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  iround(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  iround(Matrix&& m);

// round to nearest integer, rounding halfway cases away from zero
// not available for complex numbers
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::integer_or_object<S>::type
                            iround(const S& x);

// for each element in the matrix m of size M x N call the itrunc function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  itrunc(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  itrunc(Matrix&& m);

// round to nearest integer towards zero and cast to integer; 
// not available for complex numbers
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::integer_or_object<S>::type
                            itrunc(const S& x);

// different name of the itrunc function
inline Matrix               ifix(const Matrix& m)   { return itrunc(m); };
inline Matrix               ifix(Matrix&& m)        { return itrunc(std::move(m)); };

// different name of the itrunc function
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::integer_or_object<S>::type
                            ifix(const S& x);

//---------------------------------------------------------------------------
//                 VALUE DECOMPOSITION FUNCTIONS
//---------------------------------------------------------------------------

// for each element in the matrix m of size M x N call the sign function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  sign(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  sign(Matrix&& m);

// the sign of x (-1 for negative, 1 for positive, 0 for zero, NaN if x is
// NaN); for complex arguments return x / abs(x), i.e. x normalized to the
// unit circle
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::promote_scalar<S>::type
                            sign(const S& x);

// for each element in the matrix m of size M x N call the isign function; 
// return an M x N matrix 
MATCL_MATREP_EXPORT Matrix  isign(const Matrix& m);
MATCL_MATREP_EXPORT Matrix  isign(Matrix&& m);

// the sign of x returned as integer (-1 for negative, 1 for positive, 
// 0 for plus or minus zero); not available for complex numbers
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::integer_or_object<S>::type
                            isign(const S& x);

// for each element in the matrix m of size M x N call the isign function; 
// return an M x N matrix with integer values
MATCL_MATREP_EXPORT Matrix  signbit(const Matrix& x);
MATCL_MATREP_EXPORT Matrix  signbit(Matrix&& x);

// check whether the sign of x is negative; not available for complex numbers
// this is redeclaration of function defined in matcl-scalar
template<class S, class Enable>
typename md::bool_or_object<S>::type
                            signbit(const S& x);

// decompose given floating point value arg into a normalized fraction and an 
// integral power of two, i.e. x is represented as x = d x 2 ^ exp, where 
// 0.5 <= |d| < 1; if x is zero, then d = 0 and exp = 0; if is Inf or NaN, then
// d = x; not available for complex numbers
// this is redeclaration of function defined in matcl-scalar
template<class S1, class Enable>
typename md::real_unify_types_promote<S1,Float>::type
                            frexp(const S1& x, Integer& exp);

// decompose given floating point value x into integral and fractional parts, 
// each having the same type and sign as x, the integral part is stored in 
// int_part argument and fractional part is returned; not available for complex
// numbers
// this is redeclaration of function defined in matcl-scalar
template<class S1, class Enable, class Return>
Return                      modf(const S1& x, Return& int_part);

//---------------------------------------------------------------------------
//                      COMBINATORICS
//---------------------------------------------------------------------------

// return the (2 * n)-th Bernoulli number B_{2n}, n >= 0
// this is redeclaration of function defined in matcl-scalar
Real                        bernoulli_b2n(Integer n);
Float                       fbernoulli_b2n(Integer n);

// return the (2 * n)-th Bernoulli number B_{2n}, n >= 0; select single or double
// precision using the template argument
// this is redeclaration of function defined in matcl-scalar
template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
                            bernoulli_b2n(Integer n);

// return largest value, such that bernoulli_b2n does not overflow
// this is redeclaration of function defined in matcl-scalar
Integer                     max_bernoulli_b2n();
Integer                     fmax_bernoulli_b2n();

// return largest value, such that bernoulli_b2n<T> does not overflow; select
// single or double precision using the template argument
// this is redeclaration of function defined in matcl-scalar
template<class S1, class Enable>
Integer                     max_bernoulli_b2n();

// return n-th prime numbers based on table lookup, where 0 <= n < N, and 
// N is given by the function prime_max_index()
// this is redeclaration of function defined in matcl-scalar
size_t                      prime(Integer n);

// return the number of precomputed prime numbers
// this is redeclaration of function defined in matcl-scalar
Integer                     prime_max_count();

// return the factorial i!, i.e. prod_{k=1}^i k
// this is redeclaration of function defined in matcl-scalar
Real                        factorial(Integer i);
Float                       ffactorial(Integer i);

// return the factorial i!, i.e. prod_{k=1}^i k; select single or double precision
// using the template argument
// this is redeclaration of function defined in matcl-scalar
template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
                            factorial(Integer i);

// return double factorial i!!:
//     i!! = i * (i-2) * ... * 1   for i odd
//     i!! = i * (i-2) * ... * 2   for i even
//     i!! = 1                     for i = 0, -1
// this is redeclaration of function defined in matcl-scalar
Real                        double_factorial(Integer i);
Float                       fdouble_factorial(Integer i);

// return double factorial i!!; select single or double precision using the 
// template argument
// this is redeclaration of function defined in matcl-scalar
template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
                            double_factorial(Integer i);

// return the rising factorial of x and i:
//     rising_factorial(x, i)  = x*(x+1)*...*(x+i-1)
// both x and i can be positive or negative; not defined for complex values
// this is redeclaration of function defined in matcl-scalar
template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
                            rising_factorial(const S1& x, Integer i);

// return the falling factorial of x and i:
//     falling_factorial(x, i) = x*(x-1)*...*(x-i+1)
// x can be positive or negative, i must be positive
// not defined for complex values
// this is redeclaration of function defined in matcl-scalar
template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
                            falling_factorial(const S1& x, Integer i);

// returns the binomial coefficient:
//     binomial_coefficient(n,k)   = n! / k! / (n-k)!
// requires k <= n, and k,n >= 0
// this is redeclaration of function defined in matcl-scalar
Real                        binomial_coefficient(Integer n, Integer k);
Float                       fbinomial_coefficient(Integer n, Integer k);

// returns the binomial coefficient; select single or double precision using the 
// template argument
// this is redeclaration of function defined in matcl-scalar
template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
                            binomial_coefficient(Integer n, Integer k);

//---------------------------------------------------------------------------
//                      FUNCTORS
//---------------------------------------------------------------------------

// evaluate function for all elements in a matrix. The function is given by
// a template scalar_function that must implement function eval:
// 
//         Ret<T> eval(const T& arg) const
// 
// for every scalar types T (i.e. Integer, Float, Real, Float_complex, Complex, Object)
// this function eval can be a template function, or may be given by appropriate 
// overload set; return type Ret<T> can be different for different
// scalar type T, but must be one of the scalar types; function f is always evaluated
// at zero argument and may not throw for zero argument
template<class scalar_function>
Matrix                      eval_scalar_func(const Matrix& A, const scalar_function& f);

template<class scalar_function>
Matrix                      eval_scalar_func(Matrix&& A, const scalar_function& f);

// in this version function f need not be defined for zero argument; test_function
// is evaluated at zero argument and must return value of the same type as function f
// evaluated at this value; returned value is used for determining structure of resulting
// matrix (only information whether return is zero or nonzero is used) and type info
// of stored elements (if return type is Object)
template<class scalar_function, class test_function>
Matrix                      eval_scalar_func(const Matrix& A, const scalar_function& f,
                                         const test_function& t);

template<class scalar_function, class test_function>
Matrix                      eval_scalar_func(Matrix&& A, const scalar_function& f,
                                         const test_function& t);

};

#include "matcl-matrep/details/func_unary.inl"
