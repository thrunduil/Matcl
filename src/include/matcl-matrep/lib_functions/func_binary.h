/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-scalar/lib_functions/func_binary.h"

//TODO: simd version of min, max, cmp

namespace matcl
{

namespace md    = matcl::details;
namespace mrd   = matcl::raw::details;

//---------------------------------------------------------------------------
//                      ARITHMETIC FUNCTIONS
//---------------------------------------------------------------------------

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the plus function; return an M x N matrix 
MATCL_MATFUNC_EXPORT Matrix plus(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix plus(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix plus(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix plus(Matrix&& A, Matrix&& B);

// addition x + y
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
                            plus(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the plus function; return an M x N matrix
inline Matrix               operator+(const Matrix& A, const Matrix& B)    {return plus(A,B); };
inline Matrix               operator+(Matrix&& A, const Matrix& B)         {return plus(std::move(A),B); };
inline Matrix               operator+(const Matrix& A, Matrix&& B)         {return plus(A,std::move(B)); };
inline Matrix               operator+(Matrix&& A, Matrix&& B)              {return plus(std::move(A),std::move(B)); };

// addition x + y; different name of the plus function
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
                            operator+(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the minus function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix minus(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix minus(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix minus(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix minus(Matrix&& A, Matrix&& B);

// subtraction x - y
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
                            minus(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the minus function; return an M x N matrix
inline Matrix               operator-(const Matrix& A, const Matrix& B)    {return minus(A,B); };
inline Matrix               operator-(Matrix&& A, const Matrix& B)         {return minus(std::move(A),B); };
inline Matrix               operator-(const Matrix& A, Matrix&& B)         {return minus(A,std::move(B)); };
inline Matrix               operator-(Matrix&& A, Matrix&& B)              {return minus(std::move(A),std::move(B)); };

// subtraction x - y; different name of the minus function
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
                            operator-(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the mul function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix mul(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix mul(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix mul(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix mul(Matrix&& A, Matrix&& B);

// multiplication x * y
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
                            mul(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the div function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix div(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix div(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix div(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix div(Matrix&& A, Matrix&& B);

// division x / y; integer division gives Real result; if x and y are
// integer values and y is zero, then result is zero
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
                            div(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the div_0 function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix div_0(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix div_0(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix div_0(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix div_0(Matrix&& A, Matrix&& B);

// division x / y; integer division gives Real result; 0/0 division gives 0
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
                            div_0(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the div_1 function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix div_1(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix div_1(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix div_1(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix div_1(Matrix&& A, Matrix&& B);

// division x / y; integer division gives Real result; 0/0 division gives 1
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
                            div_1(const S1& x, const S2& y);

// element by element division; always return a floating point values
inline Matrix               operator/(const Matrix& A, const Matrix& B)    {return div(A,B); };
inline Matrix               operator/(Matrix&& A, const Matrix& B)         {return div(std::move(A),B); };
inline Matrix               operator/(const Matrix& A, Matrix&& B)         {return div(A,std::move(B)); };
inline Matrix               operator/(Matrix&& A, Matrix&& B)              {return div(std::move(A),std::move(B)); };

// division x / y; for build-in types default divition operator is called
// (and integer division gives Integer result), otherwise equivalent to 
// div function
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
                            operator/(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the idiv function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix idiv(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix idiv(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix idiv(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix idiv(Matrix&& A, Matrix&& B);

// division x / y; if one of argument has floating point type then result
// is equivalent to div function; however if both arguments are integers,
// then integer division is performed; if x and y are integer values and
// y is zero, then result is zero
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
                            idiv(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the mod function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix mod(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix mod(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix mod(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix mod(Matrix&& A, Matrix&& B);

// return is x - n * y where n = floor(x / y) if y != 0 and 0 otherwise;
// is mod(x,y) is nonzero, then has the same sign as y; if signs of x and
// y are the same, then mod(x,y) = rem(x,y); not available for complex 
// arguments
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::real_unify_types_promote<S1,S2>::type
                            mod(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the rem function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix rem(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix rem(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix rem(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix rem(Matrix&& A, Matrix&& B);

// return is x - n * y where n = trunc(x / y) if y != 0 and 0 otherwise
// is rem(x,y) is nonzero, then has the same sign as x; if signs of x and
// y are the same, then mod(x,y) = rem(x,y); not available for complex 
// arguments
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::real_unify_types_promote<S1,S2>::type
                            rem(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the atan2 function; return an M x N matrix with
// real values
MATCL_MATFUNC_EXPORT Matrix atan2(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix atan2(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix atan2(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix atan2(Matrix&& A, Matrix&& B);

// four quadrant arctangent of x and y; -pi <= atan2(x, y) <= pi;
// not available for complex arguments
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::real_unify_types2_promote<S1,S2,Float>::type
                            atan2(const S1& x,const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the hypot function; return an M x N matrix with
// real values
MATCL_MATFUNC_EXPORT Matrix hypot(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix hypot(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix hypot(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix hypot(Matrix&& A, Matrix&& B);

// calculate sqrt(|x|^2 + |y|^2) avoiding owerflow and underflow at
// intermediate stages of computation
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::real_unify_types2_promote<S1,S2,Float>::type
                            hypot(const S1& x,const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the pow function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix pow(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix pow(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix pow(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix pow(Matrix&& A, Matrix&& B);

// calculate power x^y
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename mrd::pow_return<S1,S2>::type
                            pow(const S1& x, const S2& y);

// return a complex conversion of the pow function
MATCL_MATFUNC_EXPORT Matrix pow_c(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix pow_c(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix pow_c(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix pow_c(Matrix&& A, Matrix&& B);

// calculate power x^y; both arguments are converted to complex values
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename mrd::pow_c_return<S1,S2>::type
                            pow_c(const S1& x, const S2& y);

// return x^y - 1; this function gives more accurate result for x and y
// close to zero, than the pow function; not defined for complex arguments
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
                            powm1(const S1& x, const S2& y);

//---------------------------------------------------------------------------
//                      LOGICAL FUNCTIONS
//---------------------------------------------------------------------------
// logical and, i.e. x && y; both x and y are converted to boolean value
// first
inline bool                 op_and(const Matrix& x, const Matrix& y)    { return (bool)x && (bool)y; };

// logical and, i.e. x && y; both x and y are converted to boolean value
// first
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
bool                        op_and(const S1& x, const S2& y);

// logical or, i.e. x || y; both x and y are converted to boolean value
// first
inline bool                 op_or(const Matrix& x, const Matrix& y)     { return (bool)x || (bool)y; };

// logical or, i.e. x || y; both x and y are converted to boolean value
// first
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
bool                        op_or(const S1& x, const S2& y);

// logical symmetric difference; both x and y are converted to boolean
// value first
inline bool                 op_xor(const Matrix& x, const Matrix& y)    { return (bool)x ^ (bool)y; };

// logical symmetric difference; both x and y are converted to boolean
// value first
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
bool                        op_xor(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the elem_and function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix elem_and(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix elem_and(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix elem_and(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix elem_and(Matrix&& A, Matrix&& B);

// element-wise and, i.e. x & y; for scalar types equivalent to x && y;
// for objects elem_and function is called
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            elem_and(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the op_and function; return an M x N matrix
inline Matrix               operator&(const Matrix& A, const Matrix& B) {return elem_and(A,B); };
inline Matrix               operator&(Matrix&& A, const Matrix& B)      {return elem_and(std::move(A),B); };
inline Matrix               operator&(const Matrix& A, Matrix&& B)      {return elem_and(A,std::move(B)); };
inline Matrix               operator&(Matrix&& A, Matrix&& B)           {return elem_and(std::move(A),std::move(B)); };

// element-wise and, i.e. x & y; for scalar types equivalent to x && y;
// for objects elem_and function is called
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            operator&(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the op_or function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix elem_or(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix elem_or(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix elem_or(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix elem_or(Matrix&& A, Matrix&& B);

// element-wise or, i.e. x | y; for scalar types equivalent to x || y;
// for objects elem_or function is called
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            elem_or(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the op_or function; return an M x N matrix
inline Matrix               operator|(const Matrix& A, const Matrix& B) {return elem_or(A,B); };
inline Matrix               operator|(Matrix&& A, const Matrix& B)      {return elem_or(std::move(A),B); };
inline Matrix               operator|(const Matrix& A, Matrix&& B)      {return elem_or(A,std::move(B)); };
inline Matrix               operator|(Matrix&& A, Matrix&& B)           {return elem_or(std::move(A),std::move(B)); };

// element-wise or, i.e. x | y; for scalar types equivalent to x || y;
// for objects elem_or function is called
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            operator|(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the op_xor function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix elem_xor(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix elem_xor(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix elem_xor(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix elem_xor(Matrix&& A, Matrix&& B);

// element-wise xor, i.e. x ^ y; for scalar types equivalent to 
// (bool)x ^ (bool)y; for objects elem_xor function is called
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            elem_xor(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the op_xor function; return an M x N matrix
inline Matrix               operator^(const Matrix& A, const Matrix& B) {return elem_xor(A,B); };
inline Matrix               operator^(Matrix&& A, const Matrix& B)      {return elem_xor(std::move(A),B); };
inline Matrix               operator^(const Matrix& A, Matrix&& B)      {return elem_xor(A,std::move(B)); };
inline Matrix               operator^(Matrix&& A, Matrix&& B)           {return elem_xor(std::move(A),std::move(B)); };

//---------------------------------------------------------------------------
//                      COMPARISON FUNCTIONS
//---------------------------------------------------------------------------

// for complex numbers the lexicographic ordering is defined

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the eeq function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix eeq(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix eeq(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix eeq(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix eeq(Matrix&& A, Matrix&& B);

// equality comparison; x == y
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            eeq(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the eeq function; return an M x N matrix with
// integer values
inline Matrix               operator==(const Matrix& A, const Matrix& B){return eeq(A,B); };
inline Matrix               operator==(Matrix&& A, const Matrix& B)     {return eeq(std::move(A),B); };
inline Matrix               operator==(const Matrix& A, Matrix&& B)     {return eeq(A,std::move(B)); };
inline Matrix               operator==(Matrix&& A, Matrix&& B)          {return eeq(std::move(A),std::move(B)); };

// equality comparison; x == y; different name for the eeq function
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            operator==(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the neq function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix neq(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix neq(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix neq(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix neq(Matrix&& A, Matrix&& B);

// inequality comparison; x != y
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            neq(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the neq function; return an M x N matrix
inline Matrix               operator!=(const Matrix& A, const Matrix& B){return neq(A,B); };
inline Matrix               operator!=(Matrix&& A, const Matrix& B)     {return neq(std::move(A),B); };
inline Matrix               operator!=(const Matrix& A, Matrix&& B)     {return neq(A,std::move(B)); };
inline Matrix               operator!=(Matrix&& A, Matrix&& B)          {return neq(std::move(A),std::move(B)); };

// inequality comparison; x != y; differet name for the neq function
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            operator!=(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the geq function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix geq(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix geq(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix geq(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix geq(Matrix&& A, Matrix&& B);

// greater or equal comparison; x >= y
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            geq(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the geq function; return an M x N matrix
inline Matrix               operator>=(const Matrix& A, const Matrix& B){return geq(A,B); };
inline Matrix               operator>=(Matrix&& A, const Matrix& B)     {return geq(std::move(A),B); };
inline Matrix               operator>=(const Matrix& A, Matrix&& B)     {return geq(A,std::move(B)); };
inline Matrix               operator>=(Matrix&& A, Matrix&& B)          {return geq(std::move(A),std::move(B)); };

// greater or equal comparison; x >= y; different name for the geq function
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            operator>=(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the gt function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix gt(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix gt(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix gt(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix gt(Matrix&& A, Matrix&& B);

// greater than comparison; x > y
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            gt(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the gt function; return an M x N matrix
inline Matrix               operator>(const Matrix& A, const Matrix& B) {return gt(A,B); };
inline Matrix               operator>(Matrix&& A, const Matrix& B)      {return gt(std::move(A),B); };
inline Matrix               operator>(const Matrix& A, Matrix&& B)      {return gt(A,std::move(B)); };
inline Matrix               operator>(Matrix&& A, Matrix&& B)           {return gt(std::move(A),std::move(B)); };

// greater than comparison; x > y; different name for the gt function
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            operator>(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the leq function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix leq(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix leq(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix leq(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix leq(Matrix&& A, Matrix&& B);

// less or equal comparison; x <= y
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            leq(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the leq function; return an M x N matrix
inline Matrix               operator<=(const Matrix& A, const Matrix& B){return leq(A,B); };
inline Matrix               operator<=(Matrix&& A, const Matrix& B)     {return leq(std::move(A),B); };
inline Matrix               operator<=(const Matrix& A, Matrix&& B)     {return leq(A,std::move(B)); };
inline Matrix               operator<=(Matrix&& A, Matrix&& B)          {return leq(std::move(A),std::move(B)); };

// less or equal comparison; x <= y; different name for the leq function
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            operator<=(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the lt function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix lt(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix lt(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix lt(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix lt(Matrix&& A, Matrix&& B);

// less than comparison; x < y
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            lt(const S1& x, const S2& y);

// less than comparison; x < y; different name for the lt function
inline Matrix               operator<(const Matrix& A, const Matrix& B) {return lt(A,B); };
inline Matrix               operator<(Matrix&& A, const Matrix& B)      {return lt(std::move(A),B); };
inline Matrix               operator<(const Matrix& A, Matrix&& B)      {return lt(A,std::move(B)); };
inline Matrix               operator<(Matrix&& A, Matrix&& B)           {return lt(std::move(A),std::move(B)); };

// less than comparison; x < y; different name for the lt function
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            operator<(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the eeq_nan function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix eeq_nan(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix eeq_nan(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix eeq_nan(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix eeq_nan(Matrix&& A, Matrix&& B);

// equality comparison; x == y; nan values are considered as equal
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            eeq_nan(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the neq_nan function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix neq_nan(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix neq_nan(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix neq_nan(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix neq_nan(Matrix&& A, Matrix&& B);

// inequality comparison; x != y; nan values are considered as equal
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
                            neq_nan(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the max function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix max(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix max(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix max(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix max(Matrix&& A, Matrix&& B);

// return maximum value of x and y
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
                            max(const S1& x, const S2& y);

// for each element of M x N matrix A (or 1x1 matrix) and M x N matrix B
// (or 1x1 matrix) call the min function; return an M x N matrix
MATCL_MATFUNC_EXPORT Matrix min(const Matrix& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix min(Matrix&& A, const Matrix& B);
MATCL_MATFUNC_EXPORT Matrix min(const Matrix& A, Matrix&& B);
MATCL_MATFUNC_EXPORT Matrix min(Matrix&& A, Matrix&& B);

// return minimum value of x and y
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
                            min(const S1& x, const S2& y);

//---------------------------------------------------------------------------
//                      MISCELLANEOUS FUNCTIONS
//---------------------------------------------------------------------------

// returns a value with the magnitude of x and the sign of y;
// not available for complex arguments
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::real_unify_types_promote<S1,Float>::type
                            copysign(const S1& x, const S2& y);

// return max(x - y, 0) or NaN if x or y is NaN; not available for complex
// numbers
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2, Float>::type
                            fdim(const S1& x, const S2& y);

// return the next representable value after x in the direction of y;
// not defined for complex arguments
// redeclaration of function defined in matcl-scalar
template<class S1, class S2, class Enable>
typename md::real_unify_types_promote<S1,Float>::type
                            nextafter(const S1& x, const S2& y);

// return the next representable value after x in the direction +INF;
// not available for complex arguments
// redeclaration of function defined in matcl-scalar
template<class S1, class Enable>
typename md::real_unify_types_promote<S1,Float>::type
                            nextabove(const S1& x);

// return the next representable value after x in the direction -INF;
// not available for complex arguments
// redeclaration of function defined in matcl-scalar
template<class S1, class Enable>
typename md::real_unify_types_promote<S1,Float>::type
                            nextbelow(const S1& x);

//TODO: fast versions?
// return a * b + c; only one rounding at the end of computation is 
// performed
// redeclaration of function defined in matcl-scalar
inline Real                 fma(Real a, Real b, Real c);
inline Float                fma(Float a, Float b, Float c);
inline Object               fma(const Object& a, const Object& b, const Object& c);

// return a * b - c; only one rounding at the end of computation is 
// performed
// redeclaration of function defined in matcl-scalar
inline Real                 fms(Real a, Real b, Real c);
inline Float                fms(Float a, Float b, Float c);
inline Object               fms(const Object& a, const Object& b, const Object& c);

// compute a * b + c * d with high accuracy using Kahan's algorithm
// redeclaration of function defined in matcl-scalar
inline Real                 dot2_ac(Real a, Real b, Real c, Real d);
inline Float                dot2_ac(Float a, Float b, Float c, Float d);
inline Object               dot2_ac(const Object& a, const Object& b, const Object& c, 
                                    const Object& d);

//---------------------------------------------------------------------------
//                              FUNCTORS
//---------------------------------------------------------------------------

// evaluate function element-by-element; both matrices are converted to dense matrices
// The function is given by a template bin_function
// that must implement function eval:
// 
//         Ret<T1,T2> eval(const T1& arg1, const T2& arg2) const
// 
// for every scalar types T1, T2 (i.e. Integer, Float, Real, Float_complex, Complex, Object)
// this function eval can be a template function, or may be given by appropriate 
// overload set; return type Ret<T1,T2> can be different for different
// scalar types T1, T2, but must be one of the scalar types; function f is always evaluated
// at zero arguments (A, B, or both) and may not throw for zero arguments
template<class bin_function>
Matrix                  eval_binary_func(const Matrix& A, const Matrix& B, const bin_function& f);
template<class bin_function>
Matrix                  eval_binary_func(Matrix&& A, const Matrix& B, const bin_function& f);
template<class bin_function>
Matrix                  eval_binary_func(const Matrix& A, Matrix&& B, const bin_function& f);
template<class bin_function>
Matrix                  eval_binary_func(Matrix&& A, Matrix&& B, const bin_function& f);

// in this version function f need not be defined for zero argument; test_function
// is evaluated at zero argument and must return value of the same type as function f
// evaluated at this value
template<class bin_function, class test_function>
Matrix                  eval_binary_func(const Matrix& A, const Matrix& B, const bin_function& f,
                                const test_function& t);
template<class bin_function, class test_function>
Matrix                  eval_binary_func(Matrix&& A, const Matrix& B, const bin_function& f,
                                const test_function& t);
template<class bin_function, class test_function>
Matrix                  eval_binary_func(const Matrix& A, Matrix&& B, const bin_function& f,
                                const test_function& t);
template<class bin_function, class test_function>
Matrix                  eval_binary_func(Matrix&& A, Matrix&& B, const bin_function& f,
                                const test_function& t);

};

#include "matcl-matrep/details/func_binary.inl"
