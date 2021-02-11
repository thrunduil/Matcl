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

#include "matcl-core/matrix/enums.h"
#include "matcl-dynamic/result_of.h"

#include "matcl-dynamic/function_name.h"
#include "matcl-dynamic/predefined_functions_names.h"
#include "matcl-dynamic/details/result_of_cast.h"

namespace matcl { namespace dynamic
{

//--------------------------------------------------------------------
//               predefined functions for object_type
//--------------------------------------------------------------------
// a unary function func is enabled for object_type<T> if a function func(a)
// can be found by unqualified call in the namespace matcl for arguments of 
// type T unless otherwise stated
// similarly a binary function func is enabled for arguments of type 
// func(object_type<T>, object_type<S>), func(object_type<T>, S), func(T, object_type<S>)
// if a function func(a,b) can be found by unqualified call in the namespace 
// matcl for arguments of type T and S unless otherwise stated

// cast an object of type From to object of type To; 
// this function is enabled if a function convert_scalar<To>(From)
// can be found by unqualified call in the namespace matcl or the type
// To is constructible from an argument of type From (constructor is preferred
// over convert_scalar function)
template<class To, class From, 
    class Enable = typename result_of::result_of_cast<To,From>::type_object>
object_type<To>        cast(const object_type<From>& x);

// convert an object of type From to object of type To; 
// this function is enabled if the type To is constructible from an
// argument of type From
template<class To, class From>
object_type<To>         convert(const object_type<From>& x, 
                                typename details::enable_convert<From, To, void*>::type = 0);

// cast to boolean value; this function is enabled if a function 
// cast_bool(T) can be found by unqualified call in the namespace
// matcl
template<class T, class Enable = typename result_of::result_of_cast_bool<T>::type_object>
bool                    cast_bool(const object_type<T>& x);

// logical negation; should be equivalent to !cast_bool(a);
// this function is enabled if a function !T
// can be found by unqualified call in the namespace matcl
template<class T, class Enable = typename result_of::result_of_op_not<T>::type_object>
bool                    operator!(const object_type<T>& x);

// unary minus
// this function is enabled if a function operator-(T) can be found
// by unqualified call in the namespace matcl
template<class T>
typename result_of::result_of_uminus<T>::type_object
                        operator-(const object_type<T>& x);

// real part of a complex number
template<class T>
typename result_of::result_of_real<T>::type_object
                        real(const object_type<T>& x);

// imaginary part of a complex number
template<class T>
typename result_of::result_of_imag<T>::type_object
                        imag(const object_type<T>& x);

// check if value is zero; this function is enabled for all types T
template<class T>
bool                    is_zero(const object_type<T>& x);

// check if value is one; this function is enabled for all types T
template<class T>
bool                    is_one(const object_type<T>& x);

// equality comparison
// this function is enabled if a function T == S can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Res = typename result_of::result_of_eeq<T,S>::type_object>
Res                     operator==(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Res = typename result_of::result_of_eeq<T,S>::type_object>
Res                     operator==(const object_type<T>& x, const S& y);
template<class T, class S, class Res = typename result_of::result_of_eeq<T,S>::type_object>
Res                     operator==(const T& x, const object_type<S>& y);

// inequality comparison
// this function is enabled if a function T != S can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Res = typename result_of::result_of_neq<T,S>::type_object>
Res                     operator!=(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Res = typename result_of::result_of_neq<T,S>::type_object>
Res                     operator!=(const object_type<T>& x, const S& y);
template<class T, class S, class Res = typename result_of::result_of_neq<T,S>::type_object>
Res                     operator!=(const T& x, const object_type<S>& y);

// greater or equal comparison
// this function is enabled if a function T >= S can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Res = typename result_of::result_of_geq<T,S>::type_object>
Res                     operator>=(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Res = typename result_of::result_of_geq<T,S>::type_object>
Res                     operator>=(const object_type<T>& x, const S& y);
template<class T, class S, class Res = typename result_of::result_of_geq<T,S>::type_object>
Res                     operator>=(const T& x, const object_type<S>& y);

// less or equal comparison
// this function is enabled if a function T <= S can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Res = typename result_of::result_of_leq<T,S>::type_object>
Res                     operator<=(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Res = typename result_of::result_of_leq<T,S>::type_object>
Res                     operator<=(const object_type<T>& x, const S& y);
template<class T, class S, class Res = typename result_of::result_of_leq<T,S>::type_object>
Res                     operator<=(const T& x, const object_type<S>& y);

// greater than comparison
// this function is enabled if a function T > S can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Res = typename result_of::result_of_gt<T,S>::type_object>
Res                     operator>(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Res = typename result_of::result_of_gt<T,S>::type_object>
Res                     operator>(const object_type<T>& x, const S& y);
template<class T, class S, class Res = typename result_of::result_of_gt<T,S>::type_object>
Res                     operator>(const T& x, const object_type<S>& y);

// less than comparison
// this function is enabled if a function T < S can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Res = typename result_of::result_of_lt<T,S>::type_object>
Res                     operator<(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Res = typename result_of::result_of_lt<T,S>::type_object>
Res                     operator<(const object_type<T>& x, const S& y);
template<class T, class S, class Res = typename result_of::result_of_lt<T,S>::type_object>
Res                     operator<(const T& x, const object_type<S>& y);

// binary plus
// this function is enabled if a function T + S can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_plus<T,S>::type_object>
Ret                     operator+(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_plus<T,S>::type_object>
Ret                     operator+(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_plus<T,S>::type_object>
Ret                     operator+(const T& x, const object_type<S>& y);

// binary minus
// this function is enabled if a function T - S can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_minus<T,S>::type_object>
Ret                     operator-(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_minus<T,S>::type_object>
Ret                     operator-(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_minus<T,S>::type_object>
Ret                     operator-(const T& x, const object_type<S>& y);

// multiplication
// this function is enabled if a function T * S can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_op_mul<T,S>::type_object>
Ret                     operator*(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_op_mul<T,S>::type_object>
Ret                     operator*(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_op_mul<T,S>::type_object>
Ret                     operator*(const T& x, const object_type<S>& y);

// division
// this function is enabled if a function T / S can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_div<T,S>::type_object>
Ret                     operator/(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_div<T,S>::type_object>
Ret                     operator/(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_div<T,S>::type_object>
Ret                     operator/(const T& x, const object_type<S>& y);

// integer division
template<class T, class S, class Ret = typename result_of::result_of_idiv<T,S>::type_object>
Ret                     idiv(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_idiv<T,S>::type_object>
Ret                     idiv(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_idiv<T,S>::type_object>
Ret                     idiv(const T& x, const object_type<S>& y);

};};

#include "matcl-dynamic/details/typed_object_functions.inl"
