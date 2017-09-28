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

#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/general/fwd_decls.h"
#include "matcl-dynamic/details/utils.h"

#include <type_traits>

namespace matcl { namespace details
{

template<class T> struct is_object
{   
    static const bool value = matcl::dynamic::details::is_object<T>::value
                            || std::is_base_of<Object, T>::value; 
};

template<class T> struct is_typed_object
{   
    static const bool value = matcl::dynamic::details::is_typed_object<T>::value;
};

template<bool val>
            struct promote_object                       {};
template<>  struct promote_object<true>                 { using type = Object;            };

template<class T>
            struct promote_scalar                       :promote_object<is_object<T>::value>{ };
template<>  struct promote_scalar<bool>                 { using type = Integer;           };
template<>  struct promote_scalar<signed char>          { using type = Integer;           };
template<>  struct promote_scalar<unsigned char>        { using type = Integer;           };
template<>  struct promote_scalar<signed short>         { using type = Integer;           };
template<>  struct promote_scalar<unsigned short>       { using type = Integer;           };
template<>  struct promote_scalar<signed int>           { using type = Integer;           };
template<>  struct promote_scalar<unsigned int>         { using type = Integer;           };
template<>  struct promote_scalar<signed long>          { using type = Integer;           };
template<>  struct promote_scalar<unsigned long>        { using type = Integer;           };
template<>  struct promote_scalar<signed long long>     { using type = Integer;           };
template<>  struct promote_scalar<unsigned long long>   { using type = Integer;           };
template<>  struct promote_scalar<Real>                 { using type = Real;              };
template<>  struct promote_scalar<Float>                { using type = Float;             };
template<>  struct promote_scalar<long double>          { using type = Real;              };
template<>  struct promote_scalar<Float_complex>        { using type = Float_complex;     };
template<>  struct promote_scalar<Complex>              { using type = Complex;           };
template<>  struct promote_scalar<Object>               { using type = Object;            };

template<bool cond,class T1,class T2> 
struct select_if                { using type = T1 ;};

template<class T1,class T2> 
struct select_if<false,T1,T2>   { using type = T2 ;};

//allow use always false conditions in static_assert 
template<class T>
struct dependent_false          { static const bool value = false; };

template<class... T>
struct dependent_false_var      { static const bool value = false; };

template<class T, T Val>
struct dependent_value_false    { static const bool value = false; };

};};