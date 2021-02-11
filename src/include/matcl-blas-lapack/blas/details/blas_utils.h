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

#include "matcl-blas-lapack/blas/config_blas.h"

namespace matcl { namespace lapack { namespace details
{

//-----------------------------------------------------------------
//                      ENABLER
//-----------------------------------------------------------------
template<bool cond, class type> 
struct blas_enable_if_impl {};

template<class T>
struct blas_enable_if_impl<true,T>
{
    using type = T;
};

template<class cond, class lazy_type>
struct blas_enable_if : public blas_enable_if_impl<cond::value,lazy_type>{};

//-----------------------------------------------------------------
//                          TEST TYPE
//-----------------------------------------------------------------

template<class T> struct is_valid_type              {   static const bool value = false; };
template<>        struct is_valid_type<s_type>      {   static const bool value = true; };
template<>        struct is_valid_type<d_type>      {   static const bool value = true; };
template<>        struct is_valid_type<c_type>      {   static const bool value = true; };
template<>        struct is_valid_type<z_type>      {   static const bool value = true; };

template<class T1, class T2> 
struct is_valid_type2
{   
    static const bool value = is_valid_type<T1>::value && is_valid_type<T2>::value; 
};

template<class T> struct is_valid_type_re           {   static const bool value = false; };
template<>        struct is_valid_type_re<s_type>   {   static const bool value = true; };
template<>        struct is_valid_type_re<d_type>   {   static const bool value = true; };

template<class T> struct is_valid_type_im           {   static const bool value = false; };
template<>        struct is_valid_type_im<z_type>   {   static const bool value = true; };
template<>        struct is_valid_type_im<c_type>   {   static const bool value = true; };

template<class T> struct complex_type               {};
template<>        struct complex_type<s_type>       {   using type = c_type; };
template<>        struct complex_type<d_type>       {   using type = z_type; };
template<>        struct complex_type<c_type>       {   using type = c_type; };
template<>        struct complex_type<z_type>       {   using type = z_type; };

template<class T> struct real_type                  {};
template<>        struct real_type<s_type>          {   using type = s_type; };
template<>        struct real_type<d_type>          {   using type = d_type; };
template<>        struct real_type<c_type>          {   using type = s_type; };
template<>        struct real_type<z_type>          {   using type = d_type; };

template<class T1, class T2> struct is_equal        {   static const bool value = false; };
template<class T1>           struct is_equal<T1,T1> {   static const bool value = true; };

template<class T>	struct is_complex				{ static const bool value = false;};
template<>	        struct is_complex<c_type>		{ static const bool value = true;};
template<>	        struct is_complex<z_type>		{ static const bool value = true;};

template<class V1,class V2>
struct is_real_type1 : public is_equal<typename complex_type<V1>::type,V2> {};

template<class V1,class V2>
struct is_real_type2 : public is_equal<V1,typename complex_type<V2>::type> {};

enum lapack_code { LAPACK_C,LAPACK_D,LAPACK_S,LAPACK_Z };

template<class T>	struct lapack_type			    {};
template<>			struct lapack_type<d_type>	    { static const lapack_code value = LAPACK_D;};
template<>			struct lapack_type<s_type>	    { static const lapack_code value = LAPACK_S;};
template<>			struct lapack_type<z_type>      { static const lapack_code value = LAPACK_Z;};
template<>			struct lapack_type<c_type>      { static const lapack_code value = LAPACK_C;};

//-----------------------------------------------------------------
//                  CONDITIONS
//-----------------------------------------------------------------

template<class V1, class V2>
struct equal_or_real1
{
    static const bool value =  is_valid_type<V1>::value
                            && is_valid_type<V2>::value
                            && ( is_equal<V1,V2>::value || is_real_type1<V1,V2>::value);
};

template<class ret, class V1, class V2>
struct enable_if_equal_or_real1 : public blas_enable_if<equal_or_real1<V1,V2>,ret> {};

template<class ret, class V>
struct enable_if_valid : public blas_enable_if<is_valid_type<V>,ret> {};

template<class ret, class V>
struct enable_if_valid_real : public blas_enable_if<is_valid_type_re<V>,ret> {};

template<class ret, class V>
struct enable_if_valid_complex : public blas_enable_if<is_valid_type_im<V>,ret> {};

template<class ret, class V1, class V2>
struct enable_if_valid2 : public blas_enable_if<is_valid_type2<V1,V2>,ret> {};

};};};

