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
#include "matcl-core/details/mpl.h"
#include "matcl-core/details/isa.h"
#include "matcl-scalar/details/utils.h"
#include "matcl-core/details/val_struct_codes.h"
#include "matcl-core/general/type_traits.h"

namespace matcl { namespace details
{

template<class ret_type,bool>
struct enable_if_scalar_impl
{};
template<class ret_type>
struct enable_if_scalar_impl<ret_type,true>
{
    using type = ret_type;
};

template<class T, class Real_T>
struct enable_if_val_complex : 
    public enable_if
            <	(std::is_same<Real_T,Real>::value || std::is_same<Real_T,Float>::value) 
                    && (std::is_same<T,Complex>::value || std::is_same<T,Float_complex>::value),
                const void*
            >
{};

template<class M, class ret_type>
struct enable_if_scalar : public enable_if_scalar_impl<ret_type,is_scalar<M>::value>
{};

//typed objects not allowed
template<class M, class ret_type>
struct enable_if_scalar_ntobj 
    : public enable_if_scalar_impl<ret_type,is_scalar<M>::value && !is_typed_object<M>::value>
{};

template<class M, class ret_type>
struct enable_if_matcl_scalar : public enable_if_scalar_impl<ret_type,is_matcl_scalar<M>::value>
{};

template<class M, bool Cond, class ret_type = void>
struct enable_if_external_scalar
    : public enable_if_scalar_impl<ret_type, Cond && is_external_scalar<M>::value>
{};

template<class S1, class S2, bool Cond, class ret_type = void>
struct enable_if_external_scalar2
    : public enable_if_scalar_impl
            < ret_type, Cond &&
                    (
                        is_external_scalar<S1>::value && is_any_scalar<S2>::value
                        ||is_external_scalar<S2>::value && is_any_scalar<S1>::value
                    )
            >
{};

template<class T, class Ret>
struct enable_if_float
    : public std::enable_if<std::is_same<T,Real>:: value || std::is_same<T,Float>:: value
                            || std::is_same<T,Complex>:: value || std::is_same<T,Float_complex>:: value,
                            Ret>
{};

template<class T, class Ret>
struct enable_if_complex
    : public std::enable_if<std::is_same<T,Complex>:: value || std::is_same<T,Float_complex>:: value, Ret>
{};

template<class T, class Ret>
struct enable_if_real
    : public std::enable_if<std::is_same<T,Real>:: value || std::is_same<T,Float>:: value
                            || std::is_same<T,Integer>:: value, Ret>
{};

template<class T, class Ret>
struct enable_if_real_float
    : public std::enable_if<std::is_same<T,Real>:: value || std::is_same<T,Float>:: value, Ret>
{};

template<class T, class Ret>
struct enable_if_not_object
    : public std::enable_if<is_object<typename std::decay<T>::type>:: value == false, Ret>
{};

template<class T, class Ret>
struct enable_if_matcl_scalar_not_object
    : public std::enable_if<is_matcl_scalar<typename std::decay<T>::type>:: value == true
                            && is_object<typename std::decay<T>::type>:: value == false
                            && std::is_same<T, typename std::decay<T>::type> :: value == true
                    , Ret>
{};

template<class T, class Ret>
struct enable_if_object
    : public std::enable_if<is_object<typename std::decay<T>::type>:: value == true, Ret>
{};

template<class M1,class M2, class ret_type>
struct enable_if_scalar2
    : public enable_if_scalar_impl
            <	ret_type,
                is_scalar<M1>::value && is_scalar<M2>::value
            >
{};

template<class M1,class M2, class ret_type>
struct enable_if_scalar2_ntobj
    : public enable_if_scalar_impl
            <	ret_type,
                is_scalar<M1>::value && is_scalar<M2>::value
                && is_typed_object<M1>::value == false
                && is_typed_object<M2>::value == false
            >
{};

template<class To, class From, class Ret>
struct enable_if_matcl_scalars2
    : public enable_if_scalar_impl
            <	Ret,
                is_matcl_scalar<To>::value && is_scalar<From>::value
            >
{};

template<class M1,class M2, class M3, class ret_type>
struct enable_if_scalar3
    : public enable_if_scalar_impl
            <	ret_type,
                is_scalar<M1>::value && is_scalar<M2>::value && is_scalar<M3>::value
            >
{};

template<class M1,class M2, class M3, class M4, class ret_type>
struct enable_if_scalar4
    : public enable_if_scalar_impl
            <	ret_type,
                is_scalar<M1>::value && is_scalar<M2>::value && is_scalar<M3>::value
                    && is_scalar<M4>::value
            >
{};

};};
