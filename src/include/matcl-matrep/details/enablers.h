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

#include "matcl-scalar/details/enablers.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-matrep/details/mpl.h"
#include "matcl-matrep/details/isa.h"
#include "matcl-matrep/details/utils.h"
#include "matcl-core/details/val_struct_codes.h"
#include "matcl-core/general/type_traits.h"

namespace matcl { namespace details
{

template<class ret_type,bool>
struct enable_if_matrix_impl
{};

template<class ret_type>
struct enable_if_matrix_impl<ret_type,true>
{
    using type = ret_type;
};

template<class T0, bool scalar, class T = typename std::decay<T0>::type>
struct enable_if_is_conv_to_mat :
    public enable_if
            <
                scalar == true && is_scalar<T>::value && !is_typed_object<T>::value
                    || scalar == false && is_matrix<T>::value,
                const void*
            >
{};

template<class M,class ret_type>
struct enable_if_matrix : public enable_if_matrix_impl<ret_type,is_matrix<M>::value>
{};

template<class M,class ret_type>
struct enable_if_matrix_or_scalar 
    : public enable_if_matrix_impl<ret_type,is_matrix<M>::value || is_scalar<M>::value>
{};

template<class T, class Ret>
struct enable_convertible_to_matrix
    : public std::enable_if<is_convertible_to_matrix<typename std::decay<T>::type>::value, Ret>
{};

};};

