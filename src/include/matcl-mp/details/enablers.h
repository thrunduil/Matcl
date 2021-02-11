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

#include "matcl-mp/details/fwd_decls.h"
#include "matcl-mp/details/initializer.h"
#include "matcl-core/details/isa.h"
#include "matcl-mp/details/utils.h"

namespace matcl { namespace mp { namespace details
{

namespace md = matcl::details;

template<class T>
struct is_any_scalar
{
    static const bool is_matcl_scal = matcl::details::is_scalar<T>::value;
    static const bool is_mp_scal    = md::is_mp_scalar<T>::value;
    static const bool is_object     = matcl::details::is_object<T>::value;

    static const bool value         = (is_matcl_scal == true && is_object == false)
                                    || is_mp_scal == true;
};

template<class T1, class Ret = void>
struct enable_mp
    : std::enable_if<md::is_mp_scalar<T1>::value,Ret>
{};

template<class T1, class T2, class Ret>
struct enable_mp_bin
    : std::enable_if
        <   is_any_scalar<T1>::value && is_any_scalar<T2>::value
                && (md::is_mp_scalar<T1>::value || md::is_mp_scalar<T2>::value),
            Ret
        >
{};

//unify types and convert resulting rational to float
template<class T1, class T2>
struct real_unify_types_rat_to_float
{
    using type0 = typename md::unify_types<T1, T2>::type;
    using type  = typename md::select_if
                  <
                        std::is_same<type0,mp_rational>::value,
                        mp_float,
                        typename md::real_type<type0>::type
                  > ::type;
};

};};}
