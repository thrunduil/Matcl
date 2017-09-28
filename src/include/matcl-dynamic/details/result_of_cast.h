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

#include "matcl-core/details/result_of_impl.h"

//result of cast
namespace matcl { namespace impl
{
    struct cast_functor
    {
        template<class To, class From>
        auto operator()(const To&, const From&) 
            -> decltype(convert_scalar<To, From>(std::declval<From>()));
        auto operator()(...) -> false_type;
    };
}};

namespace matcl { namespace result_of
{
    template<class To, class From>
    struct has_cast
    {
        using type              = typename std::result_of<matcl::impl::cast_functor (To,From)>::type;
        static const bool value = (std::is_same<type, matcl::impl::false_type>::value == false);
    };
    template<class To, class From, 
            bool Has_cast = has_cast<To,From>::value,
            bool Has_cons = matcl::dynamic::details::is_convertible_any<From, To>::value>
    struct result_of_cast;

    template<class To, class From>
    struct result_of_cast<To, From, true, false>
    {
        static const bool use_cast  = true;
        using type          = typename has_cast<To,From>::type;
        using type_object   = dynamic::object_type<typename std::decay<type>::type>;
    };

    template<class To, class From>
    struct result_of_cast<To, From, false, true>
    {
        static const bool use_cast  = false;
        using type          = To;
        using type_object   = dynamic::object_type<typename std::decay<type>::type>;
    };

    template<class To, class From>
    struct result_of_cast<To, From, true, true>
    {
        static const bool use_cast  = false;
        using type          = To;
        using type_object   = dynamic::object_type<typename std::decay<type>::type>;
    };

    template<class To, class From>
    struct result_of_cast<To, From, false, false>
    {};

}}
