/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019 - 2021
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

#include "mkgen/mkgen.h"
#include "mkgen/expression/expressions.h"
#include "matcl-core/details/mpl.h"
#include "mkgen/details/utils/mpl.h"

#include <iostream>

namespace mk    = matcl::mkgen;
namespace mkl   = mk::list;
namespace mkd   = mk::details;

template <int X>
using int_t     = std::integral_constant<int, X>;

template<class V1, class V2>
struct cmp_ints;

template<int V1, int V2>
struct cmp_ints<int_t<V1>, int_t<V2>>
{
    static const bool value = (V1 <= V2);
};

struct comp_int
{
    template<class V1, class V2>
    static constexpr bool value = cmp_ints<V1, V2>::value;

    //template<class V1, class V2>
    //static constexpr bool value = mkd::less_types<V1, V2>::value;
};

static_assert(std::is_same<mkl::list<int_t<1>, int_t<2>, int_t<3>>,
                        typename mkd::merge_sort<comp_int, mkl::list<int_t<1>, int_t<2>, int_t<3>>>::type
                    >::value, "Failed a unit test");

static_assert(std::is_same<mkl::list<int_t<1>, int_t<2>, int_t<3>>,
                        typename mkd::merge_sort<comp_int, mkl::list<int_t<2>, int_t<1>, int_t<3>>>::type
                    >::value, "Failed a unit test");

static_assert(std::is_same<mkl::list<int_t<1>, int_t<2>, int_t<3>>,
                        typename mkd::merge_sort<comp_int, mkl::list<int_t<3>, int_t<2>, int_t<1>>>::type
                    >::value, "Failed a unit test");

template<int ... v>
struct int_list;

template<class L>
struct make_int_list
{};

template<class I>
struct tint_to_int;

template<int V>
struct tint_to_int<int_t<V>>
{
    static const int value = V;
};

template<class ... Elems>
struct make_int_list<mkl::list<Elems...>>
{
    using type = int_list<tint_to_int<Elems>::value ...>;
};

template<class T>
struct seq_to_list;

template<int ... V>
struct seq_to_list<std::index_sequence<V ...>>
{
    using type = mkl::list< int_t<V>...>;
};

using list_K    = typename seq_to_list<std::make_index_sequence<500>>::type;

static_assert(std::is_same<mkl::list<int_t<1>, int_t<2>, int_t<3>>,
                        typename mkd::merge_sort<comp_int, mkl::list<int_t<1>, int_t<2>, int_t<3>>>::type
                    >::value);

static_assert(std::is_same<mkl::list<int_t<1>, int_t<2>, int_t<3>>,
                        typename mkd::merge_sort<comp_int, mkl::list<int_t<2>, int_t<1>, int_t<3>>>::type
                    >::value);

static_assert(std::is_same<mkl::list<int_t<1>, int_t<2>, int_t<3>>,
                        typename mkd::merge_sort<comp_int, mkl::list<int_t<3>, int_t<2>, int_t<1>>>::type
                    >::value);

static_assert(std::is_same<list_K,
                        typename mkd::merge_sort<comp_int, list_K>::type>::value);
