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

#pragma once

#include "mkgen/mkgen_fwd.h"
#include "mkgen/details/utils/less_types_impl.h"
#include "mkgen/details/utils/merge_sort_impl.h"
#include "mkgen/details/utils/random_impl.h"

namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              utils
//----------------------------------------------------------------------------------

//allow use always false conditions in static_assert 
// already defined in "matcl-core/details/mpl.h"
//template<class T>
//struct dependent_false        { static const bool value = false; };

template<class... T>
struct dependent_false_var      { static const bool value = false; };

template<class T, T Val>
struct dependent_value_false    { static const bool value = false; };

// hide type
template<class T>
struct lazy_type
{
    using type = T;
};

//----------------------------------------------------------------------------------
//                              type comparison
//----------------------------------------------------------------------------------
// return true if T1 is ordered less, than T2,
// ordering is based on type names
template<class T1, class T2>
struct less_types
{
    static const bool value = less_impl::less_impl<T1, T2>();
};

//----------------------------------------------------------------------------------
//                              merge_sort
//----------------------------------------------------------------------------------

// sort list of types TL = list<...> according to comparison function
// Compare; return sorted list of types;
// Compare must implement:
//
//      template<class T1, class T2>
//      static constexpr bool value = [impl]
//
// which should return true if T1 < T2
template<class Compare, class TL>
struct merge_sort
{
    using type      = typename merge_sort_impl :: merge_sort<Compare, TL>::type;
};

//----------------------------------------------------------------------------------
//                              random number
//----------------------------------------------------------------------------------

// random number generator before use must be initialized with this function
//      using RG = initialize<Engine>::type
template<class Engine>
struct initialize;

// class storing random state
template<class Engine>
struct random_state
{
    // check if engine is initialized
    static_assert(random_impl::is_initialized<Engine>::value, "initialize not called");

    // get current random number
    using value_type                = decltype(Engine::value);
    static const value_type value   = Engine::value;

    // generate next random number; next has type random_state<Engine>
    using next                      = typename random_impl::make_next<Engine>::type;
};

// constant integer generated based on current time
constexpr int get_time_seed()   { return random_impl::time_seed; };

// linear congruential engine
// generate random numbers according to
//  v(i) = (a * v(i-1) + c) % m
//  v(0) = seed
// If m == 0, then m = 2^64
// default parameters are taken from Knuth’s LCG MMIX algorithm, which
// performs well according to 
//      Bhattacharjee, Kamalika, Krishnendu Maity, and Sukanta Das. 
//      "A Search for Good Pseudo-random Number Generators: Survey and Empirical
//      Studies." arXiv preprint arXiv:1811.04035 (2018).
template<uint64_t seed = get_time_seed(),
         uint64_t a = 6364136223846793005ULL,
         uint64_t c = 1442695040888963407ULL,
         uint64_t m = 0ULL>
struct LCE;

}}}

#include "mkgen/details/utils/random_impl_engines.h"