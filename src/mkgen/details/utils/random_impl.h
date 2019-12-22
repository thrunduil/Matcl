/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019
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

// modified version of metarnd.h
// https://github.com/cr-lupin/metarand

#pragma once

#include <cstdint>
#include "matcl-core/details/mpl.h"

//----------------------------------------------------------------------------------
//                              forward declarations
//----------------------------------------------------------------------------------

namespace matcl { namespace mkgen { namespace details
{

template<class Engine>
struct random_state;

}}};

namespace matcl { namespace mkgen { namespace details { namespace random_impl
{

template<class Engine>
struct initialized_engine;

template<class Engine>
struct init;

template<class Eng>
struct make_next;

template<class Engine>
struct eval;

}}}}

namespace matcl { namespace mkgen { namespace details { namespace random_impl
{

namespace mkd = matcl :: mkgen :: details;

//----------------------------------------------------------------------------------
//                              time_seed
//----------------------------------------------------------------------------------

constexpr char inits[]  = __TIME__;

constexpr int time_seed = (inits[0]-'0')*100000+(inits[1]-'0')*10000 +
                     (inits[3]-'0')*1000+(inits[4]-'0')*100+(inits[6]-'0')*10+inits[7]-'0';


//----------------------------------------------------------------------------------
//                              eval fallback
//----------------------------------------------------------------------------------

template<class Engine>
struct eval
{
    static_assert(md::dependent_false<Engine>::value, "unknown engine");
};

//----------------------------------------------------------------------------------
//                              init fallback
//----------------------------------------------------------------------------------

template<class Engine>
struct init
{
    static_assert(md::dependent_false<Engine>::value, "unknown engine");
};

//----------------------------------------------------------------------------------
//                              initialization
//----------------------------------------------------------------------------------

template<class Engine>
struct initialized_engine
{
    using value_type                = decltype(Engine::value);
    static const value_type value   = Engine::value;
};

template<class Eng>
struct is_initialized
{
    static const bool value = false;
};

template<class Eng>
struct is_initialized<initialized_engine<Eng>>
{
    static const bool value = true;
};

//----------------------------------------------------------------------------------
//                              make_next
//----------------------------------------------------------------------------------

template<class Eng>
struct make_next
{};

template<class Eng>
struct make_next<initialized_engine<Eng>>
{
    using type0  = typename random_impl::eval<Eng>::type;
    using type  = mkd::random_state<initialized_engine<type0>>;
};

}}}}

namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              initialize
//----------------------------------------------------------------------------------

template<class Engine>
struct initialize
{
    using type  = random_impl::initialized_engine<typename random_impl::init<Engine>::type>;
};

}}}