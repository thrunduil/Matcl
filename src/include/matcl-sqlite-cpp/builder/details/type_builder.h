/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

namespace matcl { namespace sql { namespace details
{

template<unsigned flags, unsigned current>
struct unprocessed_flags
{
    static const bool result = 0 != (~(current-1) & flags);
};

template<unsigned flag, bool set, class T, unsigned all_flags, class TAG>
struct flag_to_type
{
    using result = T;
};

template<unsigned flag, class T, unsigned all_flags, class TAG>
struct flag_to_type<flag,true,T,all_flags,TAG>
{
    // no 'result', to force failure in case of missing specialization
};

template<unsigned flags, unsigned current, bool, class Base>
struct choose_base_impl
{
    private:
        static const unsigned next = current << 1;
        static const bool continue_choosing = unprocessed_flags< flags, next >::result;
        using next_type = typename choose_base_impl<flags, next, continue_choosing, Base>::result;
        static const bool flag_set = 0 != (current & flags);

    public:
        using result    = typename flag_to_type<
                                current,
                                flag_set,
                                next_type,
                                flags,
                                Base
                            >::result;
};

template<unsigned flags, unsigned current, class Base>
struct choose_base_impl<flags,current,false,Base>
{
    using result = Base;
};

template<unsigned flags, class Base>
struct choose_base
{
    using result = typename details::choose_base_impl<flags, 1, 
                        unprocessed_flags<flags, 1>::result, Base>::result;
};

}}}
