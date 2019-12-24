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

#pragma once

#include "mkgen/details/mkgen_fwd.h"
#include "mkgen/details/matrix/scalar_data.h"

namespace matcl { namespace mkgen { namespace details
{

// convert rational number to constant scalar value
template<Integer N, Integer D>
struct scal_data_const_value_tag_rational
    : mk::scal_data_const_value_tag<scal_data_const_value_tag_rational<N, D>>
{
    using rational  = mkd::scal_data_rational<N, D>;

    static void print(std::ostream& os, int prior)
    {
        rational::print(os, prior);
    };

    template<class Val>
    static constexpr Val value()    { return rational::value<Val>(); }
};

// convert rational number to scalar value
template<Integer  N, Integer D>
struct scal_data_value_tag_rational
    : mk::scal_data_value_tag<scal_data_value_tag_rational<N, D>>
{
    using rational  = mkd::scal_data_rational<N, D>;

    static void print(std::ostream& os, int prior)
    {
        rational::print(os, prior);
    };

    template<class Val>
    static Val value()              { return rational::value<Val>(); }
};

// convert constant scalar value to scalar value
template<class Tag, class Val>
struct scal_data_value_tag_const
    : mk::scal_data_value_tag<scal_data_value_tag_const<Tag, Val>>
{
    using const_scal    = mkd::scal_data_const_value<Tag, Val>;

    static void print(std::ostream& os, int prior)
    {
        const_scal::print(os, prior);
    };

    template<class Val_loc>
    static Val_loc value()          { return const_scal::value<Val_loc>(); }
};

// tag representing multiplication Tag1 x Tag2
template<class Tag1, class Tag2>
struct scal_data_const_value_tag_mult 
    : mk::scal_data_const_value_tag<scal_data_const_value_tag_mult<Tag1, Tag2>>
{
    template<class Val>
    static constexpr Val value()    { return Tag1::value<Val>() * Tag2::value<Val>(); }
};

// tag representing multiplication Tag1 x Tag2
template<class Tag1, class Tag2>
struct scal_data_value_tag_mult 
    : mk::scal_data_value_tag<scal_data_value_tag_mult<Tag1, Tag2>>
{
    template<class Val>
    static Val value()  { return Tag1::value<Val>() * Tag2::value<Val>(); }
};

// tag representing addition Tag1 + Tag2
template<class Tag1, class Tag2>
struct scal_data_const_value_tag_plus 
    : mk::scal_data_const_value_tag<scal_data_const_value_tag_plus<Tag1, Tag2>>
{
    template<class Val>
    static constexpr Val value()    { return Tag1::value<Val>() + Tag2::value<Val>(); }
};

// tag representing addition Tag1 + Tag2
template<class Tag1, class Tag2>
struct scal_data_value_tag_plus 
    : mk::scal_data_value_tag<scal_data_value_tag_plus<Tag1, Tag2>>
{
    template<class Val>
    static Val value()  { return Tag1::value<Val>() + Tag2::value<Val>(); }
};

// tag representing division 1 / Tag1
template<class Tag1>
struct scal_data_const_value_tag_inv 
    : mk::scal_data_const_value_tag<scal_data_const_value_tag_inv<Tag1>>
{
    static_assert(Tag1::value<double>() != 0, "inversion of zero");

    template<class Val>
    static constexpr Val value()    { return Val(1) / Tag1::value<Val>(); }
};

// tag representing division 1 / Tag1
template<class Tag1>
struct scal_data_value_tag_inv 
    : mk::scal_data_value_tag<scal_data_value_tag_inv<Tag1>>
{
    template<class Val>
    static Val value()  { return Val(1) / Tag1::value<Val>(); }
};

}}}