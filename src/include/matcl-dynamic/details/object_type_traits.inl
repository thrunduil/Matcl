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

#include "matcl-dynamic/object_type_traits.h"

namespace matcl { namespace dynamic
{

template<class T>
void object_type_traits_default::disp(const T& t, matcl::details::printer& pr,
                    Integer elem_width, align_type at, Integer value_pos)
{ 
    return t.disp(pr,elem_width,at,value_pos); 
};

template<class T>
bool object_type_traits_default::read(std::istream& is, T& val)
{
    is >> val;
    if (is.fail() || is.bad())  return false;
    else						return true;
};

template<class T>
void object_type_traits_default::write(std::ostream& os, const T& val)
{
    os << val;
};

template<class T>
void object_type_traits_default::load_data(iarchive_impl& ar, T& data, 
                                           unsigned int version)
{
    data.serialize(ar,version);
};

template<class T>
void object_type_traits_default::save_data(oarchive_impl& ar, const T& data, 
                                           unsigned int version)
{
	data.serialize(ar,version);
};

};};
