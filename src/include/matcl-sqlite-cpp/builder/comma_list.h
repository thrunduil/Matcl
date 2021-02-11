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

#include <cassert>
#include <string>
#include <vector>

namespace matcl { namespace sql
{

template<class T>
class comma_list : public std::vector<T>
{
    private:
        using base_class = std::vector<T>;

    public:
        comma_list(){}

        comma_list(const T& t);

        comma_list<T>&  operator()(const T& t);
        std::string     to_str(std::string (T::*f)() const = &T::to_str) const;
};

template<class T>
comma_list<T>::comma_list(const T& t)
{
    base_class::push_back(t);
}

template<class T>
comma_list<T>& comma_list<T>::operator()(const T& t)
{
    base_class::push_back(t);
    return *this;
}

template<class T>
std::string comma_list<T>::to_str(std::string (T::*f)() const) const
{
    assert(!base_class::empty());

    typename base_class::const_iterator it = base_class::begin();

    std::string txt = (*it.*f)();

    for(++it; it != base_class::end(); ++it)
    {
        txt += ", ";
        txt += (*it.*f)();
    }

    return txt;
}

}}
