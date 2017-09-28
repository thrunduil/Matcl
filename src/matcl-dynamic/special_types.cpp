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

#include "matcl-dynamic/special_types.h"
#include "matcl-core/IO/archive.h"
#include "matcl-dynamic/details/object_data.h"
#include "matcl-dynamic/object_type.h"

namespace matcl { namespace dynamic
{

//-------------------------------------------------------------------
//                      unit_type
//-------------------------------------------------------------------
std::ostream& dynamic::operator<<(std::ostream& os, const unit_type& )
{
    //nothing to do
    return os;
};
std::istream& dynamic::operator>>(std::istream& is, unit_type& )
{
    //nothing to do
    return is;
};

//-------------------------------------------------------------------
//                      null_type
//-------------------------------------------------------------------
std::ostream& dynamic::operator<<(std::ostream& os, const null_type& )
{
    //nothing to do
    return os;
};
std::istream& dynamic::operator>>(std::istream& is, null_type& )
{
    //nothing to do
    return is;
};

//-------------------------------------------------------------------
//                      any_type
//-------------------------------------------------------------------
any_type::any_type()
{};

any_type::any_type(const any_type& other)
    :m_stored(other.m_stored)
{}

any_type::any_type(any_type&& other)
    :m_stored(std::move(other.m_stored))
{};

any_type::any_type(const object& other)
    :m_stored(other)
{
    //other may store any_type object, we need to extract real object
    remove_indirections();
};

any_type::any_type(object&& other)
    :m_stored(std::move(other))
{
    //other may store any_type object, we need to extract real object
    remove_indirections();
};

any_type::~any_type()
{};

any_type& any_type::operator=(const any_type& other) &
{
    m_stored.reset(other.m_stored);
    return *this;
};

any_type& any_type::operator=(any_type&& other) &
{
    m_stored.reset(std::move(other.m_stored));
    return *this;
};

any_type& any_type::operator=(const object& other) &
{
    m_stored.reset(other);
    //other may store any_type object, we need to extract real object
    remove_indirections();
    return *this;
};

any_type& any_type::operator=(object&& other) &
{
    m_stored.reset(std::move(other));
    //other may store any_type object, we need to extract real object
    remove_indirections();
    return *this;
};

any_type::any_type(const dynamic::object_type<any_type>& other)
    :m_stored(other.get().m_stored)
{};

any_type::any_type(dynamic::object_type<any_type>&& other)
    :m_stored(other.get().m_stored)
{};

any_type& any_type::operator=(const dynamic::object_type<any_type>& other) &
{
    m_stored.reset(other.get().m_stored);
    return *this;
};

any_type& any_type::operator=(dynamic::object_type<any_type>&& other) &
{
    m_stored.reset(other.get().m_stored);
    return *this;
};

void any_type::remove_indirections()
{
    if (m_stored.get_type() != predefined::type_any())
        return;

    using impl_type = const details::object_data<any_type>*;

    impl_type impl          = static_cast<impl_type>(m_stored.get_data());
    const any_type& other   = impl->get();

    this->operator=(other);
};

bool any_type::is_zero() const
{
    return m_stored.is_null();
};

any_type::operator bool() const
{
    return is_zero() == false;
}

any_type any_type::clone() const
{
    return any_type(m_stored.clone());
};

std::string any_type::to_string(printer& pr) const
{
    if (m_stored.is_null() == true)
        return "null";

    return m_stored.to_string(pr);
};

void any_type::disp(printer& pr, Integer elem_width, align_type at, Integer value_pos) const
{
    if (m_stored.is_null() == true)
    {
        pr.disp_elem(elem_width, "null", at, 0);
        return;
    };

    m_stored.disp(pr,elem_width,at,value_pos);
};

bool any_type::operator==(const any_type& other) const
{
    return m_stored.get_data() == other.m_stored.get_data();
};

bool any_type::operator!=(const any_type& other) const
{
    return m_stored.get_data() != other.m_stored.get_data();
};

void any_type::serialize(oarchive_impl& ar, unsigned int version) const
{
    const_cast<object&>(m_stored).serialize(ar,version);
};

void any_type::serialize(iarchive_impl& ar, unsigned int version)
{
    m_stored.serialize(ar,version);
};;

std::ostream& dynamic::operator<<(std::ostream& os, const any_type& data)
{
    os << data.get_stored();
    return os;
};

std::istream& dynamic::operator>>(std::istream& is, any_type& data)
{
    object stored;
    is >> stored;

    data = any_type(stored);
    return is;
};

};};
