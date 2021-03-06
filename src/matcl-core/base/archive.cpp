/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2021
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

#include "matcl-core/config.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/IO/archive.h"
#include "matcl-core/IO/serialize.h"
#include "matcl-core/memory/alloc.h"

namespace matcl 
{

iarchive::iarchive(std::istream& is)
{
    m_impl          = matcl_new<archive_type>(is);
    m_need_delete   = true;
};

iarchive::iarchive(archive_type& ar)
{
    m_impl          = &ar;
    m_need_delete   = false;
}

iarchive::~iarchive()
{
    if (m_need_delete)
        matcl_delete(m_impl);
};

iarchive& iarchive::operator>>(Integer& v)
{
    details::serialize_load(get(),v,0);
    return *this;
};
iarchive& iarchive::operator>>(Real& v)
{
    details::serialize_load(get(),v,0);
    return *this;
};
iarchive& iarchive::operator>>(Float& v)
{
    details::serialize_load(get(),v,0);
    return *this;
};
iarchive& iarchive::operator>>(Complex& v)
{
    details::serialize_load(get(),v,0);
    return *this;
};
iarchive& iarchive::operator>>(Float_complex& v)
{
    details::serialize_load(get(),v,0);
    return *this;
};

oarchive::oarchive(std::ostream& os)
{
    m_impl          = matcl_new<archive_type>(os);
    m_need_delete   = true;
};
oarchive::oarchive(archive_type& ar)
{
    m_impl          = &ar;
    m_need_delete   = false;
};
oarchive::~oarchive()
{
    if (m_need_delete == true)
        matcl_delete(m_impl); 
};

oarchive& oarchive::operator<<(Integer v)
{
    details::serialize_save(get(),v,0);
    return *this;
};
oarchive& oarchive::operator<<(Real v)
{
    details::serialize_save(get(),v,0);
    return *this;
};
oarchive& oarchive::operator<<(Float v)
{
    details::serialize_save(get(),v,0);
    return *this;
};
oarchive& oarchive::operator<<(const Complex& v)
{
    details::serialize_save(get(),v,0);
    return *this;
};
oarchive& oarchive::operator<<(const Float_complex& v)
{
    details::serialize_save(get(),v,0);
    return *this;
};

};

// explicit instantiations of functions required by boost::serialize
// library

#include <boost/archive/impl/basic_binary_iprimitive.ipp>
#include <boost/archive/impl/basic_binary_oprimitive.ipp>
#include <boost/archive/impl/basic_binary_iarchive.ipp>

using iarchive_impl = boost::archive::basic_binary_iprimitive
            <class eos::portable_iarchive,char,struct std::char_traits<char> >;

using oarchive_impl = boost::archive::basic_binary_oprimitive
            <class eos::portable_oarchive,char,struct std::char_traits<char> >;

template iarchive_impl;
template oarchive_impl;

template void MATCL_CORE_EXPORT iarchive_impl::load(std::string &s);
template void MATCL_CORE_EXPORT iarchive_impl::load(char * t);
template void MATCL_CORE_EXPORT iarchive_impl::init();

template void MATCL_CORE_EXPORT oarchive_impl::save(const std::string &s);
template void MATCL_CORE_EXPORT oarchive_impl::save(const char * t);
template void MATCL_CORE_EXPORT oarchive_impl::init();
