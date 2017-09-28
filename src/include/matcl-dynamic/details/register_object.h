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

#include "matcl-core/config.h"
#include "matcl-dynamic/details/fwd_decls.h"
#include "matcl-dynamic/type.h"

namespace matcl { namespace dynamic { namespace details
{

struct MATCL_DYN_EXPORT register_object_helper
{
    register_object_helper(const char*, Type (*creator)());

    template<class T>           
    static Type             create_object();
    static std::string      get_name(const char* name);
    Type                    get_object(const char* n);
};

template<class T> 
class mark_type
{
    public:
        Type			get() const;
        std::string		name() const;
};

template<> 
class mark_type<null_type>
{
    public:
        Type			get() const;
        std::string		name() const;
};

template<class T> 
class mark_reference_type
{
    public:
        Type			get() const;
};

template<class T>
struct register_object
{
    private:
        static details::register_object_helper  m_hook;
        friend mark_type<T>;
};

};};};

#include "matcl-dynamic/details/register_object.inl"
