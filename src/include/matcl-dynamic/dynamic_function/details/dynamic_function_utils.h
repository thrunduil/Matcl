/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

namespace matcl { namespace dynamic { namespace details
{

struct empty_type{};    

template<class T> struct get_base_type  { using type = T; };
template<> struct get_base_type<void>   { using type = empty_type; };

template<class T>
struct get_reference_type     
{ 
    using type = T&;
    static type get_value(T v)      { return v; };
};

template<class T>            
struct get_reference_type<T*> 
{ 
    using type = T&;
    static type get_value(T* v)     { return *v; };
};

};};};