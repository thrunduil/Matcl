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

#include "matcl-core/matrix/enums.h"
#include "matcl-dynamic/type.h"

namespace matcl { namespace dynamic { namespace details 
{

enum class predef_fun
{
    is_zero,
    is_one
};

class type_impl : public matcl_new_delete
{
    protected:
        using data_type             = details::object_data_base;

    private:
        std::string                 m_name;
        mutable bool                m_has_trivial_assignment;

        friend details::type_table;
        friend details::function_table;

    protected:
        type_impl(const std::string& name);

    public:
        //operations on data; defined only for class types

        template<class T>           type_impl(const mark_type<T>& derived);

        virtual function            generate_function(predef_fun fun) const = 0;
        virtual data_type*          clone(const data_type*) const = 0;
        virtual data_type*          copy(const object_data_base*) const = 0;
        virtual bool                is_zero(const data_type* d) const = 0;
        virtual bool				is_one(const data_type* d) const  = 0;
        virtual void				disp(const data_type*, matcl::details::printer& pr, Integer elem_width, 
                                         align_type at, Integer value_pos) const = 0;
        virtual data_type*          create() const = 0;
        virtual bool                has_one() const = 0;
        virtual data_type*          create_one() const = 0;

        virtual void				save(oarchive_impl& ar, unsigned int version, const data_type* A) const = 0;
        virtual data_type*          load(iarchive_impl& ar, unsigned int version) const = 0;

        virtual void				save(std::ostream& os, const data_type* data) const = 0;
        virtual data_type*	        load(std::istream& is) const = 0;

    public:
        // operations on type
        
        static const type_impl*     get(Type ti)                { return ti.m_impl; };
        const char*                 get_class_name() const      { return m_name.c_str(); };
        bool						has_trivial_assign() const  { return m_has_trivial_assignment; };
        std::string                 to_string() const           { return m_name; };

        virtual bool                is_reference() const = 0;
        virtual Type                decay() const = 0;
};

};};};
