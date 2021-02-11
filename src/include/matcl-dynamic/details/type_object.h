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

#include "matcl-dynamic/details/type_impl.h"

#pragma warning(push)
#pragma warning(disable:4505) //unreferenced local function has been removed
#pragma warning(disable:4251)	//needs to have dll-interface 

namespace matcl  { namespace dynamic { namespace details 
{

template<class T>
class type_object : public type_impl
{
    public:
        virtual function            generate_function(predef_fun fun) const override;
        virtual data_type*          clone(const data_type*) const override;
        virtual data_type*          copy(const data_type*) const override;
        virtual bool                is_zero(const data_type* d) const override;
        virtual bool				is_one(const data_type* d) const  override;
        virtual void				disp(const data_type*, matcl::details::printer& pr, Integer elem_width, 
                                         align_type at, Integer value_pos) const override;
        virtual data_type*          create() const override;
        virtual bool                has_one() const override;
        virtual data_type*          create_one() const override;

        virtual void				save(oarchive_impl& ar, unsigned int version, 
                                         const data_type* A) const override;
        virtual data_type*          load(iarchive_impl& ar, unsigned int version) const override;

        virtual void				save(std::ostream& os, const data_type* data) const override;
        virtual data_type*	        load(std::istream& is) const override;        

    public:
        type_object()              : type_impl(mark_type<T>()){};

        virtual bool                is_reference() const override   { return false; }
        virtual Type                decay() const override          { return Type(this); };
};

};};};

#pragma warning(pop)

