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

#include "matcl-dynamic/details/type_impl.h"

namespace matcl { namespace dynamic { namespace details 
{

namespace md = matcl::details;

// reference type for internal use onlu
class reference_type : type_impl
{
    private:
        Type        m_base;

        //reference type to given base type can be created only once
        reference_type(Type base);

        friend type_table;

    public:                

        // empty implementation of these functions; reference type
        // cannot have instances; this type is for internal use only
        virtual function            generate_function(predef_fun fun) const override;
        virtual data_type*          clone(const data_type*) const override;
        virtual data_type*          copy(const object_data_base*) const override;
        virtual bool                is_zero(const data_type* d) const override;
        virtual bool				is_one(const data_type* d) const  override;
        virtual std::string         to_string(const data_type*, md::printer& pr) const override;
        virtual void				disp(const data_type*, md::printer& pr, Integer elem_width, 
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
        virtual bool                is_reference() const override   { return true; }
        virtual Type                decay() const override          { return m_base; };
};

};};};
