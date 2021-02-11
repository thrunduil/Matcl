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

#include "matcl-core/IO/archive.h"
#include "matcl-matrep/general/config.h"

namespace matcl { namespace details
{

struct serialization_helper_base
{
    virtual ~serialization_helper_base(){};
};

}};

namespace matcl
{

// helper class for serializing/deserializing pointers to virtual classes
// derived from Base_class
template<class Base_class>
class serialization_helper : public details::serialization_helper_base
{
    public:
        // get serialization helper for objects of runtime type Derived;
        // one must assign the unique string identifier to this class
        template<class Derived>
        static serialization_helper* 
                            get(const std::string& unique_name);

        // call function virtual void save(oarchive& ar) const, which must
        // be defined in Derived class, where Derived is runtime type of ptr
        virtual void        save(oarchive& ar, const Base_class* ptr) const = 0;

        // call function virtual void save(std::ostream& ar) const, which must
        // be defined in Derived class, where Derived is runtime type of ptr
        virtual void        save(std::ostream& os, const Base_class* ptr) const = 0;

        // call function static Base_class* load(iarchive& ar), which must be
        // defined in Derived class, where Derived is the runtime type of saved 
        // pointer
        static Base_class*  load(iarchive& ar);

        // call function static Base_class* load(std::istream& ar), which must be
        // defined in Derived class, where Derived is the runtime type of saved 
        // pointer
        static Base_class*  load(std::istream& is);

    private:
        virtual Base_class* make_load(iarchive& ar) const = 0;
        virtual Base_class* make_load(std::istream& is) const = 0;

        static serialization_helper* 
                            get_impl(const std::string& unique_id);
};

};

#include "matcl-matrep/details/serialization_helper.inl"
