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

#include "matcl-mp/details/initializer.h"
#include "matcl-mp/mp_float.h"
#include "matcl-core/general/thread.h"
#include "matcl-core/memory/alloc.h"

#include <map>

namespace matcl
{

// class representing a unique identifier of values stored in cache
// each instance of this class is a different identifier
class unique_id
{
    private:
        int         m_code;
        std::string m_name;

    public: 
        // create new unique identifier; optional_name is associated
        // name used only for debuging purposes
        unique_id(const std::string& optional_name);

        // comparison operators
        bool        operator<(const unique_id& other) const;
        bool        operator>(const unique_id& other) const;
        bool        operator<=(const unique_id& other) const;
        bool        operator>=(const unique_id& other) const;
        bool        operator==(const unique_id& other) const;
        bool        operator!=(const unique_id& other) const;

        // get optional name
        const std::string&
                    to_string() const;
};

// a thread local class that allow for caching computed values
// this class is thread save
class cache : matcl_new_delete, global_object
{
    private:
        using value_float   = std::pair<precision, mp_float>;
        using map_float     = std::map<unique_id, value_float>;

    private:
        map_float       m_map_float;

    public:
        // check if value with given id and precision at least 
        // prec is stored in cache
        static bool     has(unique_id id, precision prec);

        // get value with given id and precision at least prec;
        // throw exception if such value does not exis
        static mp_float get(unique_id id, precision prec);

        // add value with given id, precision prec to cache; if value
        // with greater precision is already stored, then stored value
        // is not updated
        static void     set(unique_id id, precision prec, const mp_float& val);

        // remove stored value with given id from cache
        static void     remove(unique_id id);

        // remove all values from cache; it also clears caches of 
        // implementation libraries (i.e. gmp, mpfr)
        static void     clear();

    private:
        static cache*   get();

        virtual void    clear_global() override;
        virtual void    close_global() override;
};

};
