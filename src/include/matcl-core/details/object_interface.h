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

#include "matcl-core/general/fwd_decls.h"
#include "matcl-core/matrix/enums.h"

namespace matcl { namespace details
{

class MATCL_CORE_EXPORT const_object_interface
{
    private:
        const dynamic::object*  m_object;

    public:
        const_object_interface(const dynamic::object* obj);

        bool        is_zero() const;

        void        disp_elem(printer& p, Integer w, align_type at, 
                        Integer value_pos) const;

        // stream formatting should not be changed
        void        write(std::ostream& os);
};

class MATCL_CORE_EXPORT nonconst_object_interface
{
    private:
        dynamic::object*  m_object;

    public:
        nonconst_object_interface(dynamic::object* obj);

        void        read(std::istream& is);
};

class MATCL_CORE_EXPORT const_type_interface
{
    private:
        const dynamic::Type*    m_type;

    public:
        const_type_interface(const dynamic::Type* ty);

        void        write(std::ostream& os);
};

class MATCL_CORE_EXPORT nonconst_type_interface
{
    private:
        dynamic::Type*  m_type;

    public:
        nonconst_type_interface(dynamic::Type* ty);

        void        read(std::istream& is);
};

class object_interface_impl
{
    public:
        virtual bool    is_zero(const dynamic::object& v) const = 0;
        virtual void    disp_elem(const dynamic::object& v, printer& p, Integer w, 
                            align_type at, Integer value_pos) const = 0;
        virtual void    read(dynamic::object& v, std::istream& is) const = 0;

        // stream formatting should not be changed
        virtual void    write(const dynamic::object& v, std::ostream& os) const = 0;

        virtual void    read_type(dynamic::Type& ty, std::istream& is) const = 0;
        virtual void    write_type(const dynamic::Type& ty, std::ostream& os) const = 0;
};

MATCL_CORE_EXPORT 
void set_object_intrface(object_interface_impl* oi);

}}
