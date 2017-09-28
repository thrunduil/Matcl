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

#include "matcl-dynamic/details/register_object.h"
#include "matcl-core/general/thread.h"

#pragma warning(push)
#pragma warning(disable:4505) //unreferenced local function has been removed
#pragma warning(disable:4251)	//needs to have dll-interface 

namespace matcl  { namespace dynamic { namespace details 
{

//---------------------------------------------------------------------------
//                  OBJECT DATA
//---------------------------------------------------------------------------
class MATCL_DYN_EXPORT object_data_base
{
    private:
        using along         = matcl::atomic<long>;

    private:
        mutable along	    m_refcount;

    public:     

        virtual ~object_data_base();

        void                increase_refcount() const;
        void                decrease_refcount();
        bool                is_unique() const;

        template<class T> 
        const T&            get_value() const;

        template<class T>
        T&                  get_value();

    protected:
        object_data_base();
};

template<class T>
class object_data : public details::object_data_base, public register_object<T>
{
    public:
        using value_type        = T;

    private:
        value_type              m_data;

    private:
        object_data();
        object_data(const T& val);
        object_data(T&& val);

    public:
        static object_data*     create(const T& val);
        static object_data*     create(T&& val);

        const value_type&       get() const;
        value_type&             get();

        void*					operator new(size_t size);
        void					operator delete(void* ptr, size_t size);
};

};};};

#pragma warning(pop)

#include "matcl-dynamic/details/object_data.inl"
