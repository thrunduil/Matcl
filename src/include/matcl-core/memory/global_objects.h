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

#include "matcl-core/config.h"

namespace matcl
{

class MATCL_CORE_EXPORT global_object
{
    protected:
        global_object();

    private:
        virtual void    clear_global() = 0;
        virtual void    close_global() = 0;

        friend class global_objects_impl;
};

class MATCL_CORE_EXPORT global_objects
{
    private:
        static void     open();
        static void     close();

        static void     register_global(global_object* gl);

        friend class global_object;
        friend class matcl_initializer;
};

class MATCL_CORE_EXPORT matcl_initializer
{
    public:
        matcl_initializer();
        ~matcl_initializer();

    private:
        template<class T>
        static void     open_global();

        template<class T>
        static void     close_global();

        friend class global_objects_impl;

        void            open();
        void            close();
};

static matcl_initializer matcl_init;

};

#include "matcl-core/details/global_objects.inl"