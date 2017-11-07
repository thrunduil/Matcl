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
#include "matcl-core/memory/memory.h"
#include "matcl-core/details/global_objects.inl"
#include "matcl-core/details/leak_detector.h"
#include "matcl-core/IO/out_stream_initializer.h"

#include <iostream>

namespace matcl
{

//------------------------------------------------------------
//                      matcl_initializer
//------------------------------------------------------------

// nifty counter
static int g_counter = 0;

void matcl_initializer::open()
{
    global_objects::open();
};

void matcl_initializer::close()
{
    global_objects::close();
};

matcl_initializer::matcl_initializer()
{
    if (g_counter == 0)
        open();

    ++g_counter;
}

matcl_initializer::~matcl_initializer()
{
    --g_counter;

    if (g_counter == 0)   
        close();
};

//------------------------------------------------------------
//                      global_objects_impl
//------------------------------------------------------------
class global_objects_impl
{
    private:
        using global_vec    = std::vector<global_object*>;

    private:
        global_vec  m_globals;

    public:
        global_objects_impl();
        ~global_objects_impl();

        void        open();
        void        close();

    public:
        void        clear_globals();
        void        close_globals();
        void        register_global(global_object* gl);
};

global_objects_impl::global_objects_impl()
{
};

global_objects_impl::~global_objects_impl()
{
};

void global_objects_impl::open()
{
    matcl_initializer::open_global<details::out_stream_initializer>();
    matcl_initializer::open_global<details::leak_detector>();
};

void global_objects_impl::close()
{
    matcl::free_caches();
    clear_globals();

    close_globals();

    #if MATCL_DEBUG_MEMORY
        details::leak_detector::report_leaks(std::cout);
    #endif
        
    matcl_initializer::close_global<details::leak_detector>();
    matcl_initializer::close_global<details::out_stream_initializer>();
};

void global_objects_impl::register_global(global_object* gl)
{
    m_globals.push_back(gl);
};

void global_objects_impl::clear_globals()
{
    for (auto gl : m_globals)
        gl->clear_global();
};

void global_objects_impl::close_globals()
{
    for (auto gl : m_globals)
        gl->close_global();
};

//------------------------------------------------------------
//                      global_object
//------------------------------------------------------------

global_object::global_object()
{
    global_objects::register_global(this);
};

//------------------------------------------------------------
//                      global_objects
//------------------------------------------------------------

static global_objects_impl* g_instance  = nullptr;

void global_objects::open()
{
    g_instance = new global_objects_impl();
    g_instance->open();
};

void global_objects::close()
{
    g_instance->close();
    delete g_instance;
};

void global_objects::register_global(global_object* gl)
{
    g_instance->register_global(gl);
};

};