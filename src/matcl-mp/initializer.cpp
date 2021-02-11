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

#include "matcl-mp/details/initializer.h"
#include "utils/impl_types.h"
#include "matcl-core/memory/alloc.h"

namespace matcl 
{

static void open()
{
    using allocator = default_allocator<true>;

    //set allocator
    mp_set_memory_functions 
    ( 
        &allocator::malloc,
        &allocator::realloc,
        &allocator::free
    );
};

static void close()
{
};

// nifty counter
static int g_counter = 0;

matcl_mp_initializer::matcl_mp_initializer()
{
    if (g_counter == 0)
        open();

    ++g_counter;
};

matcl_mp_initializer::~matcl_mp_initializer()
{
    --g_counter;

    if (g_counter == 0)   
        close();
}

};