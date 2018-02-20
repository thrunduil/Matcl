/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2018
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

#include "matcl-matrep/details/debug.h"

#ifdef MULTI_THREADED
	#include <atomic>
	using ATOMIC_LONG = std::atomic<long>;
#else
	using ATOMIC_LONG = long;
#endif

namespace matcl { namespace details
{

static ATOMIC_LONG n_containers;

void container_created()
{
    ++n_containers;
};

void container_destroyed()
{
    --n_containers;
};

int get_n_containers()
{
    return n_containers;
};

}}
