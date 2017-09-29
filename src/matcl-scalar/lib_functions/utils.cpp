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

#include "matcl-scalar/lib_functions/utils.h"
#include "matcl-core/general/thread.h"

namespace matcl
{

MATCL_THREAD_LOCAL static details::timer_base global_timer = {0, 0};

void tic()
{
    global_timer.tic();
};

Real toc()
{
    return global_timer.toc();
};

std::string	tocstr()
{
    return global_timer.tocstr();
};

void tocdisp()
{
    return global_timer.tocdisp();
};

};