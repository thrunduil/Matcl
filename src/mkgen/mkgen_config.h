/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019
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

namespace matcl { namespace mkgen 
{

//TODO
#define VEC_ALIGN   32

// compiler directive to mark nonaliasing
#define restricted  __restrict

// compiler directive to force inline expansion
#define force_inline        __forceinline
//#define force_inline        inline

// force inlining on different level
#define inline_lev_1        force_inline
#define inline_lev_root     force_inline
#define inline_force        force_inline
#define inline_initializer  force_inline
#define inline_expr         force_inline
#define inline_expr_split   force_inline
#define inline_loop         force_inline
#define inline_no           

}}