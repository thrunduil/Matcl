/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

typedef float                   s_type_wr;
typedef double                  d_type_wr;
typedef struct {float r, i;}    c_type_wr;
typedef struct {double r, i;}   z_type_wr;
typedef /*long*/ int            i_type_wr;
typedef /*long*/ int            l_type_wr;
typedef l_type_wr               (*sel_fun_wr)(...);
typedef /*long*/ i_type_wr      ftn_len_wr;
