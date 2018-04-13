/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

namespace matcl { namespace details
{

class ksp_solver_impl;
class scotch_partit_impl;
class metis_impl;

template<class V, class S>
struct schur_str;

template<class Val>
struct gschur_str;

template<class Val, class Struct>
struct gschur_sym_str;

struct reorder_diag;

class pschur_impl;

template<class V, class S>
struct dmperm_impl;

};};

typedef struct _p_KSP*     KSP;

namespace matcl
{
    class linear_operator;
    class linsolve_obj;
    class ksp_solver;
    class preconditioner;
};
