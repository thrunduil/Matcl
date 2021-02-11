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

#include "matcl-matrep/base/colon_info.h"
#include "matcl-matrep/lib_functions/manip.h"

namespace matcl { namespace details
{

colon_info::~colon_info()
{
    if (m_ci)
        delete m_ci;

    if (m_ri)
        delete m_ri;
};

void colon_info::set_ci(const raw::integer_dense& m) 
{
    if (m_ci == nullptr)
        m_ci    = new raw::integer_dense(ti::ti_int());

    m_ci->assign_to_fresh(m.make_explicit()); 
};

void colon_info::set_ri(const raw::integer_dense& m) 
{
    if (m_ri == nullptr)
        m_ri    = new raw::integer_dense(ti::ti_int());

    m_ri->assign_to_fresh(m.make_explicit()); 
};

void colon_info::set_ri_2(const Mat_I& r, const Mat_I& c)
{
    if (m_ri == nullptr)
        m_ri    = new raw::integer_dense(ti::ti_int());
    if (m_ci == nullptr)
        m_ci    = new raw::integer_dense(ti::ti_int());

    m_ri->assign_to_fresh(r.make_explicit()); 
    m_ci->assign_to_fresh(c.make_explicit()); 
};

};};
