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

#include "matcl-mp-obj/mp_object.h"
#include "matcl-mp/details/utils.h"
#include "matcl-core/details/IO/printer.h"

namespace matcl { namespace dynamic
{

namespace md    = matcl::details;
namespace mmd   = matcl::mp::details;

void object_type_traits<mp_int>::disp(const mp_int& t, md::printer& pr, Integer elem_width,
                            align_type at, Integer value_pos)
{
    pr.disp_elem(elem_width, t.to_string(), at, value_pos);
};

void object_type_traits<mp_float>::disp(const mp_float& t, md::printer& pr, Integer elem_width,
                            align_type at, Integer value_pos)
{
    Integer prec    = pr.get_precision();

    std::ostringstream os;
    os << t.to_string(precision(prec)) << " (" << t.get_precision().get() << ")";
    pr.disp_elem(elem_width, os.str(), at, value_pos);
};

void object_type_traits<mp_complex>::disp(const mp_complex& t, md::printer& pr, Integer elem_width,
                            align_type at, Integer value_pos)
{
    std::ostringstream os;
    os << t.to_string(precision(5)) << " (" << t.get_precision().get() << ")";
    pr.disp_elem(elem_width, os.str(), at, value_pos);
};

void object_type_traits<mp_rational>::disp(const mp_rational& t, md::printer& pr, Integer elem_width,
                            align_type at, Integer value_pos)
{
    pr.disp_elem(elem_width, t.to_string(), at, value_pos);
};

}}
