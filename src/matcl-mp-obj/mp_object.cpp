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

namespace matcl { namespace dynamic
{

namespace md    = matcl::details;
namespace mmd   = matcl::mp::details;

std::string object_type_traits<mp_int>::to_string(const mp_int& t, printer& pr)
{ 
    (void)pr;
    return t.to_string();
};

void object_type_traits<mp_int>::disp(const mp_int& t, printer& pr, Integer elem_width,
                            align_type at, Integer value_pos)
{
    pr.disp_elem(elem_width, to_string(t, pr), at, value_pos);
};


std::string object_type_traits<mp_float>::to_string(const mp_float& t, printer& pr)
{ 
    Integer prec    = pr.get_precision();
    return t.to_string(precision(prec));

};

void object_type_traits<mp_float>::disp(const mp_float& t, printer& pr, Integer elem_width,
                            align_type at, Integer value_pos)
{
    std::ostringstream os;
    os << to_string(t, pr) << " (" << t.get_precision().get() << ")";
    pr.disp_elem(elem_width, os.str(), at, value_pos);
};

std::string object_type_traits<mp_complex>::to_string(const mp_complex& t, printer& pr)
{ 
    Integer prec    = pr.get_precision();
    return t.to_string(precision(prec));
};

void object_type_traits<mp_complex>::disp(const mp_complex& t, printer& pr, Integer elem_width,
                            align_type at, Integer value_pos)
{
    std::ostringstream os;
    os << t.to_string(precision(5)) << " (" << t.get_precision().get() << ")";
    pr.disp_elem(elem_width, os.str(), at, value_pos);
};

std::string object_type_traits<mp_rational>::to_string(const mp_rational& t, printer& pr)
{
    (void)pr;
    return t.to_string();
};

void object_type_traits<mp_rational>::disp(const mp_rational& t, printer& pr, Integer elem_width,
                            align_type at, Integer value_pos)
{
    pr.disp_elem(elem_width, to_string(t, pr), at, value_pos);
};

}}
