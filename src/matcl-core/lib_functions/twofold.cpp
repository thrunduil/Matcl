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

#include "matcl-core/float/twofold.h"
#include "matcl-core/details/scalfunc_real.h"
#include "matcl-core/IO/scalar_io.h"

namespace matcl
{

std::ostream& matcl::operator<<(std::ostream& os, const twofold& x)
{
    os << "{";
    
    details::saveload_scalar_helper::eval_save(os, x.value);
    os << ", ";
    
    details::saveload_scalar_helper::eval_save(os, x.error);
    os << "}";

    return os;
};

std::istream& matcl::operator>>(std::istream& is, twofold& v)
{
    char c  = 0;

    // consume whitespaces
    while (is)
    {
        is.get(c);

        if (c != ' ' && c != '\t'  && c != '\n')
            break;
    }

    v = twofold(std::numeric_limits<double>::quiet_NaN());

    if (is.good() == false)
        return is;

    double val, err;

    if (c != '{')
    {
        // this is an error
        is.setstate(std::ios::failbit);
        return is;
    }

    char sep    = ' ';
    char fin    = ' ';

    bool ok     = true;
    
    ok &= details::saveload_scalar_helper::eval_load(is, val);

    is >> sep;

    ok &= details::saveload_scalar_helper::eval_load(is, err);

    is >> fin;

    if (sep != ',' || fin != '}')
    {
        // this is an error
        is.setstate(std::ios::failbit);
        return is;
    }

    if (ok == false)
    {
        is.setstate(std::ios::failbit);
        return is;
    };

    v   = twofold(val, err);
    return is;
};


}