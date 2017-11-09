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

//-----------------------------------------------------------------------
//                      ERROR RELATED FUNCTIONS
//-----------------------------------------------------------------------

double matcl::eps(const twofold& x)
{
    namespace mrds = matcl::raw::details::scal_func;

    double e = mrds::eps(x.value) * matcl::constants::eps();
    return e;
};

double matcl::float_distance(const twofold& x, const twofold& y)
{
    namespace mrds = matcl::raw::details::scal_func;

    if (x.value == y.value && x.error == y.error)
        return 0.0;

    double d_err_1;
    double d_err_2;

    double d1   = mrds::float_distance(x.value, y.value);    
    double eps  = matcl::constants::eps();

    if (d1 > 1.0)
        return d1 / eps;

    double d    = d1/eps;

    if (x.error == 0.0)
    {
        d_err_1 = 0;
    }
    else
    {
        double e1   = mrds::eps(x.value);
        d_err_1     = mrds::float_distance(e1, e1 + mrds::abs(x.error));

        if (x.error < 0.0)
            d_err_1 = -d_err_1;
    }

    if (y.error == 0.0)
    {
        d_err_2 = 0;
    }
    else
    {
        double e2   = mrds::eps(y.value);
        d_err_2     = mrds::float_distance(e2, e2 + mrds::abs(y.error));

        if (y.error < 0.0)
            d_err_2 = -d_err_2;
    }    

    if (d1 == 0.0)
        return mrds::abs(d_err_1 - d_err_2);

    if (x.value < y.value)
        d       = d + (d_err_2 - d_err_1);
    else
        d       = d + (d_err_1 - d_err_2);

    return mrds::abs(d);
};

//-----------------------------------------------------------------------
//                      IO FUNCTIONS
//-----------------------------------------------------------------------

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