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

#pragma once

#include "matcl-core/float/twofold.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/details/scalfunc_real.h"
#include "matcl-core/IO/scalar_io.h"

namespace matcl
{

//-----------------------------------------------------------------------
//                      ERROR RELATED FUNCTIONS
//-----------------------------------------------------------------------
template<class T>
T details::func_float_distance<T, false>::eval(const twofold<T>& x, const twofold<T>& y)
{
    namespace mrds = matcl::raw::details::scal_func;

    if (x.value == y.value && x.error == y.error)
        return 0.0;

    T d_err_1;
    T d_err_2;

    T d1   = mrds::float_distance(x.value, y.value);    
    T eps  = matcl::constants::eps<T>() * T(0.5);

    if (d1 > 1.0)
        return d1 / eps;

    T d    = d1/eps;

    if (x.error == 0.0)
    {
        d_err_1 = 0;
    }
    else
    {
        T e1        = mrds::eps(x.value) * T(0.5);
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
        T e2        = mrds::eps(y.value) * T(0.5);
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

template<class T>
void details::func_save<T, false>::eval(std::ostream& os, const twofold<T>& x)
{
    os << "{";
    
    details::saveload_scalar_helper::eval_print(os, x.value);
    os << ", ";
    
    details::saveload_scalar_helper::eval_print(os, x.error);
    os << "}";
};

template<class T>
void details::func_load<T, false>::eval(std::istream& is, twofold<T>& v)
{
    char c  = 0;

    // consume whitespaces
    while (is)
    {
        is.get(c);

        if (c != ' ' && c != '\t'  && c != '\n')
            break;
    }

    v = twofold<T>(std::numeric_limits<T>::quiet_NaN());

    if (is.good() == false)
        return;

    T val, err;

    if (c != '{')
    {
        // this is an error
        is.setstate(std::ios::failbit);
        return;
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
        return;
    }

    if (ok == false)
    {
        is.setstate(std::ios::failbit);
        return;
    };

    v   = twofold<T>(val, err);
};

template MATCL_CORE_EXPORT std::ostream& matcl::operator<<(std::ostream& os, const twofold<double>& x);
template MATCL_CORE_EXPORT std::istream& matcl::operator>>(std::istream& is, twofold<double>& v);
template MATCL_CORE_EXPORT std::ostream& matcl::operator<<(std::ostream& os, const twofold<float>& x);
template MATCL_CORE_EXPORT std::istream& matcl::operator>>(std::istream& is, twofold<float>& v);

template MATCL_CORE_EXPORT double float_distance(const twofold<double>& x, const twofold<double>& y);
template MATCL_CORE_EXPORT float float_distance(const twofold<float>& x, const twofold<float>& y);

}