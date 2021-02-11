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

#include "matcl-mp/details/initializer.h"
#include "matcl-mp/mp_float.h"
#include <algorithm>

namespace matcl { namespace mp { namespace details
{

inline precision result_prec(precision req_prec, precision arg_prec)
{
    if (req_prec == 0)
        return arg_prec;
    else
        return req_prec;
}

inline precision result_prec(precision req_prec, precision a1_prec, precision a2_prec)
{
    if (req_prec == 0)
        return precision(std::max(a1_prec, a2_prec));
    else
        return req_prec;
}

inline precision result_prec(precision req_prec, precision a1_prec, precision a2_prec, 
                             precision a3_prec)
{
    if (req_prec == 0)
        return precision(std::max(std::max(a1_prec, a2_prec), a3_prec));
    else
        return req_prec;
}

inline precision get_precision(const mp_float& a)   { return a.get_precision(); };
inline precision get_precision(const mp_complex& a) { return a.get_precision(); };
inline precision get_precision(Integer)             { return mp_float::get_default_precision(); };
inline precision get_precision(Float)               { return mp_float::get_default_precision(); };
inline precision get_precision(Real)                { return mp_float::get_default_precision(); };
inline precision get_precision(Float_complex)       { return mp_float::get_default_precision(); };
inline precision get_precision(Complex)             { return mp_float::get_default_precision(); };

};};};