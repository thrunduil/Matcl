/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-matrep/details/extract_type2_switch.h"

namespace matcl { namespace details
{

inline const char* get_trans_code(trans_type trans)
{
    switch(trans)
    {
        case trans_type::no_trans:
            return "No transpose";
        case trans_type::trans:
            return "Transpose";
        case trans_type::conj_trans:
            return "Conj transpose";
        default:
        {
            matcl_assert(0,"invalid case");
            throw error::error_general("invalid case");
        }
    };
};

template<class T> struct correct_type
{
    public:
        using type  = const T&;

        static const T& eval(const T& val)  { return val; };

    private:
        static const T& eval(T&& val);
};
template<> struct correct_type<Integer>
{
    using type  = Real;
    static Real eval(Integer val)           { return val; };
};

}};