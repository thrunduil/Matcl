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

#include "matcl-matrep/general/config.h"
#include "matcl-dynamic/object_type_traits.h"

namespace matcl { namespace dynamic
{

template<>
struct MATCL_MATREP_EXPORT object_type_traits<Matrix> : object_type_traits_default
{
    using T = Matrix;
    static const bool has_one       = true;
    static const bool is_clonable   = true;

    static bool         is_zero(const T& t);
    static bool         is_one(const T& t);

    static T            make_one(const T*);

    static void         disp(const T& t, matcl::details::printer& pr, Integer elem_width,
                             align_type at, Integer value_pos);

    // save and load to stream
    static bool         read(std::istream&, T& t);
    static void         write(std::ostream&, const T& t);

    static void	        load_data(iarchive_impl& ar, T& ret, unsigned int version);
    static void	        save_data(oarchive_impl& ar, const T& val, unsigned int version);
};

};};
