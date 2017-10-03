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

#include "matcl-dynamic/object_type_traits.h"
#include "matcl-mp/mp_int.h"
#include "matcl-mp/mp_float.h"
#include "matcl-mp/mp_complex.h"
#include "matcl-mp/mp_rational.h"
#include "matcl-mp/func_binary.h"
#include "matcl-mp/func_unary.h"

namespace matcl { namespace dynamic
{

template<>
struct object_type_traits<mp_int> : object_type_traits_default
{
    using T = mp_int;
    static const bool has_one = true;

    static bool is_zero(const mp_int& t)    { return matcl::is_zero(t);};
    static bool is_one(const mp_int& t)     { return matcl::is_one(t);};
    static T    make_one(const T*)          { return T(1); };

    MATCL_MP_OBJ_EXPORT
    static void         disp(const mp_int& t, md::printer& pr, Integer elem_width,
                             align_type at, Integer value_pos);

};

template<>
struct object_type_traits<mp_float> : object_type_traits_default
{
    using T = mp_float;
    static const bool has_one = true;

    static bool is_zero(const mp_float& t)  { return matcl::is_zero(t);};
    static bool is_one(const mp_float& t)   { return matcl::is_one(t);};
    static T    make_one(const T*)          { return T(1); };

    MATCL_MP_OBJ_EXPORT
    static void         disp(const mp_float& t, md::printer& pr, Integer elem_width,
                             align_type at, Integer value_pos);
};

template<>
struct object_type_traits<mp_complex> : object_type_traits_default
{
    using T = mp_complex;
    static const bool has_one = true;

    static bool is_zero(const mp_complex& t){ return matcl::is_zero(t);};
    static bool is_one(const mp_complex& t) { return matcl::is_one(t);};
    static T    make_one(const T*)          { return T(1); };

    MATCL_MP_OBJ_EXPORT
    static void         disp(const mp_complex& t, md::printer& pr, Integer elem_width,
                             align_type at, Integer value_pos);
};

template<>
struct object_type_traits<mp_rational> : object_type_traits_default
{
    using T = mp_rational;
    static const bool has_one = true;

    static bool is_zero(const mp_rational& t)   { return matcl::is_zero(t);};
    static bool is_one(const mp_rational& t)    { return matcl::is_one(t);};
    static T    make_one(const T*)              { return T(1); };

    MATCL_MP_OBJ_EXPORT
    static void         disp(const mp_rational& t, md::printer& pr, Integer elem_width,
                             align_type at, Integer value_pos);
};

};};
