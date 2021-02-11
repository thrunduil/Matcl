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

#include "matcl-scalar/details/utils.h"
#include "matcl-matrep/details/mpl.h"
#include "matcl-core/details/type_codes.h"
#include "matcl-core/details/val_struct_codes.h"
#include "matcl-core/details/utils.h"

namespace matcl { namespace details
{

namespace md = matcl::details;

template<class S1,class S2,class T>
struct object_or_type
{
    using type  = typename details::select_if
                <
                    std::is_same<S1,Object>::value || std::is_same<S2,Object>::value,
                    Object,
                    T
                >:: type;
};

template<class S1,class S2,class T>
struct lazy_object_or_type
{
    using type  = typename details::select_if
                <
                    std::is_same<S1,Object>::value || std::is_same<S2,Object>::value,
                    Object,
                    typename T::type
                >::type;
};

//convert compile time value to runtime value
template<class T, T val>
struct runtime_value
{
    static T eval() { return val; };
};

template<class Val>
struct get_type_name{};

template<> struct get_type_name<Integer>        { static std::string eval() { return "Integer"; } };
template<> struct get_type_name<Real>           { static std::string eval() { return "Real"; } };
template<> struct get_type_name<Float>          { static std::string eval() { return "Float"; } };
template<> struct get_type_name<Complex>        { static std::string eval() { return "Complex"; } };
template<> struct get_type_name<Float_complex>  { static std::string eval() { return "Float_complex"; } };
template<> struct get_type_name<Object>         { static std::string eval() { return "Object"; } };

struct MATCL_MATREP_EXPORT trans_manip
{
    static trans_type_ext   link_trans(trans_type_ext t1, trans_type_ext t2);
    static trans_type_ext   convert_trans(trans_type t1);

    // if t == trans_type_ext::conj then return trans_type::no_trans and set conj = true
    // otherwise return the same trans type and set conj = false
    static trans_type       convert_trans(trans_type_ext t1, bool& conj);
};

};};