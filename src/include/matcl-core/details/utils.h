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

#include "matcl-core/details/mpl.h"
#include "matcl-core/general/fwd_decls.h"

namespace matcl { namespace details
{

template<class T> struct complex_type           {};
template<> struct complex_type<Integer>         { using type = Complex;         };
template<> struct complex_type<Real>            { using type = Complex;         };
template<> struct complex_type<Float>           { using type = Float_complex;   };
template<> struct complex_type<Complex>         { using type = Complex;         };
template<> struct complex_type<Float_complex>   { using type = Float_complex;   };
template<> struct complex_type<Object>          { using type = Object;          };

template<class T> struct real_type              {};
template<> struct real_type<Integer>            { using type = Integer; };
template<> struct real_type<Real>               { using type = Real;    };
template<> struct real_type<Float>              { using type = Float;   };
template<> struct real_type<Complex>            { using type = Real;    };
template<> struct real_type<Float_complex>      { using type = Float;   };
template<> struct real_type<Object>             { using type = Object;  };

template<class T>   struct is_complex                   { static const bool value = false;};
template<>          struct is_complex<Complex>          { static const bool value = true; };
template<>          struct is_complex<Float_complex>    { static const bool value = true; };

template<class T> struct real_type_int_real         {};
template<> struct real_type_int_real<Integer>       { using type = Real;    };
template<> struct real_type_int_real<Real>          { using type = Real;    };
template<> struct real_type_int_real<Float>         { using type = Float;   };
template<> struct real_type_int_real<Complex>       { using type = Real;    };
template<> struct real_type_int_real<Float_complex> { using type = Float;   };
template<> struct real_type_int_real<Object>        { using type = Object;  };

template<class Val1, class Val> 
struct unify_real_types;

template<> struct unify_real_types<Integer,Integer> { using type = Integer; };
template<> struct unify_real_types<Integer,Float>   { using type = Real; };
template<> struct unify_real_types<Integer,Real>    { using type = Real; };
template<> struct unify_real_types<Integer,Object>  { using type = Object; };
template<> struct unify_real_types<Float,Integer>   { using type = Real; };
template<> struct unify_real_types<Float,Float>     { using type = Float; };
template<> struct unify_real_types<Float,Real>      { using type = Real; };
template<> struct unify_real_types<Float,Object>    { using type = Object; };
template<> struct unify_real_types<Real,Integer>    { using type = Real; };
template<> struct unify_real_types<Real,Float>      { using type = Real; };
template<> struct unify_real_types<Real,Real>       { using type = Real; };
template<> struct unify_real_types<Real,Object>     { using type = Object; };
template<> struct unify_real_types<Object,Integer>  { using type = Object; };
template<> struct unify_real_types<Object,Float>    { using type = Object; };
template<> struct unify_real_types<Object,Real>     { using type = Object; };
template<> struct unify_real_types<Object,Object>   { using type = Object; };

template<bool Is_complex, class Value>
struct make_complex_type                        { using type = Value; };

template<class Value>
struct make_complex_type<true,Value>            { using type = typename complex_type<Value>::type; };

template<class Value_1, class Value_2>
struct unify_types
{
    using real_1            = typename real_type<Value_1>::type;
    using real_2            = typename real_type<Value_2>::type;
    using ret_real          = typename unify_real_types<real_1, real_2>::type;
    static const bool isc   = is_complex<Value_1>::value || is_complex<Value_2>::value;

    using type              = typename make_complex_type<isc,ret_real>::type;
};

template<class Value_1>
struct unify_types<Value_1,Value_1>
{
    using type              = Value_1;
};

template<class Value_1, class Value_2>
struct unify_types<dynamic::object_type<Value_1>, Value_2>
{
    using real_1            = typename real_type<Value_1>::type;
    using real_2            = typename real_type<Value_2>::type;
    using ret_real          = typename unify_real_types<real_1, real_2>::type;
    static const bool isc   = is_complex<Value_1>::value || is_complex<Value_2>::value;

    using type0             = typename make_complex_type<isc,ret_real>::type;
    using type              = dynamic::object_type<type0>;
};

template<class Value_1, class Value_2>
struct unify_types<dynamic::object_type<Value_1>, dynamic::object_type<Value_2>>
{
    using real_1            = typename real_type<Value_1>::type;
    using real_2            = typename real_type<Value_2>::type;
    using ret_real          = typename unify_real_types<real_1, real_2>::type;
    static const bool isc   = is_complex<Value_1>::value || is_complex<Value_2>::value;

    using type0             = typename make_complex_type<isc,ret_real>::type;
    using type              = dynamic::object_type<type0>;
};

template<class Value_1, class Value_2>
struct unify_types<Value_1, dynamic::object_type<Value_2>>
{
    using real_1            = typename real_type<Value_1>::type;
    using real_2            = typename real_type<Value_2>::type;
    using ret_real          = typename unify_real_types<real_1, real_2>::type;
    static const bool isc   = is_complex<Value_1>::value || is_complex<Value_2>::value;

    using type0             = typename make_complex_type<isc,ret_real>::type;
    using type              = dynamic::object_type<type0>;
};

};};