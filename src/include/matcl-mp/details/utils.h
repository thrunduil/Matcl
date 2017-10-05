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

#include "matcl-mp/details/fwd_decls.h"
#include "matcl-mp/details/initializer.h"
#include "matcl-core/details/utils.h"
#include "matcl-core/general/type_traits.h"

namespace matcl
{

// specialization of templates from "matcl-core/general/type_traits.h"

template<> struct is_external_scalar<mp_int>        { static const bool value = true; };
template<> struct is_external_scalar<mp_float>      { static const bool value = true; };
template<> struct is_external_scalar<mp_complex>    { static const bool value = true; };
template<> struct is_external_scalar<mp_rational>   { static const bool value = true; };

template<> struct make_complex_type<mp_int>         { using type = mp_complex; };
template<> struct make_complex_type<mp_float>       { using type = mp_complex; };
template<> struct make_complex_type<mp_complex>     { using type = mp_complex; };
template<> struct make_complex_type<mp_rational>    { using type = mp_complex; };

};

namespace matcl { namespace details
{

template<class T> 
            struct is_mp_scalar                 { static const bool value = false; };
template<>  struct is_mp_scalar<mp_int>         { static const bool value = true; };
template<>  struct is_mp_scalar<mp_float>       { static const bool value = true; };
template<>  struct is_mp_scalar<mp_complex>     { static const bool value = true; };
template<>  struct is_mp_scalar<mp_rational>    { static const bool value = true; };

template<> struct real_type<mp_int>             { using type = mp_int;      };
template<> struct real_type<mp_float>           { using type = mp_float;    };
template<> struct real_type<mp_rational>        { using type = mp_rational; };
template<> struct real_type<mp_complex>         { using type = mp_float;    };

template<> struct is_complex<mp_complex>        { static const bool value = true; };

template<> struct real_type_int_real<mp_int>    { using type = mp_float;    };
template<> struct real_type_int_real<mp_float>  { using type = mp_float;    };
template<> struct real_type_int_real<mp_rational>{ using type = mp_float;   };
template<> struct real_type_int_real<mp_complex>{ using type = mp_complex;  };

template<> struct complex_type<mp_int>          { using type = mp_complex;  };
template<> struct complex_type<mp_float>        { using type = mp_complex;  };
template<> struct complex_type<mp_rational>     { using type = mp_complex;  };
template<> struct complex_type<mp_complex>      { using type = mp_complex;  };

template<> struct unify_real_types<mp_int,Integer>      { using type = mp_int; };
template<> struct unify_real_types<mp_int,Float>        { using type = mp_float; };
template<> struct unify_real_types<mp_int,Real>         { using type = mp_float; };
template<> struct unify_real_types<mp_int,mp_int>       { using type = mp_int; };
template<> struct unify_real_types<mp_int,mp_float>     { using type = mp_float; };
template<> struct unify_real_types<mp_int,mp_rational>  { using type = mp_rational; };

template<> struct unify_real_types<mp_rational,Integer>      { using type = mp_rational; };
template<> struct unify_real_types<mp_rational,Float>        { using type = mp_float; };
template<> struct unify_real_types<mp_rational,Real>         { using type = mp_float; };
template<> struct unify_real_types<mp_rational,mp_int>       { using type = mp_rational; };
template<> struct unify_real_types<mp_rational,mp_float>     { using type = mp_float; };
template<> struct unify_real_types<mp_rational,mp_rational>  { using type = mp_rational; };

template<> struct unify_real_types<mp_float,Integer>        { using type = mp_float; };
template<> struct unify_real_types<mp_float,Float>          { using type = mp_float; };
template<> struct unify_real_types<mp_float,Real>           { using type = mp_float; };
template<> struct unify_real_types<mp_float,mp_int>         { using type = mp_float; };
template<> struct unify_real_types<mp_float,mp_float>       { using type = mp_float; };
template<> struct unify_real_types<mp_float,mp_rational>    { using type = mp_float; };

template<> struct unify_real_types<Integer, mp_int>         { using type = mp_int; };
template<> struct unify_real_types<Float, mp_int>           { using type = mp_float; };
template<> struct unify_real_types<Real, mp_int>            { using type = mp_float; };
template<> struct unify_real_types<Integer, mp_rational>    { using type = mp_rational; };
template<> struct unify_real_types<Float, mp_rational>      { using type = mp_float; };
template<> struct unify_real_types<Real, mp_rational>       { using type = mp_float; };
template<> struct unify_real_types<Integer, mp_float>       { using type = mp_float; };
template<> struct unify_real_types<Float, mp_float>         { using type = mp_float; };
template<> struct unify_real_types<Real, mp_float>          { using type = mp_float; };

template<> struct promote_scalar<mp_int>                    { using type = mp_int; };
template<> struct promote_scalar<mp_float>                  { using type = mp_float; };
template<> struct promote_scalar<mp_rational>               { using type = mp_rational; };
template<> struct promote_scalar<mp_complex>                { using type = mp_complex; };

};}

namespace matcl { namespace mp { namespace details
{

template<class T> struct promote_floats                     { using type = T; };
template<>        struct promote_floats<Float>              { using type = mp_float; };
template<>        struct promote_floats<Real>               { using type = mp_float; };
template<>        struct promote_floats<Float_complex>      { using type = mp_complex; };
template<>        struct promote_floats<Complex>            { using type = mp_complex; };

}}}