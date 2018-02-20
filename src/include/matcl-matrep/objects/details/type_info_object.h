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

#include "matcl-core/matrix/enums.h"
#include "matcl-dynamic/matcl_function_names.h"
#include "matcl-dynamic/type.h"

#include "matcl-matrep/general/config.h"

namespace matcl { namespace ti
{

struct ti_empty {};

template<class T>           struct make_ti                  {};
template<>                  struct make_ti<Integer>         {using type = ti_empty;};
template<>                  struct make_ti<Real>            {using type = ti_empty;};
template<>                  struct make_ti<Float>           {using type = ti_empty;};
template<>                  struct make_ti<Complex>         {using type = ti_empty;};
template<>                  struct make_ti<Float_complex>   {using type = ti_empty;};
template<>                  struct make_ti<Object>          {using type = dynamic::Type;};

struct MATCL_MATREP_EXPORT predefined
{
    using Type              = dynamic::Type;

    static Type             get_ti_int();
    static Type             get_ti_real();
    static Type             get_ti_float();
    static Type             get_ti_complex();
    static Type             get_ti_float_complex();
};

template<class T> 
class MATCL_MATREP_EXPORT ti_type : public make_ti<T>::type 
{
    private:
        using base_type     = typename make_ti<T>::type;
        using Type          = dynamic::Type;

    public:
        ti_type(base_type)  {};
        ti_type()           {};
};

inline bool operator==(ti_empty, const dynamic::Type&)  { return false;};
inline bool operator==(const dynamic::Type&, ti_empty)  { return false;};
inline bool operator==(ti_empty, ti_empty )             { return true;};
inline bool operator!=(const dynamic::Type&, ti_empty)  { return true;};
inline bool operator!=(ti_empty, const dynamic::Type&)  { return true;};
inline bool operator!=(ti_empty, ti_empty )             { return false;};

template<> 
struct MATCL_MATREP_EXPORT ti_type<Object> : public make_ti<Object>::type 
{
    private:
        using base_type     = make_ti<Object>::type;

    public:
        ti_type() : base_type(){};
        ti_type(base_type base) :base_type(base){};
};

using ti_int        = ti_type<Integer>;
using ti_float      = ti_type<Float>;
using ti_real       = ti_type<Real>;
using ti_float_compl= ti_type<Float_complex>;
using ti_compl      = ti_type<Complex>;
using ti_object     = ti_type<Object>;

template<class Val> ti_object   ti_object_type();
template<> inline   ti_object   ti_object_type<Integer>()       {return predefined::get_ti_int(); };
template<> inline   ti_object   ti_object_type<Float>()         {return predefined::get_ti_float(); };
template<> inline   ti_object   ti_object_type<Real>()          {return predefined::get_ti_real(); };
template<> inline   ti_object   ti_object_type<Complex>()       {return predefined::get_ti_complex(); };
template<> inline   ti_object   ti_object_type<Float_complex>() {return predefined::get_ti_float_complex(); };

// object value code is not allowed
MATCL_MATREP_EXPORT ti_object  ti_object_type(value_code vc);

//----------------------------------------------------------------------
//                          get_ti_type
//----------------------------------------------------------------------
template<class T>   struct get_ti_type                  {  };
template<>          struct get_ti_type<Integer>         { using type = ti_type<Integer>; };
template<>          struct get_ti_type<Float>           { using type = ti_type<Float>; };
template<>          struct get_ti_type<Real>            { using type = ti_type<Real>; };
template<>          struct get_ti_type<Float_complex>   { using type = ti_type<Float_complex>; };
template<>          struct get_ti_type<Complex>         { using type = ti_type<Complex>; };
template<>          struct get_ti_type<Object>          { using type = ti_type<Object>; };

template<class V, class S>  
struct get_ti_type<matcl::raw::Matrix<V,S> >            { using type = typename get_ti_type<V>::type; };

//----------------------------------------------------------------------
//                          operations
//----------------------------------------------------------------------
template<class TI2>
bool has_trivial_assignment(ti_object lsh, TI2 rhs);

template<class TI_ret, class TI1, class TI2>
TI_ret unify_ti(TI1 t1, TI2 t2);

template<class TI_ret, class TI1, class TI2>
TI_ret unify_ti_assign(TI1 t1, TI2 t2);

template<class T> 
MATCL_MATREP_EXPORT auto get_ti(const T& m) -> typename get_ti_type<T>::type;

template<class T>
MATCL_MATREP_EXPORT auto get_object_ti() -> typename get_ti_type<T>::type;

template<class TI_ret, class TI1, class TI2>
TI_ret get_return_ti(const matcl::dynamic::function_name& fn, TI1 t1, TI2 t2);

template<class TI_ret, class TI1>
TI_ret get_return_ti(const matcl::dynamic::function_name& fn, TI1 t1);

template<class val_type>
MATCL_MATREP_EXPORT ti_object convert_ti_object(ti_empty t);

template<class val_type>
inline ti_object            convert_ti_object(ti_object t)  { return t; };

MATCL_MATREP_EXPORT std::ostream&  operator<<(std::ostream& os, ti_object);
MATCL_MATREP_EXPORT std::istream&  operator>>(std::istream& is, ti_object&);

};};

#include "matcl-matrep/objects/details/type_info_object.inl"
