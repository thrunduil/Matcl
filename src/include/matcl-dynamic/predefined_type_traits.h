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
#include "matcl-dynamic/special_types.h"
#include "matcl-core/details/complex_details.h"

namespace matcl { namespace dynamic
{

// object_type_traits for predefined types

template<>
struct MATCL_DYN_EXPORT object_type_traits<bool> : object_type_traits_default
{
    using T = bool;

    static bool         is_one(bool)    { return false;};

    static std::string  to_string(bool t, matcl::details::printer& pr);
    static void         disp(bool t, matcl::details::printer& pr, Integer elem_width,
                            align_type at, Integer);

    static void         load_data(iarchive_impl& ar, T& val, 
                            unsigned int version);
    static void         save_data(oarchive_impl& ar, const T& val, 
                            unsigned int version);
};

template<>
struct MATCL_DYN_EXPORT object_type_traits<Integer> : object_type_traits_default
{
    using T = Integer;

    static const bool has_one = true;

    static bool         is_one(const Integer& t)    { return t == 1;};
    static T            make_one(const T*)          { return T(1); };

    static std::string  to_string(const Integer& t, matcl::details::printer& pr);
    static void         disp(const T& t, matcl::details::printer& pr, Integer elem_width, 
                            align_type at, Integer);

    static void         load_data(iarchive_impl& ar, T& val, 
                            unsigned int version);
    static void	        save_data(oarchive_impl& ar, const T& val, 
                            unsigned int version);
};

template<>
struct MATCL_DYN_EXPORT object_type_traits<Real> : object_type_traits_default
{
    using T = Real;
    static const bool has_one = true;

    static bool         is_one(const T& t)      { return t == 1.;};
    static T            make_one(const T*)      { return T(1); };

    static std::string  to_string(const Real& t, matcl::details::printer& pr);
    static void         disp(const T& t, matcl::details::printer& pr, Integer elem_width, 
                            align_type at, Integer);

    static void         load_data(iarchive_impl& ar, T& val, 
                            unsigned int version);
    static void	        save_data(oarchive_impl& ar, const T& val, 
                            unsigned int version);
};

template<>
struct MATCL_DYN_EXPORT object_type_traits<Float> : object_type_traits_default
{
    using T = Float;
    static const bool has_one = true;

    static bool         is_one(const T& t)      { return t == 1.;};
    static T            make_one(const T*)      { return T(1); };

    static std::string  to_string(const Float& t, matcl::details::printer& pr);
    static void         disp(const T& t, matcl::details::printer& pr, Integer elem_width, 
                            align_type at, Integer);

    static void         load_data(iarchive_impl& ar, T& val, 
                            unsigned int version);
    static void	        save_data(oarchive_impl& ar, const T& val, 
                            unsigned int version);
};

template<>
struct MATCL_DYN_EXPORT object_type_traits<Complex> : object_type_traits_default
{
    using T = Complex;
    static const bool has_one = true;

    static bool         is_one(const T& t)      { return t == Complex(1.);};
    static T            make_one(const T*)      { return T(1); };

    static std::string  to_string(const Complex& t, matcl::details::printer& pr);
    static void         disp(const T& t, matcl::details::printer& pr, Integer elem_width, 
                            align_type at, Integer);

    static bool         read(std::istream&, T& t);
    static void         write(std::ostream&, const T& t);

    static void         load_data(iarchive_impl& ar, T& val, 
                            unsigned int version);
    static void         save_data(oarchive_impl& ar, const T& val, 
                            unsigned int version);
};

template<>
struct MATCL_DYN_EXPORT object_type_traits<Float_complex> : object_type_traits_default
{
    using T = Float_complex;
    static const bool has_one = true;

    static bool         is_one(const T& t)      { return t == Float_complex(1.f);};
    static T            make_one(const T*)      { return T(1); };

    static std::string  to_string(const Float_complex& t, matcl::details::printer& pr);
    static void         disp(const T& t, matcl::details::printer& pr, Integer elem_width, 
                            align_type at, Integer);

    static bool         read(std::istream&, T& t);
    static void         write(std::ostream&, const T& t);

    static void         load_data(iarchive_impl& ar, T& val, 
                            unsigned int version);
    static void         save_data(oarchive_impl& ar, const T& val, 
                            unsigned int version);
};

template<>
struct MATCL_DYN_EXPORT object_type_traits<std::string> : object_type_traits_default
{
    using T = std::string;

    static bool         is_one(const std::string&)  { return false;};

    static std::string  to_string(const std::string& t, matcl::details::printer&)
                                { return t;};
    static void         disp(const T& t, matcl::details::printer& pr, Integer elem_width, 
                            align_type at, Integer);

    static void         load_data(iarchive_impl& ar, T& val,
                            unsigned int version);
    static void         save_data(oarchive_impl& ar, const T& val, 
                            unsigned int version);
};

template<>
struct MATCL_DYN_EXPORT object_type_traits<unit_type> : object_type_traits_default
{
    using T = unit_type;

    static bool         is_one(const T&)                { return false;};
    static bool         is_zero(const T&)               { return true; };
    static std::string  to_string(const T&, matcl::details::printer&)
                                                        { return "unit";};

    static void         disp(const T& t, matcl::details::printer& pr, Integer elem_width, 
                            align_type at, Integer);

    static void         load_data(iarchive_impl& ar, T& val, unsigned int version);
    static void         save_data(oarchive_impl& ar, const T& val,
                            unsigned int version);
};

template<>
struct MATCL_DYN_EXPORT object_type_traits<null_type> : object_type_traits_default
{
    using T = null_type;

    static bool         is_one(const T&)                { return false;};
    static bool         is_zero(const T&)               { return true; };
    static std::string  to_string(const T&, matcl::details::printer&)
                                                        { return "null";};

    static void         disp(const T& t, matcl::details::printer& pr, Integer elem_width, 
                            align_type at, Integer);

    static void         load_data(iarchive_impl& ar, T& val, 
                            unsigned int version);
    static void         save_data(oarchive_impl& ar, const T& val, 
                            unsigned int version);
};

template<>
struct MATCL_DYN_EXPORT object_type_traits<any_type> : object_type_traits_default
{
    using T                         = any_type;
    static const bool is_clonable   = true;

    static bool         is_one(const T&)    { return false;};
    static bool         is_zero(const T&);
};

template<>
struct MATCL_DYN_EXPORT object_type_traits<Type> : object_type_traits_default
{
    using T                         = Type;
    static const bool is_clonable   = false;

    static bool         is_one(const T&)    { return false;};
    static bool         is_zero(const T& t) { return t == Type(); };

    static std::string  to_string(const T&, matcl::details::printer&);

    static void         disp(const T& t, matcl::details::printer& pr, Integer elem_width, 
                            align_type at, Integer);

    static void         load_data(iarchive_impl& ar, T& val, 
                            unsigned int version);
    static void         save_data(oarchive_impl& ar, const T& val, 
                            unsigned int version);
};

template<>
struct MATCL_DYN_EXPORT object_type_traits<object> : object_type_traits_default
{
    using T                         = object;
    static const bool is_clonable   = true;

    static bool         is_zero(const T& t) { return !t; };
};

};};
