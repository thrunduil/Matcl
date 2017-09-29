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

#include "matcl-mp/mp_int.h"
#include "matcl-mp/mp_float.h"
#include "matcl-mp/mp_rational.h"
#include "matcl-mp/mp_complex.h"

#pragma warning(push)
#pragma warning(disable:4127)// conditional expression is constant
#pragma warning(disable:4146)// unary minus operator applied to unsigned 
#pragma warning(disable:4800)// 'int': forcing value to bool 'true
#pragma warning(disable:4244)// conversion from 'mp_limb_t' to 'unsigned long', possible loss of data

    #include "gmpxx.h"
    #include <boost/multiprecision/gmp.hpp>
    #include <boost/multiprecision/mpfr.hpp>

#pragma warning(pop)

namespace matcl { namespace mp { namespace details
{

namespace bm = boost::multiprecision;

using impl_int          = bm::mpz_int;
using impl_float        = mpfr_t;
using impl_rat          = bm::mpq_rational;

using gmp_impl_rat      = __gmp_expr<mpq_t, mpq_t>;
using gmp_impl_int      = __gmp_expr<mpz_t, mpz_t>;

template<class T>
struct get_ptr
{};

template<>
struct get_ptr<mp_int>
{
    using impl_type = impl_int;

    template<Integer Size>
    static impl_type* eval(char (&data)[Size])
    {
        static_assert(Size == sizeof(impl_type), "invalid array size");
        return reinterpret_cast<impl_type*>(data);
    };
    template<Integer Size>
    static const impl_type* eval(const char (&data)[Size])
    {
        static_assert(Size == sizeof(impl_type), "invalid array size");
        return reinterpret_cast<const impl_type*>(data);
    };
};

template<>
struct get_ptr<mp_float>
{
    using impl_type = impl_float;

    template<Integer Size>
    static impl_type* eval(char (&data)[Size])
    {
        static_assert(Size == sizeof(impl_type), "invalid array size");
        return reinterpret_cast<impl_type*>(data);
    };
    template<Integer Size>
    static const impl_type* eval(const char (&data)[Size])
    {
        static_assert(Size == sizeof(impl_type), "invalid array size");
        return reinterpret_cast<const impl_type*>(data);
    };
};

template<>
struct get_ptr<mp_rational>
{
    using impl_type = impl_rat;

    template<Integer Size>
    static impl_type* eval(char (&data)[Size])
    {
        static_assert(Size == sizeof(impl_type), "invalid array size");
        return reinterpret_cast<impl_type*>(data);
    };
    template<Integer Size>
    static const impl_type* eval(const char (&data)[Size])
    {
        static_assert(Size == sizeof(impl_type), "invalid array size");
        return reinterpret_cast<const impl_type*>(data);
    };
};

inline Integer          impl_value(Integer val)         { return val; };
inline Real             impl_value(Real val)            { return val; };
inline const Complex&   impl_value(const Complex&  val) { return val; };

inline const impl_int&      impl_value(const mp_int& val)       { return *get_ptr<mp_int>::eval(val.m_data); }
inline const impl_float&    impl_value(const mp_float& val)     { return *get_ptr<mp_float>::eval(val.m_data); }
inline const impl_rat&      impl_value(const mp_rational& val)  { return *get_ptr<mp_rational>::eval(val.m_data); }

inline impl_int&            impl_value(mp_int& val)             { return *get_ptr<mp_int>::eval(val.m_data); }
inline impl_float&          impl_value(mp_float& val)           { return *get_ptr<mp_float>::eval(val.m_data); }
inline impl_rat&            impl_value(mp_rational& val)        { return *get_ptr<mp_rational>::eval(val.m_data); }

}};};
