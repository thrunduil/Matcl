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

#include "matcl-mp/matcl_mp.h"
#include "matcl-scalar/lib_functions/func_unary.h"

#include <limits>
#include <algorithm>

namespace matcl { namespace seq_helpers
{

template<class T>
struct max_value
{
    static T value(int prec)
    {
        (void)prec;
        return std::numeric_limits<T>::max();
    }
};

template<>
struct max_value<mp_float>
{
    static mp_float value(int prec)
    {
        return constants::mp_max_real(precision(prec));
    }
};

template<class T>
struct min_value
{
    static T value(int prec)
    {
        (void)prec;
        return std::numeric_limits<T>::min();
    }
};

template<>
struct min_value<mp_float>
{
    static mp_float value(int prec)
    {
        return constants::mp_min_real(precision(prec));
    }
};

template<class T>
struct epsilon
{
    static T value(int prec)
    {
        (void)prec;
        return std::numeric_limits<T>::epsilon();
    }
};

template<>
struct epsilon<mp_float>
{
    static mp_float value(int prec)
    {
        return constants::mp_eps(precision(prec));
    }
};

template<class T>
struct inf_value
{
    static T value(int prec)
    {
        (void)prec;
        return std::numeric_limits<T>::infinity();
    }
};

template<>
struct inf_value<mp_float>
{
    static mp_float value(int prec)
    {
        return constants::mp_inf(precision(prec));
    }
};

template<class T>
struct nan_value
{
    static T value(int prec)
    {
        (void)prec;
        return std::numeric_limits<T>::quiet_NaN();
    }
};

template<>
struct nan_value<mp_float>
{
    static mp_float value(int prec)
    {
        return constants::mp_nan(precision(prec));
    }
};

template<class T>
struct prepare_value
{
    static T eval(const T& v, int prec)
    {
        (void)prec;
        return v;
    }
};

template<>
struct prepare_value<mp_float>
{
    static mp_float eval(const mp_float& v, int prec)
    {
        return mp_float(v, precision(prec));
    }
};

template<class T>
inline
T sqr(const T& x)
{
    return x * x;
};

// function assumes, that x, y are exact; result * epsilon/2 is
// the absolute error
template<class Float_type>
Float_type  error_minus(const Float_type& x, const Float_type& y,
                const Float_type& abs_xmy);

template<class Float_type>
Float_type  error_plus(const Float_type& x, const Float_type& y,
                const Float_type& abs_xpy);

template<class Float>
struct cast_double
{
    static double eval(const Float& f)
    {
        return (double)f;
    }
};

template<>
struct cast_double<mp_float>
{
    static double eval(const mp_float& f)
    {
        return f.cast_float();
    }
};

}}

#include "matcl-scalar/details/sequence/accelerators_helpers.inl"