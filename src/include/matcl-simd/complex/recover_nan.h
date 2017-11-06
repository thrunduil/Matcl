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

#include "matcl-simd/config.h"
#include "matcl-core/details/complex_details.h"

namespace matcl { namespace simd { namespace details
{

template<class T>
struct recover_nan_mul
{};

template<class T>
struct recover_nan_div
{};

template<class T>
struct recover_nan_div_rc
{};

template<>
struct recover_nan_mul<double>
{
    using complex_type  = simd_double_complex;
    using matcl_complex = Complex;

    static complex_type eval(const complex_type& x, const complex_type& y,
                                    double r_re, double r_im)
    {
        using impl_type = md::mul_helper<Complex, Complex>;

        if (is_nan(r_re) == false || is_nan(r_im) == false)
            return complex_type(r_re, r_im);

        matcl_complex mx(real(x), imag(x));
        matcl_complex my(real(y), imag(y));

        matcl_complex ret = impl_type::recover_nan(mx, my, r_re, r_im);
        return complex_type(real(ret), imag(ret));
    };
};

template<>
struct recover_nan_mul<float>
{
    using complex_type  = simd_single_complex;
    using matcl_complex = Float_complex;

    static complex_type eval(const complex_type& x, const complex_type& y,
                                    float r_re, float r_im)
    {
        using impl_type = md::mul_helper<Float_complex, Float_complex>;

        if (is_nan(r_re) == false || is_nan(r_im) == false)
            return complex_type(r_re, r_im);

        matcl_complex mx(real(x), imag(x));
        matcl_complex my(real(y), imag(y));

        matcl_complex ret = impl_type::recover_nan(mx, my, r_re, r_im);
        return complex_type(real(ret), imag(ret));

    };
};


template<>
struct recover_nan_div<double>
{
    using complex_type  = simd_double_complex;
    using matcl_complex = Complex;

    static complex_type eval(const complex_type& x, const complex_type& y)
    {
        using impl_type = md::div_helper<Complex, Complex>;

        matcl_complex mx(real(x), imag(x));
        matcl_complex my(real(y), imag(y));

        matcl_complex ret = impl_type::eval(mx, my);
        return complex_type(real(ret), imag(ret));
    };
};

template<>
struct recover_nan_div<float>
{
    using complex_type  = simd_single_complex;
    using matcl_complex = Float_complex;

    static complex_type eval(const complex_type& x, const complex_type& y)
    {
        using impl_type = md::div_helper<Float_complex, Float_complex>;

        matcl_complex mx(real(x), imag(x));
        matcl_complex my(real(y), imag(y));

        matcl_complex ret = impl_type::eval(mx, my);
        return complex_type(real(ret), imag(ret));

    };
};


template<>
struct recover_nan_div_rc<double>
{
    using complex_type  = simd_double_complex;
    using matcl_complex = Complex;

    static complex_type eval(double x, const complex_type& y)
    {
        using impl_type = md::div_helper<double, Complex>;

        matcl_complex my(real(y), imag(y));

        matcl_complex ret = impl_type::eval(x, my);
        return complex_type(real(ret), imag(ret));
    };
};

template<>
struct recover_nan_div_rc<float>
{
    using complex_type  = simd_single_complex;
    using matcl_complex = Float_complex;

    static complex_type eval(float x, const complex_type& y)
    {
        using impl_type = md::div_helper<float, Float_complex>;

        matcl_complex my(real(y), imag(y));

        matcl_complex ret = impl_type::eval(x, my);
        return complex_type(real(ret), imag(ret));

    };
};

}}}