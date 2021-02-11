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

#include "matcl-specfunc/func/raw/complex_func.h"
#include "matcl-internals/func/converter.h"
#include "matcl-matrep/lib_functions/func_unary.h"

namespace matcl { namespace details 
{

// Abramowitz and Stegun: Eq. (7.1.14) gives this continued fraction for erfc(z)
// The continued fraction is true providing real(z) > 0.
static std::complex<double> erf_continued_fraction(const std::complex<double>& z)
{
    const double inv_sqrtpi = 1./::sqrt(constants::pi());
    const double eps        = constants::eps();
    double a                = 0.;

    std::complex<double> f(z), C(z), D(0.), delta;    

    do 
    {
        a       += .5;
        D       = z + a*D;
        C       = z + a/C;

        if (D.real() == 0. && D.imag() == 0.)
            D = constants::min_real();

        D       = 1./D;
        delta   = C*D;
        f       *= delta;
    } 
    while (std::abs(1.-delta) > eps);

    f           = 1./f;
    f           *= std::exp(-z*z)*inv_sqrtpi;

    return f;
}

// Abramawitz and Stegun: Eq. (7.1.5) gives a series for erf(z) valid
// for all z, faster convergence for small abs(z)
static std::complex<double> erf_series(const std::complex<double>& z)
{
    const double inv_2_sqrtpi   = 2./::sqrt(constants::pi());
    const double min_real       = constants::min_real();

    const std::complex<double> z2(z*z);

    std::complex<double> sum(0.), term(z);    

    for (int n = 0; n < 3 || std::abs(term) > std::abs(sum)*min_real; ++n)
    {
        sum     += term/static_cast<double>(2*n + 1);
        term    *= -z2/static_cast<double>(n + 1);
    }

    return sum * inv_2_sqrtpi;
}

// Rybicki method adapted from Numerical Recipes
static std::complex<double> erf_rybicki(const std::complex<double>& z)
{
    const double h  = .2;
    int n0          = 2*static_cast<int>(z.imag()/(2.*h) + .5);

    std::complex<double> z0(0., static_cast<double>(n0)*h);
    std::complex<double> zp(z - z0);
    std::complex<double> sum(0., 0.);

    for(int np = -35; np <= 35; np += 2) 
    {
        std::complex<double> t(zp.real(), zp.imag()-static_cast<double>(np)*h);
        std::complex<double> b(std::exp(t*t) / static_cast<double>(np + n0));

        sum     += b; 
    }
    
    sum         *= 2.0 * std::exp(-z*z) / constants::pi();

    return std::complex<double>(-sum.imag(), sum.real());
}

Float_complex scal_func::impl_erf(const Float_complex &z)
{
    Complex zc  = mr::converter<Complex,Float_complex>::eval(z);
    Complex ret = scal_func::impl_erf(zc);
    return mr::converter<Float_complex,Complex>::eval(ret);
};

Complex scal_func::impl_erf(const Complex &x)
{
    using std_complex   = std::complex<double>;

    std_complex return_value;
    std_complex z(x.value);

    const std_complex nan_result(constants::nan(), constants::nan());

    const double re         = real(x);
    const double im         = imag(x);
    const double arg_norm   = std::norm(z);

    // If argument is small, then Taylor expansion is used:
    if(arg_norm < 4.) 
        return Complex(erf_series(z));

    // Cutoff of calculations if Inf's occur:
    if(mrd::scal_func::isinf(x))
    {
        if(mrd::scal_func::isinf(re) && im == 0.)
        {
            // Results for (+/-Inf,0) are hardcoded, just in case that algorithm fails there.
            if(re > 0.) 
                return Complex(1., 0.);
            else
                return Complex(-1., 0.);
        }
        else 
        {
            return Complex(nan_result);
        };
    }

    // Cutoff of calculations if NaNs occur:
    if(mrd::scal_func::isnan(arg_norm))
        return Complex(nan_result);

    // If z is close to imaginary axis, Rybicki method is used:
    if(mrd::abs_helper<Real>::eval(re) < .5)
        return Complex(erf_rybicki(z));

    // Otherwise continued fraction expansion is used:

    if(re > 0.) 
        return_value =  1. - erf_continued_fraction( z);
    else
        return_value = -1. + erf_continued_fraction(-z);

    return Complex(return_value);
}

Float_complex scal_func::impl_erfc(const Float_complex &z)
{
    Complex zc  = mr::converter<Complex,Float_complex>::eval(z);
    Complex ret = scal_func::impl_erfc(zc);
    return mr::converter<Float_complex,Complex>::eval(ret);
};
Complex scal_func::impl_erfc(const Complex &x)
{
    using std_complex   = std::complex<double>;

    std_complex return_value;
    std_complex z(x.value);

    const std_complex nan_result(constants::nan(), constants::nan());

    const double re         = real(x);
    const double im         = imag(x);
    const double arg_norm   = std::norm(z);

    // If argument is small, then Taylor expansion is used:
    if(arg_norm < 4.) 
        return Complex(1.0 - erf_series(z));

    // Cutoff of calculations if Inf's occur:
    if(mrd::scal_func::isinf(x))
    {
        if(mrd::scal_func::isinf(re) && im == 0.)
        {
            // Results for (+/-Inf,0) are hardcoded, just in case that algorithm fails there.
            if(re > 0.) 
                return Complex(0.0, 0.);
            else
                return Complex(2.0, 0.);
        }
        else 
        {
            return Complex(nan_result);
        };
    }

    // Cutoff of calculations if NaNs occur:
    if(mrd::scal_func::isnan(arg_norm))
        return Complex(nan_result);

    // If z is close to imaginary axis, Rybicki method is used:
    if(mrd::abs_helper<Real>::eval(re) < .5)
        return Complex(1.0 - erf_rybicki(z));

    // Otherwise continued fraction expansion is used:

    if(re > 0.) 
        return_value =  erf_continued_fraction( z);
    else
        return_value = 2. - erf_continued_fraction(-z);

    return Complex(return_value);
}

}}