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
#include "matcl-matrep/lib_functions/func_unary.h"

namespace matcl { namespace details 
{

Float_complex scal_func::impl_gammaln(const Float_complex &z)
{   
    Complex zc  = Complex(z);
    Complex ret = scal_func::impl_gammaln(zc);
    return Float_complex(ret);
};

Complex scal_func::impl_gammaln(const Complex &z)
{   
    //based on cgamma(x,y,kf) function by E. Cojocaru, January 2009
    //available on Matlab file exchange, bsd licence

    using std_complex = ::std::complex<double>;
    std_complex g, x0y;

    std_complex xy  = z.value;
    double re       = real(z);
    
    // Bernoulli's numbers B2n divided by [2*n(2*n-1)], n = 1,2,...
    const static double coefs[] = {8.333333333333333e-2, -2.777777777777778e-3,
                                   7.936507936507937e-4, -5.952380952380952e-4,
                                   8.417508417508418e-4, -1.917526917526918e-3,
                                   6.410256410256410e-3, -2.955065359477124e-2,
                                   1.796443723688307e-1, -1.39243221690590 };

    if(std::trunc(re) == re  && re <= 0. && imag(z) == 0.)
    {
        // For negative real integers gamma = INF
        return Complex(constants::inf());
    }
    else if(mrd::scal_func::isnan(z))
    {
        return Complex(constants::nan(), constants::nan());
    };

    bool re_lt_0    = re < 0.;

    if(re_lt_0) 
        xy = -xy;

    double x0       = xy.real();
    bool re_lt_7    = x0 < 7.;
    int na          = 0;

    if(re_lt_7)
    {
        na          = static_cast<int>(7.-x0);
        x0          += static_cast<double>(na);
    }

    double pi   = constants::pi();

    x0y         = std_complex(x0, xy.imag());
    g           = (x0y-.5) * std::log(x0y) - x0y + .5 * std::log(2. * pi);

    for(int k = 0; k < 10; ++k)
        g       += coefs[k] * std::pow(x0y, -1-2*k);

    if(re_lt_7)
    {
        std_complex g1(0.);

        for(int j = 0; j < na; ++j)
        {
            g1  += log(xy + static_cast<double>(j));
        }

        g       -= g1;
    }

    if(re_lt_0)
    {
        g = std::log(pi) - std::log(- std::sin(pi * xy) * xy) - g;
    }

    return Complex(g);
}

Float_complex scal_func::impl_gamma(const Float_complex &z)
{
    Complex zc  = Complex(z);
    Complex ret = scal_func::impl_gamma(zc);
    return Float_complex(ret);
};

Complex scal_func::impl_gamma(const Complex &z)
{
    if(real(z) > 171.) // Biggest representable double ~= 171!
    {
        if(imag(z) != 0.)
            return Complex(constants::inf(), constants::inf());
        else
            return Complex(constants::inf());
    };

    return mrd::scal_func::exp(impl_gammaln(z));
}

}}