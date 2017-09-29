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

#include "rand_scalars.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/lib_functions/scalar_gen.h"

namespace matcl { namespace test
{

void rand_scalars::make(std::vector<Scalar>& scal, Integer n, bool fractions)
{
    for (Integer i = 0; i < n; ++i)
        scal.push_back(rand_scalar(fractions));
}

Scalar rand_scalars::rand_scalar(bool fractions)
{
    int type = abs(std::rand()) % 5;
    switch (type)
    {
        case 0:     return rand_int();
        case 1:     return rand_float(fractions);
        case 2:     return rand_real(fractions);
        case 3:     return rand_compl(fractions);
        case 4:     return rand_fcompl(fractions);
        default:    return rand_int();
    };
};

Integer rand_scalars::rand_int()
{
    int type = abs(std::rand()) % 10;
    switch (type)
    {
        case 0:     return Integer(0);
        case 1:     return Integer(1);
        case 2:     return Integer(-1);
        case 3:     return Integer(10);
        case 4:     return Integer(-11);
        default:    return Integer(matcl::irand() % 100);
    };
}
Float rand_scalars::rand_float(bool fractions)
{
    int n       = (fractions ? 15 : 0);
    int type    = abs(std::rand()) % (10 + n);
    switch (type)
    {
        case 0:     return Float(0.0f);
        case 1:     return Float(1.0f);
        case 2:     return Float(-1.0f);
        case 3:     return Float(constants::f_nan());
        case 4:     return Float(constants::f_inf());
        case 5:     return Float(-constants::f_inf());
        case 6:     return Float(+0.0f);
        case 7:     return Float(-0.0f);
        case 8:     return Float(7.0f);
        case 9:     return Float(-11.0f);
        case 10:    return Float(7.1f);
        case 11:    return Float(-1.11e4f);
        case 12:    return Float(1.45e-3f);
        case 13:    return Float(-1.45e-3f);
        case 14:    return Float(0.5f);
        default:    return Float(matcl::frandn());
    };
};
Real rand_scalars::rand_real(bool fractions)
{
    int n       = (fractions ? 15 : 0);
    int type    = abs(std::rand()) % (10 + n);
    switch (type)
    {
        case 0:     return Real(0.0);
        case 1:     return Real(1.0);
        case 2:     return Real(-1.0);
        case 3:     return Real(constants::nan());
        case 4:     return Real(constants::inf());
        case 5:     return Real(-constants::inf());
        case 6:     return Real(+0.0);
        case 7:     return Real(-0.0);
        case 8:     return Real(10.0);
        case 9:     return Real(-11.0);
        case 10:    return Float(7.1);
        case 11:    return Float(-1.11e4);
        case 12:    return Float(1.45e-3);
        case 13:    return Float(-1.45e-4);
        case 14:    return Float(-0.5);
        default:    return Real(matcl::randn());
    };
};
Complex rand_scalars::rand_compl(bool fractions)
{
    Complex c(rand_real(fractions), rand_real(fractions));
    return c;
};
Float_complex rand_scalars::rand_fcompl(bool fractions)
{
    Float_complex c(rand_float(fractions), rand_float(fractions));
    return c;
};

//--------------------------------------------------------------
//                  rand_scalars
//--------------------------------------------------------------
void rand_scalars_ext::make(std::vector<Scalar_ext>& scal, Integer n, bool fractions)
{
    for (Integer i = 0; i < n; ++i)
        scal.push_back(rand_scalar(fractions));
}

Scalar_ext rand_scalars_ext::rand_scalar(bool fractions)
{
    int type = abs(std::rand()) % 9;
    switch (type)
    {
        case 0:     return rand_scalars::rand_int();
        case 1:     return rand_scalars::rand_float(fractions);
        case 2:     return rand_scalars::rand_real(fractions);
        case 3:     return rand_scalars::rand_compl(fractions);
        case 4:     return rand_scalars::rand_fcompl(fractions);
        case 5:     return rand_mpint(fractions);
        case 6:     return rand_mpfloat(fractions);
        case 7:     return rand_mpcompl(fractions);
        case 8:     return rand_mprat(fractions);
        default:    return rand_scalars::rand_int();
    };
};

mp_int rand_scalars_ext::rand_mpint(bool fractions)
{
    (void)fractions;
    return mp_int(rand_scalars::rand_int());
}
mp_float rand_scalars_ext::rand_mpfloat(bool fractions)
{
    return mp_float(rand_scalars::rand_float(fractions));
}
mp_complex rand_scalars_ext::rand_mpcompl(bool fractions)
{
    return mp_complex(rand_scalars::rand_compl(fractions));
}
mp_rational rand_scalars_ext::rand_mprat(bool fractions)
{
    if (fractions == false)
        return mp_rational(rand_scalars::rand_int());

    Integer val1 = rand_scalars::rand_int();
    Integer val2 = rand_scalars::rand_int();

    if (val2 == 0)
        val2    = 1;

    return mp_rational(val1, val2);
}

}};
