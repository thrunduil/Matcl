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

#include "matcl-core/details/scalfunc_real.h"
#include "matcl-core/details/scalfunc_complex.h"
#include "matcl-core/error/exception_classes.h"

#pragma warning(push)
#pragma warning(disable:4127)   //conditional expression is constant

#include <boost/math/special_functions.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/complex.hpp>

#pragma warning(pop)

namespace matcl { namespace raw { namespace details
{
namespace mrd = matcl::raw::details;

namespace scal_func
{
    namespace ignore_errors 
    {
        using boost__ignore_policy
            = ::boost::math::policies::policy 
              <
                ::boost::math::policies::domain_error<::boost::math::policies::ignore_error>,
                ::boost::math::policies::pole_error<::boost::math::policies::ignore_error>,
                ::boost::math::policies::overflow_error<::boost::math::policies::ignore_error>,
                ::boost::math::policies::underflow_error<::boost::math::policies::ignore_error>,
                ::boost::math::policies::rounding_error<::boost::math::policies::ignore_error>,
                ::boost::math::policies::denorm_error<::boost::math::policies::ignore_error>,
                ::boost::math::policies::indeterminate_result_error<::boost::math::policies::ignore_error>,
                ::boost::math::policies::evaluation_error<::boost::math::policies::ignore_error>
             >;
        BOOST_MATH_DECLARE_SPECIAL_FUNCTIONS(boost__ignore_policy)
    }
};

double scal_func::sqrt1pm1(double arg)
{
    return ignore_errors::sqrt1pm1(arg);
}
float scal_func::sqrt1pm1(float arg)
{
    return ignore_errors::sqrt1pm1(arg);
}
Complex scal_func::sqrt1pm1(const Complex& arg)
{
    return md::minus_c(scal_func::sqrt(md::plus_c(1.0,arg)), 1.0);
}
Float_complex scal_func::sqrt1pm1(const Float_complex& arg)
{
    return md::minus_c(scal_func::sqrt(md::plus_c(1.0f,arg)), 1.0f);
}

//--------------------------------------------------------------------
Complex scal_func::cbrt(const Complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::cbrt(real(x));

    throw error::function_not_defined_for_complex("cbrt");
}
Float_complex scal_func::cbrt(const Float_complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::cbrt(real(x));

    throw error::function_not_defined_for_complex("cbrt");
}

//--------------------------------------------------------------------
bool scal_func::signbit(const Complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::signbit(real(x));

    throw error::function_not_defined_for_complex("signbit");
}
bool scal_func::signbit(const Float_complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::signbit(real(x));

    throw error::function_not_defined_for_complex("signbit");
}

//--------------------------------------------------------------------
Integer scal_func::ifloor(const Complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::ifloor(real(x));

    throw error::function_not_defined_for_complex("ifloor");
}
Integer scal_func::ifloor(const Float_complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::ifloor(real(x));

    throw error::function_not_defined_for_complex("ifloor");
}

//--------------------------------------------------------------------
Integer scal_func::iceil(const Complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::iceil(real(x));

    throw error::function_not_defined_for_complex("iceil");
}
Integer scal_func::iceil(const Float_complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::iceil(real(x));

    throw error::function_not_defined_for_complex("iceil");
}

//--------------------------------------------------------------------
Integer scal_func::itrunc(const Complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::itrunc(real(x));

    throw error::function_not_defined_for_complex("itrunc");
}
Integer scal_func::itrunc(const Float_complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::itrunc(real(x));

    throw error::function_not_defined_for_complex("itrunc");
}

//--------------------------------------------------------------------
Integer scal_func::iround(const Complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::iround(real(x));

    throw error::function_not_defined_for_complex("iround");
}
Integer scal_func::iround(const Float_complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::iround(real(x));

    throw error::function_not_defined_for_complex("iround");
}

//--------------------------------------------------------------------
Integer scal_func::isign(const Complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::isign(real(x));

    throw error::function_not_defined_for_complex("isign");
}
Integer scal_func::isign(const Float_complex& x)
{
    if (is_zero(imag(x)))
        return scal_func::isign(real(x));

    throw error::function_not_defined_for_complex("isign");
}

//--------------------------------------------------------------------
Real scal_func::frexp(const Complex& x, Integer& n)
{
    if (is_zero(imag(x)))
        return scal_func::frexp(real(x), n);

    throw error::function_not_defined_for_complex("frexp");
}
Float scal_func::frexp(const Float_complex& x, Integer& n)
{
    if (is_zero(imag(x)))
        return scal_func::frexp(real(x), n);

    throw error::function_not_defined_for_complex("frexp");
}

//--------------------------------------------------------------------
Real scal_func::modf(const Complex& x, Real& n)
{
    if (is_zero(imag(x)))
        return scal_func::modf(real(x), n);

    throw error::function_not_defined_for_complex("modf");
}
Float scal_func::modf(const Float_complex& x, Float& n)
{
    if (is_zero(imag(x)))
        return scal_func::modf(real(x), n);

    throw error::function_not_defined_for_complex("modf");
}

//--------------------------------------------------------------------
double scal_func::asinh(double x)
{
    return ignore_errors::asinh(x);
}
float scal_func::asinh(float x)
{
    return ignore_errors::asinh(x);
}

//--------------------------------------------------------------------
double scal_func::acosh(double x)
{
    return ignore_errors::acosh(x);
}
float scal_func::acosh(float x)
{
    return ignore_errors::acosh(x);
}

//--------------------------------------------------------------------
double scal_func::atanh(double x)
{
    return ignore_errors::atanh(x);
}
float scal_func::atanh(float x)
{
    return ignore_errors::atanh(x);
}

//--------------------------------------------------------------------
Complex scal_func::asin(const Complex &cm)
{
    // No `asin' for complex in ignore_errors, therefore:
    return Complex(::boost::math::asin(cm.value));
};
Float_complex scal_func::asin(const Float_complex &cm)
{
    // No `asin' for complex in ignore_errors, therefore:
    return Float_complex(::boost::math::asin(cm.value));
};

//--------------------------------------------------------------------
Complex scal_func::acos(const Complex &cm)
{
    // No `asin' for complex in ignore_errors, therefore:
    return Complex(::boost::math::acos(cm.value));
};
Float_complex scal_func::acos(const Float_complex &cm)
{
    // No `asin' for complex in ignore_errors, therefore:
    return Float_complex(::boost::math::acos(cm.value));
};

//--------------------------------------------------------------------
Complex scal_func::atan(const Complex &cm)
{
    // No `atan' for complex in ignore_errors, therefore:
    return Complex(::boost::math::atan(cm.value));
};
Float_complex scal_func::atan(const Float_complex &cm)
{
    // No `atan' for complex in ignore_errors, therefore:
    return Float_complex(::boost::math::atan(cm.value));
};

//--------------------------------------------------------------------
Complex scal_func::asinh(const Complex &cm)
{
    // No `asinh' for complex in ignore_errors, therefore:
    return Complex(::boost::math::asinh(cm.value));
}
Float_complex scal_func::asinh(const Float_complex &cm)
{
    // No `asinh' for complex in ignore_errors, therefore:
    return Float_complex(::boost::math::asinh(cm.value));
}

//--------------------------------------------------------------------
Complex scal_func::atanh(const Complex &cm)
{
    // No `atanh' for complex in ignore_errors, therefore:
    return Complex(::boost::math::atanh(cm.value));
}
Float_complex scal_func::atanh(const Float_complex &cm)
{
    // No `atanh' for complex in ignore_errors, therefore:
    return Float_complex(::boost::math::atanh(cm.value));
}

//--------------------------------------------------------------------
Complex scal_func::acosh(const Complex& cm)
{
    return Complex(::boost::math::acosh(cm.value));
};
Float_complex scal_func::acosh(const Float_complex& cm)
{
    return Float_complex(::boost::math::acosh(cm.value));
};

}}};
