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

#include "test_scal_accuracy.h"
#include "functions.h"

namespace matcl { namespace details
{

template<class Ty>
struct real_type<std::complex<Ty>>
{
    using type = Ty;
};

}};

namespace matcl { namespace test
{

template<class Ty>
Real min_normal();

template<>
Real min_normal<double>()
{
    return std::numeric_limits<double>::min();
};

template<>
Real min_normal<float>()
{
    return std::numeric_limits<float>::min();
};

template<>
Real min_normal<Complex>()
{
    return std::numeric_limits<double>::min();
};

template<>
Real min_normal<Float_complex>()
{
    return std::numeric_limits<float>::min();
};

template<class Type, class Func>
struct value_allowed
{
    static bool test(Real v)
    {
        (void)v;
        return true;
    };
};

template<>
precision scal_accuracy_tester::prec_type<Real>()
{
    return precision(std::numeric_limits<double>::digits);
}

template<>
precision scal_accuracy_tester::prec_type<Float>()
{
    return precision(std::numeric_limits<float>::digits);
}

template<>
precision scal_accuracy_tester::prec_type<Float_complex>()
{
    return precision(std::numeric_limits<float>::digits);
}

template<>
precision scal_accuracy_tester::prec_type<std::complex<float>>()
{
    return precision(std::numeric_limits<float>::digits);
}

template<>
precision scal_accuracy_tester::prec_type<std::complex<double>>()
{
    return precision(std::numeric_limits<double>::digits);
}

template<>
precision scal_accuracy_tester::prec_type<Complex>()
{
    return precision(std::numeric_limits<double>::digits);
}

template<class Type, class Func>
struct value_allowed_real
{
    static bool test(Real v)
    {
        (void)v;
        return true;
    };
};

template<class Type, class Func>
struct value_allowed_imag
{
    static bool test(Real v)
    {
        (void)v;
        return true;
    };
};

template<class Type>
struct value_allowed<Type, Func_sin>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed<Type, Func_cos>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed<Type, Func_tan>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed<Type, Func_cot>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed<Type, Func_sec>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed<Type, Func_csc>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_imag<Type, Func_exp>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_imag<Type, Func_expm1>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_real<Type, Func_sin>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_real<Type, Func_cos>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_real<Type, Func_tan>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_real<Type, Func_cot>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_real<Type, Func_sec>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_real<Type, Func_csc>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_imag<Type, Func_sinh>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_imag<Type, Func_cosh>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_imag<Type, Func_tanh>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_imag<Type, Func_coth>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_imag<Type, Func_sech>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct value_allowed_imag<Type, Func_csch>
{
    static bool test(Real v)
    {
        // only allow |v| < 2*pi
        if (v > 6.283185307179586 || v < -6.283185307179586)
            return false;
        else
            return true;
    };
};

template<class Type>
struct type_name;

template<>
struct type_name<Float>
{
    static std::string value()
    {
        return "Float";
    }
};

template<>
struct type_name<Real>
{
    static std::string value()
    {
        return "Real";
    }
};

template<class Type>
struct complex_type_name;

template<>
struct complex_type_name<Float>
{
    static std::string value()
    {
        return "Float_complex";
    }
};

template<>
struct complex_type_name<Real>
{
    static std::string value()
    {
        return "Complex";
    }
};

template<class Ty>
struct std_complex_type;

template<>
struct std_complex_type<Float_complex>
{
    using type = std::complex<float>;
};

template<>
struct std_complex_type<Complex>
{
    using type = std::complex<double>;
};

std::complex<float> make_std_complex(const Float_complex& v)
{
    return v.value;
}
std::complex<double> make_std_complex(const Complex& v)
{
    return v.value;
}

}};
