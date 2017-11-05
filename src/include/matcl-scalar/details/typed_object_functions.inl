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

#include "matcl-scalar/objects/typed_object_functions.h"

#pragma warning(push)
#pragma warning(disable:4800) //forcing value to bool 'true' or 'false' (performance warning)

namespace matcl { namespace unqualified_call
{

//evaluators;
//evaluations must be performed in matcl namespace

template<class Ret, class T>
struct eval_conj
{
    static Ret eval(const T& a)
    {
        return Ret(conj(a));
    }
};

template<class Ret, class T>
struct eval_arg
{
    static Ret eval(const T& a)
    {
        return Ret(arg(a));
    }
};

template<class Ret, class T>
struct eval_abs
{
    static Ret eval(const T& a)
    {
        return Ret(abs(a));
    }
};

template<class Ret, class T>
struct eval_abs2
{
    static Ret eval(const T& a)
    {
        return Ret(abs2(a));
    }
};

template<class Ret, class T>
struct eval_inv
{
    static Ret eval(const T& a)
    {
        return Ret( inv(a) );
    }
};

template<class Ret, class T>
struct eval_invs
{
    static Ret eval(const T& a)
    {
        return Ret( invs(a) );
    }
};

template<class Ret, class T, class S>
struct eval_div_0
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(div_0(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_div_1
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(div_1(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_pow
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(pow(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_powm1
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(powm1(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_pow_c
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(pow_c(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_eeq_nan
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(eeq_nan(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_neq_nan
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(neq_nan(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_elem_mul
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(mul(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_atan2
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(atan2(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_mod
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(mod(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_rem
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(rem(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_hypot
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(hypot(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_copysign
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(copysign(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_nextafter
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(nextafter(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_float_distance
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(float_distance(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_min
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(min(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_max
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(max(a,b));
    }
};

template<class Ret, class T, class S>
struct eval_op_and
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(a && b);
    }
};

template<class Ret, class T, class S>
struct eval_op_or
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(a || b);
    }
};

template<class Ret, class T, class S>
struct eval_op_xor
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(op_xor(a, b));
    }
};

template<class Ret, class T, class S>
struct eval_elem_or
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(elem_or(a, b));
    }
};

template<class Ret, class T, class S>
struct eval_elem_and
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(elem_and(a, b));
    }
};

template<class Ret, class T, class S>
struct eval_elem_xor
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(elem_xor(a, b));
    }
};

template<class Ret, class T>
struct eval_isnan
{
    static Ret eval(const T& a)
    {
        return Ret(is_nan(a));
    };
};

template<class Ret, class T>
struct eval_isinf
{
    static Ret eval(const T& a)
    {
        return Ret(is_inf(a));
    };
};

template<class Ret, class T>
struct eval_isfinite
{
    static Ret eval(const T& a)
    {
        return Ret(is_finite(a));
    };
};

template<class Ret, class T>
struct eval_isregular
{
    static Ret eval(const T& a)
    {
        return Ret(is_regular(a));
    };
};

template<class Ret, class T>
struct eval_isnormal
{
    static Ret eval(const T& a)
    {
        return Ret(is_normal(a));
    };
};

template<class Ret, class T>
struct eval_isint
{
    static Ret eval(const T& a)
    {
        return Ret(is_int(a));
    };
};

template<class Ret, class T>
struct eval_isreal
{
    static Ret eval(const T& a)
    {
        return Ret(is_real(a));
    };
};

template<class Ret, class T>
struct eval_neg
{
    static Ret eval(const T& a)
    {
        return Ret(neg(a));
    };
};

template<class Ret, class T>
struct eval_is_true
{
    static Ret eval(const T& a)
    {
        return Ret(is_true(a));
    };
};

template<class T>
struct eval_fpclassify
{
    static fp_type eval(const T& a)
    {
        Integer res = (Integer)fpclassify(a);
        return mrd::scal_func::int_to_fptype(res);
    };
};

template<class Ret, class T>
struct eval_nextabove
{
    static Ret eval(const T& a)
    {
        return Ret(nextabove(a));
    }
};

template<class Ret, class T>
struct eval_nextbelow
{
    static Ret eval(const T& a)
    {
        return Ret(nextbelow(a));
    }
};

template<class Ret, class T>
struct eval_eps
{
    static Ret eval(const T& a)
    {
        return Ret(eps(a));
    }
};

template<class Ret, class T>
struct eval_sqrt
{
    static Ret eval(const T& a)
    {
        return Ret(sqrt(a));
    }
};

template<class Ret, class T>
struct eval_cbrt
{
    static Ret eval(const T& a)
    {
        return Ret(cbrt(a));
    }
};

template<class Ret, class T>
struct eval_sqrt_c
{
    static Ret eval(const T& a)
    {
        return Ret(sqrt_c(a));
    }
};

template<class Ret, class T>
struct eval_sqrt1pm1
{
    static Ret eval(const T& a)
    {
        return Ret(sqrt1pm1(a));
    }
};

template<class Ret, class T>
struct eval_sqrt1pm1_c
{
    static Ret eval(const T& a)
    {
        return Ret(sqrt1pm1_c(a));
    }
};

template<class Ret, class T>
struct eval_exp
{
    static Ret eval(const T& a)
    {
        return Ret(exp(a));
    }
};

template<class Ret, class T>
struct eval_expm1
{
    static Ret eval(const T& a)
    {
        return Ret(expm1(a));
    }
};

template<class Ret, class T>
struct eval_expi
{
    static Ret eval(const T& a)
    {
        return Ret(expi(a));
    }
};

template<class Ret, class T>
struct eval_exp2
{
    static Ret eval(const T& a)
    {
        return Ret(exp2(a));
    }
};

template<class Ret, class T>
struct eval_exp10
{
    static Ret eval(const T& a)
    {
        return Ret(exp10(a));
    }
};

template<class Ret, class T>
struct eval_log
{
    static Ret eval(const T& a)
    {
        return Ret(log(a));
    }
};

template<class Ret, class T>
struct eval_log_c
{
    static Ret eval(const T& a)
    {
        return Ret(log_c(a));
    }
};

template<class Ret, class T>
struct eval_log1p
{
    static Ret eval(const T& a)
    {
        return Ret(log1p(a));
    }
};

template<class Ret, class T>
struct eval_log1p_c
{
    static Ret eval(const T& a)
    {
        return Ret(log1p_c(a));
    }
};

template<class Ret, class T>
struct eval_log2
{
    static Ret eval(const T& a)
    {
        return Ret(log2(a));
    }
};

template<class Ret, class T>
struct eval_log2_c
{
    static Ret eval(const T& a)
    {
        return Ret(log2_c(a));
    }
};

template<class Ret, class T>
struct eval_log10
{
    static Ret eval(const T& a)
    {
        return Ret(log10(a));
    }
};

template<class Ret, class T>
struct eval_log10_c
{
    static Ret eval(const T& a)
    {
        return Ret(log10_c(a));
    }
};

template<class Ret, class T>
struct eval_sin
{
    static Ret eval(const T& a)
    {
        return Ret(sin(a));
    }
};

template<class Ret, class T>
struct eval_cos
{
    static Ret eval(const T& a)
    {
        return Ret(cos(a));
    }
};

template<class Ret, class T>
struct eval_tan
{
    static Ret eval(const T& a)
    {
        return Ret(tan(a));
    }
};

template<class Ret, class T>
struct eval_cot
{
    static Ret eval(const T& a)
    {
        return Ret(cot(a));
    }
};

template<class Ret, class T>
struct eval_sec
{
    static Ret eval(const T& a)
    {
        return Ret(sec(a));
    }
};

template<class Ret, class T>
struct eval_csc
{
    static Ret eval(const T& a)
    {
        return Ret(csc(a));
    }
};

template<class Ret, class T>
struct eval_sinh
{
    static Ret eval(const T& a)
    {
        return Ret(sinh(a));
    }
};

template<class Ret, class T>
struct eval_cosh
{
    static Ret eval(const T& a)
    {
        return Ret(cosh(a));
    }
};

template<class Ret, class T>
struct eval_tanh
{
    static Ret eval(const T& a)
    {
        return Ret(tanh(a));
    }
};

template<class Ret, class T>
struct eval_coth
{
    static Ret eval(const T& a)
    {
        return Ret(coth(a));
    }
};

template<class Ret, class T>
struct eval_sech
{
    static Ret eval(const T& a)
    {
        return Ret(sech(a));
    }
};

template<class Ret, class T>
struct eval_csch
{
    static Ret eval(const T& a)
    {
        return Ret(csch(a));
    }
};

template<class Ret, class T>
struct eval_asin
{
    static Ret eval(const T& a)
    {
        return Ret(asin(a));
    }
};

template<class Ret, class T>
struct eval_asin_c
{
    static Ret eval(const T& a)
    {
        return Ret(asin_c(a));
    }
};

template<class Ret, class T>
struct eval_acos
{
    static Ret eval(const T& a)
    {
        return Ret(acos(a));
    }
};

template<class Ret, class T>
struct eval_acos_c
{
    static Ret eval(const T& a)
    {
        return Ret(acos_c(a));
    }
};

template<class Ret, class T>
struct eval_atan
{
    static Ret eval(const T& a)
    {
        return Ret(atan(a));
    }
};

template<class Ret, class T>
struct eval_acot
{
    static Ret eval(const T& a)
    {
        return Ret(acot(a));
    }
};

template<class Ret, class T>
struct eval_asec
{
    static Ret eval(const T& a)
    {
        return Ret(asec(a));
    }
};

template<class Ret, class T>
struct eval_asec_c
{
    static Ret eval(const T& a)
    {
        return Ret(asec_c(a));
    }
};

template<class Ret, class T>
struct eval_acsc
{
    static Ret eval(const T& a)
    {
        return Ret(acsc(a));
    }
};

template<class Ret, class T>
struct eval_acsc_c
{
    static Ret eval(const T& a)
    {
        return Ret(acsc_c(a));
    }
};

template<class Ret, class T>
struct eval_asinh
{
    static Ret eval(const T& a)
    {
        return Ret(asinh(a));
    }
};

template<class Ret, class T>
struct eval_acosh
{
    static Ret eval(const T& a)
    {
        return Ret(acosh(a));
    }
};

template<class Ret, class T>
struct eval_acosh_c
{
    static Ret eval(const T& a)
    {
        return Ret(acosh_c(a));
    }
};

template<class Ret, class T>
struct eval_atanh
{
    static Ret eval(const T& a)
    {
        return Ret(atanh(a));
    }
};

template<class Ret, class T>
struct eval_acoth
{
    static Ret eval(const T& a)
    {
        return Ret(acoth(a));
    }
};

template<class Ret, class T>
struct eval_acoth_c
{
    static Ret eval(const T& a)
    {
        return Ret(acoth_c(a));
    }
};

template<class Ret, class T>
struct eval_asech
{
    static Ret eval(const T& a)
    {
        return Ret(asech(a));
    }
};

template<class Ret, class T>
struct eval_asech_c
{
    static Ret eval(const T& a)
    {
        return Ret(asech_c(a));
    }
};

template<class Ret, class T>
struct eval_acsch
{
    static Ret eval(const T& a)
    {
        return Ret(acsch(a));
    }
};

template<class Ret, class T>
struct eval_atanh_c
{
    static Ret eval(const T& a)
    {
        return Ret(atanh_c(a));
    }
};

template<class Ret, class T>
struct eval_floor
{
    static Ret eval(const T& a)
    {
        return Ret(floor(a));
    }
};

template<class Ret, class T>
struct eval_ceil
{
    static Ret eval(const T& a)
    {
        return Ret(ceil(a));
    }
};

template<class Ret, class T>
struct eval_trunc
{
    static Ret eval(const T& a)
    {
        return Ret(trunc(a));
    }
};

template<class Ret, class T>
struct eval_round
{
    static Ret eval(const T& a)
    {
        return Ret(round(a));
    }
};

template<class Ret, class T>
struct eval_ifloor
{
    static Ret eval(const T& a)
    {
        return Ret(ifloor(a));
    }
};

template<class Ret, class T>
struct eval_iceil
{
    static Ret eval(const T& a)
    {
        return Ret(iceil(a));
    }
};

template<class Ret, class T>
struct eval_iround
{
    static Ret eval(const T& a)
    {
        return Ret(iround(a));
    }
};

template<class Ret, class T>
struct eval_itrunc
{
    static Ret eval(const T& a)
    {
        return Ret(itrunc(a));
    }
};

template<class Ret, class T>
struct eval_signbit
{
    static Ret eval(const T& a)
    {
        return Ret(signbit(a));
    }
};

template<class Ret, class T>
struct eval_sign
{
    static Ret eval(const T& a)
    {
        return Ret(sign(a));
    }
};

template<class Ret, class T>
struct eval_isign
{
    static Ret eval(const T& a)
    {
        return Ret(isign(a));
    }
};

template<class Ret, class T>
struct eval_ldexp
{
    static Ret eval(const T& a, Integer n)
    {
        return Ret(ldexp(a, n));
    }
};

template<class Ret, class T>
struct eval_scalbn
{
    static Ret eval(const T& a, Integer n)
    {
        return Ret(scalbn(a, n));
    }
};

template<class Ret, class T, class S>
struct eval_frexp
{
    static Ret eval(const T& a, S& n)
    {
        return Ret(frexp(a, n));
    }
};

template<class Ret, class T, class S>
struct eval_modf
{
    static Ret eval(const T& a, S& n)
    {
        return Ret(modf(a, n));
    }
};

template<class Ret, class T, class S>
struct eval_fdim
{
    static Ret eval(const T& a, const S& n)
    {
        return Ret(fdim(a, n));
    }
};

template<class Ret, class T>
struct eval_logb
{
    static Ret eval(const T& a)
    {
        return Ret(logb(a));
    }
};

template<class Ret, class T>
struct eval_ilogb
{
    static Ret eval(const T& a)
    {
        return Ret(ilogb(a));
    }
};

template<class Ret, class T>
struct eval_rising_factorial
{
    static Ret eval(const T& a, Integer i)
    {
        return Ret(rising_factorial(a,i));
    };
};

template<class Ret, class T>
struct eval_falling_factorial
{
    static Ret eval(const T& a, Integer i)
    {
        return Ret(falling_factorial(a,i));
    };
};

template<class Ret, class T>
struct eval_bernoulli_b2n
{
    static Ret eval(Integer i)
    {
        return Ret(bernoulli_b2n<T>(i));
    };
};

template<class Ret, class T>
struct eval_max_bernoulli_b2n
{
    static Ret eval()
    {
        return Ret(max_bernoulli_b2n<T>());
    };
};

template<class Ret, class T>
struct eval_factorial
{
    static Ret eval(Integer i)
    {
        return Ret(factorial<T>(i));
    };
};

template<class Ret, class T>
struct eval_double_factorial
{
    static Ret eval(Integer i)
    {
        return Ret(double_factorial<T>(i));
    };
};

template<class Ret, class T>
struct eval_binomial_coefficient
{
    static Ret eval(Integer n, Integer k)
    {
        return Ret(binomial_coefficient<T>(n,k));
    };
};

}};

namespace matcl { namespace dynamic
{

template<class T>
typename result_of::result_of_arg<T>::type_object
dynamic::arg(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_arg<T>::type_object;
    return matcl::unqualified_call::eval_arg<return_type, T>::eval(a.get());
};

template<class T>
typename result_of::result_of_arg<T>::type_object
dynamic::angle(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_arg<T>::type_object;
    return matcl::unqualified_call::eval_arg<return_type, T>::eval(a.get());
};

template<class T>
typename result_of::result_of_conj<T>::type_object
dynamic::conj(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_conj<T>::type_object;
    return matcl::unqualified_call::eval_conj<return_type, T>::eval(a.get());
};

template<class T>
typename result_of::result_of_abs<T>::type_object
dynamic::abs(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_abs<T>::type_object;
    return matcl::unqualified_call::eval_abs<return_type, T>::eval(a.get());
};

template<class T>
typename result_of::result_of_abs2<T>::type_object
dynamic::abs2(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_abs2<T>::type_object;
    return matcl::unqualified_call::eval_abs2<return_type, T>::eval(a.get());
};

template<class T, class Ret>
Ret dynamic::is_nan(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_isnan<Ret, T>::eval(a.get());
};

template<class T, class Ret>
Ret dynamic::is_finite(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_isfinite<Ret,T>::eval(a.get());
};

template<class T, class Ret>
Ret dynamic::is_inf(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_isinf<Ret,T>::eval(a.get());
};

template<class T, class Ret>
Ret dynamic::is_regular(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_isregular<Ret,T>::eval(a.get());
};

template<class T, class Ret>
Ret dynamic::is_normal(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_isnormal<Ret,T>::eval(a.get());
};

template<class T, class Ret>
Ret dynamic::is_int(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_isint<Ret,T>::eval(a.get());
};

template<class T, class Ret>
Ret dynamic::is_real(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_isreal<Ret,T>::eval(a.get());
};

template<class T, class Ret>
Ret dynamic::neg(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_neg<Ret,T>::eval(a.get());
};

template<class T, class Ret>
Ret dynamic::is_true(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_is_true<Ret,T>::eval(a.get());
};

// cast to boolean value; equivalent to cast_bool
template<class T, class Res = typename result_of::result_of_op_true<T>::type_object>
Res                     is_true(const object_type<T>& x);

template<class T, class Enable>
fp_type dynamic::fpclassify(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_fpclassify<T>::eval(a.get());    
};

template<class T>
typename result_of::result_of_nextabove<T>::type_object
dynamic::nextabove(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_nextabove<T>::type_object;
    return matcl::unqualified_call::eval_nextabove<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_nextbelow<T>::type_object
dynamic::nextbelow(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_nextbelow<T>::type_object;
    return matcl::unqualified_call::eval_nextbelow<return_type, T>::eval(a.get());
}

template<class T, class Ret>
Ret dynamic::signbit(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_signbit<Ret,T>::eval(a.get());
};

template<class T, class Ret>
Ret dynamic::isign(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_isign<Ret,T>::eval(a.get());
};

template<class T, class Ret>
Ret dynamic::sign(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_sign<Ret, T>::eval(a.get());
};

template<class T, class Ret>
Ret dynamic::ldexp(const object_type<T>& a, Integer exp)
{
    return matcl::unqualified_call::eval_ldexp<Ret,T>::eval(a.get(), exp);
}

template<class T, class Ret>
Ret dynamic::ldexp(const object_type<T>& a, const object_type<Integer>& exp)
{
    return matcl::unqualified_call::eval_ldexp<Ret,T>::eval(a.get(), exp.get);
}

template<class T, class Ret>
Ret dynamic::scalbn(const object_type<T>& a, Integer exp)
{
    return matcl::unqualified_call::eval_scalbn<Ret,T>::eval(a.get(), exp);
}

template<class T, class Ret>
Ret dynamic::scalbn(const object_type<T>& a, const object_type<Integer>& exp)
{
    return matcl::unqualified_call::eval_scalbn<Ret,T>::eval(a.get(), exp.get);
}

template<class T, class S, class Ret>
Ret dynamic::frexp(const object_type<T>& a, object_type<S>& exp)
{
    return matcl::unqualified_call::eval_frexp<Ret,T, S>::eval(a.get(), exp.get_unique());
}

template<class T, class S, class Ret>
Ret dynamic::frexp(const object_type<T>& a, S& exp)
{
    return matcl::unqualified_call::eval_frexp<Ret,T, S>::eval(a.get(), exp);
}

template<class T, class S, class Ret>
Ret dynamic::frexp(const T& a, object_type<S>& exp)
{
    return matcl::unqualified_call::eval_frexp<Ret,T, S>::eval(a, exp.get_unique());
}

template<class T, class S, class Ret>
Ret dynamic::modf(const object_type<T>& a, object_type<S>& exp)
{
    return matcl::unqualified_call::eval_modf<Ret,T, S>::eval(a.get(), exp.get_unique());
}

template<class T, class S, class Ret>
Ret dynamic::modf(const object_type<T>& a, S& exp)
{
    return matcl::unqualified_call::eval_modf<Ret,T, S>::eval(a.get(), exp);
}

template<class T, class S, class Ret>
Ret dynamic::modf(const T& a, object_type<S>& exp)
{
    return matcl::unqualified_call::eval_modf<Ret,T, S>::eval(a, exp.get_unique());
}

template<class T, class S, class Ret>
Ret dynamic::fdim(const object_type<T>& a, const object_type<S>& exp)
{
    return matcl::unqualified_call::eval_fdim<Ret,T, S>::eval(a.get(), exp.get());
}

template<class T, class S, class Ret>
Ret dynamic::fdim(const object_type<T>& a, const S& exp)
{
    return matcl::unqualified_call::eval_fdim<Ret,T, S>::eval(a.get(), exp);
}

template<class T, class S, class Ret>
Ret dynamic::fdim(const T& a, const object_type<S>& exp)
{
    return matcl::unqualified_call::eval_fdim<Ret,T, S>::eval(a, exp.get());
}

template<class T, class Ret>
Ret dynamic::logb(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_logb<Ret,T>::eval(a.get());
}

template<class T, class Ret>
Ret dynamic::ilogb(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_ilogb<Ret,T>::eval(a.get());
};

template<class T>
typename result_of::result_of_eps<T>::type_object
dynamic::eps(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_eps<T>::type_object;
    return matcl::unqualified_call::eval_eps<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_sqrt<T>::type_object
dynamic::sqrt(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_sqrt<T>::type_object;
    return matcl::unqualified_call::eval_sqrt<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_cbrt<T>::type_object
dynamic::cbrt(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_cbrt<T>::type_object;
    return matcl::unqualified_call::eval_cbrt<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_sqrt_c<T>::type_object
dynamic::sqrt_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_sqrt_c<T>::type_object;
    return matcl::unqualified_call::eval_sqrt_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_sqrt1pm1<T>::type_object
dynamic::sqrt1pm1(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_sqrt1pm1<T>::type_object;
    return matcl::unqualified_call::eval_sqrt1pm1<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_sqrt1pm1_c<T>::type_object
dynamic::sqrt1pm1_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_sqrt1pm1_c<T>::type_object;
    return matcl::unqualified_call::eval_sqrt1pm1_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_exp<T>::type_object
dynamic::exp(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_exp<T>::type_object;
    return matcl::unqualified_call::eval_exp<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_expm1<T>::type_object
dynamic::expm1(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_expm1<T>::type_object;
    return matcl::unqualified_call::eval_expm1<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_expi<T>::type_object
dynamic::expi(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_expi<T>::type_object;
    return matcl::unqualified_call::eval_expi<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_exp2<T>::type_object
dynamic::exp2(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_exp2<T>::type_object;
    return matcl::unqualified_call::eval_exp2<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_exp10<T>::type_object
dynamic::exp10(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_exp10<T>::type_object;
    return matcl::unqualified_call::eval_exp10<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_log<T>::type_object
dynamic::log(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_log<T>::type_object;
    return matcl::unqualified_call::eval_log<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_log_c<T>::type_object
dynamic::log_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_log_c<T>::type_object;
    return matcl::unqualified_call::eval_log_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_log1p<T>::type_object
dynamic::log1p(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_log1p<T>::type_object;
    return matcl::unqualified_call::eval_log1p<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_log1p_c<T>::type_object
dynamic::log1p_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_log1p_c<T>::type_object;
    return matcl::unqualified_call::eval_log1p_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_log2<T>::type_object
dynamic::log2(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_log2<T>::type_object;
    return matcl::unqualified_call::eval_log2<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_log2_c<T>::type_object
dynamic::log2_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_log2_c<T>::type_object;
    return matcl::unqualified_call::eval_log2_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_log10<T>::type_object
dynamic::log10(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_log10<T>::type_object;
    return matcl::unqualified_call::eval_log10<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_log10_c<T>::type_object
dynamic::log10_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_log10_c<T>::type_object;
    return matcl::unqualified_call::eval_log10_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_sin<T>::type_object
dynamic::sin(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_sin<T>::type_object;
    return matcl::unqualified_call::eval_sin<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_cos<T>::type_object
dynamic::cos(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_cos<T>::type_object;
    return matcl::unqualified_call::eval_cos<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_tan<T>::type_object
dynamic::tan(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_tan<T>::type_object;
    return matcl::unqualified_call::eval_tan<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_cot<T>::type_object
dynamic::cot(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_cot<T>::type_object;
    return matcl::unqualified_call::eval_cot<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_sec<T>::type_object
dynamic::sec(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_sec<T>::type_object;
    return matcl::unqualified_call::eval_sec<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_csc<T>::type_object
dynamic::csc(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_csc<T>::type_object;
    return matcl::unqualified_call::eval_csc<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_sinh<T>::type_object
dynamic::sinh(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_sinh<T>::type_object;
    return matcl::unqualified_call::eval_sinh<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_cosh<T>::type_object
dynamic::cosh(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_cosh<T>::type_object;
    return matcl::unqualified_call::eval_cosh<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_tanh<T>::type_object
dynamic::tanh(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_tanh<T>::type_object;
    return matcl::unqualified_call::eval_tanh<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_coth<T>::type_object
dynamic::coth(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_coth<T>::type_object;
    return matcl::unqualified_call::eval_coth<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_sech<T>::type_object
dynamic::sech(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_sech<T>::type_object;
    return matcl::unqualified_call::eval_sech<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_csch<T>::type_object
dynamic::csch(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_csch<T>::type_object;
    return matcl::unqualified_call::eval_csch<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_asin<T>::type_object
dynamic::asin(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_asin<T>::type_object;
    return matcl::unqualified_call::eval_asin<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_asin_c<T>::type_object
dynamic::asin_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_asin_c<T>::type_object;
    return matcl::unqualified_call::eval_asin_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_acos<T>::type_object
dynamic::acos(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_acos<T>::type_object;
    return matcl::unqualified_call::eval_acos<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_acos_c<T>::type_object
dynamic::acos_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_acos_c<T>::type_object;
    return matcl::unqualified_call::eval_acos_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_atan<T>::type_object
dynamic::atan(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_atan<T>::type_object;
    return matcl::unqualified_call::eval_atan<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_acot<T>::type_object
dynamic::acot(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_acot<T>::type_object;
    return matcl::unqualified_call::eval_acot<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_asec<T>::type_object
dynamic::asec(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_asec<T>::type_object;
    return matcl::unqualified_call::eval_asec<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_asec_c<T>::type_object
dynamic::asec_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_asec_c<T>::type_object;
    return matcl::unqualified_call::eval_asec_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_acsc<T>::type_object
dynamic::acsc(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_acsc<T>::type_object;
    return matcl::unqualified_call::eval_acsc<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_acsc_c<T>::type_object
dynamic::acsc_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_acsc_c<T>::type_object;
    return matcl::unqualified_call::eval_acsc_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_asinh<T>::type_object
dynamic::asinh(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_asinh<T>::type_object;
    return matcl::unqualified_call::eval_asinh<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_acosh<T>::type_object
dynamic::acosh(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_acosh<T>::type_object;
    return matcl::unqualified_call::eval_acosh<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_acosh_c<T>::type_object
dynamic::acosh_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_acosh_c<T>::type_object;
    return matcl::unqualified_call::eval_acosh_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_atanh<T>::type_object
dynamic::atanh(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_atanh<T>::type_object;
    return matcl::unqualified_call::eval_atanh<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_atanh_c<T>::type_object
dynamic::atanh_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_atanh_c<T>::type_object;
    return matcl::unqualified_call::eval_atanh_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_acoth<T>::type_object
dynamic::acoth(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_acoth<T>::type_object;
    return matcl::unqualified_call::eval_acoth<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_acoth_c<T>::type_object
dynamic::acoth_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_acoth_c<T>::type_object;
    return matcl::unqualified_call::eval_acoth_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_asech<T>::type_object
dynamic::asech(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_asech<T>::type_object;
    return matcl::unqualified_call::eval_asech<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_asech_c<T>::type_object
dynamic::asech_c(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_asech_c<T>::type_object;
    return matcl::unqualified_call::eval_asech_c<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_acsch<T>::type_object
dynamic::acsch(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_acsch<T>::type_object;
    return matcl::unqualified_call::eval_acsch<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_floor<T>::type_object
dynamic::floor(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_floor<T>::type_object;
    return matcl::unqualified_call::eval_floor<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_ceil<T>::type_object
dynamic::ceil(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_ceil<T>::type_object;
    return matcl::unqualified_call::eval_ceil<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_round<T>::type_object
dynamic::round(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_round<T>::type_object;
    return matcl::unqualified_call::eval_round<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_trunc<T>::type_object
dynamic::trunc(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_trunc<T>::type_object;
    return matcl::unqualified_call::eval_trunc<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_ifloor<T>::type_object
dynamic::ifloor(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_ifloor<T>::type_object;
    return matcl::unqualified_call::eval_ifloor<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_iceil<T>::type_object
dynamic::iceil(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_iceil<T>::type_object;
    return matcl::unqualified_call::eval_iceil<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_iround<T>::type_object
dynamic::iround(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_iround<T>::type_object;
    return matcl::unqualified_call::eval_iround<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_itrunc<T>::type_object
dynamic::itrunc(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_itrunc<T>::type_object;
    return matcl::unqualified_call::eval_itrunc<return_type, T>::eval(a.get());
}

template<class T>
typename result_of::result_of_inv<T>::type_object
dynamic::inv(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_inv<T>::type_object;
    return matcl::unqualified_call::eval_inv<return_type, T>::eval(a.get());
};

template<class T>
typename result_of::result_of_invs<T>::type_object
dynamic::invs(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_invs<T>::type_object;
    return matcl::unqualified_call::eval_invs<return_type, T>::eval(a.get());
};

template<class T, class S, class Ret>
Ret dynamic::div(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_div<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::div(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_div<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::div(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_div<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::div_0(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_div_0<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::div_0(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_div_0<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::div_0(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_div_0<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::div_1(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_div_1<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::div_1(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_div_1<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::div_1(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_div_1<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::pow(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_pow<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::pow(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_pow<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::pow(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_pow<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::pow_c(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_pow_c<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::pow_c(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_pow_c<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::pow_c(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_pow_c<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::powm1(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_powm1<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::powm1(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_powm1<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::powm1(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_powm1<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::atan2(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_atan2<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::atan2(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_atan2<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::atan2(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_atan2<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::mod(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_mod<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::mod(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_mod<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::mod(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_mod<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::rem(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_rem<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::rem(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_rem<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::rem(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_rem<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::hypot(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_hypot<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::hypot(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_hypot<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::hypot(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_hypot<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::copysign(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_copysign<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::copysign(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_copysign<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::copysign(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_copysign<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::nextafter(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_nextafter<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::nextafter(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_nextafter<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::nextafter(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_nextafter<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::float_distance(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_float_distance<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::float_distance(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_float_distance<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::float_distance(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_float_distance<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::plus(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_plus<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::plus(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_plus<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::plus(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_plus<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::minus(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_minus<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::minus(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_minus<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::minus(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_minus<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::mul(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_elem_mul<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::mul(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_elem_mul<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::mul(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_elem_mul<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::mmul(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_mul<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::mmul(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_op_mul<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::mmul(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_mul<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::kron(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_mul<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::kron(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_op_mul<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::kron(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_mul<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::eeq(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_eeq<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::eeq(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_eeq<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::eeq(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_eeq<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::eeq_nan(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_eeq_nan<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::eeq_nan(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_eeq_nan<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::eeq_nan(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_eeq_nan<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::neq(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_neq<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::neq(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_neq<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::neq(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_neq<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::neq_nan(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_neq_nan<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::neq_nan(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_neq_nan<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::neq_nan(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_neq_nan<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::geq(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_geq<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::geq(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_geq<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::geq(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_geq<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::leq(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_leq<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::leq(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_leq<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::leq(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_leq<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::gt(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_gt<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::gt(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_gt<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::gt(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_gt<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::lt(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_lt<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::lt(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_lt<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::lt(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_lt<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::min(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_min<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class Ret>
Ret dynamic::min(const object_type<T>& a, const object_type<T>& b)
{
    return matcl::unqualified_call::eval_min<Ret,T,T>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::min(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_min<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::min(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_min<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::max(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_max<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class Ret>
Ret dynamic::max(const object_type<T>& a, const object_type<T>& b)
{
    return matcl::unqualified_call::eval_max<Ret,T,T>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::max(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_max<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::max(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_max<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
bool dynamic::operator&&(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_and<bool,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
bool dynamic::operator&&(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_op_and<bool,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
bool dynamic::operator&&(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_and<bool,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
bool dynamic::operator||(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_or<bool,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
bool dynamic::operator||(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_op_or<bool,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
bool dynamic::operator||(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_or<bool,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
bool dynamic::op_and(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_and<bool,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
bool dynamic::op_and(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_op_and<bool,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
bool dynamic::op_and(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_and<bool,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
bool dynamic::op_or(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_or<bool,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
bool dynamic::op_or(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_op_or<bool,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
bool dynamic::op_or(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_or<bool,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
bool dynamic::op_xor(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_xor<bool,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
bool dynamic::op_xor(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_op_xor<bool,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
bool dynamic::op_xor(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_xor<bool,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::elem_and(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_elem_and<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::elem_and(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_elem_and<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::elem_and(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_elem_and<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::elem_or(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_elem_or<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::elem_or(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_elem_or<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::elem_or(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_elem_or<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::elem_xor(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_elem_xor<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::elem_xor(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_elem_xor<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::elem_xor(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_elem_xor<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::operator&(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_elem_or<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::operator&(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_elem_or<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::operator&(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_elem_or<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::operator|(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_elem_or<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::operator|(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_elem_or<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::operator|(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_elem_or<Ret,T,S>::eval(a,b.get());
}

template<class T, class Ret>
Ret dynamic::rising_factorial(const object_type<T>& x, Integer i)
{
    return matcl::unqualified_call::eval_rising_factorial<Ret,T>::eval(x.get(), i);
}

template<class T, class Ret>
Ret dynamic::falling_factorial(const object_type<T>& x, Integer i)
{
    return matcl::unqualified_call::eval_falling_factorial<Ret,T>::eval(x.get(), i);
}

template<class T, class TB, class Ret>
Ret dynamic::bernoulli_b2n(Integer n)
{
    return matcl::unqualified_call::eval_bernoulli_b2n<Ret,TB>::eval(n);
}

template<class T, class TB, class Ret>
Integer dynamic::max_bernoulli_b2n()
{
    return matcl::unqualified_call::eval_max_bernoulli_b2n<Integer,TB>::eval();
}

template<class T, class TB, class Ret>
Ret dynamic::factorial(Integer n)
{
    return matcl::unqualified_call::eval_factorial<Ret,TB>::eval(n);
}

template<class T, class TB, class Ret>
Ret dynamic::double_factorial(Integer n)
{
    return matcl::unqualified_call::eval_double_factorial<Ret,TB>::eval(n);
}

template<class T, class TB, class Ret>
Ret dynamic::binomial_coefficient(Integer n, Integer k)
{
    return matcl::unqualified_call::eval_binomial_coefficient<Ret,TB>::eval(n, k);
}

};};

#pragma warning(pop)