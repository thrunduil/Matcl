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

#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl
{

namespace mrd   = matcl::raw::details;
namespace mr    = matcl::raw;

template<class S, class Enable>
typename md::bool_or_object<S>::type
matcl::is_finite(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::isfinite_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::bool_or_object<S>::type
matcl::is_nan(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::isnan_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::bool_or_object<S>::type
matcl::is_inf(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::isinf_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::bool_or_object<S>::type
matcl::is_regular(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::isregular_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::bool_or_object<S>::type
matcl::is_normal(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::isnormal_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::bool_or_object<S>::type
matcl::is_int(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::isint_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::bool_or_object<S>::type
matcl::is_real(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::isreal_helper<SP>::eval(m);
};

template<class S, class Enable>
bool matcl::is_zero(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::iszero_helper<SP>::eval(m);
};

template<class S, class Enable>
bool matcl::is_one(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::isone_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::real_type_promote<S>::type imag(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::imag_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::real_type_promote<S>::type real(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::real_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::real_unify_types_promote<S,Float>::type eps(const S& x)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::eps_helper<SP>::eval(SP(x));
};

template<class S, class Enable>
typename md::real_type_promote<S>::type abs(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::abs_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::real_type_promote<S>::type abs2(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::abs2_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::real_type_promote_unify<S,Float>::type arg(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::arg_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::real_type_promote_unify<S,Float>::type angle(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::arg_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::promote_scalar<S>::type matcl::conj(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::conj_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
matcl::sqrt(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::sqrt_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type 
sqrt_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::sqrt_c_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type 
matcl::sqrt1pm1(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::sqrt1pm1_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type 
matcl::sqrt1pm1_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::sqrt1pm1_c_helper<SP>::eval(SP(m));
};


template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type matcl::cbrt(const S& m)
{
    static_assert(md::is_complex<S>::value == false, "complex scalars not allowed");

    using SP = typename md::promote_scalar<S>::type;
    return mrd::cbrt_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type  matcl::expm1(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::expm1_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type  matcl::expi(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::expi_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
matcl::log(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::log_helper<SP>::eval(SP(m));
};
template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type 
log_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::log_c_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type 
matcl::log1p(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::log1p_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type 
log1p_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::log1p_c_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type exp(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::exp_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
matcl::log2(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::log2_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
log2_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::log2_c_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
matcl::log10(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::log10_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
log10_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::log10_c_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::promote_scalar<S>::type matcl::floor(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::floor_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::promote_scalar<S>::type matcl::ceil(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::ceil_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::promote_scalar<S>::type matcl::round(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::round_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::promote_scalar<S>::type matcl::trunc(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::trunc_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::integer_or_object<S>::type matcl::ifloor(const S& m)
{
    static_assert(md::is_complex<S>::value == false, "complex scalars not allowed");
    using SP = typename md::promote_scalar<S>::type;
    return mrd::ifloor_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::integer_or_object<S>::type matcl::iceil(const S& m)
{
    static_assert(md::is_complex<S>::value == false, "complex scalars not allowed");
    using SP = typename md::promote_scalar<S>::type;
    return mrd::iceil_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::integer_or_object<S>::type matcl::iround(const S& m)
{
    static_assert(md::is_complex<S>::value == false, "complex scalars not allowed");
    using SP = typename md::promote_scalar<S>::type;
    return mrd::iround_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::integer_or_object<S>::type matcl::itrunc(const S& m)
{
    static_assert(md::is_complex<S>::value == false, "complex scalars not allowed");
    using SP = typename md::promote_scalar<S>::type;
    return mrd::itrunc_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::promote_scalar<S>::type matcl::sign(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::sign_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::integer_or_object<S>::type
matcl::isign(const S& m)
{
    static_assert(md::is_complex<S>::value == false, "complex scalars not allowed");
    using SP = typename md::promote_scalar<S>::type;
    return mrd::isign_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type sin(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::sin_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type cos(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::cos_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type tan(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::tan_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type cot(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::cot_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type sec(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::sec_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type csc(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::csc_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
matcl::asin(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::asin_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
asin_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::asin_c_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
matcl::acos(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::acos_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
acos_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::acos_c_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
matcl::asec(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::asec_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
asec_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::asec_c_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
matcl::acsc(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::acsc_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
acsc_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::acsc_c_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type atan(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::atan_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type acot(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::acot_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type sinh(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::sinh_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type cosh(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::cosh_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type tanh(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::tanh_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type coth(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::coth_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type sech(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::sech_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type csch(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::csch_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type asinh(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::asinh_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
matcl::acosh(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::acosh_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
acosh_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::acosh_c_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
matcl::atanh(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::atanh_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
atanh_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::atanh_c_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
matcl::acoth(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::acoth_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
acoth_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::acoth_c_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type
matcl::asech(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;

    return mrd::asech_helper<SP>::eval(SP(m));
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float_complex>::type
asech_c(const S& m)
{
    using SP = typename md::unify_types_promote<S,Float>::type;
    return mrd::asech_c_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type acsch(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::acsch_helper<SP>::eval(m);
};

template<class S, class Enable>
bool matcl::cast_bool(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::cast_bool_helper<SP>::eval(m);
};

template<class S, class Enable>
bool matcl::op_not(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::op_not_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::bool_or_object<S>::type
matcl::neg(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::op_neg_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::bool_or_object<S>::type
matcl::is_true(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::op_true_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::promote_scalar<S>::type matcl::uminus(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::uminus_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S, Float>::type
matcl::invs(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::invs_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S, Float>::type
matcl::inv(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::inv_helper<SP>::eval(m);
};

template<class S, class Enable>
bool matcl::is_scalar_true(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::is_true_helper<SP>::eval(m);
}

template<class S, class Enable>
bool matcl::is_scalar_false(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::is_false_helper<SP>::eval(m);
}

template<class S, class Enable>
typename md::bool_or_object<S>::type
matcl::signbit(const S& m)
{
    static_assert(md::is_complex<S>::value == false, "complex scalars not allowed");
    using SP = typename md::promote_scalar<S>::type;
    return mrd::signbit_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::real_unify_types_promote<S,Float>::type 
matcl::logb(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::logb_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::integer_or_object<S>::type matcl::ilogb(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::ilogb_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type matcl::exp2(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::exp2_helper<SP>::eval(m);
};

template<class S, class Enable>
typename md::unify_types_promote<S,Float>::type matcl::exp10(const S& m)
{
    using SP = typename md::promote_scalar<S>::type;
    return mrd::exp10_helper<SP>::eval(m);
};

template<class S, class Enable>
fp_type matcl::fpclassify(const S& m)
{
    static_assert(md::is_complex<S>::value == false, "complex scalars not allowed");

    using SP    = typename md::promote_scalar<S>::type;
    fp_type val = mrd::fpclassify_helper<SP>::eval(m);
    return val;
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::ldexp(const S1& x, Integer exp)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return mrd::ldexp_helper<S>::eval(S(x), exp);
}

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::scalbn(const S1& x, Integer n)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return mrd::scalbn_helper<S>::eval(S(x), n);
}

template<class S1, class Enable>
typename md::real_unify_types_promote<S1,Float>::type
matcl::frexp(const S1& x, Integer& n)
{
    static_assert(md::is_complex<S1>::value == false, "complex scalars not allowed");

    using S = typename md::unify_types_promote<S1,Float>::type;
    return mrd::frexp_helper<S>::eval(S(x), n);
}

template<class S1, class Enable, class Return>
Return matcl::modf(const S1& x, Return& int_part)
{
    static_assert(md::is_complex<S1>::value == false, "complex scalars not allowed");
    using S = typename md::unify_types_promote<S1,Float>::type;
    return mrd::modf_helper<S>::eval(S(x), int_part);
}

};
