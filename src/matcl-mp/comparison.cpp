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

#include "matcl-mp/details/comparison.h"
#include "utils/impl_types.h"
#include "matcl-core/matrix/complex_type.h"
#include "matcl-core/details/scalfunc_real.h"
#include "matcl-mp/func_unary.h"

namespace matcl
{

static fp_type fpclassify(Real val)
{
    int code = std::fpclassify(val);

    switch(code)
    {
        case FP_INFINITE:   return fp_type::fp_infinite;
        case FP_NAN:        return fp_type::fp_nan;
        case FP_ZERO:       return fp_type::fp_zero;
        case FP_SUBNORMAL:  return fp_type::fp_subnormal;
        case FP_NORMAL:     return fp_type::fp_normal;
        default:
            return fp_type::fp_unknown;
    }
};

};

namespace matcl { namespace mp { namespace details
{

namespace bm = boost::multiprecision;
namespace mmd = matcl::mp::details;
namespace mrd = matcl::raw::details;

//---------------------------------------------------------------
//                  COMPARE
//---------------------------------------------------------------
static cmp_type get_cmp_type(int cmp)
{
    if (cmp < 0)
        return cmp_type::less;
    else if (cmp == 0)
        return cmp_type::equal;
    else
        return cmp_type::greater;
};

inline bool is_gt(cmp_type cmp)
{
    return cmp == cmp_type::greater;
}
inline bool is_lt(cmp_type cmp)
{
    return cmp == cmp_type::less;
}
inline bool is_eq(cmp_type cmp)
{
    return cmp == cmp_type::equal;
}
static cmp_type compare(const mp_int& val1, Real val2)
{
    if (mrd::scal_func::isnan(val2) == true)
        return cmp_type::not_ordered;

    int cmp = mpz_cmp_d(mmd::impl_value(val1).backend().data(), val2);
    return get_cmp_type(cmp);
}

static cmp_type compare(const mp_rational& val1, Real val2)
{
    fp_type fp = matcl::fpclassify(val2);
    if (fp == fp_type::fp_nan)
    {
        return cmp_type::not_ordered;
    }
    else if (fp == fp_type::fp_infinite)
    {
        if (val2 > 0.0)
            return cmp_type::less;
        else
            return cmp_type::greater;
    }
    else
    {
        int cmp = mmd::impl_value(val1).compare(mmd::impl_value(val2));
        return get_cmp_type(cmp);
    }
}
static cmp_type compare(const mp_float& val1, Real val2)
{
    mpfr_clear_erangeflag();
    int cmp = mpfr_cmp_d(mmd::impl_value(val1), val2);
    if (mpfr_erangeflag_p() != 0)
        return cmp_type::not_ordered;
    else
        return get_cmp_type(cmp);
}

static cmp_type compare(const mp_float& val1, Integer val2)
{
    mpfr_clear_erangeflag();
    int cmp = mpfr_cmp_si(mmd::impl_value(val1), val2); 
    if (mpfr_erangeflag_p() != 0)
        return cmp_type::not_ordered;
    else
        return get_cmp_type(cmp);
}
static cmp_type compare(const mp_float& val1, const mp_int& val2)
{
    mpfr_clear_erangeflag();
    int cmp = mpfr_cmp_z(mmd::impl_value(val1), mmd::impl_value(val2).backend().data());
    if (mpfr_erangeflag_p() != 0)
        return cmp_type::not_ordered;
    else
        return get_cmp_type(cmp);
}
static cmp_type compare(const mp_float& val1, const mp_rational& val2)
{
    mpfr_clear_erangeflag();
    int cmp = mpfr_cmp_q(mmd::impl_value(val1), mmd::impl_value(val2).backend().data());
    if (mpfr_erangeflag_p() != 0)
        return cmp_type::not_ordered;
    else
        return get_cmp_type(cmp);
}
cmp_type mmd::compare(const mp_float& val1, const mp_float& val2)
{
    mpfr_clear_erangeflag();
    int cmp = mpfr_cmp(mmd::impl_value(val1), mmd::impl_value(val2));
    if (mpfr_erangeflag_p() != 0)
        return cmp_type::not_ordered;
    else
        return get_cmp_type(cmp);
}
cmp_type mmd::compare_abs(const mp_float& val1, const mp_float& val2)
{
    mpfr_clear_erangeflag();
    int cmp = mpfr_cmpabs(mmd::impl_value(val1), mmd::impl_value(val2));
    if (mpfr_erangeflag_p() != 0)
        return cmp_type::not_ordered;
    else
        return get_cmp_type(cmp);
}

//---------------------------------------------------------------
//                  EQUALITY COMPARISON
//---------------------------------------------------------------
bool mmd::eeq(const mp_int& val1, const mp_int& val2)
{
    return bm::operator==(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::eeq(const mp_int& val1, Integer val2)
{
    return bm::operator==(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::eeq(Integer val1, const mp_int& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_int& val1, Real val2)
{
    if (mrd::scal_func::isnan(val2) == true)
        return false;

    return mpz_cmp_d(mmd::impl_value(val1).backend().data(), val2) == 0;
}
bool mmd::eeq(Real val1, const mp_int& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_int& val1, const Complex& val2)
{
    return matcl::imag(val2) == 0.0 
            && mmd::eeq(val1, matcl::real(val2));
}
bool mmd::eeq(const Complex& val1, const mp_int& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_int& val1, const Float_complex& val20)
{
    Complex val2 = val20;
    return mmd::eeq(val1,val2);
}
bool mmd::eeq(const Float_complex& val1, const mp_int& val2)
{
    return mmd::eeq(val2,val1);
}

bool mmd::eeq(const mp_rational& val1, const mp_rational& val2)
{
    return bm::operator==(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::eeq(const mp_rational& val1, Integer val2)
{
    return bm::operator==(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::eeq(Integer val1, const mp_rational& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_rational& val1, Real val2)
{
    if (std::isfinite(val2) == false)
        return false;

    return mmd::eeq(val1, mp_rational(val2));
}
bool mmd::eeq(Real val1, const mp_rational& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_rational& val1, const Complex& val2)
{
    return matcl::imag(val2) == 0.0 && mmd::eeq(val1, matcl::real(val2));
}
bool mmd::eeq(const Complex& val1, const mp_rational& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_rational& val1, const Float_complex& val20)
{
    Complex val2 = val20;
    return mmd::eeq(val1, val2);
}
bool mmd::eeq(const Float_complex& val1, const mp_rational& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_rational& val1, const mp_int& val2)
{
    return bm::operator==(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::eeq(const mp_int& val1, const mp_rational& val2)
{
    return mmd::eeq(val2,val1);
}

bool mmd::eeq(const mp_float& val1, const mp_float& val2)
{
    return mpfr_equal_p(mmd::impl_value(val1), mmd::impl_value(val2)) != 0;
}
bool mmd::eeq(const mp_float& val1, Integer val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_si(mmd::impl_value(val1), val2) == 0 
                && mpfr_erangeflag_p() == 0;
}
bool mmd::eeq(Integer val1, const mp_float& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_float& val1, Real val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_d(mmd::impl_value(val1), val2) == 0 
                && mpfr_erangeflag_p() == 0;
}
bool mmd::eeq(Real val1, const mp_float& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_float& val1, const Complex& val2)
{
    return matcl::imag(val2) == 0.0 && mmd::eeq(val1, matcl::real(val2));
}
bool mmd::eeq(const Complex& val1, const mp_float& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_float& val1, const Float_complex& val20)
{
    Complex val2 = val20;
    return mmd::eeq(val1, val2);
}
bool mmd::eeq(const Float_complex& val1, const mp_float& val2)
{
    return mmd::eeq(val2,val1);
}

bool mmd::eeq(const mp_float& val1, const mp_int& val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_z(mmd::impl_value(val1), mmd::impl_value(val2).backend().data()) == 0
                && mpfr_erangeflag_p() == 0;
}
bool mmd::eeq(const mp_int& val1, const mp_float& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_float& val1, const mp_rational& val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_q(mmd::impl_value(val1), mmd::impl_value(val2).backend().data()) == 0
                && mpfr_erangeflag_p() == 0;
}
bool mmd::eeq(const mp_rational& val1, const mp_float& val2)
{
    return mmd::eeq(val2,val1);
}

bool mmd::eeq(const mp_complex& val1, const mp_complex& val2)
{
    return mmd::eeq(val1.imag(), val2.imag()) && mmd::eeq(val1.real(), val2.real());
}
bool mmd::eeq(const mp_complex& val1, Integer val2)
{
    return mmd::eeq(val1.imag(), 0.0) && mmd::eeq(val1.real(), val2);
}
bool mmd::eeq(Integer val1, const mp_complex& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_complex& val1, Real val2)
{
    return mmd::eeq(val1.imag(), 0.0) && mmd::eeq(val1.real(), val2);
}
bool mmd::eeq(Real val1, const mp_complex& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_complex& val1, const Complex& val2)
{
    return mmd::eeq(val1.imag(), matcl::imag(val2)) && mmd::eeq(val1.real(), matcl::real(val2));
}
bool mmd::eeq(const Complex& val1, const mp_complex& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_complex& val1, const Float_complex& val20)
{
    Complex val2 = val20;
    return mmd::eeq(val1, val2);
}
bool mmd::eeq(const Float_complex& val1, const mp_complex& val2)
{
    return mmd::eeq(val2,val1);
}

bool mmd::eeq(const mp_complex& val1, const mp_int& val2)
{
    return mmd::eeq(val1.imag(), 0.0) && mmd::eeq(val1.real(), val2);
}
bool mmd::eeq(const mp_int& val1, const mp_complex& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_complex& val1, const mp_float& val2)
{
    return mmd::eeq(val1.imag(), 0.0) && mmd::eeq(val1.real(), val2);
}
bool mmd::eeq(const mp_float& val1, const mp_complex& val2)
{
    return mmd::eeq(val2,val1);
}
bool mmd::eeq(const mp_complex& val1, const mp_rational& val2)
{
    return mmd::eeq(val1.imag(), 0.0) && mmd::eeq(val1.real(), val2);
}
bool mmd::eeq(const mp_rational& val1, const mp_complex& val2)
{
    return mmd::eeq(val2,val1);
}

//---------------------------------------------------------------
//                  EQUALITY COMPARISON, NaN == NaN
//---------------------------------------------------------------
bool mmd::eeq_nan(const mp_int& val1, const mp_int& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_int& val1, Integer val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(Integer val1, const mp_int& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_int& val1, Float val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_int& val1, Real val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(Float val1, const mp_int& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(Real val1, const mp_int& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_int& val1, const Complex& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const Complex& val1, const mp_int& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_int& val1, const Float_complex& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const Float_complex& val1, const mp_int& val2)
{
    return mmd::eeq(val1, val2);
}

bool mmd::eeq_nan(const mp_rational& val1, const mp_rational& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_rational& val1, Integer val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(Integer val1, const mp_rational& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_rational& val1, Float val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_rational& val1, Real val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(Float val1, const mp_rational& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(Real val1, const mp_rational& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_rational& val1, const Complex& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const Complex& val1, const mp_rational& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_rational& val1, const Float_complex& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const Float_complex& val1, const mp_rational& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_rational& val1, const mp_int& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_int& val1, const mp_rational& val2)
{
    return mmd::eeq(val1, val2);
}

bool mmd::eeq_nan(const mp_float& val1, const mp_float& val2)
{
    if (is_nan(val1) == true)
    {
        if (is_nan(val2) == true)
            return true;
        else
            return false;
    }
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_float& val1, Integer val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(Integer val1, const mp_float& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_float& val1, Float val2)
{
    return eeq_nan(val1, (Real)val2);
}
bool mmd::eeq_nan(const mp_float& val1, Real val2)
{
    if (is_nan(val1) == true)
    {
        if (mrd::scal_func::isnan(val2) == true)
            return true;
        else
            return false;
    }
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(Float val1, const mp_float& val2)
{
    return eeq_nan((Real)val1, val2);
}
bool mmd::eeq_nan(Real val1, const mp_float& val2)
{
    if (mrd::scal_func::isnan(val1) == true)
    {
        if (is_nan(val2) == true)
            return true;
        else
            return false;
    }
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_float& val1, const Complex& val2)
{
    return matcl::imag(val2) == 0.0 && eeq_nan(val1, matcl::real(val2));
}
bool mmd::eeq_nan(const Complex& val1, const mp_float& val2)
{
    return eeq_nan(val2, val1);
}
bool mmd::eeq_nan(const mp_float& val1, const Float_complex& val2)
{
    return matcl::imag(val2) == 0.0 && eeq_nan(val1, matcl::real(val2));
}
bool mmd::eeq_nan(const Float_complex& val1, const mp_float& val2)
{
    return eeq_nan(val2, val1);
}

bool mmd::eeq_nan(const mp_float& val1, const mp_int& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_int& val1, const mp_float& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_float& val1, const mp_rational& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_rational& val1, const mp_float& val2)
{
    return mmd::eeq(val1, val2);
}

bool mmd::eeq_nan(const mp_complex& val1, const mp_complex& val2)
{
    return eeq_nan(val1.imag(), val2.imag()) && eeq_nan(val1.real(), val2.real());
}
bool mmd::eeq_nan(const mp_complex& val1, Integer val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(Integer val1, const mp_complex& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_complex& val1, Float val2)
{
    return mmd::eeq(val1.imag(), 0.0) && eeq_nan(val1.real(), val2);
}
bool mmd::eeq_nan(const mp_complex& val1, Real val2)
{
    return mmd::eeq(val1.imag(), 0.0) && eeq_nan(val1.real(), val2);
}
bool mmd::eeq_nan(Float val1, const mp_complex& val2)
{
    return eeq_nan(val2,val1);
}
bool mmd::eeq_nan(Real val1, const mp_complex& val2)
{
    return eeq_nan(val2,val1);
}
bool mmd::eeq_nan(const mp_complex& val1, const Complex& val2)
{
    return eeq_nan(val1.imag(), matcl::imag(val2))
        && eeq_nan(val1.real(), matcl::real(val2));
}
bool mmd::eeq_nan(const Complex& val1, const mp_complex& val2)
{
    return eeq_nan(val2,val1);
}
bool mmd::eeq_nan(const mp_complex& val1, const Float_complex& val20)
{
    Complex val2 = val20;
    return eeq_nan(val1, val2);
}
bool mmd::eeq_nan(const Float_complex& val1, const mp_complex& val2)
{
    return eeq_nan(val2,val1);
}

bool mmd::eeq_nan(const mp_complex& val1, const mp_int& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_int& val1, const mp_complex& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_complex& val1, const mp_float& val2)
{
    return mmd::eeq(val1.imag(), 0.0) && eeq_nan(val1.real(), val2);
}
bool mmd::eeq_nan(const mp_float& val1, const mp_complex& val2)
{
    return eeq_nan(val2,val1);
}
bool mmd::eeq_nan(const mp_complex& val1, const mp_rational& val2)
{
    return mmd::eeq(val1, val2);
}
bool mmd::eeq_nan(const mp_rational& val1, const mp_complex& val2)
{
    return mmd::eeq(val1, val2);
}

//---------------------------------------------------------------
//                  INEQUALITY COMPARISON
//---------------------------------------------------------------
bool mmd::neq(const mp_int& val1, const mp_int& val2)
{
    return bm::operator!=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::neq(const mp_int& val1, Integer val2)
{
    return bm::operator!=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::neq(Integer val1, const mp_int& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_int& val1, Real val2)
{
    if (mrd::scal_func::isnan(val2) == true)
        return true;

    return mpz_cmp_d(mmd::impl_value(val1).backend().data(), val2) != 0;
}
bool mmd::neq(Real val1, const mp_int& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_int& val1, const Complex& val2)
{
    return matcl::imag(val2) != 0.0 || mmd::neq(val1, matcl::real(val2));
}
bool mmd::neq(const Complex& val1, const mp_int& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_int& val1, const Float_complex& val2)
{
    return mmd::neq(val1, (Complex)val2);
}
bool mmd::neq(const Float_complex& val1, const mp_int& val2)
{
    return mmd::neq(val2, val1);
}

bool mmd::neq(const mp_rational& val1, const mp_rational& val2)
{
    return bm::operator!=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::neq(const mp_rational& val1, Integer val2)
{
    return bm::operator!=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::neq(Integer val1, const mp_rational& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_rational& val1, Real val2)
{
    if (mrd::scal_func::finite(val2) == false)
        return true;
    else
        return mmd::neq(val1, mp_rational(val2));
}
bool mmd::neq(Real val1, const mp_rational& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_rational& val1, const Complex& val2)
{
    return matcl::imag(val2) != 0.0 || mmd::neq(val1, matcl::real(val2));
}
bool mmd::neq(const Complex& val1, const mp_rational& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_rational& val1, const Float_complex& val2)
{
    return mmd::neq(val1, (Complex)val2);
}
bool mmd::neq(const Float_complex& val1, const mp_rational& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_rational& val1, const mp_int& val2)
{
    return bm::operator!=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::neq(const mp_int& val1, const mp_rational& val2)
{
    return mmd::neq(val2, val1);
}

bool mmd::neq(const mp_float& val1, const mp_float& val2)
{
    return !mmd::eeq(val1, val2);
};
bool mmd::neq(const mp_float& val1, Integer val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_si(mmd::impl_value(val1), val2) != 0
            || mpfr_erangeflag_p() != 0;
}
bool mmd::neq(Integer val1, const mp_float& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_float& val1, Real val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_d(mmd::impl_value(val1), val2) != 0
                || mpfr_erangeflag_p() != 0;
}
bool mmd::neq(Real val1, const mp_float& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_float& val1, const Complex& val2)
{
    return matcl::imag(val2) != 0.0 || mmd::neq(val1, matcl::real(val2));
}
bool mmd::neq(const Complex& val1, const mp_float& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_float& val1, const Float_complex& val2)
{
    return mmd::neq(val1, (Complex)val2);
}
bool mmd::neq(const Float_complex& val1, const mp_float& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_float& val1, const mp_int& val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_z(mmd::impl_value(val1), mmd::impl_value(val2).backend().data()) != 0
            || mpfr_erangeflag_p() != 0;
}
bool mmd::neq(const mp_int& val1, const mp_float& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_float& val1, const mp_rational& val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_q(mmd::impl_value(val1), mmd::impl_value(val2).backend().data()) != 0
            || mpfr_erangeflag_p() != 0;
}
bool mmd::neq(const mp_rational& val1, const mp_float& val2)
{
    return mmd::neq(val2, val1);
};

bool mmd::neq(const mp_complex& val1, const mp_complex& val2)
{
    return mmd::neq(val1.imag(), val2.imag()) || mmd::neq(val1.real(), val2.real());
}
bool mmd::neq(const mp_complex& val1, Integer val2)
{
    return mmd::neq(val1.imag(), 0.0) || mmd::neq(val1.real(), val2);
}
bool mmd::neq(Integer val1, const mp_complex& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_complex& val1, Real val2)
{
    return mmd::neq(val1.imag(), 0.0) || mmd::neq(val1.real(), val2);
}
bool mmd::neq(Real val1, const mp_complex& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_complex& val1, const Complex& val2)
{
    return mmd::neq(val1.imag(), matcl::imag(val2)) 
            || mmd::neq(val1.real(), matcl::real(val2));
}
bool mmd::neq(const Complex& val1, const mp_complex& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_complex& val1, const Float_complex& val2)
{
    return mmd::neq(val1, (Complex)val2);
}
bool mmd::neq(const Float_complex& val1, const mp_complex& val2)
{
    return mmd::neq(val2, val1);
}

bool mmd::neq(const mp_complex& val1, const mp_int& val2)
{
    return mmd::neq(val1.imag(), 0.0) || mmd::neq(val1.real(), val2);
}
bool mmd::neq(const mp_int& val1, const mp_complex& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_complex& val1, const mp_float& val2)
{
    return mmd::neq(val1.imag(), 0.0) || mmd::neq(val1.real(), val2);
}
bool mmd::neq(const mp_float& val1, const mp_complex& val2)
{
    return mmd::neq(val2, val1);
}
bool mmd::neq(const mp_complex& val1, const mp_rational& val2)
{
    return mmd::neq(val1.imag(), 0.0) || mmd::neq(val1.real(), val2);
}
bool mmd::neq(const mp_rational& val1, const mp_complex& val2)
{
    return mmd::neq(val2, val1);
}

//---------------------------------------------------------------
//                  INEQUALITY COMPARISON NAN
//---------------------------------------------------------------
bool mmd::neq_nan(const mp_int& val1, const mp_int& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_int& val1, Integer val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(Integer val1, const mp_int& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_int& val1, Float val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_int& val1, Real val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(Float val1, const mp_int& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(Real val1, const mp_int& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_int& val1, const Complex& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const Complex& val1, const mp_int& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_int& val1, const Float_complex& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const Float_complex& val1, const mp_int& val2)
{
    return mmd::neq(val1, val2);
}

bool mmd::neq_nan(const mp_rational& val1, const mp_rational& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_rational& val1, Integer val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(Integer val1, const mp_rational& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_rational& val1, Float val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_rational& val1, Real val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(Float val1, const mp_rational& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(Real val1, const mp_rational& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_rational& val1, const Complex& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const Complex& val1, const mp_rational& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_rational& val1, const Float_complex& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const Float_complex& val1, const mp_rational& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_rational& val1, const mp_int& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_int& val1, const mp_rational& val2)
{
    return mmd::neq(val1, val2);
}

bool mmd::neq_nan(const mp_float& val1, const mp_float& val2)
{
    return !eeq_nan(val1, val2);
};
bool mmd::neq_nan(const mp_float& val1, Integer val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(Integer val1, const mp_float& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_float& val1, Float val2)
{
    return !eeq_nan(val1, val2);
}
bool mmd::neq_nan(const mp_float& val1, Real val2)
{
    return !eeq_nan(val1, val2);
}
bool mmd::neq_nan(Float val1, const mp_float& val2)
{
    return neq_nan(val2, val1);
}
bool mmd::neq_nan(Real val1, const mp_float& val2)
{
    return neq_nan(val2, val1);
}
bool mmd::neq_nan(const mp_float& val1, const Complex& val2)
{
    return matcl::imag(val2) != 0.0 || neq_nan(val1, matcl::real(val2));
}
bool mmd::neq_nan(const Complex& val1, const mp_float& val2)
{
    return neq_nan(val2, val1);
}
bool mmd::neq_nan(const mp_float& val1, const Float_complex& val2)
{
    return neq_nan(val1, (Complex)val2);
}
bool mmd::neq_nan(const Float_complex& val1, const mp_float& val2)
{
    return neq_nan(val2, val1);
}
bool mmd::neq_nan(const mp_float& val1, const mp_int& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_int& val1, const mp_float& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_float& val1, const mp_rational& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_rational& val1, const mp_float& val2)
{
    return mmd::neq(val1, val2);
};

bool mmd::neq_nan(const mp_complex& val1, const mp_complex& val2)
{
    return neq_nan(val1.imag(), val2.imag()) || neq_nan(val1.real(), val2.real());
}
bool mmd::neq_nan(const mp_complex& val1, Integer val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(Integer val1, const mp_complex& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_complex& val1, Float val2)
{
    return mmd::neq(val1.imag(), 0.0) || neq_nan(val1.real(), val2);
}
bool mmd::neq_nan(const mp_complex& val1, Real val2)
{
    return mmd::neq(val1.imag(), 0.0) || neq_nan(val1.real(), val2);
}
bool mmd::neq_nan(Float val1, const mp_complex& val2)
{
    return neq_nan(val2, val1);
}
bool mmd::neq_nan(Real val1, const mp_complex& val2)
{
    return neq_nan(val2, val1);
}
bool mmd::neq_nan(const mp_complex& val1, const Complex& val2)
{
    return neq_nan(val1.imag(), matcl::imag(val2)) 
        || neq_nan(val1.real(), matcl::real(val2));
}
bool mmd::neq_nan(const Complex& val1, const mp_complex& val2)
{
    return neq_nan(val2, val1);
}
bool mmd::neq_nan(const mp_complex& val1, const Float_complex& val2)
{
    return neq_nan(val1, (Complex)val2);
}
bool mmd::neq_nan(const Float_complex& val1, const mp_complex& val2)
{
    return neq_nan(val2, val1);
}

bool mmd::neq_nan(const mp_complex& val1, const mp_int& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_int& val1, const mp_complex& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_complex& val1, const mp_float& val2)
{
    return mmd::neq(val1.imag(), 0.0) || neq_nan(val1.real(), val2);
}
bool mmd::neq_nan(const mp_float& val1, const mp_complex& val2)
{
    return neq_nan(val2, val1);
}
bool mmd::neq_nan(const mp_complex& val1, const mp_rational& val2)
{
    return mmd::neq(val1, val2);
}
bool mmd::neq_nan(const mp_rational& val1, const mp_complex& val2)
{
    return mmd::neq(val1, val2);
}

//---------------------------------------------------------------
//                  GREATER OR EQUAL COMPARISON
//---------------------------------------------------------------
bool mmd::geq(const mp_int& val1, const mp_int& val2)
{
    return bm::operator>=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::geq(const mp_int& val1, Integer val2)
{
    return bm::operator>=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::geq(Integer val1, const mp_int& val2)
{
    return mmd::leq(val2, val1);
}
bool mmd::geq(const mp_int& val1, Real val2)
{
    if (mrd::scal_func::isnan(val2) == true)
        return false;

    return mpz_cmp_d(mmd::impl_value(val1).backend().data(), val2) >= 0;
}
bool mmd::geq(Real val1, const mp_int& val2)
{
    return mmd::leq(val2, val1);
}
bool mmd::geq(const mp_int& val1, const Complex& val2)
{
    cmp_type cmp = compare(val1, matcl::real(val2));
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return 0.0 >= matcl::imag(val2);
    else
        return false;
}
bool mmd::geq(const Complex& val1, const mp_int& val2)
{
    return mmd::leq(val2, val1);
}
bool mmd::geq(const mp_int& val1, const Float_complex& val2)
{
    return mmd::geq(val1, Complex(val2));
}
bool mmd::geq(const Float_complex& val1, const mp_int& val2)
{
    return mmd::leq(val2, val1);
}
bool mmd::geq(const mp_rational& val1, const mp_rational& val2)
{
    return bm::operator>=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::geq(const mp_rational& val1, Integer val2)
{
    return bm::operator>=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::geq(Integer val1, const mp_rational& val2)
{
    return mmd::leq(val2, val1);
}
bool mmd::geq(const mp_rational& val1, Real val2)
{
    fp_type fp = matcl::fpclassify(val2);
    if (fp == fp_type::fp_nan)
    {
        return false;
    }
    else if (fp == fp_type::fp_infinite)
    {
        if (val2 > 0.0)
            return false;
        else
            return true;
    }
    else
    {
        return mmd::geq(val1, mp_rational(val2));
    }
}
bool mmd::geq(Real val1, const mp_rational& val2)
{
    return mmd::leq(val2, val1);
}
bool mmd::geq(const mp_rational& val1, const Complex& val2)
{
    cmp_type cmp = compare(val1, matcl::real(val2));
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return 0.0 >= matcl::imag(val2);
    else
        return false;
}
bool mmd::geq(const Complex& val1, const mp_rational& val2)
{
    return mmd::leq(val2, val1);
}
bool mmd::geq(const mp_rational& val1, const Float_complex& val2)
{
    return mmd::geq(val1, Complex(val2));
}
bool mmd::geq(const Float_complex& val1, const mp_rational& val2)
{
    return mmd::leq(val2, val1);
}
bool mmd::geq(const mp_rational& val1, const mp_int& val2)
{
    return bm::operator>=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::geq(const mp_int& val1, const mp_rational& val2)
{
    return mmd::leq(val2, val1);
}

bool mmd::geq(const mp_float& val1, const mp_float& val2)
{
    return mpfr_greaterequal_p(mmd::impl_value(val1), mmd::impl_value(val2)) != 0;
};
bool mmd::geq(const mp_float& val1, Integer val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_si(mmd::impl_value(val1), val2) >= 0 
                && mpfr_erangeflag_p() == 0;
}
bool mmd::geq(Integer val1, const mp_float& val2)
{
    return mmd::leq(val2, val1);
}
bool mmd::geq(const mp_float& val1, Real val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_d(mmd::impl_value(val1), val2) >= 0 
                && mpfr_erangeflag_p() == 0;
}
bool mmd::geq(Real val1, const mp_float& val2)
{
    return mmd::leq(val2, val1);
}
bool mmd::geq(const mp_float& val1, const Complex& val2)
{
    cmp_type cmp = compare(val1, matcl::real(val2));
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return 0.0 >= matcl::imag(val2);
    else
        return false;
}
bool mmd::geq(const Complex& val1, const mp_float& val2)
{
    return mmd::leq(val2,val1);
}
bool mmd::geq(const mp_float& val1, const Float_complex& val2)
{
    return mmd::geq(val1, Complex(val2));
}
bool mmd::geq(const Float_complex& val1, const mp_float& val2)
{
    return mmd::leq(val2,val1);
}
bool mmd::geq(const mp_float& val1, const mp_int& val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_z(mmd::impl_value(val1), mmd::impl_value(val2).backend().data()) >= 0
                && mpfr_erangeflag_p() == 0;
}
bool mmd::geq(const mp_int& val1, const mp_float& val2)
{
    return mmd::leq(val2,val1);
}
bool mmd::geq(const mp_float& val1, const mp_rational& val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_q(mmd::impl_value(val1), mmd::impl_value(val2).backend().data()) >= 0
                && mpfr_erangeflag_p() == 0;
}
bool mmd::geq(const mp_rational& val1, const mp_float& val2)
{
    return mmd::leq(val2,val1);
}

bool mmd::geq(const mp_complex& val1, const mp_complex& val2)
{
    cmp_type cmp = compare(matcl::real(val1), matcl::real(val2));
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::geq(matcl::imag(val1), matcl::imag(val2));
    else
        return false;
}
bool mmd::geq(const mp_complex& val1, Integer val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::geq(matcl::imag(val1), 0);
    else
        return false;
}
bool mmd::geq(Integer val1, const mp_complex& val2)
{
    return mmd::leq(val2,val1);
}
bool mmd::geq(const mp_complex& val1, Real val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::geq(matcl::imag(val1), 0.0);
    else
        return false;         
}
bool mmd::geq(Real val1, const mp_complex& val2)
{
    return mmd::leq(val2,val1);
}
bool mmd::geq(const mp_complex& val1, const Complex& val2)
{
    cmp_type cmp = compare(matcl::real(val1), matcl::real(val2));
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::geq(matcl::imag(val1), matcl::imag(val2));
    else
        return false;         
}
bool mmd::geq(const Complex& val1, const mp_complex& val2)
{
    return mmd::leq(val2,val1);
}
bool mmd::geq(const mp_complex& val1, const Float_complex& val2)
{
    return mmd::geq(val1, Complex(val2));
}
bool mmd::geq(const Float_complex& val1, const mp_complex& val2)
{
    return mmd::leq(val2,val1);
}
bool mmd::geq(const mp_complex& val1, const mp_int& val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::geq(matcl::imag(val1), 0);
    else
        return false;         
}
bool mmd::geq(const mp_int& val1, const mp_complex& val2)
{
    return mmd::leq(val2,val1);
}
bool mmd::geq(const mp_complex& val1, const mp_float& val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::geq(matcl::imag(val1), 0.0);
    else
        return false;         
}
bool mmd::geq(const mp_float& val1, const mp_complex& val2)
{
    return mmd::leq(val2, val1);
}
bool mmd::geq(const mp_complex& val1, const mp_rational& val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::geq(matcl::imag(val1), 0.0);
    else
        return false;         
}
bool mmd::geq(const mp_rational& val1, const mp_complex& val2)
{
    return mmd::leq(val2, val1);
}

//---------------------------------------------------------------
//                  LESS OR EQUAL COMPARISON
//---------------------------------------------------------------
bool mmd::leq(const mp_int& val1, const mp_int& val2)
{
    return bm::operator<=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::leq(const mp_int& val1, Integer val2)
{
    return bm::operator<=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::leq(Integer val1, const mp_int& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_int& val1, Real val2)
{
    if (mrd::scal_func::isnan(val2) == true)
        return false;

    return mpz_cmp_d(mmd::impl_value(val1).backend().data(), val2) <= 0;
}
bool mmd::leq(Real val1, const mp_int& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_int& val1, const Complex& val2)
{
    cmp_type cmp = compare(val1, matcl::real(val2));
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return 0.0 <= matcl::imag(val2);
    else
        return false;
}
bool mmd::leq(const Complex& val1, const mp_int& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_int& val1, const Float_complex& val2)
{
    return mmd::leq(val1, Complex(val2));
}
bool mmd::leq(const Float_complex& val1, const mp_int& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_rational& val1, const mp_rational& val2)
{
    return bm::operator<=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::leq(const mp_rational& val1, Integer val2)
{
    return bm::operator<=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::leq(Integer val1, const mp_rational& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_rational& val1, Real val2)
{
    fp_type fp = matcl::fpclassify(val2);
    if (fp == fp_type::fp_nan)
    {
        return false;
    }
    else if (fp == fp_type::fp_infinite)
    {
        if (val2 > 0.0)
            return true;
        else
            return false;
    }
    else
    {
        return mmd::leq(val1, mp_rational(val2));
    }
}
bool mmd::leq(Real val1, const mp_rational& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_rational& val1, const Complex& val2)
{
    cmp_type cmp = compare(val1, matcl::real(val2));
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return 0.0 <= matcl::imag(val2);
    else
        return false;         
}
bool mmd::leq(const Complex& val1, const mp_rational& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_rational& val1, const Float_complex& val2)
{
    return mmd::leq(val1, Complex(val2));
}
bool mmd::leq(const Float_complex& val1, const mp_rational& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_rational& val1, const mp_int& val2)
{
    return bm::operator<=(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::leq(const mp_int& val1, const mp_rational& val2)
{
    return mmd::geq(val2, val1);
}

bool mmd::leq(const mp_float& val1, const mp_float& val2)
{
    return mpfr_lessequal_p(mmd::impl_value(val1), mmd::impl_value(val2)) != 0;
};
bool mmd::leq(const mp_float& val1, Integer val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_si(mmd::impl_value(val1), val2) <= 0 
                && mpfr_erangeflag_p() == 0;
}
bool mmd::leq(Integer val1, const mp_float& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_float& val1, Real val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_d(mmd::impl_value(val1), val2) <= 0 
                && mpfr_erangeflag_p() == 0;
}
bool mmd::leq(Real val1, const mp_float& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_float& val1, const Complex& val2)
{
    cmp_type cmp = compare(val1, matcl::real(val2));
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return 0.0 <= matcl::imag(val2);
    else
        return false;         
}
bool mmd::leq(const Complex& val1, const mp_float& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_float& val1, const Float_complex& val2)
{
    return mmd::leq(val1, Complex(val2));
}
bool mmd::leq(const Float_complex& val1, const mp_float& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_float& val1, const mp_int& val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_z(mmd::impl_value(val1), mmd::impl_value(val2).backend().data()) <= 0
                && mpfr_erangeflag_p() == 0;
}
bool mmd::leq(const mp_int& val1, const mp_float& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_float& val1, const mp_rational& val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_q(mmd::impl_value(val1), mmd::impl_value(val2).backend().data()) <= 0
                && mpfr_erangeflag_p() == 0;
}
bool mmd::leq(const mp_rational& val1, const mp_float& val2)
{
    return mmd::geq(val2, val1);
}

bool mmd::leq(const mp_complex& val1, const mp_complex& val2)
{
    cmp_type cmp = compare(matcl::real(val1), matcl::real(val2));
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::leq(matcl::imag(val1), matcl::imag(val2));
    else
        return false;         
};
bool mmd::leq(const mp_complex& val1, Integer val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::leq(matcl::imag(val1), 0);
    else
        return false;         
}
bool mmd::leq(Integer val1, const mp_complex& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_complex& val1, Real val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::leq(matcl::imag(val1), 0.0);
    else
        return false;         
}
bool mmd::leq(Real val1, const mp_complex& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_complex& val1, const Complex& val2)
{
    cmp_type cmp = compare(matcl::real(val1), matcl::real(val2));
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::leq(matcl::imag(val1), matcl::imag(val2));
    else
        return false;         
}
bool mmd::leq(const Complex& val1, const mp_complex& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_complex& val1, const Float_complex& val2)
{
    return mmd::leq(val1, Complex(val2));
}
bool mmd::leq(const Float_complex& val1, const mp_complex& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_complex& val1, const mp_int& val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::leq(matcl::imag(val1), 0);
    else
        return false;         
}
bool mmd::leq(const mp_int& val1, const mp_complex& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_complex& val1, const mp_float& val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::leq(matcl::imag(val1), 0.0);
    else
        return false;         
}
bool mmd::leq(const mp_float& val1, const mp_complex& val2)
{
    return mmd::geq(val2, val1);
}
bool mmd::leq(const mp_complex& val1, const mp_rational& val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::leq(matcl::imag(val1), 0.0);
    else
        return false;         
}
bool mmd::leq(const mp_rational& val1, const mp_complex& val2)
{
    return mmd::geq(val2, val1);
}

//---------------------------------------------------------------
//                  GREATER THAN COMPARISON
//---------------------------------------------------------------
bool mmd::gt(const mp_int& val1, const mp_int& val2)
{
    return bm::operator>(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::gt(const mp_int& val1, Integer val2)
{
    return bm::operator>(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::gt(Integer val1, const mp_int& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_int& val1, Real val2)
{
    if (mrd::scal_func::isnan(val2) == true)
        return false;

    return mpz_cmp_d(mmd::impl_value(val1).backend().data(), val2) > 0;
}
bool mmd::gt(Real val1, const mp_int& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_int& val1, const Complex& val2)
{
    cmp_type cmp = compare(val1, matcl::real(val2));
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return 0.0 > matcl::imag(val2);
    else
        return false;         
}
bool mmd::gt(const Complex& val1, const mp_int& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_int& val1, const Float_complex& val2)
{
    return mmd::gt(val1, Complex(val2));
}
bool mmd::gt(const Float_complex& val1, const mp_int& val2)
{
    return mmd::lt(val2, val1);
}

bool mmd::gt(const mp_rational& val1, const mp_rational& val2)
{
    return bm::operator>(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::gt(const mp_rational& val1, Integer val2)
{
    return bm::operator>(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::gt(Integer val1, const mp_rational& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_rational& val1, Real val2)
{
    fp_type fp = matcl::fpclassify(val2);
    if (fp == fp_type::fp_nan)
    {
        return false;
    }
    else if (fp == fp_type::fp_infinite)
    {
        if (val2 > 0.0)
            return false;
        else
            return true;
    }
    else
    {
        return mmd::gt(val1, mp_rational(val2));
    }
}
bool mmd::gt(Real val1, const mp_rational& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_rational& val1, const Complex& val2)
{
    cmp_type cmp = compare(val1, matcl::real(val2));
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return 0.0 > matcl::imag(val2);
    else
        return false;         
}
bool mmd::gt(const Complex& val1, const mp_rational& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_rational& val1, const Float_complex& val2)
{
    return mmd::gt(val1, Complex(val2));
}
bool mmd::gt(const Float_complex& val1, const mp_rational& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_rational& val1, const mp_int& val2)
{
    return bm::operator>(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::gt(const mp_int& val1, const mp_rational& val2)
{
    return mmd::lt(val2, val1);
}

bool mmd::gt(const mp_float& val1, const mp_float& val2)
{
    return mpfr_greater_p(mmd::impl_value(val1), mmd::impl_value(val2)) != 0;
};
bool mmd::gt(const mp_float& val1, Integer val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_si(mmd::impl_value(val1), val2) > 0 
                && mpfr_erangeflag_p() == 0;
}
bool mmd::gt(Integer val1, const mp_float& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_float& val1, Real val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_d(mmd::impl_value(val1), val2) > 0 
                && mpfr_erangeflag_p() == 0;
}
bool mmd::gt(Real val1, const mp_float& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_float& val1, const Complex& val2)
{
    cmp_type cmp = compare(val1, matcl::real(val2));
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return 0.0 > matcl::imag(val2);
    else
        return false;         
}
bool mmd::gt(const Complex& val1, const mp_float& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_float& val1, const Float_complex& val2)
{
    return mmd::gt(val1, Complex(val2));
}
bool mmd::gt(const Float_complex& val1, const mp_float& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_float& val1, const mp_int& val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_z(mmd::impl_value(val1), mmd::impl_value(val2).backend().data()) > 0
                && mpfr_erangeflag_p() == 0;
}
bool mmd::gt(const mp_int& val1, const mp_float& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_float& val1, const mp_rational& val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_q(mmd::impl_value(val1), mmd::impl_value(val2).backend().data()) > 0
                && mpfr_erangeflag_p() == 0;
}
bool mmd::gt(const mp_rational& val1, const mp_float& val2)
{
    return mmd::lt(val2, val1);
}

bool mmd::gt(const mp_complex& val1, const mp_complex& val2)
{
    cmp_type cmp = compare(matcl::real(val1), matcl::real(val2));
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::gt(matcl::imag(val1), matcl::imag(val2));
    else
        return false;         
};
bool mmd::gt(const mp_complex& val1, Integer val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::gt(matcl::imag(val1), 0);
    else
        return false;         
}
bool mmd::gt(Integer val1, const mp_complex& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_complex& val1, Real val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::gt(matcl::imag(val1), 0.0);
    else
        return false;         
}
bool mmd::gt(Real val1, const mp_complex& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_complex& val1, const Complex& val2)
{
    cmp_type cmp = compare(matcl::real(val1), matcl::real(val2));
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::gt(matcl::imag(val1), matcl::imag(val2));
    else
        return false;         
}
bool mmd::gt(const Complex& val1, const mp_complex& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_complex& val1, const Float_complex& val2)
{
    return mmd::gt(val1, Complex(val2));
}
bool mmd::gt(const Float_complex& val1, const mp_complex& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_complex& val1, const mp_int& val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::gt(matcl::imag(val1), 0);
    else
        return false;         
}
bool mmd::gt(const mp_int& val1, const mp_complex& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_complex& val1, const mp_float& val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::gt(matcl::imag(val1), 0.0);
    else
        return false;         
}
bool mmd::gt(const mp_float& val1, const mp_complex& val2)
{
    return mmd::lt(val2, val1);
}
bool mmd::gt(const mp_complex& val1, const mp_rational& val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_gt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::gt(matcl::imag(val1), 0.0);
    else
        return false;         
}
bool mmd::gt(const mp_rational& val1, const mp_complex& val2)
{
    return mmd::lt(val2, val1);
}

//---------------------------------------------------------------
//                  LESS THAN COMPARISON
//---------------------------------------------------------------
bool mmd::lt(const mp_int& val1, const mp_int& val2)
{
    return bm::operator<(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::lt(const mp_int& val1, Integer val2)
{
    return bm::operator<(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::lt(Integer val1, const mp_int& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_int& val1, Real val2)
{
    if (mrd::scal_func::isnan(val2) == true)
        return false;

    return mpz_cmp_d(mmd::impl_value(val1).backend().data(), val2) < 0;
}
bool mmd::lt(Real val1, const mp_int& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_int& val1, const Complex& val2)
{
    cmp_type cmp = compare(val1, matcl::real(val2));
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return 0.0 < matcl::imag(val2);
    else
        return false;         
}
bool mmd::lt(const Complex& val1, const mp_int& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_int& val1, const Float_complex& val2)
{
    return mmd::lt(val1, Complex(val2));
}
bool mmd::lt(const Float_complex& val1, const mp_int& val2)
{
    return mmd::gt(val2, val1);
}

bool mmd::lt(const mp_rational& val1, const mp_rational& val2)
{
    return bm::operator<(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::lt(const mp_rational& val1, Integer val2)   
{
    return bm::operator<(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::lt(Integer val1, const mp_rational& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_rational& val1, Real val2)
{
    fp_type fp = matcl::fpclassify(val2);
    if (fp == fp_type::fp_nan)
    {
        return false;
    }
    else if (fp == fp_type::fp_infinite)
    {
        if (val2 > 0.0)
            return true;
        else
            return false;
    }
    else
    {
        return mmd::lt(val1, mp_rational(val2));
    }
}
bool mmd::lt(Real val1, const mp_rational& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_rational& val1, const Complex& val2)
{
    cmp_type cmp = compare(val1, matcl::real(val2));
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return 0.0 < matcl::imag(val2);
    else
        return false;         
}
bool mmd::lt(const Complex& val1, const mp_rational& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_rational& val1, const Float_complex& val2)
{
    return mmd::lt(val1, Complex(val2));
}
bool mmd::lt(const Float_complex& val1, const mp_rational& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_rational& val1, const mp_int& val2)
{
    return bm::operator<(mmd::impl_value(val1), mmd::impl_value(val2));
}
bool mmd::lt(const mp_int& val1, const mp_rational& val2)
{
    return mmd::gt(val2, val1);
}

bool mmd::lt(const mp_float& val1, const mp_float& val2)
{
    return mpfr_less_p(mmd::impl_value(val1), mmd::impl_value(val2)) != 0;
};
bool mmd::lt(const mp_float& val1, Integer val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_si(mmd::impl_value(val1), val2) < 0 
                && mpfr_erangeflag_p() == 0;
};
bool mmd::lt(Integer val1, const mp_float& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_float& val1, Real val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_d(mmd::impl_value(val1), val2) < 0 
                && mpfr_erangeflag_p() == 0;
}
bool mmd::lt(Real val1, const mp_float& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_float& val1, const Complex& val2)
{
    cmp_type cmp = compare(val1, matcl::real(val2));
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return 0.0 < matcl::imag(val2);
    else
        return false;         
}
bool mmd::lt(const Complex& val1, const mp_float& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_float& val1, const Float_complex& val2)
{
    return mmd::lt(val1, Complex(val2));
}
bool mmd::lt(const Float_complex& val1, const mp_float& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_float& val1, const mp_int& val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_z(mmd::impl_value(val1), mmd::impl_value(val2).backend().data()) < 0
                && mpfr_erangeflag_p() == 0;
}
bool mmd::lt(const mp_int& val1, const mp_float& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_float& val1, const mp_rational& val2)
{
    mpfr_clear_erangeflag();
    return mpfr_cmp_q(mmd::impl_value(val1), mmd::impl_value(val2).backend().data()) < 0
                && mpfr_erangeflag_p() == 0;
}
bool mmd::lt(const mp_rational& val1, const mp_float& val2)
{
    return mmd::gt(val2, val1);
}

bool mmd::lt(const mp_complex& val1, const mp_complex& val2)
{
    cmp_type cmp = compare(matcl::real(val1), matcl::real(val2));
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::lt(matcl::imag(val1), matcl::imag(val2));
    else
        return false;         
};
bool mmd::lt(const mp_complex& val1, Integer val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::lt(matcl::imag(val1), 0);
    else
        return false;         
}
bool mmd::lt(Integer val1, const mp_complex& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_complex& val1, Real val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::lt(matcl::imag(val1), 0.0);
    else
        return false;         
}
bool mmd::lt(Real val1, const mp_complex& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_complex& val1, const Complex& val2)
{
    cmp_type cmp = compare(matcl::real(val1), matcl::real(val2));
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::lt(matcl::imag(val1), matcl::imag(val2));
    else
        return false;         
}
bool mmd::lt(const Complex& val1, const mp_complex& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_complex& val1, const Float_complex& val2)
{
    return mmd::lt(val1, Complex(val2));
}
bool mmd::lt(const Float_complex& val1, const mp_complex& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_complex& val1, const mp_int& val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::lt(matcl::imag(val1), 0);
    else
        return false;         
}
bool mmd::lt(const mp_int& val1, const mp_complex& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_complex& val1, const mp_float& val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::lt(matcl::imag(val1), 0.0);
    else
        return false;         
}
bool mmd::lt(const mp_float& val1, const mp_complex& val2)
{
    return mmd::gt(val2, val1);
}
bool mmd::lt(const mp_complex& val1, const mp_rational& val2)
{
    cmp_type cmp = compare(matcl::real(val1), val2);
    if (is_lt(cmp) == true)
        return true;
    else if (is_eq(cmp) == true)
        return mmd::lt(matcl::imag(val1), 0.0);
    else
        return false;         
}
bool mmd::lt(const mp_rational& val1, const mp_complex& val2)
{
    return mmd::gt(val2, val1);
}

};}};
