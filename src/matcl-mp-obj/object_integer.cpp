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

#include "matcl-mp-obj/mp_object.h"
#include "matcl-dynamic/details/register_function_macro.h"
#include "matcl-dynamic/matcl_function_names.h"
#include "matcl-scalar/lib_functions/func_forwarding.h"

namespace matcl 
{

namespace mdy = matcl::dynamic;
namespace mdyf = matcl::dynamic::functions;

//--------------------------------------------------------------
//              functions
//--------------------------------------------------------------
#define MATCL_REGISTER_SCALAR_FUNC_SIMP(fn)                             \
MATCL_REGISTER_SCALAR_FUNC(f_##fn, fn, mp_int, mdyf::fn)

MATCL_REGISTER_SCALAR_FUNC_SIMP(is_nan)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_inf)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_regular)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_finite)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_normal)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_int)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_real)
MATCL_REGISTER_SCALAR_FUNC_SIMP(fpclassify)
MATCL_REGISTER_SCALAR_FUNC_SIMP(real)
MATCL_REGISTER_SCALAR_FUNC_SIMP(imag)
MATCL_REGISTER_SCALAR_FUNC_SIMP(abs)
MATCL_REGISTER_SCALAR_FUNC_SIMP(abs2)
MATCL_REGISTER_SCALAR_FUNC_SIMP(conj)
MATCL_REGISTER_SCALAR_FUNC_SIMP(arg)
MATCL_REGISTER_SCALAR_FUNC_SIMP(signbit)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sign)
MATCL_REGISTER_SCALAR_FUNC_SIMP(isign)
MATCL_REGISTER_SCALAR_FUNC_SIMP(nextabove)
MATCL_REGISTER_SCALAR_FUNC_SIMP(nextbelow)
MATCL_REGISTER_SCALAR_FUNC_SIMP(eps)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sqrt)
MATCL_REGISTER_SCALAR_FUNC_SIMP(cbrt)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sqrt_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(exp)
MATCL_REGISTER_SCALAR_FUNC_SIMP(expm1)
MATCL_REGISTER_SCALAR_FUNC_SIMP(expi)
MATCL_REGISTER_SCALAR_FUNC_SIMP(exp2)
MATCL_REGISTER_SCALAR_FUNC_SIMP(exp10)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log1p)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log1p_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log2)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log2_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log10)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log10_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(logb)
MATCL_REGISTER_SCALAR_FUNC_SIMP(ilogb)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sin)
MATCL_REGISTER_SCALAR_FUNC_SIMP(cos)
MATCL_REGISTER_SCALAR_FUNC_SIMP(tan)
MATCL_REGISTER_SCALAR_FUNC_SIMP(cot)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sec)
MATCL_REGISTER_SCALAR_FUNC_SIMP(csc)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sinh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(cosh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(tanh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(coth)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sech)
MATCL_REGISTER_SCALAR_FUNC_SIMP(csch)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asin)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asin_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acos)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acos_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(atan)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asinh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acosh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acosh_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(atanh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(atanh_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(inv)
MATCL_REGISTER_SCALAR_FUNC_SIMP(invs)
MATCL_REGISTER_SCALAR_FUNC_SIMP(floor)
MATCL_REGISTER_SCALAR_FUNC_SIMP(ceil)
MATCL_REGISTER_SCALAR_FUNC_SIMP(round)
MATCL_REGISTER_SCALAR_FUNC_SIMP(trunc)
MATCL_REGISTER_SCALAR_FUNC_SIMP(ifloor)
MATCL_REGISTER_SCALAR_FUNC_SIMP(iceil)
MATCL_REGISTER_SCALAR_FUNC_SIMP(iround)
MATCL_REGISTER_SCALAR_FUNC_SIMP(itrunc)

MATCL_REGISTER_SCALAR_FUNC(uminus, operator-, mp_int, mdyf::op_uminus)
MATCL_REGISTER_SCALAR_FUNC(op_not, !, mp_int, mdyf::op_not)
MATCL_REGISTER_SCALAR_FUNC(op_bool, (bool), mp_int, mdyf::op_bool)

MATCL_REGISTER_BIN_FUNC(ldexp1, matcl::ldexp, mp_int, Integer, mdyf::ldexp)
MATCL_REGISTER_BIN_FUNC(frexp1, matcl::frexp, mp_int, Integer&, mdyf::frexp)
MATCL_REGISTER_BIN_FUNC(modf1, matcl::modf, mp_int, mp_float&, mdyf::modf)
MATCL_REGISTER_BIN_FUNC(f_mod, matcl::mod, mp_int, mp_int, mdyf::mod)
MATCL_REGISTER_BIN_FUNC(f_rem, matcl::rem, mp_int, mp_int, mdyf::rem)
MATCL_REGISTER_BIN_FUNC(f_dim, matcl::fdim, mp_int, mp_int, mdyf::fdim)

MATCL_REGISTER_TEMPL_1_FUNC(f_factorial,matcl::factorial,mp_int,Integer,mdyf::factorial)
MATCL_REGISTER_TEMPL_1_FUNC(f_double_factorial,double_factorial,mp_int,Integer,mdyf::double_factorial)
MATCL_REGISTER_TEMPL_2_FUNC(f_binomial_coefficient,binomial_coefficient,mp_int,Integer,Integer,
                            mdyf::binomial_coefficient)

//--------------------------------------------------------------
//              conversion
//--------------------------------------------------------------
//promotions

//equivalent
struct conv_int_big : mdy::register_convert<conv_int_big, mdy::convert_equivalent>
{
	static MP_int eval(const OInteger& val)    { return MP_int(val.get()); };
};

//decay

//explicit
struct conv_mpint_int : mdy::register_convert<conv_mpint_int, mdy::convert_explicit>
{
	static OInteger eval(const MP_int& val) { return OInteger(val.get().cast_int()); };
};

//cast
struct conv_float_mpint : mdy::register_convert<conv_float_mpint, mdy::convert_cast>
{
	static MP_int eval(const OFloat& val)   { return MP_int(Integer(val.get())); };
};
struct conv_real_mpint : mdy::register_convert<conv_real_mpint, mdy::convert_cast>
{
	static MP_int eval(const OReal& val)    { return MP_int(Integer(val.get())); };
};
struct conv_mpreal_mpint : mdy::register_convert<conv_mpreal_mpint, mdy::convert_cast>
{
	static MP_int eval(const MP_float& val) { return MP_int(val.get().cast_mp_int()); };
};
struct conv_mprat_mpint : mdy::register_convert<conv_mprat_mpint, mdy::convert_cast>
{
	static MP_int eval(const MP_rational& val) { return MP_int(val.get().cast_mp_int()); };
};

struct conv_fcompl_mpint : mdy::register_convert<conv_fcompl_mpint, mdy::convert_cast>
{
	static MP_int eval(const OFloat_complex& val)
    { 
        return MP_int(convert_scalar<mp_int>(val.get())); 
    };
};
struct conv_compl_mpint : mdy::register_convert<conv_compl_mpint, mdy::convert_cast>
{
	static MP_int eval(const OComplex& val) 
    { 
        return MP_int(convert_scalar<mp_int>(val.get())); 
    };
};
struct conv_mc_mpint : mdy::register_convert<conv_mc_mpint, mdy::convert_cast>
{
	static MP_int eval(const MP_complex& val) 
    { 
        return MP_int(convert_scalar<mp_int>(val.get())); 
    };
};

};