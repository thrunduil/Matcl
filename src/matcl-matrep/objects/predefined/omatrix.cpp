/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-scalar/objects/object_functions.h"
#include "matcl-dynamic/register_function.h"
#include "matcl-dynamic/object_type.h"
#include "matcl-matrep/objects/omatrix.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/IO/matrix_io.h"
#include "matcl-dynamic/details/register_function_macro.h"
#include "matcl-core/details/IO/printer.h"

namespace matcl { namespace dynamic
{

// save and load to stream
bool object_type_traits<Matrix>::read(std::istream& is, T& t)
{
    is >> t;

    if (is.fail() || is.bad())  return false;
    else						return true;
}

void object_type_traits<Matrix>::write(std::ostream& os, const T& t)
{
    os << t;
}

void object_type_traits<Matrix>::load_data(iarchive_impl& ar, T& ret, unsigned int version)
{
    (void)version;
    iarchive iar(ar);
    matcl::load(iar, ret);
};

void object_type_traits<Matrix>::save_data(oarchive_impl& ar, const T& val, unsigned int version)
{
    (void)version;
    oarchive oar(ar);
    save(oar, val);
}

bool object_type_traits<Matrix>::is_zero(const T& t)
{ 
    return t.get_matrix_code() == mat_code::integer_scalar
            && t.get_scalar<Integer>() == 0;
};

bool object_type_traits<Matrix>::is_one(const T& t)
{ 
    return t.get_matrix_code() == mat_code::integer_scalar
            && t.get_scalar<Integer>() == 1;
};

Matrix object_type_traits<Matrix>::make_one(const T*)
{ 
    return Matrix(1);
};

static std::string get_struct_name(struct_code sc)
{
    switch(sc)
    {
        case struct_code::struct_dense:     return "dense";
        case struct_code::struct_sparse:    return "sparse";
        case struct_code::struct_banded:    return "band";
        case struct_code::struct_scalar:    return "scalar";

        //impossible case
        default:                            return "unknown";
    };
};

static std::string get_value_name(value_code vc, dynamic::Type type)
{
    switch (vc)
    {
        case value_code::v_integer:         return "Integer";
        case value_code::v_float:           return "Float";  
        case value_code::v_real:            return "Real";
        case value_code::v_float_complex:   return "Float_complex";
        case value_code::v_complex:         return "Complex";
        case value_code::v_object:          return type.to_string();

        //impossible case
        default:                            return "unknown";
    };
};

static std::string to_string(const Matrix& t)
{
    value_code vc   = t.get_value_code();
    struct_code sc  = t.get_struct_code();
    Integer r       = t.rows();
    Integer c       = t.cols();

    std::string str_name    = get_struct_name(sc);
    std::string val_name    = get_value_name(vc, t.get_type());    

    std::ostringstream buf;
    buf << str_name << " " << val_name << " (" << r << 'x' << c << ")";

    return buf.str();
};

void object_type_traits<Matrix>::disp(const T& t, md::printer& pr, Integer elem_width,
                            align_type at, Integer value_pos)
{
    //disp short header only; printing matrix object is implemented
    //differently in disp function for scalar; this function can be
    //called only if this is a matrix element

    pr.disp_elem(elem_width, to_string(t), at, value_pos);
};

}};

namespace matcl 
{

namespace mdd = matcl::dynamic::details;
namespace mdyf = matcl::dynamic::functions;

//---------------------------------------------------------------
//              special functions
//---------------------------------------------------------------
struct conv_any_mat : dynamic::register_convert<conv_any_mat, dynamic::convert_explicit>
{
    //conversion from any object to matrix; this conversion must be explicit
    //otherwise infinite recurrsions can appear; notice that for a Matrix
    //conversion from any object is implicit
	static OMatrix eval(const dynamic::Any& val)
    {
        return OMatrix(Matrix(val.get().get_stored()));
    };
};

//OMatrix is a unifier
struct unif_mat : dynamic::register_unifier<unif_mat>
{
    static OMatrix eval() {return OMatrix();}
};

//---------------------------------------------------------------
//              unary functions
//---------------------------------------------------------------
#define MATCL_REGISTER_SCALAR_FUNC_SIMP(fn)                             \
MATCL_REGISTER_SCALAR_FUNC(f_##fn, matcl::fn, Matrix, mdyf::fn)

//unary functions from dynamic
MATCL_REGISTER_SCALAR_FUNC_SIMP(real)
MATCL_REGISTER_SCALAR_FUNC_SIMP(imag)

MATCL_REGISTER_SCALAR_FUNC(uminus, -, Matrix, mdyf::op_uminus)
MATCL_REGISTER_SCALAR_FUNC(op_not, !, Matrix, mdyf::op_not)
MATCL_REGISTER_SCALAR_FUNC(op_bool, (bool), Matrix, mdyf::op_bool)

//unary functions from mmlib

MATCL_REGISTER_SCALAR_FUNC(op_neg, neg, Matrix, mdyf::op_neg)
MATCL_REGISTER_SCALAR_FUNC(op_true, is_true, Matrix, mdyf::op_true)

MATCL_REGISTER_SCALAR_FUNC_SIMP(invs)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_nan)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_inf)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_regular)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_finite)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_normal)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_int)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_real)
MATCL_REGISTER_SCALAR_FUNC_SIMP(abs)
MATCL_REGISTER_SCALAR_FUNC_SIMP(abs2)
MATCL_REGISTER_SCALAR_FUNC_SIMP(conj)
MATCL_REGISTER_SCALAR_FUNC_SIMP(arg)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sign)
MATCL_REGISTER_SCALAR_FUNC_SIMP(isign)
MATCL_REGISTER_SCALAR_FUNC_SIMP(eps)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sqrt)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sqrt_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(cbrt)
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
MATCL_REGISTER_SCALAR_FUNC_SIMP(acot)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asec)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asec_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acsc)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acsc_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asinh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acosh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acosh_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(atanh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(atanh_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acoth)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acoth_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asech)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asech_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acsch)
MATCL_REGISTER_SCALAR_FUNC_SIMP(floor)
MATCL_REGISTER_SCALAR_FUNC_SIMP(ceil)
MATCL_REGISTER_SCALAR_FUNC_SIMP(round)
MATCL_REGISTER_SCALAR_FUNC_SIMP(trunc)
MATCL_REGISTER_SCALAR_FUNC_SIMP(ifloor)
MATCL_REGISTER_SCALAR_FUNC_SIMP(iceil)
MATCL_REGISTER_SCALAR_FUNC_SIMP(iround)
MATCL_REGISTER_SCALAR_FUNC_SIMP(itrunc)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sqrt1pm1)

//---------------------------------------------------------------
//              bin functions
//---------------------------------------------------------------
#define REGISTER_BIN_FUNC(fname, ret, func, func_name)                              \
MATCL_REGISTER_BIN_FUNC(fname##1, func, Matrix, Matrix, func_name)                  \
                                                                                    \
struct fname##2 : dynamic::register_function_template_return<fname##2, func_name>   \
{                                                                                   \
	static ret eval(const OMatrix& obj1, const dynamic::Template& obj2)             \
	{                                                                               \
		return ret(func(obj1.get(), obj2.get()));                                   \
	};                                                                              \
    static dynamic::Type eval_return(int n_template, const dynamic::Type* templ,    \
                int n_arg, const dynamic::Type* args)                               \
    {                                                                               \
        (void)n_template; (void)templ; (void)n_arg;                                 \
        if (args[1] == dynamic::Any::get_static_type())                             \
            return dynamic::Type();                                                 \
        else                                                                        \
            return OMatrix::get_static_type();                                      \
    };                                                                              \
};                                                                                  \
struct fname##3 : dynamic::register_function_template_return<fname##3, func_name>   \
{                                                                                   \
	static ret eval(const dynamic::Template& obj1, const OMatrix& obj2)             \
	{                                                                               \
		return ret(func(obj1.get(), obj2.get()));                                   \
	};                                                                              \
    static dynamic::Type eval_return(int n_template, const dynamic::Type* templ,    \
                int n_arg, const dynamic::Type* args)                               \
    {                                                                               \
        (void)n_template; (void)templ; (void)n_arg;                                 \
        if (args[0] == dynamic::Any::get_static_type())                             \
            return dynamic::Type();                                                 \
        else                                                                        \
            return OMatrix::get_static_type();                                      \
    };                                                                              \
};                                                                                  \

REGISTER_BIN_FUNC(f_plus, OMatrix, operator+, mdyf::op_plus)
REGISTER_BIN_FUNC(f_minus, OMatrix, operator-, mdyf::op_minus)
REGISTER_BIN_FUNC(f_mul, OMatrix, mul, mdyf::elem_mul)
REGISTER_BIN_FUNC(f_div, OMatrix, operator/, mdyf::op_div)
REGISTER_BIN_FUNC(f_idiv, OMatrix, idiv, mdyf::idiv)
REGISTER_BIN_FUNC(f_div_0, OMatrix, div_0, mdyf::div_0)
REGISTER_BIN_FUNC(f_div_1, OMatrix, div_1, mdyf::div_1)
REGISTER_BIN_FUNC(f_atan2, OMatrix, atan2, mdyf::atan2)
REGISTER_BIN_FUNC(f_hypot, OMatrix, hypot, mdyf::hypot)
REGISTER_BIN_FUNC(f_mod, OMatrix, mod, mdyf::mod)
REGISTER_BIN_FUNC(f_rem, OMatrix, rem, mdyf::rem)
REGISTER_BIN_FUNC(f_pow, OMatrix, pow, mdyf::pow)
REGISTER_BIN_FUNC(f_pow_c, OMatrix, pow_c, mdyf::pow_c)
REGISTER_BIN_FUNC(f_elem_and, OMatrix, elem_and, mdyf::elem_and)
REGISTER_BIN_FUNC(f_elem_or, OMatrix, elem_or, mdyf::elem_or)
REGISTER_BIN_FUNC(f_elem_xor, OMatrix, elem_xor, mdyf::elem_xor)
REGISTER_BIN_FUNC(f_eeq, OMatrix, eeq, mdyf::op_eeq)
REGISTER_BIN_FUNC(f_eeq_nan, OMatrix, eeq_nan, mdyf::eeq_nan)
REGISTER_BIN_FUNC(f_neq, OMatrix, neq, mdyf::op_neq)
REGISTER_BIN_FUNC(f_neq_nan, OMatrix, neq_nan, mdyf::neq_nan)
REGISTER_BIN_FUNC(f_leq, OMatrix, leq, mdyf::op_leq)
REGISTER_BIN_FUNC(f_geq, OMatrix, geq, mdyf::op_geq)
REGISTER_BIN_FUNC(f_gt, OMatrix, gt, mdyf::op_gt)
REGISTER_BIN_FUNC(f_lt, OMatrix, lt, mdyf::op_lt)
REGISTER_BIN_FUNC(f_max, OMatrix, max, mdyf::max)
REGISTER_BIN_FUNC(f_min, OMatrix, min, mdyf::min)
REGISTER_BIN_FUNC(f_op_and, dynamic::OBool, op_and, mdyf::op_and)
REGISTER_BIN_FUNC(f_op_or, dynamic::OBool, op_or, mdyf::op_or)
REGISTER_BIN_FUNC(f_op_xor, dynamic::OBool, op_xor, mdyf::op_xor)

//not defined for a matrix:
//nextabove, nextbelow, signbit, fpclassify, ldexp, frexp, modf, inv, powm1, fdim
//copysign, nextafter, fma, dot2_ac

};
