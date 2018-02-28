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

//unary functions from matcl

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

};
