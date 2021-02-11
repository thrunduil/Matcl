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

#include "test_gmp.h"
#include "matcl-mp/matcl_mp.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-scalar/objects/object_functions.h"
#include "matcl-scalar/lib_functions/scalar_gen.h"

#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-scalar/lib_functions/func_matrix.h"
#include "matcl-scalar/lib_functions/func_forwarding.h"
#include "matcl-scalar/IO/scalar_io.h"

#include "rand_scalars.h"
#include "eval_cons.h"
#include "utils.h"
#include <iostream>

#pragma warning(push)
#pragma warning(disable:4127) // conditional expression is constant

namespace matcl { namespace test
{

namespace mdy = matcl::dynamic;

template<class Derived>
struct eval_operator_arithm_obj : eval_scalars<eval_operator_arithm_obj<Derived>>
{
    Integer code;
    bool    m_div;
    bool    complex_allowed;

    eval_operator_arithm_obj(Integer c, bool div, bool complex_allowed_) 
        : code(c), m_div(div), complex_allowed(complex_allowed_)
    {};

    template<class T1, class T2>
    auto eval_op(const T1& a1, const T2& a2) -> decltype(Derived::eval(std::declval<T1>(), std::declval<T2>()))
    {
        return Derived::eval(a1,a2);
    };

    template<class T1, class T2>
    bool different(const T1& a, const T2& b)
    {
        if (a == b)
            return false;

        if ((bool)is_nan(a) == true && (bool)is_nan(b) == true)
            return false;

        auto dif    = abs(a - b);
        auto tol    = eps(b);

        if (dif <= tol)
            return false;

        out_stream << a << " " << b << "\n";
        out_stream << dif << " " << tol << "\n";
        return true;
    };

    template<class T1, class T2>
    double eval_scal_func(const T1& s1, const T2& s2)
    {
        //out_stream << code << "\n";
        if (code == -1)
            disp("break");

        if ((md::is_complex<T1>::value || md::is_complex<T2>::value) && complex_allowed == false)
            return 0.0;

        using res_type  = decltype(eval_op(s1,s2));

        bool int_res    = std::is_same<res_type, Integer>::value
                        ||std::is_same<res_type, mp_int>::value
                        || std::is_same<res_type, mp_rational>::value;

        if (m_div == true)
        {
            if (is_zero(s2) && int_res == true)
            {
                return 0.;
            }
        };

        auto res_1      = eval_op(s1,s2);

        using type_1    = object_type<T1>;
        using type_2    = object_type<T2>;
        using type_res  = object_type<decltype(res_1)>;

        type_1 o1(s1);
        type_2 o2(s2);

        auto res_2      = eval_op(o1, o2);
        auto res_21     = eval_op(s1, o2);
        auto res_22     = eval_op(o1, s2);

        matcl::dynamic::object oo1      = matcl::dynamic::object(o1);
        matcl::dynamic::object oo2      = matcl::dynamic::object(o2);

        auto res_3      = eval_op(oo1, oo2);

        double out      = 0;
        if (different(res_1, res_2))
            out         += 1;
        if (different(res_1, res_21))
            out         += 1;
        if (different(res_1, res_22))
            out         += 1;
        if (different(matcl::dynamic::object(type_res(res_1)), res_3))
            out         += 1;
        if (res_3.get_type() != res_2.get_static_type())
            out         += 1;
        if (res_2.get_static_type() != res_21.get_static_type())
            out         += 1;
        if (res_2.get_static_type() != res_22.get_static_type())
            out         += 1;
        if (type_res::get_static_type() != res_2.get_static_type())
            out         += 1;
        if (different(matcl::dynamic::object(type_res(res_1)) , res_3))
            out         += 1;

        if (out != 0.0)
        {
            out_stream << s1 << " " << s2 << "\n";
            out_stream << code << "\n";
        }
        return out;
    };
};

struct eval_add_obj : eval_operator_arithm_obj<eval_add_obj>
{
    eval_add_obj(Integer c) : eval_operator_arithm_obj(c,false, true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(std::declval<T1>() + std::declval<T2>())
    {
        return a1 + a2;
    };

    static auto eval(const Integer& a1, const Float& a2) -> Real
    {
        return a1 + a2;
    };
    static auto eval(const Float& a1, const Integer & a2) -> Real
    {
        return a1 + a2;
    };

};

struct eval_sub_obj : eval_operator_arithm_obj<eval_sub_obj>
{
    eval_sub_obj(Integer c) : eval_operator_arithm_obj(c,false, true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(std::declval<T1>() - std::declval<T2>())
    {
        return a1 - a2;
    };

    static auto eval(const Integer& a1, const Float& a2) -> Real
    {
        return a1 - a2;
    };
    static auto eval(const Float& a1, const Integer & a2) -> Real
    {
        return a1 - a2;
    };
};

struct eval_mul_obj : eval_operator_arithm_obj<eval_mul_obj>
{
    eval_mul_obj(Integer c) : eval_operator_arithm_obj(c,false, true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(mul(a1,a2))
    {
        return mul(a1, a2);
    };
};

struct eval_mmul_obj : eval_operator_arithm_obj<eval_mmul_obj>
{
    eval_mmul_obj(Integer c) : eval_operator_arithm_obj(c,false, true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(std::declval<T1>() * std::declval<T2>())
    {
        return a1 * a2;
    };

    static auto eval(const Integer& a1, const Float& a2) -> Real
    {
        return a1 * a2;
    };
    static auto eval(const Float& a1, const Integer & a2) -> Real
    {
        return a1 * a2;
    };
};

struct eval_div_obj : eval_operator_arithm_obj<eval_div_obj>
{
    eval_div_obj(Integer c) : eval_operator_arithm_obj(c,true, true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(std::declval<T1>() / std::declval<T2>())
    {
        return a1 / a2;
    };

    static auto eval(const Integer& a1, const Integer& a2) -> Real
    {
        return Real(a1) / Real(a2);
    };

    static auto eval(const Integer& a1, const Float& a2) -> Real
    {
        return Real(a1) / Real(a2);
    };
    static auto eval(const Float& a1, const Integer & a2) -> Real
    {
        return Real(a1) / Real(a2);
    };
};
struct eval_div2_obj : eval_operator_arithm_obj<eval_div2_obj>
{
    eval_div2_obj(Integer c) : eval_operator_arithm_obj(c,true, true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(div(std::declval<T1>() , std::declval<T2>()))
    {
        return div(a1 , a2);
    };
};
struct eval_div_0_obj : eval_operator_arithm_obj<eval_div_0_obj>
{
    eval_div_0_obj(Integer c) : eval_operator_arithm_obj(c,true, true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(div_0(std::declval<T1>() , std::declval<T2>()))
    {
        return div_0(a1 , a2);
    };
};
struct eval_div_1_obj : eval_operator_arithm_obj<eval_div_1_obj>
{
    eval_div_1_obj(Integer c) : eval_operator_arithm_obj(c,true, true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(div_1(std::declval<T1>() , std::declval<T2>()))
    {
        return div_1(a1 , a2);
    };
};

struct eval_idiv_obj : eval_operator_arithm_obj<eval_idiv_obj>
{
    eval_idiv_obj(Integer c) : eval_operator_arithm_obj(c,true, true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(idiv(std::declval<T1>() , std::declval<T2>()))
    {
        return idiv(a1 , a2);
    };
};
struct eval_pow_obj : eval_operator_arithm_obj<eval_pow_obj>
{
    eval_pow_obj(Integer c) : eval_operator_arithm_obj(c,false, true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(pow(std::declval<T1>() , std::declval<T2>()))
    {
        return pow(a1 , a2);
    };
};
struct eval_pow_c_obj : eval_operator_arithm_obj<eval_pow_c_obj>
{
    eval_pow_c_obj(Integer c) : eval_operator_arithm_obj(c,false, true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(pow_c(std::declval<T1>() , std::declval<T2>()))
    {
        return pow_c(a1 , a2);
    };
};

double gmp_tester_bin::test_add_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_add_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_sub_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_sub_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_mul_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_mul_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_mmul_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_mmul_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_div_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_div_obj test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_div2_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_div2_obj test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_div_0_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_div_0_obj test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_div_1_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_div_1_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_idiv_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_idiv_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_pow_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_pow_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_pow_c_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_pow_c_obj test(code);
    return test.make(s1, s2);
};

struct eval_copysign_obj : eval_operator_arithm_obj<eval_copysign_obj>
{
    eval_copysign_obj(Integer c) : eval_operator_arithm_obj(c,false,false){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> decltype(copysign(std::declval<T1>(), std::declval<T2>()))
    {
        return copysign(a1, a2);
    };

    template<class T1, class T2, class Enable = typename enable_complex_nobj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> Real
    {
        (void)a1;
        (void)a2;
        return Real(0);
    };

    template<class T1, class T2, class Enable = typename enable_complex_obj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> OReal
    {
        (void)a1;
        (void)a2;
        return OReal(0);
    };
};

struct eval_nextafter_obj : eval_operator_arithm_obj<eval_nextafter_obj>
{
    eval_nextafter_obj(Integer c) : eval_operator_arithm_obj(c,false,false){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> decltype(nextafter(std::declval<T1>(), std::declval<T2>()))
    {
        return nextafter(a1, a2);
    };

    template<class T1, class T2, class Enable = typename enable_complex_nobj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> Real
    {
        (void)a1;
        (void)a2;
        return Real(0);
    };

    template<class T1, class T2, class Enable = typename enable_complex_obj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> OReal
    {
        (void)a1;
        (void)a2;
        return OReal(0);
    };
};

struct eval_float_distance_obj : eval_operator_arithm_obj<eval_float_distance_obj>
{
    eval_float_distance_obj(Integer c) : eval_operator_arithm_obj(c,false,false){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) 
        -> decltype(float_distance(std::declval<T1>(), std::declval<T2>()))
    {
        return float_distance(a1, a2);
    };

    template<class T1, class T2, class Enable = typename enable_complex_nobj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> Real
    {
        (void)a1;
        (void)a2;
        return Real(0);
    };

    template<class T1, class T2, class Enable = typename enable_complex_obj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> OReal
    {
        (void)a1;
        (void)a2;
        return OReal(0);
    };
};

struct eval_hypot_obj : eval_operator_arithm_obj<eval_hypot_obj>
{
    eval_hypot_obj(Integer c) : eval_operator_arithm_obj(c,false, true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(hypot(std::declval<T1>(), std::declval<T2>()))
    {
        return hypot(a1, a2);
    };
};

struct eval_atan2_obj : eval_operator_arithm_obj<eval_atan2_obj>
{
    eval_atan2_obj(Integer c) : eval_operator_arithm_obj(c,false,false){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> decltype(atan2(std::declval<T1>(), std::declval<T2>()))
    {
        return atan2(a1, a2);
    };

    template<class T1, class T2, class Enable = typename enable_complex_nobj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> Real
    {
        (void)a1;
        (void)a2;
        return Real(0);
    };

    template<class T1, class T2, class Enable = typename enable_complex_obj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> OReal
    {
        (void)a1;
        (void)a2;
        return OReal(0);
    };

};

struct eval_mod_obj : eval_operator_arithm_obj<eval_mod_obj>
{
    eval_mod_obj(Integer c) : eval_operator_arithm_obj(c,false,false){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> decltype(mod(std::declval<T1>(), std::declval<T2>()))
    {
        return mod(a1, a2);
    };

    template<class T1, class T2, class Enable = typename enable_complex_nobj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> Real
    {
        (void)a1;
        (void)a2;
        return Real(0);
    };

    template<class T1, class T2, class Enable = typename enable_complex_obj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> OReal
    {
        (void)a1;
        (void)a2;
        return OReal(0);
    };
};

struct eval_rem_obj : eval_operator_arithm_obj<eval_rem_obj>
{
    eval_rem_obj(Integer c) : eval_operator_arithm_obj(c,false,false){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> decltype(rem(std::declval<T1>(), std::declval<T2>()))
    {
        return rem(a1, a2);
    };

    template<class T1, class T2, class Enable = typename enable_complex_nobj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> Real
    {
        (void)a1;
        (void)a2;
        return Real(0);
    };

    template<class T1, class T2, class Enable = typename enable_complex_obj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> OReal
    {
        (void)a1;
        (void)a2;
        return OReal(0);
    };
};

struct eval_max_obj : eval_operator_arithm_obj<eval_max_obj>
{
    eval_max_obj(Integer c) : eval_operator_arithm_obj(c,false,true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(max(std::declval<T1>(), std::declval<T2>()))
    {
        return max(a1, a2);
    };
};
struct eval_min_obj : eval_operator_arithm_obj<eval_min_obj>
{
    eval_min_obj(Integer c) : eval_operator_arithm_obj(c,false,true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(min(std::declval<T1>(), std::declval<T2>()))
    {
        return min(a1, a2);
    };
};

struct eval_elem_or_obj : eval_operator_arithm_obj<eval_elem_or_obj>
{
    eval_elem_or_obj(Integer c) : eval_operator_arithm_obj(c,false,true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(elem_or(std::declval<T1>(), std::declval<T2>()))
    {
        return elem_or(a1, a2);
    };
};
struct eval_elem_xor_obj : eval_operator_arithm_obj<eval_elem_xor_obj>
{
    eval_elem_xor_obj(Integer c) : eval_operator_arithm_obj(c,false,true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(elem_xor(std::declval<T1>(), std::declval<T2>()))
    {
        return elem_xor(a1, a2);
    };
};
struct eval_elem_and_obj : eval_operator_arithm_obj<eval_elem_and_obj>
{
    eval_elem_and_obj(Integer c) : eval_operator_arithm_obj(c,false,true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(elem_and(std::declval<T1>(), std::declval<T2>()))
    {
        return elem_and(a1, a2);
    };
};

struct eval_fdim_obj : eval_operator_arithm_obj<eval_fdim_obj>
{
    eval_fdim_obj(Integer c) : eval_operator_arithm_obj(c,false,false){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> decltype(fdim(std::declval<T1>(), std::declval<T2>()))
    {
        return fdim(a1, a2);
    };

    template<class T1, class T2, class Enable = typename enable_complex_nobj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> Real
    {
        (void)a1;
        (void)a2;
        return Real(0);
    };

    template<class T1, class T2, class Enable = typename enable_complex_obj<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> OReal
    {
        (void)a1;
        (void)a2;
        return OReal(0);
    };
};

double gmp_tester_bin::test_copysign_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_copysign_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_nextafter_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_nextafter_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_float_distance_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_float_distance_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_hypot_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_hypot_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_atan2_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_atan2_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_mod_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_mod_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_rem_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_rem_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_max_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_max_obj test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_min_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_min_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_elem_or_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_elem_or_obj test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_elem_xor_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_elem_xor_obj test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_elem_and_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_elem_and_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_fdim_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_fdim_obj test(code);
    return test.make(s1, s2);
};

}};

#pragma warning(pop)