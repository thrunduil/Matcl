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

#include "test_gmp.h"
#include "matcl-mp/matcl_mp.h"
#include "matcl-core/IO/archive.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-scalar/objects/object_functions.h"
#include "matcl-scalar/lib_functions/scalar_gen.h"

#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-scalar/lib_functions/func_matrix.h"
#include "matcl-scalar/lib_functions/func_forwarding.h"

#include "rand_scalars.h"
#include "eval_cons.h"
#include "utils.h"
#include <iostream>

namespace matcl { namespace test
{

namespace mdy = matcl::dynamic;

template<class Derived>
struct eval_operator_arithm : eval_scalars<eval_operator_arithm<Derived>>
{
    Integer code;
    bool    is_div;
    bool    is_idiv;

    eval_operator_arithm(Integer c_, bool is_div_, bool is_idiv_) 
        : code(c_), is_div(is_div_), is_idiv(is_idiv_){};

    template<class T1, class T2>
    auto eval_op(const T1& a1, const T2& a2) -> decltype(Derived::eval(std::declval<T1>(), std::declval<T2>()))
    {
        return Derived::eval(a1,a2);
    };

    template<class T1, class T2>
    bool different(const T1& a, const T2& b)
    {
        if (is_div == false)
            return a != b;

        if (a == b)
            return false;

        mp_complex v1(a);
        mp_complex v2(b, get_prec<T2>());

        mp_float v11(abs(a));
        mp_float v21(abs(b), get_prec<T2>());

        mp_float dif    = abs(v1 - v2);
        mp_float tol    = 2 * (eps(v11) + eps(v21));

        if (dif < tol)
            return false;

        if (is_nan(v1) && is_nan(v2))
            return false;

        if (code == -1)
        {
            std::cout << a << " " << b << "\n";
            std::cout << v1 << " " << v2 << "\n";
            std::cout << dif << " " << tol << "\n";
        };

        return true;
    };

    template<class T1, class T2>
    double eval_scal_func(const T1& s1, const T2& s2)
    {
        //std::cout << code << "\n";
        if (code == -1)
            std::cout << "break\n";

        auto res        = eval_op(s1,s2);
        bool is_real_1  = matcl::imag(s1) == 0;
        bool is_real_2  = matcl::imag(s2) == 0;
        bool is_fin_1   = matcl::is_finite(s1) && is_real_1;
        bool is_fin_2   = matcl::is_finite(s2) && is_real_2;
        bool is_fin     = matcl::is_finite(res);
        bool is_zero_1  = (s1 == 0) || is_fin_1 == false;
        bool is_zero_2  = (s2 == 0) || is_fin_2 == false;
        bool is_int_1   = std::is_same<T1,Integer>::value;
        bool is_int_2   = std::is_same<T2,Integer>::value;

        mp_int      v1_1    = convert_scalar<mp_int>(s1);
        mp_float    v1_2    = convert_scalar<mp_float>(s1);
        mp_rational v1_3    = convert_scalar<mp_rational>(s1);
        mp_complex  v1_4    = convert_scalar<mp_complex>(s1);

        mp_int v2_1         = convert_scalar<mp_int>(s2);
        mp_float v2_2       = convert_scalar<mp_float>(s2);
        mp_rational v2_3    = convert_scalar<mp_rational>(s2);
        mp_complex v2_4     = convert_scalar<mp_complex>(s2);

        double out = 0;

        if ((is_div == true && is_zero_2 == true && is_int_2 == true) == false
            && (is_idiv == true && is_int_2 == true && is_int_1 == false) == false)
        {
            if (different(eval_op(v1_1, s2), res) && is_fin_1 == true)
                out     += 1;
        }
        if (different(eval_op(v1_2, s2), res) && is_real_1 == true
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        }
        if ((is_div == true && is_zero_2 == true && is_int_2 == true) == false
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            if (different(eval_op(v1_3, s2), res) && is_fin_1 == true)
                out     += 1;
        };
        if (different(eval_op(v1_4, s2), res)
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        };

        if ((is_div == true && is_zero_2 == true) == false
            && (is_idiv == true && is_int_2 == false && is_int_1 == true) == false)
        {
            if (different(eval_op(s1, v2_1), res) && is_fin_2 == true)
                out     += 1;
        }
        if (different(eval_op(s1, v2_2), res) && is_real_2 == true
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        }
        if ((is_div == true && is_zero_2 == true) == false
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            if (different(eval_op(s1, v2_3), res) && is_fin_2 == true)
                out     += 1;
        }
        if (different(eval_op(s1, v2_4), res)
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        };

        if ((is_div == true && is_zero_2 == true) == false
            && (is_idiv == true && (is_int_2 == false || is_int_1 == false)) == false)
        {
            if (different(eval_op(v1_1, v2_1), res) && is_fin_1 == true && is_fin_2 == true)
                out     += 1;
        };
        if (different(eval_op(v1_2, v2_1), res) && is_real_1 == true && is_fin_2 == true
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        }
        if ((is_div == true && is_zero_2 == true) == false
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            if (different(eval_op(v1_3, v2_1), res) && is_fin_1 == true && is_fin_2 == true)
                out     += 1;
        };
        if (different(eval_op(v1_4, v2_1), res) && is_fin_2 == true
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        };

        if (different(eval_op(v1_1, v2_2), res) && is_fin_1 == true && is_real_2 == true
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        }
        if (different(eval_op(v1_2, v2_2), res) && is_real_1 == true && is_real_2 == true
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        }
        if (different(eval_op(v1_3, v2_2), res) && is_fin_1 == true && is_real_2 == true
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        }
        if (different(eval_op(v1_4, v2_2), res) && is_real_2 == true
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        }

        if ((is_div == true && is_zero_2 == true) == false
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            if (different(eval_op(v1_1, v2_3), res) && is_fin_1 == true && is_fin_2 == true)
                out     += 1;
        }
        if (different(eval_op(v1_2, v2_3), res) && is_real_1 == true && is_fin_2 == true
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        }
        if ((is_div == true && is_zero_2 == true) == false
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            if (different(eval_op(v1_3, v2_3), res) && is_fin_1 == true&& is_fin_2 == true)
                out     += 1;
        }
        if (different(eval_op(v1_4, v2_3), res) && (is_fin_2 == true)
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        };

        if (different(eval_op(v1_1, v2_4), res) && is_fin_1 == true
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        }
        if (different(eval_op(v1_2, v2_4), res) && is_real_1 == true
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        }
        if (different(eval_op(v1_3, v2_4), res) && is_fin_1 == true
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        }
        if (different(eval_op(v1_4, v2_4), res)
            && (is_idiv == true && is_int_2 == true && is_int_1 == true) == false)
        {
            out     += 1;
        }

        if (is_fin == false)
            out     = 0;

        if (out != 0)
            std::cout << code << " " << s1 << " " << s2 << "\n";

        return out;
    };
};


struct eval_add : eval_operator_arithm<eval_add>
{
    eval_add(Integer c) : eval_operator_arithm(c,false,false){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(std::declval<T1>() + std::declval<T2>())
    {
        return a1 + a2;
    };
};

struct eval_sub : eval_operator_arithm<eval_sub>
{
    eval_sub(Integer c) : eval_operator_arithm(c,false,false){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(std::declval<T1>() - std::declval<T2>())
    {
        return a1 - a2;
    };
};

struct eval_mul : eval_operator_arithm<eval_mul>
{
    eval_mul(Integer c) : eval_operator_arithm(c,false,false){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(mul(a1,a2))
    {
        return mul(a1, a2);
    };
};

struct eval_mmul : eval_operator_arithm<eval_mmul>
{
    eval_mmul(Integer c) : eval_operator_arithm(c,false,false){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(std::declval<T1>() * std::declval<T2>())
    {
        return a1 * a2;
    };
};

struct eval_div : eval_operator_arithm<eval_div>
{
    eval_div(Integer c) : eval_operator_arithm(c,true,false){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(std::declval<T1>() / std::declval<T2>())
    {
        return a1 / a2;
    };

    static auto eval(Integer a1, Integer a2) -> decltype(std::declval<Real>() / std::declval<Real>())
    {
        return Real(a1) / Real(a2);
    };
    static auto eval(Integer a1, Float a2) -> decltype(std::declval<Real>() / std::declval<Real>())
    {
        return Real(a1) / Real(a2);
    };
    static auto eval(Float a1, Integer a2) -> decltype(std::declval<Real>() / std::declval<Real>())
    {
        return Real(a1) / Real(a2);
    };
};
struct eval_div2 : eval_operator_arithm<eval_div2>
{
    eval_div2(Integer c) : eval_operator_arithm(c,true,false){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(div(std::declval<T1>() , std::declval<T2>()))
    {
        return div(a1 , a2);
    };
};
struct eval_div_0 : eval_operator_arithm<eval_div_0>
{
    eval_div_0(Integer c) : eval_operator_arithm(c,true,false){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(div_0(std::declval<T1>() , std::declval<T2>()))
    {
        return div_0(a1 , a2);
    };
};
struct eval_div_1 : eval_operator_arithm<eval_div_1>
{
    eval_div_1(Integer c) : eval_operator_arithm(c,true,false){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(div_1(std::declval<T1>(), std::declval<T2>()))
    {
        return matcl::div_1(a1 , a2);
    };
};

struct eval_idiv : eval_operator_arithm<eval_idiv>
{
    eval_idiv(Integer c) : eval_operator_arithm(c,true,true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(idiv(std::declval<T1>() , std::declval<T2>()))
    {
        return idiv(a1,a2);
    };
};

double gmp_tester_bin::test_add(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_add test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_sub(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_sub test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_mul(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_mul test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_mmul(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_mmul test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_div(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_div test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_div2(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_div2 test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_div_0(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_div_0 test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_div_1(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_div_1 test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_idiv(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_idiv test(code);
    return test.make(s1, s2);
};

}};
