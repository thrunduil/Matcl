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

#include "rand_scalars.h"
#include "matcl-scalar/lib_functions/scalar_gen.h"
#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-scalar/lib_functions/func_matrix.h"
#include "matcl-scalar/lib_functions/func_forwarding.h"
#include "matcl-mp/matcl_mp.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-scalar/objects/object_functions.h"
#include "matcl-scalar/IO/scalar_io.h"

#include "utils.h"

#include <iostream>

namespace matcl { namespace test
{

namespace mdy = matcl::dynamic;

template<class Derived>
struct eval_binfunc : eval_scalars<eval_binfunc<Derived>>
{
    Integer code;
    bool    use_signed_zero;
    bool    zeroinf_allowed;
    bool    intint_different;
    bool    complex_allowed;
    bool    compare_nonfinite;
    double  mult_tol;

    eval_binfunc(Integer c_, bool use_signed_zero_, bool zeroinf_allowed_, bool intint_different_,
                 bool complex_allowed_, double mult, bool cmp_inf) 
        : code(c_), use_signed_zero(use_signed_zero_), zeroinf_allowed(zeroinf_allowed_)
        , intint_different(intint_different_), complex_allowed(complex_allowed_), mult_tol(mult)
        , compare_nonfinite(cmp_inf)
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
        if (is_nan(a) && is_nan(b))
            return false;

        mp_complex v1(a);
        mp_complex v2(b, get_prec<T2>());

        mp_float v11(abs(a));
        mp_float v21(abs(b), get_prec<T2>());

        mp_float dif    = abs(v1 - v2);
        mp_float tol    = (2 * mult_tol) * (eps(v11) + eps(v21));

        if (dif <= tol)
            return false;

        if (code == -1)
        {
            out_stream << a << " " << b << "\n";
            out_stream << v1 << " " << v2 << "\n";
            out_stream << dif << " " << tol << "\n";
        };

        return true;
    };

    template<class T1, class T2>
    double eval_scal_func(const T1& s1, const T2& s2)
    {
        //out_stream << code << "\n";
        if (code == -1)
            disp("break");

        bool is_real_1  = matcl::imag(s1) == 0;
        bool is_real_2  = matcl::imag(s2) == 0;

        bool is_compl_1 = md::is_complex<T1>::value;
        bool is_compl_2 = md::is_complex<T2>::value;

        if (complex_allowed == false && (is_compl_1 == true || is_compl_2 == true))
            return 0.0;

        auto res        = eval_op(s1,s2);
        bool test_int_1 = use_signed_zero == false || s1 != 0 || (s1 == 0 && signbit(real(s1)) == false);
        bool test_int_2 = use_signed_zero == false || s2 != 0 || (s2 == 0 && signbit(real(s2)) == false);
        bool is_fin_1   = matcl::is_finite(s1) && is_real_1 && test_int_1;
        bool is_fin_2   = matcl::is_finite(s2) && is_real_2 && test_int_2;
        bool is_fin     = matcl::is_finite(res);
        bool is_zero_1  = (s1 == 0) || is_fin_1 == false;
        bool is_zero_2  = (s2 == 0) || is_fin_2 == false;
        bool is_int_1   = std::is_same<T1,Integer>::value;
        bool is_int_2   = std::is_same<T2,Integer>::value;
        bool both_int   = is_int_1 && is_int_2;

        (void)is_zero_1;
        (void)is_zero_2;

        mp_int      v1_1    = convert_scalar<mp_int>(s1);
        mp_float    v1_2    = convert_scalar<mp_float>(s1);
        mp_rational v1_3    = convert_scalar<mp_rational>(s1);
        mp_complex  v1_4    = convert_scalar<mp_complex>(s1);

        mp_int v2_1         = convert_scalar<mp_int>(s2);
        mp_float v2_2       = convert_scalar<mp_float>(s2);
        mp_rational v2_3    = convert_scalar<mp_rational>(s2);
        mp_complex v2_4     = convert_scalar<mp_complex>(s2);

        bool zeroinf_s1 = is_zero(real(s1)) || is_finite(real(s1)) == false;
        bool zeroinf_1  = is_zero(v1_1) || is_finite(v1_1) == false;
        bool zeroinf_2  = is_zero(v1_2) || is_finite(v1_2) == false;
        bool zeroinf_3  = is_zero(v1_3) || is_finite(v1_3) == false;
        bool zeroinf_4  = is_zero(real(v1_4)) || is_finite(real(v1_4)) == false;

        bool test_intint    = (intint_different == false) || (both_int == true);
        bool test_intint_1  = (intint_different == false) || (is_int_1 == false || both_int == true);
        bool test_intint_2  = (intint_different == false) || (is_int_2 == false || both_int == true);
        bool test_f         = (intint_different == false) || (both_int == false)|| true;

        bool test_s1    = zeroinf_allowed == true || zeroinf_s1 == false;
        bool test_1     = zeroinf_allowed == true || zeroinf_1 == false;
        bool test_2     = zeroinf_allowed == true || zeroinf_2 == false;
        bool test_3     = zeroinf_allowed == true || zeroinf_3 == false;
        bool test_4     = zeroinf_allowed == true || zeroinf_4 == false;

        double out = 0;

        if (different(eval_op(v1_1, s2), res) && is_fin_1 == true && test_1 && test_intint_2)
            out     += 1;
        if (different(eval_op(v1_2, s2), res) && is_real_1 == true && test_2 && test_f)
            out     += 1;
        if (different(eval_op(v1_3, s2), res) && is_fin_1 == true && test_3 && test_f)
            out     += 1;
        if (different(eval_op(v1_4, s2), res) && test_4 && test_f && complex_allowed)
            out     += 1;

        if (different(eval_op(s1, v2_1), res) && is_fin_2 == true && test_s1 && test_intint_1)
            out     += 1;
        if (different(eval_op(s1, v2_2), res) && is_real_2 == true && test_s1 && test_f)
            out     += 1;
        if (different(eval_op(s1, v2_3), res) && is_fin_2 == true && test_s1 && test_f)
            out     += 1;
        if (different(eval_op(s1, v2_4), res) && test_s1 && test_f && complex_allowed)
            out     += 1;

        if (different(eval_op(v1_1, v2_1), res) && is_fin_1 == true && is_fin_2 == true && test_1
            && test_intint)
            out     += 1;
        if (different(eval_op(v1_2, v2_1), res) && is_real_1 == true && is_fin_2 == true && test_2 && test_f)
            out     += 1;
        if (different(eval_op(v1_3, v2_1), res) && is_fin_1 == true && is_fin_2 == true && test_3 && test_f)
            out     += 1;
        if (different(eval_op(v1_4, v2_1), res) && is_fin_2 == true && test_4 && test_f && complex_allowed)
            out     += 1;

        if (different(eval_op(v1_1, v2_2), res) && is_fin_1 == true && is_real_2 == true && test_1 && test_f)
            out     += 1;
        if (different(eval_op(v1_2, v2_2), res) && is_real_1 == true && is_real_2 == true && test_2 && test_f)
            out     += 1;
        if (different(eval_op(v1_3, v2_2), res) && is_fin_1 == true && is_real_2 == true && test_3 && test_f)
            out     += 1;
        if (different(eval_op(v1_4, v2_2), res) && is_real_2 == true && test_4 && test_f && complex_allowed)
            out     += 1;

        if (different(eval_op(v1_1, v2_3), res) && is_fin_1 == true && is_fin_2 == true && test_1 && test_f)
            out     += 1;
        if (different(eval_op(v1_2, v2_3), res) && is_real_1 == true && is_fin_2 == true && test_2 && test_f)
            out     += 1;
        if (different(eval_op(v1_3, v2_3), res) && is_fin_1 == true&& is_fin_2 == true && test_3 && test_f)
            out     += 1;
        if (different(eval_op(v1_4, v2_3), res) && (is_fin_2 == true) && test_4 && test_f && complex_allowed)
            out     += 1;

        if (different(eval_op(v1_1, v2_4), res) && is_fin_1 == true && test_1 && test_f && complex_allowed)
            out     += 1;
        if (different(eval_op(v1_2, v2_4), res) && is_real_1 == true && test_2 && test_f && complex_allowed)
            out     += 1;
        if (different(eval_op(v1_3, v2_4), res) && is_fin_1 == true && test_3 && test_f && complex_allowed)
            out     += 1;
        if (different(eval_op(v1_4, v2_4), res) && test_4 && test_f && complex_allowed)
            out     += 1;

        if (is_fin == false && compare_nonfinite == false)
            out     = 0;

        if (out != 0)
            out_stream << code << " " << s1 << " " << s2 << "\n";

        return out;
    };
};

struct eval_atan2 : eval_binfunc<eval_atan2>
{
    eval_atan2(Integer c) : eval_binfunc(c,true,true,false,false,1.0,true){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> decltype(atan2(std::declval<T1>() , std::declval<T2>()))
    {
        return atan2(a1,a2);
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
struct eval_hypot : eval_binfunc<eval_hypot>
{
    eval_hypot(Integer c) : eval_binfunc(c,false,true,false,true,1.0,true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(hypot(std::declval<T1>() , std::declval<T2>()))
    {
        return hypot(a1,a2);
    };
};
struct eval_mod : eval_binfunc<eval_mod>
{
    eval_mod(Integer c) : eval_binfunc(c,false,true,false,false,1.0,true){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> decltype(mod(std::declval<T1>() , std::declval<T2>()))
    {
        return mod(a1,a2);
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
struct eval_rem : eval_binfunc<eval_rem>
{
    eval_rem(Integer c) : eval_binfunc(c,false,true,false,false,1.0,true){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> decltype(rem(std::declval<T1>() , std::declval<T2>()))
    {
        return rem(a1,a2);
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

struct eval_max : eval_binfunc<eval_max>
{
    eval_max(Integer c) : eval_binfunc(c,false,true,false,true,1.0,true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(max(std::declval<T1>() , std::declval<T2>()))
    {
        return max(a1,a2);
    };
};
struct eval_min : eval_binfunc<eval_min>
{
    eval_min(Integer c) : eval_binfunc(c,false,true,false,true,1.0,true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> decltype(min(std::declval<T1>() , std::declval<T2>()))
    {
        return min(a1,a2);
    };
};

struct eval_fdim : eval_binfunc<eval_fdim>
{
    eval_fdim(Integer c) : eval_binfunc(c,false,true,false,false,1.0,true){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> decltype(fdim(std::declval<T1>() , std::declval<T2>()))
    {
        return fdim(a1,a2);
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

struct eval_nextafter : eval_binfunc<eval_nextafter>
{
    eval_nextafter(Integer c) : eval_binfunc(c,false,false,false,false,1.0,true){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> decltype(nextafter(std::declval<T1>() , std::declval<T2>()))
    {
        return nextafter(a1,a2);
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

struct eval_float_distance : eval_binfunc<eval_float_distance>
{
    eval_float_distance(Integer c) : eval_binfunc(c,false,false,false,false,1.0,true){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) 
        -> decltype(float_distance(std::declval<T1>() , std::declval<T2>()))
    {
        return float_distance(a1,a2);
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

struct eval_copysign : eval_binfunc<eval_copysign>
{
    eval_copysign(Integer c) : eval_binfunc(c,true,true,false,false,1.0,true){};

    template<class T1, class T2, class Enable = typename enable_not_complex<T1,T2>::type>
    static auto eval(const T1& a1, const T2& a2) -> decltype(copysign(std::declval<T1>() , std::declval<T2>()))
    {
        return copysign(a1,a2);
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
        return OReal(0);
    };
};

template<class T> struct promote_float_double                   { using type = T; };
template<>        struct promote_float_double<float>            { using type = double; };
template<>        struct promote_float_double<Float_complex>    { using type = Complex; };

template<class T1, class T2>
struct result_pow2
{
    using type_1    = typename promote_float_double<T1>::type;
    using type_2    = typename promote_float_double<T2>::type;
    using type      = decltype(pow(type_1(), type_2()));
};
template<class T1, class T2>
struct result_pow_c2
{
    using type_1    = typename promote_float_double<T1>::type;
    using type_2    = typename promote_float_double<T2>::type;
    using type      = decltype(pow_c(type_1(), type_2()));
};

struct eval_pow : eval_binfunc<eval_pow>
{
    eval_pow(Integer c) : eval_binfunc(c,true,true,false,true,100.0,true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> typename result_pow2<T1,T2>::type
    {
        using T1P   = typename result_pow2<T1,T2>::type_1;
        using T2P   = typename result_pow2<T1,T2>::type_2;
        return pow(T1P(a1), T2P(a2));
    };
};
struct eval_pow_c : eval_binfunc<eval_pow_c>
{
    eval_pow_c(Integer c) : eval_binfunc(c,true,true,false,true,100.0,true){};

    template<class T1, class T2>
    static auto eval(const T1& a1, const T2& a2) -> typename result_pow_c2<T1,T2>::type
    {
        using T1P   = typename result_pow_c2<T1,T2>::type_1;
        using T2P   = typename result_pow_c2<T1,T2>::type_2;
        return pow_c(T1P(a1), T2P(a2));
    };
};


double gmp_tester_bin::test_atan2(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_atan2 test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_hypot(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_hypot test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_mod(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_mod test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_rem(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_rem test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_min(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_min test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_max(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_max test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_fdim(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_fdim test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_nextafter(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_nextafter test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_float_distance(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_float_distance test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_copysign(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_copysign test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_pow(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_pow test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_pow_c(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_pow_c test(code);
    return test.make(s1, s2);
};

static double rand_fl()
{
    int c = abs(matcl::irand()) % 5;
    switch (c)
    {
        case 0: return 0.0;
        case 1: return 1.0;
        case 2: return -1.0;
        case 3: return 2.0;
        case 4: return -2.0;
        default:
            return 0.0;
    }
};

void gmp_tester_bin::test_pow_spec()
{
    bool res    = true;

    for (int i = 0; i < 1000; ++i)
    {
        Real a          = rand_fl();
        Real b          = rand_fl();

        Real c, d;

        if (abs(matcl::irand()) % 2 == 0)
        {
            c           = Real(matcl::irand() % 40) /2.0;
            d           = 0.0;
        }
        else
        {
            c           = rand_fl();
            d           = rand_fl();
        };

        Complex ex      = Complex(c,d);
        Complex res1    = pow(Complex(a,b), ex);
        mp_complex res2 = pow(mp_complex(a,b), ex);

        mp_float dif    = abs(res1 - res2);
        Real tol        = 100.0 * eps(res1);

        if (is_nan(res1) == true && is_nan(res2) == true)
            continue;
        if (res1 == res2)
            continue;

        if (dif > tol || is_finite(res1) == false || is_finite(res2) == false)
        {
            res         = false;
            out_stream << res1 << " " << res2 << " " << dif << "\n";
        }
    };

    if (res == true)
        out_stream << "pow special: OK" << "\n";
    else
        out_stream << "pow special: FAILED" << "\n";
};

}};
