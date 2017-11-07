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
#include "helpers.h"

namespace matcl { namespace test
{

void test_scal_accuracy(std::ostream& os)
{
    scal_accuracy_tester test;
    test.make(os);
};

void scal_accuracy_tester::make(std::ostream& os)
{
  #ifdef _DEBUG
    Integer N       = 50;    
    Integer Nc      = 10;
  #else
    Integer N       = 50*50;    
    Integer Nc      = 10*20;
  #endif
  
    max_vec max_real  = {1.0, 1.0e-1, 1.0e1, 1.0e-2, 1.0e2, 1.0e-4, 1.0e4, 1.0e-8, 1.0e8,
                        1.0e-16, 1.0e16, 1.0e-32, 1.0e32, 1.0e-64, 1.0e64, 1.0e-128, 1.0e128,
                        1.0e-256, 1.0e256};
    max_vec max_compl = {1.0, 1.0e-1, 1.0e1, 1.0e-2, 1.0e2, 1.0e-4, 1.0e4, 1.0e-8, 1.0e8,
                        1.0e-16, 1.0e16, 1.0e-32, 1.0e32, 1.0e-128, 1.0e128, 0.0};

    test_unary<Real>(os, N, max_real);
    test_unary<Float>(os, N, max_real);    

    test_unary_compl<Float>(os, Nc, max_compl);    
    test_unary_compl<Real>(os, Nc, max_compl);    
};

template<class Type>
void scal_accuracy_tester::test_unary(std::ostream& os, Integer N, const max_vec& max_v)
{

    out_stream << "type: " << type_name<Type>::value() << "\n";

    formatted_disp dm;

    dm.set_row_label("function", align_type::right, 6);
    dm.add_column("matcl ulp", align_type::left, 7);
    dm.add_column("std ulp", align_type::left, 7);

    dm.disp_header();

    test_unary_func<Type, Func_exp>(dm, N, max_v);

    test_unary_func<Type, Func_sqrt>(dm, N, max_v);
    test_unary_func<Type, Func_cbrt>(dm, N, max_v);
    test_unary_func<Type, Func_sqrt1pm1>(dm, N, max_v);    
    test_unary_func<Type, Func_exp>(dm, N, max_v);
    test_unary_func<Type, Func_exp2>(dm, N, max_v);
    test_unary_func<Type, Func_exp10>(dm, N, max_v);
    test_unary_func<Type, Func_expm1>(dm, N, max_v);
    test_unary_func<Type, Func_log>(dm, N, max_v);
    test_unary_func<Type, Func_log2>(dm, N, max_v);
    test_unary_func<Type, Func_log10>(dm, N, max_v);
    test_unary_func<Type, Func_log1p>(dm, N, max_v);

    test_unary_func<Type, Func_sin>(dm, N, max_v);
    test_unary_func<Type, Func_cos>(dm, N, max_v);
    test_unary_func<Type, Func_tan>(dm, N, max_v);
    test_unary_func<Type, Func_cot>(dm, N, max_v);
    test_unary_func<Type, Func_sec>(dm, N, max_v);
    test_unary_func<Type, Func_csc>(dm, N, max_v);

    test_unary_func<Type, Func_sinh>(dm, N, max_v);
    test_unary_func<Type, Func_cosh>(dm, N, max_v);
    test_unary_func<Type, Func_tanh>(dm, N, max_v);
    test_unary_func<Type, Func_coth>(dm, N, max_v);
    test_unary_func<Type, Func_sech>(dm, N, max_v);
    test_unary_func<Type, Func_csch>(dm, N, max_v);

    test_unary_func<Type, Func_asin>(dm, N, max_v);
    test_unary_func<Type, Func_acos>(dm, N, max_v);
    test_unary_func<Type, Func_atan>(dm, N, max_v);
    test_unary_func<Type, Func_acot>(dm, N, max_v);
    test_unary_func<Type, Func_asec>(dm, N, max_v);
    test_unary_func<Type, Func_acsc>(dm, N, max_v);

    test_unary_func<Type, Func_asinh>(dm, N, max_v);
    test_unary_func<Type, Func_acosh>(dm, N, max_v);
    test_unary_func<Type, Func_atanh>(dm, N, max_v);
    test_unary_func<Type, Func_acoth>(dm, N, max_v);
    test_unary_func<Type, Func_asech>(dm, N, max_v);
    test_unary_func<Type, Func_acsch>(dm, N, max_v);

    test_unary_func<Type, Func_abs>(dm, N, max_v);
    test_unary_func<Type, Func_abs2>(dm, N, max_v);
    test_unary_func<Type, Func_arg>(dm, N, max_v);
    test_unary_func<Type, Func_inv>(dm, N, max_v);
    test_unary_func<Type, Func_sign>(dm, N, max_v);

    os << "\n";
};

template<class Type>
void scal_accuracy_tester::test_unary_compl(std::ostream& os, Integer N, const max_vec& max_v)
{

    out_stream << "type: " << complex_type_name<Type>::value() << "\n";

    formatted_disp dm;

    // it is not possible to have small error in real/imag part
    // therefore we test error of abs value
    dm.set_row_label("function", align_type::right, 6);
    dm.add_column("matcl ulp", align_type::left, 7);
    dm.add_column("std ulp", align_type::left, 7);

    dm.disp_header();

    //
    test_unary_compl_func<Type, Func_sqrt>(dm, N, max_v);
    test_unary_compl_func<Type, Func_sqrt1pm1>(dm, N, max_v);    
    test_unary_compl_func<Type, Func_exp>(dm, N, max_v);
    test_unary_compl_func<Type, Func_expi>(dm, N, max_v);
    test_unary_compl_func<Type, Func_exp2>(dm, N, max_v);
    test_unary_compl_func<Type, Func_exp10>(dm, N, max_v);
    test_unary_compl_func<Type, Func_expm1>(dm, N, max_v);
    test_unary_compl_func<Type, Func_log>(dm, N, max_v);
    test_unary_compl_func<Type, Func_log2>(dm, N, max_v);
    test_unary_compl_func<Type, Func_log10>(dm, N, max_v);
    test_unary_compl_func<Type, Func_log1p>(dm, N, max_v);

    test_unary_compl_func<Type, Func_sin>(dm, N, max_v);
    test_unary_compl_func<Type, Func_cos>(dm, N, max_v);
    test_unary_compl_func<Type, Func_tan>(dm, N, max_v);
    test_unary_compl_func<Type, Func_cot>(dm, N, max_v);
    test_unary_compl_func<Type, Func_sec>(dm, N, max_v);
    test_unary_compl_func<Type, Func_csc>(dm, N, max_v);

    test_unary_compl_func<Type, Func_sinh>(dm, N, max_v);
    test_unary_compl_func<Type, Func_cosh>(dm, N, max_v);
    test_unary_compl_func<Type, Func_tanh>(dm, N, max_v);
    test_unary_compl_func<Type, Func_coth>(dm, N, max_v);
    test_unary_compl_func<Type, Func_sech>(dm, N, max_v);
    test_unary_compl_func<Type, Func_csch>(dm, N, max_v);

    test_unary_compl_func<Type, Func_asin>(dm, N, max_v);
    test_unary_compl_func<Type, Func_acos>(dm, N, max_v);
    test_unary_compl_func<Type, Func_atan>(dm, N, max_v);
    test_unary_compl_func<Type, Func_acot>(dm, N, max_v);
    test_unary_compl_func<Type, Func_asec>(dm, N, max_v);
    test_unary_compl_func<Type, Func_acsc>(dm, N, max_v);

    test_unary_compl_func<Type, Func_asinh>(dm, N, max_v);
    test_unary_compl_func<Type, Func_acosh>(dm, N, max_v);
    test_unary_compl_func<Type, Func_atanh>(dm, N, max_v);
    test_unary_compl_func<Type, Func_acoth>(dm, N, max_v);
    test_unary_compl_func<Type, Func_asech>(dm, N, max_v);
    test_unary_compl_func<Type, Func_acsch>(dm, N, max_v);

    test_unary_compl_func<Type, Func_abs>(dm, N, max_v);
    test_unary_compl_func<Type, Func_abs2>(dm, N, max_v);
    test_unary_compl_func<Type, Func_arg>(dm, N, max_v);
    test_unary_compl_func<Type, Func_inv>(dm, N, max_v);
    test_unary_compl_func<Type, Func_sign>(dm, N, max_v);

    os << "\n";
};

template<>
Float scal_accuracy_tester::rand_scalar<Float>(double max)
{
    return matcl::frandn() * Float(max);
};

template<>
Real scal_accuracy_tester::rand_scalar<Real>(double max)
{
    return matcl::randn() * max;
};

template<class Type>
Real scal_accuracy_tester::calc_ulp(Type res, const mp_float& res_ext)
{
    mp_float dif    = abs(res - res_ext, prec_type<Type>() + 4);
    mp_float tol    = eps(mp_float(res_ext, prec_type<Type>()));
    Real ulp        = (dif / tol).cast_float();

    return ulp;
};

template<class Type>
Real scal_accuracy_tester::calc_ulp_complex(const Type& res, const mp_complex& res_ext)
{
    using Real_ty       = typename details::real_type<Type>::type;

    precision p         = prec_type<Type>();
    mp_float abs_res    = abs(mp_complex(res, p + 4));
    mp_float abs_res_ex = abs(res_ext, p + 4);

    mp_float dif        = abs(abs_res - abs_res_ex, p + 4);
    mp_float tol        = eps(mp_float(abs_res_ex, p));
    Real ulp            = (dif / tol).cast_float();

    return ulp;
};

template<class Type, class Derived>
struct eval_scalar_func_unary
{
    Integer code;

    eval_scalar_func_unary(Integer c) 
        : code(c)
    {};

    template<class T1>
    auto eval_fun_matcl(const T1& a1) 
        -> decltype(Derived::eval_matcl(std::declval<T1>()))
    {
        return Derived::eval_matcl(a1);
    };

    template<class T1>
    auto eval_fun_std(const T1& a1) 
        -> decltype(Derived::eval_std(std::declval<T1>()))
    {
        return Derived::eval_std(a1);
    };

    template<class T1>
    auto eval_fun_mp(const T1& a1, precision p) 
        -> decltype(Derived::eval_mp(std::declval<T1>(), p))
    {
        return Derived::eval_mp(a1, p);
    };

    void make(Type s, Real& ulp_matcl, Real& ulp_std)
    {
        // double precision
        int dprec           = 10;
        precision prec      = scal_accuracy_tester::prec_type<Type>() + dprec;

        Type res_matcl      = eval_fun_matcl(s);
        Type res_std        = eval_fun_std(s);
        mp_float res_ext    = eval_fun_mp(mp_float(s, prec), prec);

        // ignore underflows
        if (abs(res_std) < min_normal<Type>() && abs(res_matcl) < min_normal<Type>())
        {
            ulp_matcl = 0.0;
            ulp_std = 0.0;
            return;
        };

        ulp_matcl           = scal_accuracy_tester::calc_ulp(res_matcl, res_ext);
        ulp_std             = scal_accuracy_tester::calc_ulp(res_std, res_ext);
    };    
};

template<class Type, class Derived>
struct eval_scalar_func_unary_complex
{
    Integer code;

    eval_scalar_func_unary_complex(Integer c) 
        : code(c)
    {};

    template<class T1>
    auto eval_fun_matcl(const T1& a1) 
        -> decltype(Derived::eval_matcl(std::declval<T1>()))
    {
        return Derived::eval_matcl(a1);
    };

    template<class T1, class Ty_std = typename std_complex_type<T1>::type>
    auto eval_fun_std(const T1& a1) 
        -> decltype(Derived::eval_std(std::declval<Ty_std>()))
    {
        return Derived::eval_std(make_std_complex(a1));
    };

    template<class T1>
    auto eval_fun_mp(const T1& a1, precision p) 
        -> decltype(Derived::eval_mp(std::declval<T1>(), p))
    {
        return Derived::eval_mp(a1, p);
    };

    void make(Type s, Real& ulp_matcl, Real& ulp_std)
    {
        // double precision
        int dprec           = 10;
        precision prec      = scal_accuracy_tester::prec_type<Type>() + dprec;

        auto res_matcl      = eval_fun_matcl(s);
        auto res_std        = eval_fun_std(s);
        mp_complex res_ext  = eval_fun_mp(mp_complex(s, prec), prec);
        Real tmp            = abs(res_ext).cast_float();

        // ignore underflows
        if (abs(res_std) < min_normal<Type>() && abs(res_matcl) < min_normal<Type>())
        {
            ulp_matcl       = 0.0;
            ulp_std         = 0.0;
            return;
        };

        ulp_matcl   = scal_accuracy_tester::calc_ulp_complex(res_matcl, res_ext);
        ulp_std     = scal_accuracy_tester::calc_ulp_complex(res_std, res_ext);
    };    
};

template<class Type, class Func>
void scal_accuracy_tester::test_scalar(Type s, Integer code, Real& ulp_matcl, Real& ulp_std)
{
    eval_scalar_func_unary<Type, Func> test1(code);
    test1.make(s, ulp_matcl, ulp_std);
};

template<class Type, class Func>
void scal_accuracy_tester::test_scalar_complex(const Type& s, Integer code, 
                 Real& ulp_matcl, Real& ulp_std)
{
    eval_scalar_func_unary_complex<Type, Func> test1(code);
    test1.make(s, ulp_matcl, ulp_std);
};

template<class Type, class Func>
void scal_accuracy_tester::test_unary_func(formatted_disp& os, Integer N, const max_vec& max_v)
{
    double matcl_ulp_max    = 0.0;
    double std_ulp_max      = 0.0;

    for (const auto& max : max_v)
        test_unary_func<Type, Func>(os, N, max, matcl_ulp_max, std_ulp_max);

    os.disp_row(Func::func_name(), matcl_ulp_max, std_ulp_max);
};

template<class Type, class Func>
void scal_accuracy_tester::test_unary_compl_func(formatted_disp& os, Integer N, const max_vec& max_v)
{
    double matcl_ulp_max    = 0.0;
    double std_ulp_max      = 0.0;

    for (const auto& max_re : max_v)
    for (const auto& max_im : max_v)
    {
        test_unary_compl_func<Type, Func>(os, N, max_re, max_im, matcl_ulp_max, std_ulp_max);
    }

    os.disp_row(Func::func_name(), matcl_ulp_max, std_ulp_max);
};

template<class Type, class Func>
void scal_accuracy_tester::test_unary_func(formatted_disp& os, Integer N, Real max, Real& matcl_ulp,
            Real& std_ulp)
{
    std::vector<Type> scalars;
    
    for (Integer i = 0; i < 10*N; ++i)
    {
        Type val = rand_scalar<Type>(max);

        if (value_allowed<Type, Func>::test(val) == false)
            continue;

        scalars.push_back(val);

        if ((Integer)scalars.size() > N)
            break;
    }

    N   = (Integer)scalars.size();

    for (Integer i = 0; i < N; ++i)
    {
        double ulp_matcl;
        double ulp_std;

        test_scalar<Type, Func>(scalars[i], i, ulp_matcl, ulp_std);

        if (matcl::is_finite(ulp_matcl) == false)
            continue;

        matcl_ulp   = std::max(matcl_ulp, ulp_matcl);
        std_ulp     = std::max(std_ulp, ulp_std);

        if (matcl_ulp > 1.0)
        {
            test_scalar<Type, Func>(scalars[i], i, ulp_matcl, ulp_std);
        };
    };
};

template<class Type, class Func>
void scal_accuracy_tester::test_unary_compl_func(formatted_disp& os, Integer N, Real max_re, 
                Real max_im, Real& matcl_ulp, Real& std_ulp)
{
    using compl_type    = typename details::complex_type<Type>::type;

    std::vector<compl_type> scalars;
    
    for (Integer i = 0; i < 10*N; ++i)
    {
        Type re = rand_scalar<Type>(max_re);
        Type im = rand_scalar<Type>(max_im);

        if (value_allowed_real<Type, Func>::test(re) == false)
            continue;

        if (value_allowed_imag<Type, Func>::test(im) == false)
            continue;

        scalars.push_back(compl_type(re, im));

        if ((Integer)scalars.size() > N)
            break;
    }

    N   = (Integer)scalars.size();

    for (Integer i = 0; i < N; ++i)
    {
        double ulp_matcl;
        double ulp_std;

        test_scalar_complex<compl_type, Func>(scalars[i], i, ulp_matcl, ulp_std);

        if (matcl::is_finite(ulp_matcl) == false)
            continue;

        matcl_ulp       = std::max(matcl_ulp, ulp_matcl);
        std_ulp         = std::max(std_ulp, ulp_std);

        if (ulp_matcl > 4)
        {
            test_scalar_complex<compl_type, Func>(scalars[i], i, ulp_matcl, ulp_std);
        };
    };
};

}};
