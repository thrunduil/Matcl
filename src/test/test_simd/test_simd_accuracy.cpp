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

#include "test_simd_accuracy.h"
#include "matcl-core/IO/output_stream.h"
#include "test_functions_accuracy.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/lib_functions/func_binary.h"

#include <algorithm>

namespace matcl { namespace test
{

void test_math_accuracy()
{
    simd_accuracy_tester test;
    test.make();
};

static bool g_rang_denormals    = false;

void set_rand_denormals(bool val)
{
    g_rang_denormals = val;
}

bool simd_accuracy_tester::rand_denormals()
{
    return g_rang_denormals;
};

void simd_accuracy_tester::make()
{
  #ifdef _DEBUG
    int N               = 100000;    
  #else
    int N               = 1000*1000;    
  #endif
  
    test_unary<double>(N);
    test_unary<float>(N);
};

template<class Type>
void simd_accuracy_tester::test_unary(int N)
{
    out_stream << "type: " << typeid(Type).name() << " accuracy test" << "\n";

    formatted_disp dm;

    dm.set_row_label("function", align_type::right, 10);
    dm.add_column("func ulp", align_type::left, 7);
    dm.add_column("ref ulp",  align_type::left, 7);

    dm.disp_header();

    namespace func = test_functions_accuracy;

    test_unary_func<Type>(dm, func::Func_exp<Type>(),N);
    test_unary_func<Type>(dm, func::Func_log<Type>(),N);
    test_unary_func<Type>(dm, func::Func_sin<Type>(),N);
    test_unary_func<Type>(dm, func::Func_cos<Type>(),N);
    test_unary_func<Type>(dm, func::Func_tan<Type>(),N);
    test_unary_func<Type>(dm, func::Func_cot<Type>(),N);

    out_stream << "\n";
};

template<class Type>
void simd_accuracy_tester::test_unary_func(formatted_disp& os, const Func<Type>& f, int N)
{
    double ulp_max_base     = 0.0;
    double ulp_max_ref      = 0.0;

    test_unary_func<Type>(os, f, N, ulp_max_base, ulp_max_ref);
    os.disp_row(f.name(), ulp_max_base, ulp_max_ref);
};

template<class Type>
void simd_accuracy_tester::test_unary_func(formatted_disp& os, const Func<Type>& f, int N, 
                                            double& ulp_base, double& ulp_ref)
{
    (void)os;

    std::vector<Type> scalars;
    f.rand_values(scalars, N);
    
    N   = (int)scalars.size();

    for (Integer i = 0; i < N; ++i)
    {
        double ulp_1;
        double ulp_2;

        test_scalar<Type>(f, scalars[i], i, ulp_1, ulp_2);

        ulp_base    = std::max(ulp_base, ulp_1);
        ulp_ref     = std::max(ulp_ref, ulp_2);

        /*
        if (ulp_1 > 1.5)
        {
            test_scalar<Type>(f, scalars[i], i, ulp_1, ulp_2);
        };
        */
    };
};

template<class Val>
struct get_precision_type{};

template<>
struct get_precision_type<double>
{
    static precision eval()
    {
        return precision::precision_double();
    }
};

template<>
struct get_precision_type<float>
{
    static precision eval()
    {
        return precision::precision_float();
    }
};

template<class Type>
struct eval_scalar_func_unary
{
    template<class Val>
    using Func      = matcl::test_accuracy_function<Val>;

    int code;

    eval_scalar_func_unary(int c) 
        : code(c)
    {};

    template<class T1>
    T1 eval_fun_base(const Func<Type>& f, const T1& a1) 
    {
        return f.eval_base(a1);
    };

    template<class T1>
    T1 eval_fun_ref(const Func<Type>& f, const T1& a1) 
    {
        return f.eval_ref(a1);
    };

    mp_float eval_fun_mp(const Func<Type>& f, const mp_float& a1) 
    {
        return f.eval_mp(a1);
    };

    void make(const Func<Type>& f, Type s, double& ulp_base, double& ulp_ref)
    {
        // double precision
        int dprec           = 20;
        precision prec      = get_precision_type<Type>::eval() + dprec;

        Type res_base       = eval_fun_base(f, s);
        Type res_ref        = eval_fun_ref(f, s);
        mp_float res_mp     = eval_fun_mp(f, mp_float(s, prec));

        ulp_base            = simd_accuracy_tester::calc_ulp(res_base, res_mp);
        ulp_ref             = simd_accuracy_tester::calc_ulp(res_ref, res_mp);
    };    
};

template<class Type>
void simd_accuracy_tester::test_scalar(const Func<Type>& f, Type s, int code, double& ulp_base, 
                                      double& ulp_ref)
{
    eval_scalar_func_unary<Type> test1(code);
    test1.make(f, s, ulp_base, ulp_ref);
};

template<class Type>
double simd_accuracy_tester::calc_ulp(Type res, const mp_float& res_mp)
{
    precision prec  = get_precision_type<Type>::eval();

    if (res == res_mp)
        return 0.0;

    if (matcl::is_nan(res) == true && matcl::is_nan(res_mp) == true)
        return 0.0;

    int prec_lost   = precision_lost_subnormal(res);
    double dist     = ulp_distance(res_mp, res, prec - prec_lost).cast_float();

    if (matcl::is_nan(dist) == true)
        dist        = constants::inf();

    return dist;
};

template<class Type>
int simd_accuracy_tester::precision_lost_subnormal(Type val)
{
    Integer exp;
    matcl::frexp(val, exp);

    Integer min_exp     = std::numeric_limits<Type>::min_exponent;

    if (exp >= min_exp)
        return 0;
    
    Integer prec_lost   = min_exp - exp;
    return prec_lost;
}

}};
