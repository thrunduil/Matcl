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

#include "test_poly.h"
#include "matcl-core/IO/logger.h"
#include "matcl-core/float/twofold.h"
#include "matcl-simd/poly/poly_eval.h"
#include "matcl-scalar/lib_functions/utils.h"
#include "matcl-scalar/lib_functions/scalar_gen.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/lib_functions/func_binary.h"

#include "matcl-core/profile/benchmark.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace matcl;

force_inline
double fnms(double x, double y, double z)
{
    return z - x * y;
};

double eval_horn(double x)
{
    double poly[] = {1.0, 1.0, 1.0/2.0, 1.0/6.0, 1.0/24.};//, 1.0/120., 1.0/720., 1.0/5040.};
                    //1./40320., 1.0/362880., 1.0/3628800., 1.0/39916800.};

    static const int size = sizeof(poly)/sizeof(double);
    //using simd_type = matcl::simd::simd<double, 128, simd::sse_tag>;
    using simd_type = double;

    simd_type xs    = simd_type(x);
    simd_type z     = simd::horner<size>(xs, poly);
    //double res    = z.first();
    double res      = z;

    return res;
};

double eval_estrin(double x)
{
    double poly[] = {1.0, 1.0, 1.0/2.0, 1.0/6.0, 1.0/24.};//, 1.0/120., 1.0/720., 1.0/5040.};
                    //1./40320., 1.0/362880., 1.0/3628800., 1.0/39916800.};

    static const int size = sizeof(poly)/sizeof(double);
    //using simd_type = matcl::simd::simd<double, 256, simd::avx_tag>;
    using simd_type = double;

    simd_type xs    = simd_type(x);
    simd_type z     = simd::estrin<size>(xs, poly);
    //double res    = z.first();
    double res      = z;

    return res;
};

struct exp_impl
{
    force_inline
    static double eval(double a0)
    {
        double hi, lo, x;
        double k = reduce(a0, hi, lo, x);
        double c = approx(x);
        c        = finalize(x, c, hi, lo);

        (void)k;
        return c;
    };

    force_inline
    static double reduce(double a0, double& hi, double& lo, double& x)
    {
        const double inv_log2 = 1.442695040888963387;
        const double log2_hi  = 0.69314718036912381649;
        //const double log2_lo  = 1.9082149292705877e-10;

        //double k    = simd::scalar_func::round_impl<double>::eval(inv_log2 * a0);
        //hi          = fnms(k, log2_hi, a0); //a0-k*L
        lo          = 0.0;
        hi          = 0.0;
        //lo          = k * log2_lo;
        //x           = hi-lo;
        x           = a0;
        //return k;
        return 0.0;
    };

    force_inline
    static double approx(double x)
    {
        // horner approximation of (exp(x) - 1 - x) / x^2
        double p        = simd::small_horner(x,  
                            5.0000000000000010e-1,
                            1.6666666666666675e-1,
                            4.1666666666624158e-2,
                            8.3333333333222153e-3,
                            1.3888888917198445e-3,
                            1.9841269886566057e-4,
                            2.4801521320070370e-5,
                            2.7557242365240200e-6,
                            2.7620076836694271e-7,
                            2.5110039032532694e-8
                                                );
        return p;
    };

    force_inline
    static double finalize(double x, double c, double hi, double lo)
    {
        (void)lo;
        (void)hi;
        double c2   = fma_a(c, x, 1.0);
        double c3   = fma_a(c2, x, 1.0);
        return c3;
    }
};

/*
struct exp_impl
{
    force_inline
    static double eval(double a0)
    {
        double hi, lo, x;
        double k = reduce(a0, hi, lo, x);
        double c = approx(x);
        c        = finalize(x, c, hi, lo);

        (void)k;
        return c;
    };

    force_inline
    static double reduce(double a0, double& hi, double& lo, double& x)
    {
        const double inv_log2 = 1.442695040888963387;
        const double log2_hi  = 0.69314718036912381649;
        const double log2_lo  = 1.9082149292705877e-10;

        double k    = simd::scalar_func::round_impl<double>::eval(inv_log2 * a0);
        hi          = fnms(k, log2_hi, a0); //a0-k*L
        lo          = k * log2_lo;
        x           = hi-lo;
        return k;
    };

    force_inline
    static double approx(double x)
    {
        double const t  = x * x;
        double p        = simd::small_horner(t, 0.16666666666666601904, 
                                                -0.0027777777777015593384,
                                                6.6137563214379343612e-05,
                                                -1.6533902205465251539e-06, 
                                                4.1381367970572384604e-08
                                                );
        return fnms(t, p , x); //x-p*t
    };

    force_inline
    static double finalize(double x, double c, double hi, double lo)
    {
        return 1.0 - (((lo - (x*c)/(2.0 - c)) - hi));
    }
};
*/
constexpr int64_t fact(int64_t v)
{
    return (v == 1) ? v : v * fact(v-1);
};

struct exp_twofold_impl
{
    using twofold   = twofold<double>;

    //static const int size       = 19;
    static const int size       = 20;
    static twofold poly[size];
    static double dum_init;

    force_inline
    static double eval(double x)
    {
        return simd::compensated_horner(x, size, poly).value;
        //return simd::compensated_horner(x, size, poly);
    };

    static void init()
    {
        poly[0]     = twofold(1.0, 0.0);
        poly[1]     = twofold(1.0, 0.0);

        precision prec  = precision(53*4);

        for (int i = 2; i < size; ++i)
        {
            mp_float f  = factorial<mp_float>(i, prec);
            mp_float fi = mp_float(1.0)/f;

            double val  = fi.cast_float();
            
            mp_float er = fi - mp_float(val);
            double err  = er.cast_float();

            poly[i]     = twofold(val, err);
        }
    };

    static double initialize()
    {
        init();
        return poly[0].value;
    }
};

twofold<double> exp_twofold_impl::poly[size];
double exp_twofold_impl::dum_init = initialize();

struct Func_exp : public benchmark_function
{
    int             m_N;
    const double*   m_arr;

    Func_exp(int n, const double* arr) 
        : m_N(n), m_arr(arr)
    {};

    void eval() override
    {
        double z    = 0.0;

        for (int i = 0; i <m_N; ++i)
        {
            z       += std::exp(m_arr[i]);        
        };

        benchmark::use_value(z);
    };
};

struct Func_exp2 : public benchmark_function
{
    int             m_N;
    const double*   m_arr;

    Func_exp2(int n, const double* arr) 
        : m_N(n), m_arr(arr)
    {};

    void eval() override
    {
        double z    = 0.0;

        for (int i = 0; i <m_N; ++i)
        {
            z       += exp_impl::eval(m_arr[i]);        
        };

        benchmark::use_value(z);
    };
};

struct Func_exp_twofold : public benchmark_function
{
    int             m_N;
    const double*   m_arr;

    Func_exp_twofold(int n, const double* arr) 
        : m_N(n), m_arr(arr)
    {};

    void eval() override
    {
        double z    = 0.0;

        for (int i = 0; i <m_N; ++i)
        {
            z       += exp_twofold_impl::eval(m_arr[i]);        
        };

        benchmark::use_value(z);
    };
};

struct Func_poly : public benchmark_function
{
    int             m_N;
    const double*   m_arr;

    Func_poly(int n, const double* arr) 
        : m_N(n), m_arr(arr)
    {};

    void eval() override
    {
        double z    = 0.0;

        for (int i = 0; i < m_N; ++i)
        {
            z        += eval_horn(m_arr[i]);
        };

        benchmark::use_value(z);
    };
};

struct Func_estrin : public benchmark_function
{
    int             m_N;
    const double*   m_arr;

    Func_estrin(int n, const double* arr) 
        : m_N(n), m_arr(arr)
    {};

    void eval() override
    {
        double z    = 0.0;

        for (int i = 0; i <m_N; ++i)
        {
            z        += eval_estrin(m_arr[i]);
        };

        benchmark::use_value(z);
    };
};

int main(int argc, const char* argv[])
{
    (void)argc;
    (void)argv;

    using log_ptr   = std::shared_ptr<std::ofstream>;

    //test::test_iaca(2.0);

    if (1)
    {
        int N       = 10000;
        int M       = 1000;

        using func_ptr  = benchmark::function_ptr;

        std::vector<double> x;
        x.reserve(N);

        for (int i = 0; i < N; ++i)
            x.push_back(matcl::randn());

        const double* ptr_x = x.data();

        benchmark b1    = benchmark(func_ptr(new Func_exp(N, ptr_x)));
        benchmark b2    = benchmark(func_ptr(new Func_poly(N, ptr_x)));
        benchmark b3    = benchmark(func_ptr(new Func_estrin(N, ptr_x)));
        benchmark b4    = benchmark(func_ptr(new Func_exp2(N, ptr_x)));
        benchmark b5    = benchmark(func_ptr(new Func_exp_twofold(N, ptr_x)));

        time_stats s1   = b1.make(M);
        time_stats s2   = b2.make(M);
        time_stats s3   = b3.make(M);
        time_stats s4   = b4.make(M);
        time_stats s5   = b5.make(M);

        std::cout << "exp : " << s1 << "\n";
        std::cout << "poly: " << s2 << "\n";
        std::cout << "estr: " << s3 << "\n";
        std::cout << "exp2: " << s4 << "\n";
        std::cout << "twof: " << s5 << "\n";
        std::cout << "r p : " << s1/s2 << "\n";
        std::cout << "r e : " << s1/s3 << "\n";
        std::cout << "r e2: " << s1/s4 << "\n";
        std::cout << "r tf: " << s1/s5 << "\n";

    };

    {
        int N       = 1000000;

        std::vector<double> x;
        x.reserve(N);

        double val_log2 = 0.693147180559945309417232121458176;
        (void)val_log2;

        for (int i = 0; i < N; ++i)
        {
            x.push_back((1.0 - 2.0*matcl::rand()) * val_log2/2);
            //x.push_back((1023.0 + matcl::rand())/1024 * val_log2);
            //x.push_back(matcl::rand());
        };

        precision prec      = precision(53*4);
        precision prec_d    = precision::precision_double();

        double max_dif = 0;
        for (int i = 0; i < N; ++i)
        {
            //x[i]          = val_log2/2;
            double res1     = exp_impl::eval(x[i]);
            mp_float res2   = exp(mp_float(x[i], prec));

            double dif      = ulp_distance(res2, mp_float(res1), prec_d).cast_float();

            if (dif > max_dif)
            {
                std::cout << dif << " " << x[i] << " " << res2 << " " << res1 - res2 << "\n";
                exp_twofold_impl::eval(x[i]);
                max_dif = dif;
            };
        };
    };

    try
    {         
        {
            std::string log_file_name   = std::string("log_test_twofold.txt");
            log_ptr log = log_ptr(new std::ofstream(log_file_name));
            set_logger(log);
        };
        
        matcl::test::test_poly_cond();
        matcl::test::test_poly();
        
        matcl::test::test_poly_dyn(true);
        matcl::test::test_poly_dyn(false);        

        std::cout << "\n";
        std::cout << "finished" << "\n";
    }
    catch(std::exception& ex)
    {
        std::cout << ex.what() << "\n";
        return 1;
    };

    return 0;
}
