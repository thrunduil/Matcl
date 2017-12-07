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
        double z        = 0.0;
        using simd_type = simd::simd<double, 128, simd::sse_tag>;

        for (int i = 0; i <m_N; ++i)
        {
            simd_type x = simd_type::broadcast(m_arr[i]);
            z       += simd::exp(x).first();        
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

int main(int argc, const char* argv[])
{
    (void)argc;
    (void)argv;

    using log_ptr   = std::shared_ptr<std::ofstream>;

    if (1)
    {
        int N       = 30000;
        int M       = 3000;

        using func_ptr  = benchmark::function_ptr;

        std::vector<double> x;
        x.reserve(N);

        for (int i = 0; i < N; ++i)
            x.push_back(matcl::randn());

        const double* ptr_x = x.data();

        benchmark b1    = benchmark(func_ptr(new Func_exp(N, ptr_x)));
        benchmark b3    = benchmark(func_ptr(new Func_exp2(N, ptr_x)));
	    //benchmark b4    = benchmark(func_ptr(new Func_exp3(N, ptr_x)));
        //benchmark b5    = benchmark(func_ptr(new Func_exp_twofold(N, ptr_x)));

        time_stats s3   = b3.make(M);
        time_stats s1   = b1.make(M);
        //time_stats s4   = b4.make(M);        		
        //time_stats s5   = b5.make(M);

        std::cout << "exp : " << s1 << "\n";
        std::cout << "exp2: " << s3 << "\n";
		//std::cout << "exp3: " << s4 << "\n";
        //std::cout << "twof: " << s5 << "\n";
        std::cout << "r e2: " << s1/s3 << "\n";
		//std::cout << "r e3: " << s1/s4 << "\n";
        //std::cout << "r tf: " << s1/s5 << "\n";

    };

    if (0)
    {
        double max_a        = double(709.78271289338399684324569237317);
        double min_a        = double(-709.08956571282405153382846025171);
        double max_a1       = std::nextafter(max_a, 1000.0);
        double max_a2       = std::nextafter(max_a, -1000.0);
        double min_a1       = std::nextafter(min_a, 1000.0);
        double min_a2       = std::nextafter(min_a, -1000.0);

        double e1           = simd::exp(std::numeric_limits<double>::quiet_NaN());
        double e2           = simd::exp(std::numeric_limits<double>::infinity());
        double e3           = simd::exp(-std::numeric_limits<double>::infinity());
        double e4           = simd::exp(710.0);
        double e5           = simd::exp(-710.0);                
        double e6           = simd::exp(708.0);
        double e7           = simd::exp(-708.0);

        double e10          = simd::exp(min_a2);
        double e11          = simd::exp(min_a);
        double e12          = simd::exp(min_a1);
        double e13          = simd::exp(max_a2);
        double e14          = simd::exp(max_a);
        double e15          = simd::exp(max_a1);

        std::cout << std::exp(min_a2) << "\n";
        std::cout << std::exp(min_a) << "\n";
        std::cout << std::exp(min_a1) << "\n";
        std::cout << std::exp(max_a2) << "\n";
        std::cout << std::exp(max_a) << "\n";
        std::cout << std::exp(max_a1) << "\n";

        std::cout << e1 << "\n";
        std::cout << e2 << "\n";
        std::cout << e3 << "\n";
        std::cout << e4 << "\n";
        std::cout << e5 << "\n";
        std::cout << e6 << "\n";
        std::cout << e7 << "\n";

        std::cout << e10 << "\n";
        std::cout << e11 << "\n";
        std::cout << e12 << "\n";
        std::cout << e13 << "\n";
        std::cout << e14 << "\n";
        std::cout << e15 << "\n";
    }

    if (1)
    {
        int N       = 1000*1000*100;

        #if (_DEBUG)
            N       = N / 10;
        #endif

        std::vector<double> x;
        x.reserve(N);

        double val_log2 = 0.693147180559945309417232121458176;
        double step     = val_log2 / 1024 / 1024;
        double start    = step * 89600000;
        step            = step / 256;
        (void)val_log2;
        (void)start;

        for (int i = 0; i < N; ++i)
        {
            //double mult     = 1.0/256.0/2.0;
            //double mult     = val_log2/2.0;
            x.push_back((1.0 - 2.0*matcl::rand()) * 700.0);
            //x.push_back((1.0 - 2.0*matcl::rand()) * mult);
            //x.push_back(-(start + i*step));
            //x.push_back(matcl::rand());
        };

        precision prec      = precision(53*4);
        precision prec_d    = precision::precision_double();

        double max_dif = 0;
        for (int i = 0; i < N; ++i)
        {
            //x[i]          = val_log2/2;
            double res1     = simd::exp(x[i]);
            mp_float res2   = exp(mp_float(x[i], prec));

            double dif      = ulp_distance(res2, mp_float(res1), prec_d).cast_float();

            if (dif > max_dif)
            {
                std::cout << dif - 1.0 << " " << i << " " << x[i] << " " << res2 << " " << res1 - res2 << "\n";
                simd::exp(x[i]);
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
