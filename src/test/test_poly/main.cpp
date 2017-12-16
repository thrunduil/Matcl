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
#include "matcl-simd/poly/poly_eval_twofold.h"
#include "matcl-scalar/lib_functions/utils.h"
#include "matcl-scalar/lib_functions/scalar_gen.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/lib_functions/func_binary.h"

#include "matcl-core/profile/benchmark.h"
#include "matcl-simd/simd_math.h"


#include <iostream>
#include <fstream>
#include <iomanip>

using namespace matcl;

constexpr int64_t fact(int64_t v)
{
    return (v == 1) ? v : v * fact(v-1);
};

force_inline
float log_impl(float x0)
{
    using simd_type     = simd::simd<float, 256, simd::avx_tag>;
    simd_type x         = simd_type(x0);
    simd_type res       = simd::log(x);

    return res.first();
}

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

template<class T>
struct Func_exp : public benchmark_function
{
    int             m_N;
    const T*        m_arr;

    Func_exp(int n, const T* arr) 
        : m_N(n), m_arr(arr)
    {};

    void eval() override
    {
        T z         = T();

        for (int i = 0; i <m_N; ++i)
        {
            z       += std::exp(m_arr[i]);        
        };

        benchmark::use_value(z);
    };
};

template<class T>
struct Func_log : public benchmark_function
{
    int             m_N;
    const T*        m_arr;

    Func_log(int n, const T* arr) 
        : m_N(n), m_arr(arr)
    {};

    void eval() override
    {
        T z         = T();

        for (int i = 0; i <m_N; ++i)
        {
            z       += std::log(m_arr[i]);        
        };

        benchmark::use_value(z);
    };
};

template<class T>
struct Func_exp2 : public benchmark_function
{
    int             m_N;
    const T*        m_arr;

    Func_exp2(int n, const T* arr) 
        : m_N(n), m_arr(arr)
    {};

    void eval() override
    {        
        using simd_type = simd::simd<T, 256, simd::avx_tag>;
        static const int vector_size = simd_type::vector_size;

        simd_type z     = simd_type::zero();
        //T z             = T();

        for (int i = 0; i < m_N; i += vector_size)
        {
            simd_type x = simd_type::load(m_arr + i);
            z           += simd::exp(x);
            //z           += simd::exp(m_arr[i]);
        };

        benchmark::use_value(z);
    };
};

template<class T>
struct Func_log2 : public benchmark_function
{
    int             m_N;
    const T*        m_arr;

    Func_log2(int n, const T* arr) 
        : m_N(n), m_arr(arr)
    {};

    void eval() override
    {        
        using simd_type = simd::simd<T, 256, simd::avx_tag>;
        static const int vector_size = simd_type::vector_size;

        simd_type z     = simd_type::zero();
        //T z             = T();

        for (int i = 0; i < m_N; i += vector_size)
        {
            simd_type x = simd_type::load(m_arr + i);
            z           += simd::log(x);
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
    
    //matcl::simd::details::exp_table_double::generate_lookup_table(std::cout, 4);

    if (0)
    {
        int N       = 50000;
        int M       = 3000;

        #if (_DEBUG)
            N       = N / 10;
            M       = M / 10;
        #endif

        using func_ptr  = benchmark::function_ptr;

        std::vector<float> x;
        x.reserve(N);

        for (int i = 0; i < N; ++i)
            x.push_back(abs(matcl::frandn()));

        const float* ptr_x = x.data();

        benchmark b1    = benchmark(func_ptr(new Func_log<float>(N, ptr_x)));
        benchmark b3    = benchmark(func_ptr(new Func_log2<float>(N, ptr_x)));

        //benchmark b1    = benchmark(func_ptr(new Func_exp<double>(N, ptr_x)));
        //benchmark b3    = benchmark(func_ptr(new Func_exp2<double>(N, ptr_x)));

	    //benchmark b4    = benchmark(func_ptr(new Func_exp3(N, ptr_x)));
        //benchmark b5    = benchmark(func_ptr(new Func_exp_twofold(N, ptr_x)));

        time_stats s3   = b3.make(M);
        time_stats s1   = b1.make(M);
        //time_stats s4   = b4.make(M);        		
        //time_stats s5   = b5.make(M);

        std::cout << "log std: " << s1 << "\n";
        std::cout << "log imp: " << s3 << "\n";
		//std::cout << "exp3: " << s4 << "\n";
        //std::cout << "twof: " << s5 << "\n";
        std::cout << "r e2: " << s1/s3 << "\n";
		//std::cout << "r e3: " << s1/s4 << "\n";
        //std::cout << "r tf: " << s1/s5 << "\n";
        std::cout << "\n";
    };

    if (1)
    {
        int N       = 1000*1000*100;

        #if (_DEBUG)
            N       = N / 10;
        #endif

        std::vector<float> x;
        x.reserve(N);

        double val_log2 = 0.693147180559945309417232121458176;
        double step     = val_log2 / 1024 / 1024;
        double start    = 0.0;
        step            = step / 256;
        (void)val_log2;
        (void)start;

        float scal_den  = std::numeric_limits<float>::denorm_min() * 1.0e4f;
        for (int i = 0; i < N; ++i)
        {
            //double mult     = 1.0/256.0/2.0;
            //double mult     = val_log2/2.0;
            //x.push_back((1.0 - 2.0*matcl::rand()) * 5.0);
            //x.push_back((1.0f - 2.0f*matcl::frand()) * 80.0f);
            //x.push_back((1.0 - 2.0*matcl::rand()) * mult);
            //x.push_back(-(start + i*step));
            x.push_back(std::abs(matcl::frand()) * scal_den);

            //x.push_back(matcl::frand() * 0.66666666f + 0.6666666666f);
        };

        precision prec      = precision(53*2);
        precision prec_d    = precision::precision_float();

        x[0]                = 9.601e-39f;

        double max_dif = 0;
        for (int i = 0; i < N; ++i)
        {
            //x[i]          = val_log2/2;
            float res1      = log_impl(x[i]);
            //double res1   = std::log(x[i]);

            //mp_float y    = mp_float(x[i], prec) - mp_float(1.0);
            //mp_float y2   = y * y;
            mp_float xm1    = mp_float(x[i], prec) - 1.0;

            //mp_float res2   = (log(mp_float(x[i], prec))/xm1 - 1.0)/ xm1;
            mp_float res2   = log(mp_float(x[i], prec));

            double dif      = ulp_distance(res2, mp_float(res1), prec_d).cast_float();

            if (dif > max_dif)
            {
                std::cout << dif - 1.0 << " " << i << " " << x[i] << " " << res2.to_string(precision(20)) << " " << res1 - res2 << "\n";
                log_impl(x[i]);
                max_dif = dif;
            };
        };
    };

    try
    {         
        {
            std::string log_file_name   = std::string("log_test_poly.txt");
            log_ptr log = log_ptr(new std::ofstream(log_file_name));
            set_logger(log);
        };
        
        matcl::test::test_poly();        
        
        matcl::test::test_poly_dyn(true);
        matcl::test::test_poly_dyn(false);        

        matcl::test::test_poly_cond();        

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
