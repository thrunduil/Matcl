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

#include "matcl-mp/matcl_mp.h"
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

#include "matcl-matrep/matcl_matrep.h"
#include <iostream>
#include <fstream>
#include <iomanip>

//#include "C:/C++/iaca_2-win64/iacaMarks.h"

using namespace matcl;
namespace ms = matcl :: simd;

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
struct Func_sin : public benchmark_function
{
    int             m_N;
    int             m_K;
    const T*        m_arr;

    Func_sin(int n, int k, const T* arr) 
        : m_N(n), m_arr(arr), m_K(k)
    {};

    void eval() override
    {
        T z         = T();

        for (int j = 0; j < m_K; ++j)
        for (int i = 0; i <m_N; ++i)
        {
            z       += std::sin(m_arr[i]);        
        };

        benchmark::use_value(z);
    };
};

template<class T>
struct Func_sin2 : public benchmark_function
{
    int             m_N;
    int             m_K;
    const T*        m_arr;

    Func_sin2(int n, int k, const T* arr) 
        : m_N(n), m_arr(arr), m_K(k)
    {};

    void eval() override
    {        
        //using simd_type = simd::simd<T, 256, simd::avx_tag>;
        //static const int vector_size = simd_type::vector_size;
        static const int vector_size = 1;

        //simd_type z     = simd_type::zero();
        T z             = T();

        for (int j = 0; j < m_K; ++j)
        for (int i = 0; i < m_N; i += vector_size)
        {
            //simd_type x = simd_type::load(m_arr + i);
            //z           += ms::sin(x);
            z           += ms::sin(m_arr[i]);
        };

        benchmark::use_value(z);
    };
};

template<class T>
struct Func_tan : public benchmark_function
{
    int             m_N;
    int             m_K;
    const T*        m_arr;

    Func_tan(int n, int k, const T* arr) 
        : m_N(n), m_arr(arr), m_K(k)
    {};

    void eval() override
    {
        T z         = T();

        for (int j = 0; j < m_K; ++j)
        for (int i = 0; i <m_N; ++i)
        {
            z       += std::tan(m_arr[i]);        
        };

        benchmark::use_value(z);
    };
};

template<class T>
struct Func_tan2 : public benchmark_function
{
    int             m_N;
    int             m_K;
    const T*        m_arr;

    Func_tan2(int n, int k, const T* arr) 
        : m_N(n), m_arr(arr), m_K(k)
    {};

    void eval() override
    {        
        //using simd_type = simd::simd<T, 256, simd::avx_tag>;
        //static const int vector_size = simd_type::vector_size;
        static const int vector_size = 1;

        //simd_type z     = simd_type::zero();
        T z             = T();

        for (int j = 0; j < m_K; ++j)
        for (int i = 0; i < m_N; i += vector_size)
        {
            //simd_type x = simd_type::load(m_arr + i);
            //z           += ms::sin(x);
            z           += ms::tan(m_arr[i]);
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

double next(double x, int i)
{
    if (i == 0)
        return x;

    double dir  = (i > 0)   ? std::numeric_limits<double>::infinity()
                            : -std::numeric_limits<double>::infinity();

for (int k = 1; k <= std::abs(i); ++k)
    x   = std::nextafter(x, dir);
    
return x;
};

float next(float x, int i)
{
if (i == 0)
    return x;

float dir  = (i > 0)   ? std::numeric_limits<float>::infinity()
                        : -std::numeric_limits<float>::infinity();

for (int k = 1; k <= std::abs(i); ++k)
    x   = std::nextafter(x, dir);
    
return x;
};

_declspec(noinline)
ms::simd<double, 256, ms::avx_tag>
test_sin(const ms::simd<double, 256, ms::avx_tag>& x)
{
    //IACA_VC64_START
    auto y = sin(x);
    //IACA_VC64_END
    return y;
};

_declspec(noinline)
void test_iaca()
{
    using simd_type     = ms::simd<double, 256, ms::avx_tag>;

    simd_type x         = simd_type::broadcast(51.44);    
    simd_type y         = test_sin(x);    

    std::cout << y.first() << "\n";
};

int main(int argc, const char* argv[])
{
    (void)argc;
    (void)argv;

    using log_ptr   = std::shared_ptr<std::ofstream>;
    
    test_iaca();

    //matcl::simd::details::exp_table_double::generate_lookup_table(std::cout, 4);

    {
        Matrix x            = randn(100,100);    
        Matrix y            = randn(100,100);

        tic();
        Matrix z            = x * y;
        tocdisp();

        disp(z);
    }

    if (0)
    {
        precision prec      = precision(1300);
        //precision prec      = precision(53);
        precision prec_d    = precision::precision_float();
        precision prec_d2   = prec_d + prec_d;
        mp_float two_by_pi  = 2.0 / constants::mp_pi(prec);
        mp_float pi2        = constants::mp_pi(prec) / 2.0;
        double eps          = std::numeric_limits<double>::epsilon();
        (void)eps;

        //std::cout << two_by_pi.to_string_hex() << "\n";
        //std::cout << pi2.to_string_hex() << "\n";

        //double vd           = std::ldexp(6381956970095103.0, 797);
        for (int p = -100; p <= 128; ++p)
        //for (int p = 0; p < 10000000; ++p)
        {
            double vd0          = (std::ldexp(1.0, p) * pi2).cast_float();            
            //double vd0          = double(p) * pi2.cast_float();
            //double vd           = std::ldexp(matcl::rand(), 50);
            //double vd           = std::ldexp(6381956970095103.0, 797);
            //mp_float vdf        = mp_float(std::ldexp(1.0 + p * eps, 95), prec) * pi2;
            //double vd           = vdf.cast_float();
            //float vd            = std::ldexp(16367173.0f, 72);

            float vd            = (float)vd0;
            mp_float xf         = mp_float(vd, prec);            
            xf                  = xf * two_by_pi;

            mp_float int_part   = round(xf);
            mp_float frac       = xf - int_part;
            frac                = frac * pi2;

            int quadrant2;
            double f2;
            quadrant2           = ms::reduce_pi2_ph_float(double(vd), 53, f2);
            double err          = 0.0;

            //int num_corr_dig    = 64 - 2;
            int num_corr_dig    = 48;

            double true_frac    = frac.cast_float();
            mp_float dif        = frac - f2;
            double true_err     = dif.cast_float();
            mp_float dif2       = dif - err;
            mp_float err2       = mp_float(err, prec);

            double ulp_1        = ulp_distance(frac, f2, prec_d).cast_float();
            double ulp_2        = ulp_distance(frac, mp_float(f2, prec) + err, precision(num_corr_dig)).cast_float();

            //std::cout << frac.to_string_binary(precision(50))<< "\n";

            //std::cout << f << " " << f2 << " " << err << " " << (dif + 1).to_string(precision(100)) << "\n";

            int q1              = uint32_t(quadrant2) % 4;
            mp_float q2         = mod(int_part, 4);

            if (ulp_1 > 0.51 || ulp_2 > 1.0 || q1 != q2)
            //if (ulp_1 > 0.5)
            {
                std::cout << "\n";
                std::cout << p << "\n";
                std::cout << ulp_1 << " " << (dif).to_string(precision(30)) << " " << true_frac << " " << f2 << "\n";
                std::cout << ulp_2 << " " << dif2.to_string(precision(30)) << " " << true_err << " " << err << "\n";
                std::cout << q1 << " " << q2.to_string(precision(10)) << "\n";

                ms::reduce_pi2_ph_float(double(vd), 53, f2);
            };
        }
    }

    if (0)
    {
        int N       = 10000;
        int K       = 20;
        int M       = 100;

        #if (_DEBUG)
            N       = N / 10;
            M       = M / 10;
        #endif

        using func_ptr  = benchmark::function_ptr;

        std::vector<float> x;
        x.reserve(N);

        //double max_x  = constants::pi_2() * std::ldexp(1.0, 800);
        //double max_x  = constants::pi_2() * std::ldexp(1.0, 10);
        float max_x     = constants::f_pi_2() * std::ldexp(1.0f, 10);

        for (int i = 0; i < N; ++i)
            x.push_back(matcl::frand() * max_x);

        const float* ptr_x = x.data();

        benchmark b1    = benchmark(func_ptr(new Func_tan<float>(N, K, ptr_x)));
        benchmark b3    = benchmark(func_ptr(new Func_tan2<float>(N, K, ptr_x)));

        time_stats s3   = b3.make(M);
        time_stats s1   = b1.make(M);

        std::cout << "tan std: " << s1 << "\n";
        std::cout << "tan imp: " << s3 << "\n";
        std::cout << "r e2: " << s1/s3 << "\n";
        std::cout << "\n";

        return 0;
    };

    if (1)
    {
        using simd_type     = ms::simd<float, 256, ms::avx_tag>;

        precision prec      = precision(53*4);
        //precision prec_d    = precision::precision_double();
        precision prec_d  = precision::precision_float();

        int max_k           = (1 << 19);
        //int max_k           = 1;

        mp_float pi2        = 0.5 * constants::mp_pi(prec);
        mp_float inv_pi2    = 1.0 / pi2;
        double max_dif      = 0.0;

        for (int k = 1; k <= max_k; ++k)
        {
            if (k % 1000000 == 0)
                std::cout << k << "\n";
            
            //mp_float xt     = pi2 * mp_float(k);
            //double x0       = xt.cast_float();

            for (int i = 0; i <= 1000000; ++i)
            {
                //double x        = next(x0, i);
                float x         = float(matcl::rand() * double(max_k) * pi2.cast_float() * 0.5);
                simd_type xs    = simd_type(x);
                float res1      = ms::tan(xs).first();

                mp_float x2     = mp_float(x, prec);
                mp_float res2   = tan(x2);

                double dif      = ulp_distance(res2, mp_float(res1), prec_d).cast_float();

                if (dif > max_dif)
                {
                    std::cout   << dif - 1.0 << " " << k << " " << i << " " << x << " " 
                                << mp_float(res1).to_string(precision(20)) << " " 
                                << res2.to_string(precision(20)) << " " << res1 - res2 << "\n";
                    ms::tan(x);

                    max_dif = dif;
                };
            };
        };
    };

    return 0;
}
