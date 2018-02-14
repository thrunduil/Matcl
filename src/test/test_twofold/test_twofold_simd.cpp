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

#include "test_twofold_simd.h"

#include "matcl-scalar/lib_functions/utils.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-scalar/IO/scalar_io.h"

#include "utils.h"
#include "test_functions.h"

#include <vector>

namespace matcl { namespace test
{

void test::test_functions_simd()
{   
    test_twofold_simd().make_unary();    
    test_twofold_simd().make_binary();    
};

void test::test_io_simd()
{
    test_twofold_simd().make_io();
};

void test::test_error_simd()
{
    test_twofold_simd().make_error();
};

int test_twofold_simd::get_size() const
{
    #ifdef _DEBUG
        return 100000;
    #else
        return 1000000;
    #endif
};

int test_twofold_simd::get_size_perf() const
{
    return 1000;
};

int test_twofold_simd::get_num_rep() const
{
    #ifdef _DEBUG
        return 100;
    #else
        return 500000;
    #endif
}

template<class T>
std::string test_twofold_simd::get_header() const
{
    std::string ret = "simd ";

    ret += typeid(T).name();
    return ret;
};

void test_twofold_simd::make_unary()
{
    test_functions<double>();
};

void test_twofold_simd::make_binary()
{
    test_functions_bin<double>();
};

void test_twofold_simd::make_io()
{
    test_functions_io<double>();
};

void test_twofold_simd::make_error()
{
    out_stream << "\n";

    make_error_type<double>();
};

template<class Float_type>
void test_twofold_simd::make_error_type()
{
    using simd_type = simd::simd<Float_type, 256, simd::avx_tag>;
    using twofold   = twofold<simd_type>;    

    simd_type e     = simd_type(matcl::constants::eps<Float_type>());
    simd_type e2    = e * e;

    bool ok         = true;

    simd_type zero  = simd_type(Float_type(0.0));
    simd_type one   = simd_type(Float_type(1.0));
    simd_type two   = simd_type(Float_type(2.0));

    for (int mult = 1; mult < 100; ++mult)
    {
        ok  &= make_error<simd_type>(twofold::normalize_fast(one, e2 - e/two), mult);    

        ok  &= make_error<simd_type>(twofold::normalize_fast(one, zero), mult);
        ok  &= make_error<simd_type>(twofold::normalize_fast(one, e2), mult);
        ok  &= make_error<simd_type>(twofold::normalize_fast(one, -e2), mult);
        ok  &= make_error<simd_type>(twofold::normalize_fast(one, e/two - e2), mult);
        ok  &= make_error<simd_type>(twofold::normalize_fast(one, e2 - e/two), mult);    

        ok  &= make_error<simd_type>(twofold::normalize_fast(-one, zero), mult);
        ok  &= make_error<simd_type>(twofold::normalize_fast(-one, e2), mult);
        ok  &= make_error<simd_type>(twofold::normalize_fast(-one, -e2), mult);
        ok  &= make_error<simd_type>(twofold::normalize_fast(-one, e2 - e/two), mult);
        ok  &= make_error<simd_type>(twofold::normalize_fast(-one, e/two - e2), mult);
    };

    if (ok == true)
        out_stream << "error simd: OK" << "\n";
    else
        out_stream << "error simd: FAIL" << "\n";
};

template<class Simd_type>
static bool neq_bool(const Simd_type& x, const Simd_type& y)
{
    return any(neq(x, y));
};

template<class Simd_type>
bool test_twofold_simd::make_error(const twofold<Simd_type>& x, int mult0)
{
    using twofold   = twofold<Simd_type>;

    Simd_type zero  = Simd_type(0.0);
    Simd_type two   = Simd_type(2.0);

    Simd_type mult  = Simd_type(mult0);
    Simd_type e     = eps(x) * mult;

    twofold x1      = x + e;
    twofold x2      = x - e;

    Simd_type d1   = float_distance(x, x1);
    Simd_type d2   = float_distance(x, x2);
    Simd_type d3   = float_distance(x1, x);
    Simd_type d4   = float_distance(x2, x);
    Simd_type d5   = float_distance(x1, x2);
    Simd_type d6   = float_distance(x2, x1);
    Simd_type d7   = float_distance(x1, x1);
    Simd_type d8   = float_distance(x2, x2);

    if (neq_bool(d1, mult) || neq_bool(d2, mult) 
        || neq_bool(d3, mult) || neq_bool(d4, mult))
    {
        return false;
    }

    if (neq_bool(d5, two * mult) || neq_bool(d6, two * mult))
        return false;

    if (neq_bool(d7, zero) || neq_bool(d8, zero))
        return false;

    return true;
};

template<class T>
void test_twofold_simd::test_functions()
{
    int N           = get_size();

    using T2        = twofold<T>;

    std::vector<T>  in_1;
    std::vector<T>  in_2;

    std::vector<T>  out_1;
    std::vector<T>  out_2;
    std::vector<T>  out_base;

    in_1.resize(N);
    in_2.resize(N);
    out_1.resize(N);
    out_2.resize(N);
    out_base.resize(N);

    T* ptr_in_1     = in_1.data();
    T* ptr_in_2     = in_2.data();
    T* ptr_out_1    = out_1.data();
    T* ptr_out_2    = out_2.data();
    T* ptr_out_base = out_base.data();

    for (int i = 0; i < N; ++i)
    {
        T2 tmp      = rand_scalar<T2>::make(true, false);
        ptr_in_1[i] = tmp.value;
        ptr_in_2[i] = tmp.error;
    }

    std::string header  = get_header<T>();

    disp(" ");
    disp(header);

    formatted_disp dm;

    dm.set_row_label("func",        align_type::right, 10);
    dm.add_column("t sse",          align_type::left, 5);
    dm.add_column("t avx",          align_type::left, 5);
    dm.add_column("rt (TF sse)",    align_type::left, 5);
    dm.add_column("rt (TF avx)",    align_type::left, 5);
    dm.add_column("ulp err",        align_type::left, 5);
    dm.add_column("status",         align_type::left, 5);

    dm.disp_header();

    test_function<T, test_functions::Func_uminus_2>
                (dm, N, ptr_in_1, ptr_in_2, ptr_out_1, ptr_out_2, ptr_out_base);
    test_function<T, test_functions::Func_sum<T>>
                (dm, N, ptr_in_1, ptr_in_2, ptr_out_1, ptr_out_2, ptr_out_base);
    test_function<T, test_functions::Func_abs>
                (dm, N, ptr_in_1, ptr_in_2, ptr_out_1, ptr_out_2, ptr_out_base);

    test_function<T, test_functions::Func_sqrt_1>
                (dm, N, ptr_in_1, ptr_in_2, ptr_out_1, ptr_out_2, ptr_out_base);
    test_function<T, test_functions::Func_sqrt_2>
                (dm, N, ptr_in_1, ptr_in_2, ptr_out_1, ptr_out_2, ptr_out_base);
    test_function<T, test_functions::Func_inv_1>
                (dm, N, ptr_in_1, ptr_in_2, ptr_out_1, ptr_out_2, ptr_out_base);
    test_function<T, test_functions::Func_inv_2>
                (dm, N, ptr_in_1, ptr_in_2, ptr_out_1, ptr_out_2, ptr_out_base);
};

template<class T>
void test_twofold_simd::test_functions_io()
{
    int N           = get_size();

    using T2        = twofold<T>;

    std::vector<T>  in_1;
    std::vector<T>  in_2;

    std::vector<T>  out_1;
    std::vector<T>  out_2;
    std::vector<T>  out_base;

    in_1.resize(N);
    in_2.resize(N);
    out_1.resize(N);
    out_2.resize(N);
    out_base.resize(N);

    T* ptr_in_1     = in_1.data();
    T* ptr_in_2     = in_2.data();
    T* ptr_out_1    = out_1.data();
    T* ptr_out_2    = out_2.data();
    T* ptr_out_base = out_base.data();

    for (int i = 0; i < N; ++i)
    {
        T2 tmp      = rand_scalar<T2>::make(true, false);
        ptr_in_1[i] = tmp.value;
        ptr_in_2[i] = tmp.error;
    }

    std::string header  = get_header<T>();

    disp(" ");
    disp(header);

    formatted_disp dm;

    dm.set_row_label("func",        align_type::right, 10);
    dm.add_column("t sse",          align_type::left, 5);
    dm.add_column("t avx",          align_type::left, 5);
    dm.add_column("rt (TF sse)",    align_type::left, 5);
    dm.add_column("rt (TF avx)",    align_type::left, 5);
    dm.add_column("ulp err",        align_type::left, 5);
    dm.add_column("status",         align_type::left, 5);

    dm.disp_header();

    test_function_io<T, test_functions::Func_save_load>
                (dm, N, ptr_in_1, ptr_in_2, ptr_out_1, ptr_out_2, ptr_out_base);
};

template<class T>
void test_twofold_simd::test_functions_bin()
{
    std::string header  = get_header<T>();

    disp(" ");
    disp(header);

    formatted_disp dm;

    dm.set_row_label("func",        align_type::right, 10);
    dm.add_column("t sse",          align_type::left, 5);
    dm.add_column("t avx",          align_type::left, 5);
    dm.add_column("rt (TF sse)",    align_type::left, 5);
    dm.add_column("rt (TF avx)",    align_type::left, 5);
    dm.add_column("ulp err",        align_type::left, 5);
    dm.add_column("status",         align_type::left, 5);

    dm.disp_header();

    int N   = get_size();

    using T2        = twofold<T>;

    std::vector<T>  in_11;
    std::vector<T>  in_12;
    std::vector<T>  in_21;
    std::vector<T>  in_22;

    std::vector<T>  out_1;
    std::vector<T>  out_2;
    std::vector<T>  out_base;

    in_11.resize(N);
    in_12.resize(N);
    in_21.resize(N);
    in_22.resize(N);
    out_1.resize(N);
    out_2.resize(N);
    out_base.resize(N);

    T* ptr_in_11    = in_11.data();
    T* ptr_in_12    = in_12.data();
    T* ptr_in_21    = in_21.data();
    T* ptr_in_22    = in_22.data();
    T* ptr_out_1    = out_1.data();
    T* ptr_out_2    = out_2.data();
    T* ptr_out_base = out_base.data();

    test_functions_sum_bin<T>
        (dm, N, ptr_in_11, ptr_in_12, ptr_in_21, ptr_in_22, ptr_out_1, ptr_out_2, ptr_out_base);
    test_functions_bin<T>
        (dm, N, ptr_in_11, ptr_in_12, ptr_in_21, ptr_in_22, ptr_out_1, ptr_out_2, ptr_out_base);
};

template<class T, class Func>
void test_twofold_simd::test_function(formatted_disp& fd, int size, 
                    const T* in_1, const T* in_2, T* out_1, T* out_2, T* out_1_base)
{
    using TS1    = simd::simd<T, 128, simd::sse_tag>;
    using TS2    = simd::simd<T, 256, simd::avx_tag>;
    using TS3    = simd::simd<T, 128, simd::sse_tag>;
    using TS4    = simd::simd<T, 256, simd::avx_tag>;

    test_function_base<T, TS1, Func>(size, 1, in_1, out_1_base);
    test_function_twofold<T, TS3, Func>(size, 1, in_1, in_2, out_1, out_2);

    double max_dist = 2.0;
    double d1;
    bool ok     = test_equal(size, out_1, out_1_base, max_dist, d1);

    int N       = get_size_perf();
    int M       = get_num_rep();

    double t1   = test_function_base<T, TS1, Func>(N, M, in_1, out_1_base);
    double t2   = test_function_base<T, TS2, Func>(N, M, in_1, out_1_base);
    double t3   = test_function_twofold<T, TS3, Func>(N, M, in_1, in_2, out_1, out_2);
    double t4   = test_function_twofold<T, TS4, Func>(N, M, in_1, in_2, out_1, out_2);

    std::string status  = (ok == true) ? "OK" : "FAIL"; 
    fd.disp_row(Func::name(), t1, t2, t3/t1, t4/t2, d1, status);
};

template<class T, class Func>
void test_twofold_simd::test_function_io(formatted_disp& fd, int size, 
                    const T* in_1, const T* in_2, T* out_1, T* out_2, T* out_1_base)
{
    using TS1    = simd::simd<T, 128, simd::sse_tag>;
    using TS2    = simd::simd<T, 256, simd::avx_tag>;
    using TS3    = simd::simd<T, 128, simd::sse_tag>;
    using TS4    = simd::simd<T, 256, simd::avx_tag>;

    test_function_base<T, TS1, Func>(size, 1, in_1, out_1_base);
    test_function_twofold<T, TS3, Func>(size, 1, in_1, in_2, out_1, out_2);

    double max_dist = 0.0;
    double d1;
    bool ok     = test_equal(size, out_1, out_1_base, max_dist, d1);

    int N       = get_size_perf();
    int M       = 10;

    double t1   = test_function_base<T, TS1, Func>(N, M, in_1, out_1_base);
    double t2   = test_function_base<T, TS2, Func>(N, M, in_1, out_1_base);
    double t3   = test_function_twofold<T, TS3, Func>(N, M, in_1, in_2, out_1, out_2);
    double t4   = test_function_twofold<T, TS4, Func>(N, M, in_1, in_2, out_1, out_2);

    std::string status  = (ok == true) ? "OK" : "FAIL"; 
    fd.disp_row(Func::name(), t1, t2, t3/t1, t4/t2, d1, status);
};

template<class T, class TS, class Func>
double test_twofold_simd::test_function_base(int size, int n_rep, const T* in_1, T* out_1)
{
    tic();

    static const int vec_size = TS::vector_size;

    volatile T val = 0;

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; i += vec_size)
        {
            TS x    = TS::load(in_1 + i, std::false_type());
            TS res  = Func::eval(x);

            res.store(out_1 + i, std::false_type());
        };

        val += out_1[0];
    };

    double t = toc();
    return t;
};

template<class T, class TS, class Func>
double test_twofold_simd::test_function_twofold(int size, int n_rep, const T* in_1, const T* in_2, 
                                                T* out_1, T* out_2)
{
    tic();

    using T2        = twofold<TS>;

    static const int vec_size = TS::vector_size;

    volatile T val = 0;

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; i += vec_size)
        {
            TS v    = TS::load(in_1 + i, std::false_type());
            TS e    = TS::load(in_2 + i, std::false_type());

            T2 x    = T2(v,e);
            T2 res  = Func::eval(x);

            res.value.store(out_1 + i, std::false_type());
            res.error.store(out_2 + i, std::false_type());
        };

        val += out_1[0];
        val += out_2[0];
    };

    double t = toc();
    return t;
};

template<class T>
void test_twofold_simd::test_functions_bin(formatted_disp& dm, int N, 
                T* in_11, T* in_12, T* in_21, T* in_22, T* out_1, T* out_2, T* out_base)
{
    using T2        = twofold<T>;

    for (int i = 0; i < N; ++i)
    {
       T2 tmp1      = rand_scalar<T2>::make(true, false);
       T2 tmp2      = rand_scalar<T2>::make(true, false);

       in_11[i]     = tmp1.value;
       in_12[i]     = tmp1.error;

       in_21[i]     = tmp2.value;
       in_22[i]     = tmp2.error;
    }

    test_function_bin<T, test_functions::Func_mult_11>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, false);
    test_function_bin<T, test_functions::Func_div_11>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, false);

    test_function_bin<T, test_functions::Func_mult_12>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, false);
    test_function_bin<T, test_functions::Func_mult_21>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, false);
    test_function_bin<T, test_functions::Func_mult_22>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, false);

    test_function_bin<T, test_functions::Func_div_12>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, false);
    test_function_bin<T, test_functions::Func_div_21>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, false);
    test_function_bin<T, test_functions::Func_div_22>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, false);
};

template<class T>
void test_twofold_simd::test_functions_sum_bin(formatted_disp& dm, int N, 
                T* in_11, T* in_12, T* in_21, T* in_22, T* out_1, T* out_2, T* out_base)
{
    using T2        = twofold<T>;

    for (int i = 0; i < N; ++i)
    {
       T2 v1        = rand_scalar<T2>::make(true, false);
       T2 v2        = rand_scalar<T2>::make(true, true);
       v2           = v2 * v1;

       if (i % 2 == 0)
       {
           in_11[i]     = v1.value;
           in_12[i]     = v1.error;
           in_21[i]     = v2.value;
           in_22[i]     = v2.error;
       }
       else
       {
           in_11[i]     = v2.value;
           in_12[i]     = v2.error;
           in_21[i]     = v1.value;
           in_22[i]     = v1.error;
       }
    }

    test_function_bin<T, test_functions::Func_plus_11>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, false);
    test_function_bin<T, test_functions::Func_plus_11_sort>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, false);
    test_function_bin<T, test_functions::Func_minus_11>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, false);
    test_function_bin<T, test_functions::Func_minus_11_sort>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, false);

    test_function_bin<T, test_functions::Func_plus_12>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, true);
    test_function_bin<T, test_functions::Func_plus_21>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, true);
    test_function_bin<T, test_functions::Func_plus_22>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, true);

    test_function_bin<T, test_functions::Func_minus_12>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, true);
    test_function_bin<T, test_functions::Func_minus_21>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, true);
    test_function_bin<T, test_functions::Func_minus_22>
        (dm, N, in_11, in_12, in_21, in_22, out_1, out_2, out_base, true);
};

template<class T, class Func>
void test_twofold_simd::test_function_bin(formatted_disp& fd, int size, 
                const T* in_11, const T* in_12, const T* in_21, const T* in_22, 
                T* out_1, T* out_2, T* out_base, bool inaccurate)
{
    using TS1    = simd::simd<T, 128, simd::sse_tag>;
    using TS2    = simd::simd<T, 256, simd::avx_tag>;
    using TS3    = simd::simd<T, 128, simd::sse_tag>;
    using TS4    = simd::simd<T, 256, simd::avx_tag>;

   test_function_bin_base<T, TS2, Func>
         (size, 1, in_11, in_21, out_base);
   test_function_bin_twofold<T, TS4, Func>
         (size, 1, in_11, in_12, in_21, in_22, out_1, out_2);

    double max_dist = 2.0;
    double d1;
    bool ok;
    
    if (inaccurate == false)
        ok = test_equal(size, out_1, out_base, max_dist, d1);
    else
        ok = test_equal_inacc(size, in_11, in_21, out_1, out_base, max_dist, d1);

    int N   = get_size_perf();
    int M   = get_num_rep();

    double t1   = test_function_bin_base<T, TS1, Func>
                        (N, M, in_11, in_21, out_base);
    double t2   = test_function_bin_base<T, TS2, Func>
                        (N, M, in_11, in_21, out_base);
    double t3   = test_function_bin_twofold<T, TS3, Func>
                        (N, M, in_11, in_12, in_21, in_22, out_1, out_2);
    double t4   = test_function_bin_twofold<T, TS4, Func>
                        (N, M, in_11, in_12, in_21, in_22, out_1, out_2);

    std::string status  = (ok == true) ? "OK" : "FAIL"; 
    fd.disp_row(Func::name(), t1, t2, t3/t1, t4/t2, d1, status);
};

template<class T, class TS, class Func>
double test_twofold_simd::test_function_bin_base(int size, int n_rep, const T* in_1, 
                                             const T* in_2, T* out_1)
{
    tic();

    static const int vec_size = TS::vector_size;

    volatile T val = 0;

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; i += vec_size)
        {
            TS x1    = TS::load(in_1 + i, std::false_type());
            TS x2    = TS::load(in_2 + i, std::false_type());
            TS res  = Func::eval(x1, x2);

            res.store(out_1 + i, std::false_type());
        };

        val += out_1[0];
    };

    double t = toc();
    return t;
};

template<class T, class TS, class Func>
double test_twofold_simd::test_function_bin_twofold(int size, int n_rep, 
                    const T* in_11, const T* in_12, const T* in_21, const T* in_22,
                    T* out_1, T* out_2)
{
    tic();

    using T2        = twofold<TS>;

    static const int vec_size = TS::vector_size;

    volatile T val = 0;

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; i += vec_size)
        {
            TS v1   = TS::load(in_11 + i, std::false_type());
            TS e1   = TS::load(in_12 + i, std::false_type());
            TS v2   = TS::load(in_21 + i, std::false_type());
            TS e2   = TS::load(in_22 + i, std::false_type());

            T2 x1   = T2(v1, e1);
            T2 x2   = T2(v2, e2);

            T2 res  = Func::eval(x1, x2);

            res.value.store(out_1 + i, std::false_type());
            res.error.store(out_2 + i, std::false_type());
        };

        val += out_1[0];
        val += out_2[0];
    };

    double t = toc();
    return t;
};

template<class T>
bool test_twofold_simd::test_equal(int size, const T* res, const T* res_gen, 
                                   double max_dist, double& dist)
{
    bool eq         = true;
    dist            = 0.0;

    double loc_dist;

    for (int i = 0; i < size; ++i)
    {
        bool tmp    = test_equal(res[i], res_gen[i], max_dist, loc_dist);
        eq          = eq && tmp;
        dist        = std::max(dist, loc_dist);
    };

    return eq;
}

template<class T>
bool test_twofold_simd::test_equal_inacc(int size, const T* in_11, const T* in_21, 
                    const T* res, const T* res_gen, double max_dist, double& dist)
{
    bool eq         = true;
    dist            = 0.0;

    double loc_dist;

    for (int i = 0; i < size; ++i)
    {
        bool tmp    = test_equal_inacc(in_11[i], in_21[i], res[i], res_gen[i], max_dist, loc_dist);
        eq          = eq && tmp;
        dist        = std::max(dist, loc_dist);
    };

    return eq;
}

template<class T>
bool test_twofold_simd::test_equal(const T& res, const T& res_gen, double max_dist, double& dist)
{
    dist = 0.0;

    if (res == res_gen)
        return true;

    // twofold does not handle INV/NAN values correctly, when res is not finite
    // then true result can still be finite
    if (is_finite(res) == false)
        return true;
    if (is_finite(res_gen) == false)
        return true;

    dist = float_distance(res_gen, res);

    if (is_nan(dist) == true)
    {
        //std::cout << dist << " " << res_mp << " " << res_gen << " " << res_mp - res_gen << "\n";
        return false;
    };

    if (dist > max_dist)
    {
        //std::cout << dist << " " << res_mp << " " << res_gen << " " << res_mp - res_gen << "\n";
        return false;
    }
    else
        return true;
}

template<class T>
bool test_twofold_simd::test_equal_inacc(const T& in_11, const T& in_21, const T& res, 
                                      const T& res_gen, double max_dist, double& dist)
{
    dist = 0.0;

    if (res == res_gen)
        return true;

    // twofold does not handle INV/NAN values correctly, when res is not finite
    // then true result can still be finite
    if (is_finite(res) == false)
        return true;
    if (is_finite(res_gen) == false)
        return true;

    double sum_abs  = std::abs(in_11) + std::abs(in_21);

    (void)sum_abs;

    dist = std::abs(res_gen - res) / matcl::eps(sum_abs);

    if (is_nan(dist) == true)
    {
        //std::cout << dist << " " << res_mp << " " << res_gen << " " << res_mp - res_gen << "\n";
        return false;
    };

    if (dist > max_dist)
    {
        //std::cout << dist << " " << res_mp << " " << res_gen << " " << res_mp - res_gen << "\n";
        return false;
    }
    else
        return true;
}

}};
