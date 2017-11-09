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

#include "test_twofold.h"
#include "matcl-mp/matcl_mp.h"
#include "matcl-scalar/IO/scalar_io.h"
#include "matcl-scalar/lib_functions/utils.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/lib_functions/func_binary.h"

#include "utils.h"
#include "test_functions.h"

#include <vector>

namespace matcl { namespace test
{

void test::test_fma()
{
    test_twofold(false).make_fma();
};

void test::test_double()
{   
    test_twofold(true).make_unary();
    test_twofold(false).make_unary();

    test_twofold(true).make_binary();
    test_twofold(false).make_binary();
};

void test::test_io()
{
    test_twofold(false).make_io();
};

void test::test_error()
{
    test_twofold(false).make_error();
};

test_twofold::test_twofold(bool normalized)
    :m_normalized(normalized)
{};

void test_twofold::make_unary()
{
    test_functions<double>();
};

void test_twofold::make_binary()
{
    test_functions_bin<double>();
};

void test_twofold::make_fma()
{
    test_functions_fma<double>();
    test_functions_fma<float>();
};

void test_twofold::make_io()
{
    test_functions_io<double>();
};

void test_twofold::make_error()
{
    twofold x   = twofold(1.0, 0.0);

    double e    = matcl::constants::eps();
    double e2   = e * e;

    bool ok     = true;

    for (int mult = 1; mult < 100; ++mult)
    {
        ok  &= make_error(twofold::normalize_fast(1.0, e2 - e/2.0), mult);    

        ok  &= make_error(twofold::normalize_fast(1.0, 0.0), mult);
        ok  &= make_error(twofold::normalize_fast(1.0, e2), mult);
        ok  &= make_error(twofold::normalize_fast(1.0, -e2), mult);
        ok  &= make_error(twofold::normalize_fast(1.0, e/2.0 - e2), mult);
        ok  &= make_error(twofold::normalize_fast(1.0, e2 - e/2.0), mult);    

        ok  &= make_error(twofold::normalize_fast(-1.0, 0.0), mult);
        ok  &= make_error(twofold::normalize_fast(-1.0, e2), mult);
        ok  &= make_error(twofold::normalize_fast(-1.0, -e2), mult);
        ok  &= make_error(twofold::normalize_fast(-1.0, e2 - e/2.0), mult);
        ok  &= make_error(twofold::normalize_fast(-1.0, e/2.0 - e2), mult);
    };

    if (ok == true)
        out_stream << "error: OK" << "\n";
    else
        out_stream << "error: FAIL" << "\n";
};

bool test_twofold::make_error(const twofold& x, int mult)
{
    double e    = eps(x) * mult;

    twofold x1  = x + e;
    twofold x2  = x - e;

    double d1   = float_distance(x, x1);
    double d2   = float_distance(x, x2);
    double d3   = float_distance(x1, x);
    double d4   = float_distance(x2, x);
    double d5   = float_distance(x1, x2);
    double d6   = float_distance(x2, x1);
    double d7   = float_distance(x1, x1);
    double d8   = float_distance(x2, x2);

    if (d1 != mult || d2 != mult || d3 != mult || d4 != mult)
        return false;

    if (d5 != 2*mult || d6 != 2*mult)
        return false;

    if (d7 != 0 || d8 != 0)
        return false;

    return true;
};

template<class T>
void test_twofold::test_functions()
{
    int N   = get_N();
    int M   = get_M();

    using T2    = twofold;
    using TMP   = mp_float;

    std::vector<T2> in;
    std::vector<T2> out;
    std::vector<TMP> out_gen;

    in.resize(N);
    out.resize(N);
    out_gen.resize(N);

    T2* ptr_in          = in.data();
    T2* ptr_out         = out.data();
    TMP* ptr_out_gen    = out_gen.data();

    for (int i = 0; i < N; ++i)
        ptr_in[i]   = rand_scalar<T2>::make(m_normalized, false);

    std::string header  = get_header<T>();

    disp(" ");
    disp(header);

    formatted_disp dm;

    dm.set_row_label("func",    align_type::right, 15);
    dm.add_column("t double",   align_type::left, 5);
    dm.add_column("t(TF)/t(D)", align_type::left, 5);
    dm.add_column("ulp error",  align_type::left, 5);
    dm.add_column("status",     align_type::left, 5);

    dm.disp_header();

    test_function<T, T2, TMP, test_functions::Func_uminus_2>
                (dm, N, M, ptr_in, ptr_out, ptr_out_gen, 0.0);
    test_function<T, T2, TMP, test_functions::Func_sum>
                (dm, N, M, ptr_in, ptr_out, ptr_out_gen, 0.0);
    test_function<T, T2, TMP, test_functions::Func_abs>
                (dm, N, M, ptr_in, ptr_out, ptr_out_gen, 0.0);

    test_function<T, T2, TMP, test_functions::Func_sqrt_1>
                (dm, N, M, ptr_in, ptr_out, ptr_out_gen, 0.5);
    test_function<T, T2, TMP, test_functions::Func_sqrt_2>
                (dm, N, M, ptr_in, ptr_out, ptr_out_gen, 2.0);
};

template<class T>
void test_twofold::test_functions_io()
{
    int N   = get_N();

    using T2    = twofold;

    std::vector<T2> in;
    std::vector<T2> out;

    in.resize(N);
    out.resize(N);

    T2* ptr_in          = in.data();
    T2* ptr_out         = out.data();

    for (int i = 0; i < N; ++i)
        ptr_in[i]   = rand_scalar<T2>::make(m_normalized, false);

    std::string header  = get_header<T>();

    disp(" ");
    disp(header);

    formatted_disp dm;

    dm.set_row_label("func",    align_type::right, 15);
    dm.add_column("t double",   align_type::left, 5);
    dm.add_column("t(TF)/t(D)", align_type::left, 5);
    dm.add_column("ulp error",  align_type::left, 5);
    dm.add_column("status",     align_type::left, 5);

    dm.disp_header();

    test_function_io<T, T2, test_functions::Func_save_load>
                (dm, N, 1, ptr_in, ptr_out, 0.0);
};

template<class T>
void test_twofold::test_functions_bin()
{
    std::string header  = get_header<T>();

    disp(" ");
    disp(header);

    formatted_disp dm;

    dm.set_row_label("func",    align_type::right, 15);
    dm.add_column("t double",   align_type::left, 5);
    dm.add_column("t(TF)/t(D)", align_type::left, 5);
    dm.add_column("ulp error",  align_type::left, 5);
    dm.add_column("status",     align_type::left, 5);

    dm.disp_header();    

    int N   = get_N();
    int M   = get_M();

    using T2    = twofold;
    using TMP   = mp_float;

    std::vector<T2> in_1;
    std::vector<T2> in_2;
    std::vector<T2> out;
    std::vector<TMP> out_gen;

    in_1.resize(N);
    in_2.resize(N);
    out.resize(N);
    out_gen.resize(N);

    T2* ptr_in_1     = in_1.data();
    T2* ptr_in_2     = in_2.data();
    T2* ptr_out      = out.data();
    TMP* ptr_out_gen = out_gen.data();

    test_functions_bin<T, T2, TMP>(dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen);
    test_functions_sum_bin<T, T2, TMP>(dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen);
};

template<class T>
void test_twofold::test_functions_fma()
{
    std::string header  = get_header<T>();

    disp(" ");
    disp(header);

    formatted_disp dm;

    dm.set_row_label("func",    align_type::right, 15);
    dm.add_column("t double",   align_type::left, 5);
    dm.add_column("t(TF)/t(D)", align_type::left, 5);
    dm.add_column("ulp error",  align_type::left, 5);
    dm.add_column("status",     align_type::left, 5);

    dm.disp_header();    

    int N   = get_N();
    int M   = get_M();

    using T2    = twofold;
    using TMP   = mp_float;

    std::vector<T> in_1;
    std::vector<T> in_2;
    std::vector<T> in_3;
    std::vector<T> out;
    std::vector<TMP> out_gen;

    in_1.resize(N);
    in_2.resize(N);
    out.resize(N);
    out_gen.resize(N);

    T* ptr_in_1     = in_1.data();
    T* ptr_in_2     = in_2.data();
    T* ptr_in_3     = in_2.data();
    T* ptr_out      = out.data();
    TMP* ptr_out_gen = out_gen.data();

    for (int i = 0; i < N; ++i)
    {
        T s1        = rand_scalar<T>::make(false);
        T s2        = rand_scalar<T>::make(false);
        T s3        = rand_scalar<T>::make(true);

        if (i % 2 == 0)
            s3      = s1 * s2 * s3;
        else
            s3      = s1 * s2 / s3;

        ptr_in_1[i] = static_cast<T>(s1);
        ptr_in_2[i] = static_cast<T>(s2);
        ptr_in_3[i] = static_cast<T>(s3);
    };

    test_function_fma<T, TMP, test_functions::Func_fma>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_in_3, ptr_out, ptr_out_gen, 0.5);
};

template<class T, class T2, class TMP>
void test_twofold::test_functions_bin(formatted_disp& dm, int N, int M, 
                T2* in_1, T2* in_2, T2* out, TMP* out_gen)
{
    T2* ptr_in_1     = in_1;
    T2* ptr_in_2     = in_2;
    T2* ptr_out      = out;
    TMP* ptr_out_gen = out_gen;

    for (int i = 0; i < N; ++i)
    {
        ptr_in_1[i] = rand_scalar<T2>::make(m_normalized, false);
        ptr_in_2[i] = rand_scalar<T2>::make(m_normalized, false);
    };

    test_function_bin<T, T2, TMP, test_functions::Func_mult_11>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 0.0);
    test_function_bin<T, T2, TMP, test_functions::Func_div_11>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 0.5);

    test_function_bin<T, T2, TMP, test_functions::Func_mult_12>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 1.5);
    test_function_bin<T, T2, TMP, test_functions::Func_mult_21>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 1.5);
    test_function_bin<T, T2, TMP, test_functions::Func_mult_22>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 4.0);

    test_function_bin<T, T2, TMP, test_functions::Func_div_12>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 4.0);
    test_function_bin<T, T2, TMP, test_functions::Func_div_21>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 4.0);
    test_function_bin<T, T2, TMP, test_functions::Func_div_22>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 4.0);

    test_function_bin<T, T2, TMP, test_functions::Func_mult_dekker>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 0.0);
};

template<class T, class T2, class TMP>
void test_twofold::test_functions_sum_bin(formatted_disp& dm, int N, int M, 
                T2* in_1, T2* in_2, T2* out, TMP* out_gen)
{
    T2* ptr_in_1     = in_1;
    T2* ptr_in_2     = in_2;
    T2* ptr_out      = out;
    TMP* ptr_out_gen = out_gen;

    for (int i = 0; i < N; ++i)
    {
        T2 v1       = rand_scalar<T2>::make(m_normalized, false);
        T2 v2       = rand_scalar<T2>::make(m_normalized, true);
        v2          = v2 * v1;

        if (i % 2 == 0)
        {
            ptr_in_1[i] = v1;
            ptr_in_2[i] = v2;
        }
        else
        {
            ptr_in_1[i] = v2;
            ptr_in_2[i] = v1;
        }
    };

    test_function_bin<T, T2, TMP, test_functions::Func_plus_11>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 0.0);
    test_function_bin<T, T2, TMP, test_functions::Func_plus_11_sort>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 0.0);
    test_function_bin<T, T2, TMP, test_functions::Func_minus_11>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 0.0);
    test_function_bin<T, T2, TMP, test_functions::Func_minus_11_sort>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 0.0);    

    test_function_bin<T, T2, TMP, test_functions::Func_plus_12>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 0.5);
    test_function_bin<T, T2, TMP, test_functions::Func_plus_21>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 0.5);
    test_function_bin<T, T2, TMP, test_functions::Func_plus_22>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 0.5);

    test_function_bin<T, T2, TMP, test_functions::Func_minus_12>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 0.5);
    test_function_bin<T, T2, TMP, test_functions::Func_minus_21>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 0.5);
    test_function_bin<T, T2, TMP, test_functions::Func_minus_22>
        (dm, N, M, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 0.5);
};

template<class T, class T2, class TMP, class Func>
void test_twofold::test_function(formatted_disp& fd, int size, int n_rep, 
                    const T2* in, T2* out, TMP* out_gen, double max_dist)
{
    double t0   = test_function_gen<T2, TMP, Func>(size, 1, in, out_gen);
    double t1   = test_function_base<T, T2, Func>(size, n_rep, in, out);
    double t2   = test_function_twofold<T2, Func>(size, n_rep, in, out);

    (void)t0;

    double d1;
    bool ok1    = test_equal(size, out, out_gen, max_dist, d1);
    bool ok2    = test_constraints_1<Func>(size, out, in);
    bool ok     = ok1 && ok2;

    std::string status  = (ok == true) ? "OK" : "FAIL"; 
    fd.disp_row(Func::name(), t1, t2/t1, d1, status);
};

template<class T, class T2, class Func>
void test_twofold::test_function_io(formatted_disp& fd, int size, int n_rep, 
                    const T2* in, T2* out, double max_dist)
{
    double t1   = test_function_base<T, T2, Func>(size, n_rep, in, out);
    double t2   = test_function_twofold<T2, Func>(size, n_rep, in, out);

    double d1;
    bool ok1    = test_equal(size, out, in, max_dist, d1);
    bool ok     = ok1;

    std::string status  = (ok == true) ? "OK" : "FAIL"; 
    fd.disp_row(Func::name(), t1, t2/t1, d1, status);
};

template<class T, class T2, class TMP, class Func>
void test_twofold::test_function_bin(formatted_disp& fd, int size, int n_rep, 
                const T2* in_1, const T2* in_2, T2* out, TMP* out_gen, double max_dist)
{
    double t0   = test_function_bin_gen<T2, TMP, Func>(size, 1, in_1, in_2, out_gen);
    double t1   = test_function_bin_base<T, T2, Func>(size, n_rep, in_1, in_2, out);
    double t2   = test_function_bin_twofold<T2, Func>(size, n_rep, in_1, in_2, out);

    (void)t0;

    double d1;
    bool ok1    = test_equal(size, out, out_gen, max_dist, d1);
    bool ok2    = test_constraints_2<Func>(size, out, in_1, in_2);
    bool ok     = ok1 && ok2;

    std::string status  = (ok == true) ? "OK" : "FAIL"; 
    fd.disp_row(Func::name(), t1, t2/t1, d1, status);
};

template<class T, class TMP, class Func>
void test_twofold::test_function_fma(formatted_disp& fd, int size, int n_rep, 
                const T* in_1, const T* in_2, const T* in_3, T* out, TMP* out_gen, double max_dist)
{
    double t0   = test_function_fma_gen<T, TMP, Func>(size, 1, in_1, in_2, in_3, out_gen);
    double t1   = test_function_fma_base<T, Func>(size, n_rep, in_1, in_2, in_3, out);
    double t2   = test_function_fma_twofold<T, Func>(size, n_rep, in_1, in_2, in_3, out);

    (void)t0;

    double d1;
    bool ok1    = test_equal_fma(size, out, out_gen, max_dist, d1);
    bool ok     = ok1;

    std::string status  = (ok == true) ? "OK" : "FAIL"; 
    fd.disp_row(Func::name(), t1, t2/t1, d1, status);
};

template<class Twofold>
struct convert_to_mp{};

template<>
struct convert_to_mp<twofold>
{
    using twofold_type  = twofold;

    static mp_float eval(const twofold_type& x)
    {
        precision prec = precision(53 * 3);

        if (matcl::is_finite(x.value) == false)
            return mp_float(x.value, prec);

        mp_float res    = plus(mp_float(x.value), mp_float(x.error), prec);
        return res;
    };
};

template<>
struct convert_to_mp<double>
{
    static mp_float eval(const double& x)
    {
        precision prec = precision(53 * 3);

        if (matcl::is_finite(x) == false)
            return mp_float(x, prec);

        mp_float res    = mp_float(x, prec);
        return res;
    };
};

template<>
struct convert_to_mp<float>
{
    static mp_float eval(const float& x)
    {
        precision prec = precision(53 * 3);

        if (matcl::is_finite(x) == false)
            return mp_float(x, prec);

        mp_float res    = mp_float(x, prec);
        return res;
    };
};

template<class T2, class TMP, class Func>
double test_twofold::test_function_gen(int size, int n_rep, const T2* in, TMP* out)
{
    tic();

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; ++i)
        {
            T2 xc       = Func::convert_arg_1(in[i]);
            TMP x       = convert_to_mp<T2>::eval(xc);
            TMP res     = Func::eval(x);
            out[i]      = res;
        };
    };

    double t = toc();
    return t;
};

template<class T, class T2, class Func>
double test_twofold::test_function_base(int size, int n_rep, const T2* in, T2* out)
{
    tic();

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; ++i)
        {
            out[i].value= Func::eval(in[i].value);
        };
    };

    double t = toc();
    return t;
};

template<class T2, class Func>
double test_twofold::test_function_twofold(int size, int n_rep, const T2* in, T2* out)
{
    tic();

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; ++i)
        {
            T2 x        = in[i];
            T2 res      = Func::eval(x);
            out[i]      = res;
        };
    };

    double t = toc();
    return t;
};

template<class T2, class TMP, class Func>
double test_twofold::test_function_bin_gen(int size, int n_rep, const T2* in1, const T2* in2, TMP* out)
{
    tic();

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; ++i)
        {
            T2 x        = Func::convert_arg_1(in1[i]);
            T2 y        = Func::convert_arg_2(in2[i]);
            TMP x1      = convert_to_mp<T2>::eval(x);
            TMP x2      = convert_to_mp<T2>::eval(y);
            TMP res     = Func::eval(x1, x2);
            out[i]      = res;
        };
    };

    double t = toc();
    return t;
};

template<class T, class T2, class Func>
double test_twofold::test_function_bin_base(int size, int n_rep, const T2* in1, const T2* in2, T2* out)
{
    tic();

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; ++i)
        {
            out[i].value = Func::eval(in1[i].value, in2[i].value);
        };
    };

    double t = toc();
    return t;
};

template<class T2, class Func>
double test_twofold::test_function_bin_twofold(int size, int n_rep, const T2* in1, const T2* in2, T2* out)
{
    tic();

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; ++i)
        {
            T2 x1       = in1[i];
            T2 x2       = in2[i];
            T2 res      = Func::eval(x1, x2);
            out[i]      = res;
        };
    };

    double t = toc();
    return t;
};

//
template<class T, class TMP, class Func>
double test_twofold::test_function_fma_gen(int size, int n_rep, const T* in1, 
                                        const T* in2, const T* in3, TMP* out)
{
    tic();

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; ++i)
        {
            T x         = in1[i];
            T y         = in2[i];
            T z         = in3[i];

            TMP x1      = convert_to_mp<T>::eval(x);
            TMP x2      = convert_to_mp<T>::eval(y);
            TMP x3      = convert_to_mp<T>::eval(z);
            TMP res     = Func::eval(x1, x2, x3);

            out[i]      = res;
        };
    };

    double t = toc();
    return t;
};

template<class T, class Func>
double test_twofold::test_function_fma_base(int size, int n_rep, const T* in1, 
                                        const T* in2, const T* in3, T* out)
{
    tic();

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; ++i)
        {
            T res   = Func::eval(in1[i], in2[i], in3[i]);
            out[i]  = res;
        };
    };

    double t = toc();
    return t;
};

template<class T, class Func>
double test_twofold::test_function_fma_twofold(int size, int n_rep, const T* in1, 
                                               const T* in2, const T* in3, T* out)
{
    tic();

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; ++i)
        {
            T x1        = in1[i];
            T x2        = in2[i];
            T x3        = in3[i];
            T res       = Func::eval_twofold(x1, x2, x3);
            out[i]      = res;
        };
    };

    double t = toc();
    return t;
};

template<class T2, class TMP>
bool test_twofold::test_equal(int size, const T2* res, const TMP* res_gen, double max_dist, double& dist)
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

template<class T2>
bool test_twofold::test_equal(int size, const T2* res, const T2* res_gen, double max_dist, double& dist)
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

template<class T, class TMP>
bool test_twofold::test_equal_fma(int size, const T* res, const TMP* res_gen, double max_dist, double& dist)
{
    bool eq         = true;
    dist            = 0.0;

    double loc_dist;

    for (int i = 0; i < size; ++i)
    {
        bool tmp    = test_equal_fma(res[i], res_gen[i], max_dist, loc_dist);
        eq          = eq && tmp;
        dist        = std::max(dist, loc_dist);
    };

    return eq;
}

template<class Func, class T2>
bool test_twofold::test_constraints_1(int size, const T2* res, const T2* in)
{
    bool ok         = true;

    for (int i = 0; i < size; ++i)
    {
        bool tmp    = test_constraints_elem_1<Func, T2>(res[i], in[i]);
        ok          = ok && tmp;
    };

    return ok;
}

template<class Func, class T2>
bool test_twofold::test_constraints_2(int size, const T2* res, const T2* in1, const T2* in2)
{
    bool ok         = true;

    for (int i = 0; i < size; ++i)
    {
        bool tmp    = test_constraints_elem_2<Func, T2>(res[i], in1[i], in2[i]);
        ok          = ok && tmp;
    };

    return ok;
}

template<class T>
struct required_precision_test{};

template<>
struct required_precision_test<twofold>
{
    static int eval()
    {
        return 2 * 52 + 1;
    }
};

template<class T>
struct required_precision_test_fma{};

template<>
struct required_precision_test_fma<double>
{
    static int eval()
    {
        return (int)precision::precision_double();
    }
};

template<>
struct required_precision_test_fma<float>
{
    static int eval()
    {
        int prec    = (int)precision::precision_float();
        return prec;
    }
};

template<class T2, class TMP>
bool test_twofold::test_equal(const T2& res, const TMP& res_gen, double max_dist, double& dist)
{
    dist = 0.0;

    precision req_p = precision(required_precision_test<T2>::eval());
    TMP res_mp      = convert_to_mp<T2>::eval(res);

    if (res_mp == res_gen)
        return true;

    // twofold does not handle INV/NAN values correctly, when res is not finite
    // then true result can still be finite
    if (is_finite(res_mp) == false)
        return true;
    
    dist = ulp_distance(res_gen, res_mp, req_p).cast_float();

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

template<class T2>
bool test_twofold::test_equal(const T2& res, const T2& res_gen, double max_dist, double& dist)
{
    dist = 0.0;

    precision req_p = precision(required_precision_test<T2>::eval());

    if (res.value == res_gen.value && res.error == res_gen.error)
        return true;

    // twofold does not handle INV/NAN values correctly, when res is not finite
    // then true result can still be finite
    if (res.is_finite() == false)
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

template<class T, class TMP>
bool test_twofold::test_equal_fma(const T& res, const TMP& res_gen, double max_dist, double& dist)
{
    dist = 0.0;

    precision req_p = precision(required_precision_test_fma<T>::eval());
    TMP res_mp      = convert_to_mp<T>::eval(res);

    if (res_mp == res_gen)
        return true;

    // twofold does not handle INV/NAN values correctly, when res is not finite
    // then true result can still be finite
    if (is_finite(res_mp) == false)
        return true;
    
    dist = ulp_distance(res_gen, res_mp, req_p).cast_float();

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

template<class Func, class T2>
bool test_twofold::test_constraints_elem_1(const T2& res, const T2& in1)
{
    {
        bool ok = is_normalized(res);

        if (ok == false)
        {   
            is_normalized(res);
            return false;
        };
    }

    (void)in1;
    return true;
}

template<class Func, class T2>
bool test_twofold::test_constraints_elem_2(const T2& res, const T2& in1, const T2& in2)
{
    {
        bool ok = is_normalized(res);

        if (ok == false)
        {
            is_normalized(res);
            return false;
        }
    }

    (void)in1;
    (void)in2;
    return true;
}

template<class T2>
bool test_twofold::is_normalized(const T2& res)
{
    // NAN/INF are not handled correctly, nonfinite values
    // need not be normalized
    if (matcl::is_finite(res.value) == false)
        return true;

    if (matcl::is_finite(res.error) == false)
        return true;

    double eps  = matcl::eps(res.value)/2;
    bool ok     = abs(res.error) <= eps;

    return ok;
};

int test_twofold::get_N() const
{
    #ifdef _DEBUG
        return 10000;
    #else
        return 1000000;
    #endif
};

int test_twofold::get_M() const
{
    #ifdef _DEBUG
        return 1;
    #else
        return 100;
    #endif
}

template<class T>
std::string test_twofold::get_header() const
{
    std::string ret;

    if (m_normalized)
        ret = "normalized ";
    else
        ret = "not normalized ";

    ret += typeid(T).name();

    return ret;
};

}};
