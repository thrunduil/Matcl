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

#include "test_simd_config.h"
#include "test_simd_compl.h"
#include "test_simd.h"
#include "utils.h"

#include "matcl-scalar/lib_functions/utils.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-scalar/IO/scalar_io.h"
#include "matcl-simd/simd.h"
#include "matcl-simd/simd_complex.h"
#include "test_functions.h"

#include <vector>

namespace matcl { namespace test
{

namespace ms = matcl::simd;

void test::test_performance_complex()
{
    test_simd_compl(false).make_binary();
    test_simd_compl(false).make_unary();
};

void test::test_values_complex()
{    
    test_simd_compl(true).make_unary();
    test_simd_compl(true).make_binary();
};

test_simd_compl::test_simd_compl(bool test_values)
    :m_instr_tag(MATCL_TEST_SIMD_TAG), m_test_values(test_values)
{};

void test_simd_compl::make_unary()
{
    test_functions<Float_complex>();    
    test_functions<Complex>();    
};

void test_simd_compl::make_binary()
{
    test_functions_bin<Complex>();
    test_functions_bin<Float_complex>();    
};

int test_simd_compl::get_size() const
{
    #ifdef _DEBUG
        int N   = 1000;
    #else
        int N   = 10000;
    #endif

    if (m_test_values == true)
        #ifdef _DEBUG
            return 100000;
        #else
            return 1000000;
        #endif
    else
        return N;
};

int test_simd_compl::get_size_perf() const
{
    return 1000;
}

int test_simd_compl::get_num_rep() const
{
    #ifdef _DEBUG
        int M   = 1000;
    #else
        int M   = 10000;
    #endif

    if (m_test_values == true)
        return M/10;
    else
        return M;
}

template<class T>
void test_simd_compl::test_functions()
{
    int N   = get_size();

    std::vector<T> in;
    std::vector<T> out;
    std::vector<T> out_gen;

    in.resize(N);
    out.resize(N);
    out_gen.resize(N);

    T* ptr_in       = in.data();
    T* ptr_out      = out.data();
    T* ptr_out_gen  = out_gen.data();

    for (int i = 0; i < N; ++i)
        ptr_in[i]   = rand_scalar<T>::make(m_test_values);

    std::string header  = m_instr_tag + " " + typeid(T).name();

    disp(" ");
    disp(header);

    formatted_disp dm;

    dm.set_row_label("func",    align_type::right, 7);
    dm.add_column("nosimd",     align_type::left, 5);
    dm.add_column("128 no",     align_type::left, 5);
    dm.add_column("256 no",     align_type::left, 5);
    dm.add_column("128 sse",    align_type::left, 5);
    dm.add_column("256 sse",    align_type::left, 5);
    dm.add_column("256 avx",    align_type::left, 5);
    dm.add_column("status",     align_type::left, 10);

    dm.disp_header();

    test_function<T, test_functions::Func_uminus>(dm, N, ptr_in, ptr_out, ptr_out_gen, true);
    test_function<T, test_functions::Func_conj>(dm, N, ptr_in, ptr_out, ptr_out_gen, true);

    test_function_block<T, test_functions::Func_any_nan>(dm, N, ptr_in, ptr_out, ptr_out_gen, true);
    test_function_block<T, test_functions::Func_reverse>(dm, N, ptr_in, ptr_out, ptr_out_gen, true);    
    test_function_block<T, test_functions::Func_cast>(dm, N, ptr_in, ptr_out, ptr_out_gen, true);    

    using TR    = typename ms::details::real_type<T>::type;

    // block methods are not exact; use positive values in order to
    // reduce roundoff errors; scaling in order to avoid overflows
    for (int i = 0; i < N; ++i)
        ptr_in[i]   = T(std::abs(real(ptr_in[i])) / TR(8), std::abs(imag(ptr_in[i])) / TR(8));
    
    test_function_block<T, test_functions::Func_hor_sum>(dm, N, ptr_in, ptr_out, ptr_out_gen, true); 
};

template<class T>
void test_simd_compl::test_functions_bin()
{
    using TR    = typename ms::details::real_type<T>::type;

    int N   = get_size();

    std::vector<T> in_1;
    std::vector<T> in_2;
    std::vector<T> out;
    std::vector<T> out_gen;

    in_1.resize(N);
    in_2.resize(N);
    out.resize(N);
    out_gen.resize(N);

    T* ptr_in_1     = in_1.data();
    T* ptr_in_2     = in_2.data();
    T* ptr_out      = out.data();
    T* ptr_out_gen  = out_gen.data();

    for (int i = 0; i < N; ++i)
    {
        ptr_in_1[i] = rand_scalar<T>::make(m_test_values);
        ptr_in_2[i] = rand_scalar<T>::make(m_test_values);
    };

    std::string header  = m_instr_tag + " " + typeid(T).name();

    disp(" ");
    disp(header);

    formatted_disp dm;

    dm.set_row_label("func",    align_type::right, 7);
    dm.add_column("nosimd",     align_type::left, 5);
    dm.add_column("128 no",     align_type::left, 5);
    dm.add_column("256 no",     align_type::left, 5);
    dm.add_column("128 sse",    align_type::left, 5);
    dm.add_column("256 sse",    align_type::left, 5);
    dm.add_column("256 avx",    align_type::left, 5);
    dm.add_column("status",     align_type::left, 10);

    dm.disp_header();

    TR* ptr_in_1r = (TR*)ptr_in_1;
    TR* ptr_in_2r = (TR*)ptr_in_2;

    test_function_bin<T, test_functions::Func_mult>
        (dm, N, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 1.0, true);
    test_function_bin<T, test_functions::Func_div>
        (dm, N, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 10.0, false);
    test_function_bin<T, test_functions::Func_plus>
        (dm, N, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 1.0, true);
    test_function_bin<T, test_functions::Func_minus>
        (dm, N, ptr_in_1, ptr_in_2, ptr_out, ptr_out_gen, 1.0, true);

    test_function_bin_RC<T, test_functions::Func_mult_RC, TR>
        (dm, N, ptr_in_1r, ptr_in_2, ptr_out, ptr_out_gen, 1.0, true);
    test_function_bin_RC<T, test_functions::Func_div_RC, TR>
        (dm, N, ptr_in_1r, ptr_in_2, ptr_out, ptr_out_gen, 10.0, false);
    test_function_bin_CR<T, test_functions::Func_mult_CR, TR>
        (dm, N, ptr_in_1, ptr_in_2r, ptr_out, ptr_out_gen, 1.0, true);
    test_function_bin_CR<T, test_functions::Func_div_CR, TR>
        (dm, N, ptr_in_1, ptr_in_2r, ptr_out, ptr_out_gen, 1.0, true);
};

template<class T, class Func>
void test_simd_compl::test_function(formatted_disp& fd, int size, 
                                     const T* in, T* out, T* out_gen, bool test_componentwise)
{
    using TR    = typename ms::details::real_type<T>::type;

    double t0, t1, t2, t3, t4, t5;
    bool v1, v2, v3, v4, v5;
    double d1, d2, d3, d4, d5;

    t0  = test_function_simd<T, ms::simd_compl<TR, 128, ms::nosimd_tag>, Func>(size, 1, in, out_gen);

    t1  = test_function_simd<T, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>(size, 1, in, out);
    v1  = test_equal(size, out, out_gen, 1.0, d1, test_componentwise);

    t2  = test_function_simd<T, simd::simd_compl<TR, 256, simd::nosimd_tag>, Func>(size, 1, in, out);
    v2  = test_equal(size, out, out_gen, 1.0, d2, test_componentwise);

    t3  = test_function_simd<T, simd::simd_compl<TR, 128, simd::sse_tag>, Func>(size, 1, in, out);    
    v3  = test_equal(size, out, out_gen, 1.0, d3, test_componentwise);

    t4  = test_function_simd<T, simd::simd_compl<TR, 256, simd::sse_tag>, Func>(size, 1, in, out);
    v4  = test_equal(size, out, out_gen, 1.0, d4, test_componentwise);

    #if MATCL_ARCHITECTURE_HAS_AVX
        t5  = test_function_simd<T, simd::simd_compl<TR, 256, simd::avx_tag>, Func>(size, 1, in, out);    
        v5  = test_equal(size, out, out_gen, 1.0, d5, test_componentwise);
    #else
        t5  = t0;
        v5  = true;
        d5  = 0;
    #endif

    bool ok     = v1 && v2 && v3 && v4 && v5;
    double d    = d1;
    d           = std::max(d, d2);
    d           = std::max(d, d3);
    d           = std::max(d, d4);
    d           = std::max(d, d5);
    std::string status = make_status(ok, d);

    int N       = get_size_perf();
    int M       = get_num_rep();

    t0  = test_function_simd<T, ms::simd_compl<TR, 128, ms::nosimd_tag>, Func>
                (N, M, in, out_gen);
    t1  = test_function_simd<T, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                (N, M, in, out);
    t2  = test_function_simd<T, simd::simd_compl<TR, 256, simd::nosimd_tag>, Func>
                (N, M, in, out);
    t3  = test_function_simd<T, simd::simd_compl<TR, 128, simd::sse_tag>, Func>
                (N, M, in, out);    
    t4  = test_function_simd<T, simd::simd_compl<TR, 256, simd::sse_tag>, Func>
                (N, M, in, out);

    #if MATCL_ARCHITECTURE_HAS_AVX
        t5  = test_function_simd<T, simd::simd_compl<TR, 256, simd::avx_tag>, Func>
                (N, M, in, out);    
    #else
        t5  = t0;
    #endif

    fd.disp_row(Func::name(), t0, t0/t1, t0/t2, t0/t3, t0/t4, t0/t5, status);
};

template<class T, class Func>
void test_simd_compl::test_function_block(formatted_disp& fd, int size, 
                                const T* in, T* out, T* out_gen, bool test_componentwise)
{
    using TR    = typename ms::details::real_type<T>::type;

    double t0, t1, t2, t3, t4;
    bool v1, v3, v4;
    double d1, d3, d4;

    t0  = test_function_simd<T, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
            (size, 1, in, out_gen);

    t1  = test_function_simd<T, simd::simd_compl<TR, 128, simd::sse_tag>, Func>
            (size, 1, in, out);    
    v1  = test_equal(size, out, out_gen, 4.0, d1, test_componentwise);

    t2  = test_function_simd<T, simd::simd_compl<TR, 256, simd::nosimd_tag>, Func>
            (size, 1, in, out_gen);

    t3  = test_function_simd<T, simd::simd_compl<TR, 256, simd::sse_tag>, Func>
            (size, 1, in, out);
    v3  = test_equal(size, out, out_gen, 4.0, d3, test_componentwise);

    #if MATCL_ARCHITECTURE_HAS_AVX
        t4  = test_function_simd<T, simd::simd_compl<TR, 256, simd::avx_tag>, Func>
                (size, 1, in, out);    
        v4  = test_equal(size, out, out_gen, 4.0, d4, test_componentwise);
    #else
        t4  = t0;
        v4  = true;
        d4  = 0;
    #endif

    bool ok = v1 && v3 && v4;
    double d    = d1;
    d           = std::max(d, d3);
    d           = std::max(d, d4);
    std::string status = make_status(ok, d);

    int N       = get_size_perf();
    int M       = get_num_rep();

    t0  = test_function_simd<T, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                (N, M, in, out_gen);
    t1  = test_function_simd<T, simd::simd_compl<TR, 128, simd::sse_tag>, Func>
                (N, M, in, out);    
    t2  = test_function_simd<T, simd::simd_compl<TR, 256, simd::nosimd_tag>, Func>
                (N, M, in, out_gen);
    t3  = test_function_simd<T, simd::simd_compl<TR, 256, simd::sse_tag>, Func>
                (N, M, in, out);

    #if MATCL_ARCHITECTURE_HAS_AVX
        t4  = test_function_simd<T, simd::simd_compl<TR, 256, simd::avx_tag>, Func>
                (N, M, in, out);    
    #else
        t4  = t0;
    #endif

    fd.disp_row(Func::name(), t0, t0/t0, t2/t2, t0/t1, t2/t3, t2/t4, status);
};

template<class T, class Func>
void test_simd_compl::test_function_bin(formatted_disp& fd, int size, 
                                    const T* in_1, const T* in_2, T* out, T* out_gen, 
                                    double max_dist, bool test_componentwise)
{
    using TR    = typename ms::details::real_type<T>::type;

    double t0, t1, t2, t3, t4, t5;
    bool v1, v2, v3, v4, v5;
    double d1, d2, d3, d4, d5;

    t0  = test_function_bin_simd<T, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                                    (size, 1, in_1, in_2, out_gen);

    t1  = test_function_bin_simd<T, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                                    (size, 1, in_1, in_2, out);
    v1  = test_equal(size, out, out_gen, max_dist, d1, test_componentwise);

    t2  = test_function_bin_simd<T, simd::simd_compl<TR, 256, simd::nosimd_tag>, Func>
                                    (size, 1, in_1, in_2, out);
    v2  = test_equal(size, out, out_gen, max_dist, d2, test_componentwise);

    t3  = test_function_bin_simd<T, simd::simd_compl<TR, 128, simd::sse_tag>, Func>
                                    (size, 1, in_1, in_2, out);    
    v3  = test_equal(size, out, out_gen, max_dist, d3, test_componentwise);

    t4  = test_function_bin_simd<T, simd::simd_compl<TR, 256, simd::sse_tag>, Func>
                                    (size, 1, in_1, in_2, out);
    v4  = test_equal(size, out, out_gen,max_dist, d4, test_componentwise);

    #if MATCL_ARCHITECTURE_HAS_AVX
        t5  = test_function_bin_simd<T, simd::simd_compl<TR, 256, simd::avx_tag>, Func>
                                    (size, 1, in_1, in_2, out);    
        v5  = test_equal(size, out, out_gen, max_dist, d5, test_componentwise);
    #else
        t5  = t0;
        v5  = true;
        d5  = 0;
    #endif

    bool ok     = v1 && v2 && v3 && v4 && v5;
    double d    = d1;
    d           = std::max(d, d2);
    d           = std::max(d, d3);
    d           = std::max(d, d4);
    d           = std::max(d, d5);

    int N       = get_size_perf();
    int M       = get_num_rep();

    t0  = test_function_bin_simd<T, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                                    (N, M, in_1, in_2, out_gen);
    t1  = test_function_bin_simd<T, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                                    (N, M, in_1, in_2, out);
    t2  = test_function_bin_simd<T, simd::simd_compl<TR, 256, simd::nosimd_tag>, Func>
                                    (N, M, in_1, in_2, out);
    t3  = test_function_bin_simd<T, simd::simd_compl<TR, 128, simd::sse_tag>, Func>
                                    (N, M, in_1, in_2, out);    
    t4  = test_function_bin_simd<T, simd::simd_compl<TR, 256, simd::sse_tag>, Func>
                                    (N, M, in_1, in_2, out);

    #if MATCL_ARCHITECTURE_HAS_AVX
        t5  = test_function_bin_simd<T, simd::simd_compl<TR, 256, simd::avx_tag>, Func>
                                    (N, M, in_1, in_2, out);    
    #else
        t5  = t0;
    #endif

    std::string status = make_status(ok, d);
    fd.disp_row(Func::name(), t0, t0/t1, t0/t2, t0/t3, t0/t4, t0/t5, status);
};

template<class T, class Func, class TR>
void test_simd_compl::test_function_bin_RC(formatted_disp& fd, int size, 
                                    const TR* in_1, const T* in_2, T* out, T* out_gen, 
                                    double max_dist, bool test_componentwise)
{
    double t0, t1, t2, t3, t4, t5;
    bool v1, v2, v3, v4, v5;
    double d1, d2, d3, d4, d5;

    t0  = test_function_bin_simd_RC<T, TR, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                                    (size, 1, in_1, in_2, out_gen);

    t1  = test_function_bin_simd_RC<T, TR, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                                    (size, 1, in_1, in_2, out);
    v1  = test_equal(size, out, out_gen, max_dist, d1, test_componentwise);

    t2  = test_function_bin_simd_RC<T, TR, simd::simd_compl<TR, 256, simd::nosimd_tag>, Func>
                                    (size, 1, in_1, in_2, out);
    v2  = test_equal(size, out, out_gen, max_dist, d2, test_componentwise);

    t3  = test_function_bin_simd_RC<T, TR, simd::simd_compl<TR, 128, simd::sse_tag>, Func>
                                    (size, 1, in_1, in_2, out);    
    v3  = test_equal(size, out, out_gen, max_dist, d3, test_componentwise);

    t4  = test_function_bin_simd_RC<T, TR, simd::simd_compl<TR, 256, simd::sse_tag>, Func>
                                    (size, 1, in_1, in_2, out);
    v4  = test_equal(size, out, out_gen,max_dist, d4, test_componentwise);

    #if MATCL_ARCHITECTURE_HAS_AVX
        t5  = test_function_bin_simd_RC<T, TR, simd::simd_compl<TR, 256, simd::avx_tag>, Func>
                                    (size, 1, in_1, in_2, out);    
        v5  = test_equal(size, out, out_gen, max_dist, d5, test_componentwise);
    #else
        t5  = t0;
        v5  = true;
        d5  = 0;
    #endif

    bool ok     = v1 && v2 && v3 && v4 && v5;
    double d    = d1;
    d           = std::max(d, d2);
    d           = std::max(d, d3);
    d           = std::max(d, d4);
    d           = std::max(d, d5);

    int N       = get_size_perf();
    int M       = get_num_rep();

    t0  = test_function_bin_simd_RC<T, TR, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                                    (N, M, in_1, in_2, out_gen);
    t1  = test_function_bin_simd_RC<T, TR, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                                    (N, M, in_1, in_2, out);
    t2  = test_function_bin_simd_RC<T, TR, simd::simd_compl<TR, 256, simd::nosimd_tag>, Func>
                                    (N, M, in_1, in_2, out);
    t3  = test_function_bin_simd_RC<T, TR, simd::simd_compl<TR, 128, simd::sse_tag>, Func>
                                    (N, M, in_1, in_2, out);    
    t4  = test_function_bin_simd_RC<T, TR, simd::simd_compl<TR, 256, simd::sse_tag>, Func>
                                    (N, M, in_1, in_2, out);

    #if MATCL_ARCHITECTURE_HAS_AVX
        t5  = test_function_bin_simd_RC<T, TR, simd::simd_compl<TR, 256, simd::avx_tag>, Func>
                                    (N, M, in_1, in_2, out);    
    #else
        t5  = t0;
    #endif

    std::string status = make_status(ok, d);
    fd.disp_row(Func::name(), t0, t0/t1, t0/t2, t0/t3, t0/t4, t0/t5, status);
};

template<class T, class Func, class TR>
void test_simd_compl::test_function_bin_CR(formatted_disp& fd, int size,
                                    const T* in_1, const TR* in_2, T* out, T* out_gen, 
                                    double max_dist, bool test_componentwise)
{
    double t0, t1, t2, t3, t4, t5;
    bool v1, v2, v3, v4, v5;
    double d1, d2, d3, d4, d5;

    t0  = test_function_bin_simd_CR<T, TR, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                                    (size, 1, in_1, in_2, out_gen);

    t1  = test_function_bin_simd_CR<T, TR, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                                    (size, 1, in_1, in_2, out);
    v1  = test_equal(size, out, out_gen, max_dist, d1, test_componentwise);

    t2  = test_function_bin_simd_CR<T, TR, simd::simd_compl<TR, 256, simd::nosimd_tag>, Func>
                                    (size, 1, in_1, in_2, out);
    v2  = test_equal(size, out, out_gen, max_dist, d2, test_componentwise);

    t3  = test_function_bin_simd_CR<T, TR, simd::simd_compl<TR, 128, simd::sse_tag>, Func>
                                    (size, 1, in_1, in_2, out);    
    v3  = test_equal(size, out, out_gen, max_dist, d3, test_componentwise);

    t4  = test_function_bin_simd_CR<T, TR, simd::simd_compl<TR, 256, simd::sse_tag>, Func>
                                    (size, 1, in_1, in_2, out);
    v4  = test_equal(size, out, out_gen,max_dist, d4, test_componentwise);

    #if MATCL_ARCHITECTURE_HAS_AVX
        t5  = test_function_bin_simd_CR<T, TR, simd::simd_compl<TR, 256, simd::avx_tag>, Func>
                                    (size, 1, in_1, in_2, out);    
        v5  = test_equal(size, out, out_gen, max_dist, d5, test_componentwise);
    #else
        t5  = t0;
        v5  = true;
        d5  = 0;
    #endif

    bool ok     = v1 && v2 && v3 && v4 && v5;
    double d    = d1;
    d           = std::max(d, d2);
    d           = std::max(d, d3);
    d           = std::max(d, d4);
    d           = std::max(d, d5);

    int N       = get_size_perf();
    int M       = get_num_rep();

    t0  = test_function_bin_simd_CR<T, TR, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                                    (N, M, in_1, in_2, out_gen);
    t1  = test_function_bin_simd_CR<T, TR, simd::simd_compl<TR, 128, simd::nosimd_tag>, Func>
                                    (N, M, in_1, in_2, out);
    t2  = test_function_bin_simd_CR<T, TR, simd::simd_compl<TR, 256, simd::nosimd_tag>, Func>
                                    (N, M, in_1, in_2, out);
    t3  = test_function_bin_simd_CR<T, TR, simd::simd_compl<TR, 128, simd::sse_tag>, Func>
                                    (N, M, in_1, in_2, out);    
    t4  = test_function_bin_simd_CR<T, TR, simd::simd_compl<TR, 256, simd::sse_tag>, Func>
                                    (N, M, in_1, in_2, out);

    #if MATCL_ARCHITECTURE_HAS_AVX
        t5  = test_function_bin_simd_CR<T, TR, simd::simd_compl<TR, 256, simd::avx_tag>, Func>
                                    (N, M, in_1, in_2, out);    
    #else
        t5  = t0;
    #endif

    std::string status = make_status(ok, d);
    fd.disp_row(Func::name(), t0, t0/t1, t0/t2, t0/t3, t0/t4, t0/t5, status);
};

template<class T, class Simd_type, class Func>
double test_simd_compl::test_function_simd(int size, int n_rep, const T* in, T* out)
{
    static const int vec_size   = Simd_type::vector_size;

    tic();

    using TR        = typename ms::details::real_type<T>::type;
    volatile TR val = 0;

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; i += vec_size)
        {
            Simd_type x     = Simd_type::load(in + i, std::false_type());
            Simd_type res   = Func::eval(x);

            res.store(out + i, std::false_type());
        };

        val = real(out[0]);
    };

    double t = toc();
    return t;
};

template<class T, class Simd_type, class Func>
double test_simd_compl::test_function_bin_simd(int size, int n_rep, const T* in1, const T* in2, T* out)
{
    static const int vec_size   = Simd_type::vector_size;

    tic();

    using TR        = typename ms::details::real_type<T>::type;
    volatile TR val = 0;

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i = 0; i < size; i += vec_size)
        {
            Simd_type x1    = Simd_type::load(in1 + i, std::false_type());
            Simd_type x2    = Simd_type::load(in2 + i, std::false_type());

            Simd_type res   = Func::eval(x1, x2);

            res.store(out + i, std::false_type());
        };

        val += real(out[0]);
    };

    double t = toc();
    return t;
};

template<class T, class TR, class Simd_type, class Func>
double test_simd_compl::test_function_bin_simd_RC(int size, int n_rep, 
                                          const TR* in1, const T* in2, T* out)
{
    using Simd_type_C           = Simd_type;
    using Simd_type_R           = typename ms::details::real_type<Simd_type>::type;
    static const int vec_size   = Simd_type_C::vector_size;

    tic();

    volatile TR val = 0;

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i1 = 0, i2 = 0; i1 < size; i1 += vec_size, i2 += 2*vec_size)
        {
            Simd_type_R x1    = Simd_type_R::load(in1 + i2, std::false_type());
            Simd_type_C x2    = Simd_type_C::load(in2 + i1, std::false_type());

            Simd_type_C res   = Func::eval(x1, x2);

            res.store(out + i1, std::false_type());
        };

        val += real(out[0]);
    };

    double t = toc();
    return t;
};

template<class T, class TR, class Simd_type, class Func>
double test_simd_compl::test_function_bin_simd_CR(int size, int n_rep, 
                                          const T* in1, const TR* in2, T* out)
{
    using Simd_type_C           = Simd_type;
    using Simd_type_R           = typename ms::details::real_type<Simd_type>::type;
    static const int vec_size   = Simd_type_C::vector_size;

    tic();

    volatile TR val = 0;

    for(int j = 0; j < n_rep; ++j)
    {
        for (int i1 = 0, i2 = 0; i1 < size; i1 += vec_size, i2 += 2*vec_size)
        {
            Simd_type_C x1    = Simd_type_C::load(in1 + i1, std::false_type());
            Simd_type_R x2    = Simd_type_R::load(in2 + i2, std::false_type());

            Simd_type_C res   = Func::eval(x1, x2);

            res.store(out + i1, std::false_type());
        };

        val += real(out[0]);
    };

    double t = toc();
    return t;
};

template<class T>
bool test_simd_compl::test_equal(const T& res, const T& res_gen, double max_dist, double& dist,
                                        bool componentwise)
{
    dist = 0.0;

    if (eeq_nan(res, res_gen) == true)
        return true;

    using TR    = typename ms::details::real_type<T>::type;

    if (componentwise == true)
    {        
        double d1, d2;
        bool b1     = test_equal_real(real(res), real(res_gen), max_dist, d1);
        bool b2     = test_equal_real(imag(res), imag(res_gen), max_dist, d2);

        dist    = std::max(d1, d2);

        return b1 && b2;
    };

    TR abs_dif  = abs(res - res_gen);
    TR abs_base = abs(res_gen);

    if (is_nan(abs_base) != is_nan(abs_dif))
        return false;

    if (is_nan(abs_base) == true)
    {
        if (is_nan(real(res)) != is_nan(real(res_gen)))
            return false;

        if (is_nan(imag(res)) != is_nan(imag(res_gen)))
            return false;

        return true;
    };

    dist    = abs_dif / eps(abs_base);

    if (dist > max_dist)
        return false;
    else
        return true;
}

template<class T>
bool test_simd_compl::test_equal_real(const T& res, const T& res_gen, double max_dist, double& dist)
{
    dist = 0.0;

    if (res == res_gen)
        return true;

    if (is_nan(res) && is_nan(res_gen))
        return true;

    dist = float_distance(res, res_gen);

    if (is_nan(dist) == true)
        return false;

    if (dist > max_dist)
        return false;
    else
        return true;
}

template<class T>
bool test_simd_compl::test_equal(int size, const T* res, const T* res_gen, 
                                        double max_dist, double& dist, bool componentwise)
{
    bool eq         = true;
    dist            = 0.0;

    double loc_dist;

    for (int i = 0; i < size; ++i)
    {
        bool tmp    = test_equal(res[i], res_gen[i], max_dist, loc_dist, componentwise);
        eq          = eq && tmp;
        dist        = std::max(dist, loc_dist);
    };

    return eq;
}

std::string test_simd_compl::make_status(bool ok, double ulp_dist)
{
    if (ok == true)
        return "OK";

    std::ostringstream os;
    os << "F" << "(" << ulp_dist << ")";

    return os.str();
}

}};
