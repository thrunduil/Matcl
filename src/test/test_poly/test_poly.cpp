/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2021
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

#pragma once

#include "test_poly.h"
#include "matcl-simd/poly/poly_eval.h"
#include "matcl-simd/poly/poly_eval_twofold.h"

#include "matcl-scalar/lib_functions/utils.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/lib_functions/scalar_gen.h"

// we are calling horner for mp_float type
#pragma warning(disable :4714) // __forceinline not inlined

#include <iomanip>

namespace matcl { namespace test
{

template<class T>
struct get_type_name{};

template<>
struct get_type_name<float>
{
    static std::string eval(int size)
    {
        std::ostringstream os;
        os << "F" << size;
        return os.str();
    }
};

template<>
struct get_type_name<double>
{
    static std::string eval(int size)
    {
        std::ostringstream os;
        os << "D" << size;
        return os.str();
    }
};

template<class T, int Bits, class Tag>
struct get_type_name<matcl::simd::simd<T,Bits, Tag>>
{
    static std::string eval(int size)
    {
        std::ostringstream os;
        os << get_type_name<T>::eval(size);
        os << " ";
        os << Bits;
        
        return os.str();
    }
};

template<class T>
T test::rand_scalar()
{
    return T(matcl::randn());
}

template<>
float test::rand_scalar<float>()
{
    // avoid overflows
    return matcl::frandn()/2.0f;
}

template<class T, class TS, class Func>
struct test_function
{
    static double eval(int N, int M, const T* x, T* out)
    {
        volatile double dum = 0.0;
        static const int vec_size   = TS::vector_size;

        tic();

        for (int j = 0; j < M; ++j)
        {
            for (int i = 0; i < N; i += vec_size)
            {
                TS xs   = TS::load(x + i, std::false_type());
                TS res  = Func::eval(xs);

                res.store(out + i, std::false_type());
            };

            dum = dum + out[0];
        };

        double t1 = toc();
        return t1;
    };
};

template<class T, class Func>
struct test_function<T, T, Func>
{
    template<class TR>
    static double eval(int N, int M, const T* x, TR* out)
    {
        volatile double dum = 0.0;

        tic();

        for (int j = 0; j < M; ++j)
        {
            for (int i = 0; i < N; ++i)
            {
                TR res  = Func::eval(x[i]);
                out[i]  = res;
            };

            dum = dum + T(0);
        };

        double t1 = toc();
        return t1;
    };
};

void test::test_poly()
{
    out_stream << "\n" << "testing polynomials with compile time size" << "\n";

    formatted_disp dm;    

    dm.set_row_label("type",    align_type::right, 10);
    dm.add_column("horn 1",     align_type::left, 5);
    dm.add_column("horn 2",     align_type::left, 5);
    dm.add_column("horn comp",  align_type::left, 5);
    dm.add_column("estrin",     align_type::left, 5);
    dm.add_column("ulp",        align_type::left, 8);
    dm.add_column("status",     align_type::left, 5);

    dm.disp_header();

    test_polynomials().make<double>(dm);
    test_polynomials().make<float>(dm);    
}

void test::test_poly_cond()
{
    formatted_disp dm;

    out_stream << "\n" << "testing poly condition numbers" << "\n";

    dm.set_row_label("type",    align_type::right, 10);
    dm.add_column("base",       align_type::left, 5);
    dm.add_column("horn cond",  align_type::left, 5);
    dm.add_column("comp cond",  align_type::left, 5);
    dm.add_column("ulp",        align_type::left, 8);
    dm.add_column("status",     align_type::left, 5);

    dm.disp_header();

    test_polynomials().make_cond<double>(dm);
    test_polynomials().make_cond<float>(dm);    
}

void test::test_poly_dyn(bool pow2)
{
    if (pow2 == true)
        disp("\nruntime size, (power of two)");
    else
        disp("\nruntime size, (general)");

    formatted_disp dm;

    dm.set_row_label("type",    align_type::right, 10);
    dm.add_column("horn dyn",   align_type::left, 5);
    dm.add_column("estrin dyn",  align_type::left, 5);
    dm.add_column("ulp",        align_type::left, 8);
    dm.add_column("status",     align_type::left, 5);

    dm.disp_header();

    if (pow2 == true)
    {
        test_polynomials().make_dynamic_pow2<double>(dm);
        test_polynomials().make_dynamic_pow2<float>(dm);    
    }
    else
    {
        test_polynomials().make_dynamic<double>(dm);
        test_polynomials().make_dynamic<float>(dm);    
    }
}

template<class T>
void test_polynomials::make_dynamic(formatted_disp& dm)
{
    using TS1   = T;
    using TS2   = matcl::simd::simd<T, 128, simd::sse_tag>;
    using TS3   = matcl::simd::simd<T, 256, simd::avx_tag>;

    int start   = 8;
    int step    = 31;

    #ifdef _DEBUG
        int end = 200;
    #else
        int end = 2000;
    #endif

    for (int i = start; i < end; i += step)
        make_dynamic<T, TS1>(dm, i);

    for (int i = start; i < end; i += step)
        make_dynamic<T, TS2>(dm, i);

    for (int i = start; i < end; i += step)
        make_dynamic<T, TS3>(dm, i);
};

template<class T>
void test_polynomials::make_dynamic_pow2(formatted_disp& dm)
{
    using TS1   = T;
    using TS2   = matcl::simd::simd<T, 128, simd::sse_tag>;
    using TS3   = matcl::simd::simd<T, 256, simd::avx_tag>;

    int start   = 8;
    int step    = 2;

    #ifdef _DEBUG
        int end = 500;
    #else
        int end = 5000;
    #endif

    for (int i = start; i < end; i *= step)
        make_dynamic<T, TS1>(dm, i);

    for (int i = start; i < end; i *= step)
        make_dynamic<T, TS2>(dm, i);

    for (int i = start; i < end; i *= step)
        make_dynamic<T, TS3>(dm, i);
};

template<class T>
void test_polynomials::make(formatted_disp& dm)
{
    using TS1   = T;
    using TS2   = matcl::simd::simd<T, 128, simd::sse_tag>;
    using TS3   = matcl::simd::simd<T, 256, simd::avx_tag>;

    make2<T, TS1, 2>(dm);
    make2<T, TS2, 2>(dm);
    make2<T, TS3, 2>(dm);

    make2<T, TS1, 3>(dm);
    make2<T, TS2, 3>(dm);
    make2<T, TS3, 3>(dm);

    make2<T, TS1, 4>(dm);
    make2<T, TS2, 4>(dm);
    make2<T, TS3, 4>(dm);

    make2<T, TS1, 5>(dm);
    make2<T, TS2, 5>(dm);
    make2<T, TS3, 5>(dm);

    make2<T, TS1, 6>(dm);
    make2<T, TS2, 6>(dm);
    make2<T, TS3, 6>(dm);

    make2<T, TS1, 7>(dm);
    make2<T, TS2, 7>(dm);
    make2<T, TS3, 7>(dm);

    make2<T, TS1, 8>(dm);
    make2<T, TS2, 8>(dm);
    make2<T, TS3, 8>(dm);

    make2<T, TS1, 9>(dm);
    make2<T, TS2, 9>(dm);
    make2<T, TS3, 9>(dm);

    make2<T, TS1, 10>(dm);
    make2<T, TS2, 10>(dm);
    make2<T, TS3, 10>(dm);

    make2<T, TS1, 11>(dm);
    make2<T, TS2, 11>(dm);
    make2<T, TS3, 11>(dm);

    make2<T, TS1, 12>(dm);
    make2<T, TS2, 12>(dm);
    make2<T, TS3, 12>(dm);

    make2<T, TS1, 13>(dm);
    make2<T, TS2, 13>(dm);
    make2<T, TS3, 13>(dm);

    make2<T, TS1, 14>(dm);
    make2<T, TS2, 14>(dm);
    make2<T, TS3, 14>(dm);

    make2<T, TS1, 15>(dm);
    make2<T, TS2, 15>(dm);
    make2<T, TS3, 15>(dm);

    make2<T, TS1, 16>(dm);
    make2<T, TS2, 16>(dm);
    make2<T, TS3, 16>(dm);
    
    make2<T, TS1, 32>(dm);
    make2<T, TS2, 32>(dm);
    make2<T, TS3, 32>(dm);

    make2<T, TS1, 64>(dm);
    make2<T, TS2, 64>(dm);
    make2<T, TS3, 64>(dm);
};

template<class T>
void test_polynomials::make_cond(formatted_disp& dm)
{
    using TS1   = T;
    using TS2   = matcl::simd::simd<T, 128, simd::sse_tag>;
    using TS3   = matcl::simd::simd<T, 256, simd::avx_tag>;

    int start   = 8;
    int step    = 2;

    #ifdef _DEBUG
        int end = 500;
    #else
        int end = 5000;
    #endif

    for (int i = start; i < end; i *= step)
        make_cond2<T, TS1>(dm, i);

    for (int i = start; i < end; i *= step)
        make_cond2<T, TS2>(dm, i);

    for (int i = start; i < end; i *= step)
        make_cond2<T, TS3>(dm, i);

};

template<class T, class TS, int Size>
void test_polynomials::make2(formatted_disp& dm)
{   
    int N           = 160;

    #ifdef _DEBUG
        int M       = 1000/Size;
    #else
        int M       = 1000*1000*100/N/Size;
    #endif

    using TMP       = mp_float;

    std::vector<T>  x;
    std::vector<TMP> out_1;
    std::vector<T> out_2;

    x.resize(N);
    out_1.resize(N);
    out_2.resize(N);

    T* ptr_x        = x.data();
    TMP* ptr_out_1  = out_1.data();
    T* ptr_out_2    = out_2.data();

    for (int i = 0; i < N; ++i)
        ptr_x[i]    = rand_scalar<T>();

    using Poly  = random_polynomial<T, Size>;
    Poly::rand();

    test_function<T,T, Func_horn_mp<Poly>>::eval(N, 1, ptr_x, ptr_out_1);

    double gamma    = 1.0;
    double gamma2   = 0.5;

    double t2   = test_function<T,TS, Func_horn_1<T, Poly>>
                    ::eval(N, M, ptr_x, ptr_out_2);

    double d2;
    bool ok2    = test_equal_cond<T, Poly>(N, ptr_x, ptr_out_2, ptr_out_1, gamma, d2);

    double t3   = test_function<T,TS, Func_horn_2<T, Poly>>
                    ::eval(N, M, ptr_x, ptr_out_2);
    double d3;
    bool ok3    = test_equal_cond<T, Poly>(N, ptr_x, ptr_out_2, ptr_out_1, gamma, d3);

    double t4   = test_function<T,TS, Func_horn_twofold<T, Poly>>
                        ::eval(N, M, ptr_x, ptr_out_2);
    double d4;
    bool ok4    = test_equal_ulp<T, Poly>(N, ptr_x, ptr_out_2, ptr_out_1, gamma2, d4);    

    double t5   = test_function<T,TS, Func_estrin<T, Poly>>
                    ::eval(N, M, ptr_x, ptr_out_2);

    double d5;
    bool ok5    = test_equal_cond<T, Poly>(N, ptr_x, ptr_out_2, ptr_out_1, 20 * gamma, d5);

    std::string type    = get_type_name<TS>::eval(Size);

    bool ok             = ok2 && ok3 && ok4 && ok5;
    double err          = std::max(d2, d3);
    err                 = std::max(err, d4);
    err                 = std::max(err, d5);

    std::string status  = (ok == true) ? "OK" : "FAIL"; 
    dm.disp_row(type, t2, t3/t2, t4/t2, t5/t2, err, status);

    if (ok == false)
    {
        std::cout << d2 << " " << d3 << " " << d4 << " " << d5 << "\n";
    }
};

template<class T, class TS>
void test_polynomials::make_cond2(formatted_disp& dm, int size)
{   
    int N           = 160;

    #ifdef _DEBUG
        int M       = 1;
    #else
        int M       = 1000*1000*100/N/size;
    #endif

    using TMP       = mp_float;

    std::vector<T>  x;
    std::vector<TMP> out_1;
    std::vector<T> out_2;

    x.resize(N);
    out_1.resize(N);
    out_2.resize(N);

    T* ptr_x        = x.data();
    TMP* ptr_out_1  = out_1.data();
    T* ptr_out_2    = out_2.data();

    for (int i = 0; i < N; ++i)
        ptr_x[i]    = rand_scalar<T>();

    static const int max_size   = 5000;
    using Poly  = random_polynomial_dyn<T, max_size>;
    Poly::rand(size);

    test_function<T,T, Func_horn_mp<Poly>>::eval(N, 1, ptr_x, ptr_out_1);

    double t2   = test_function<T,TS, Func_horn_2<T, Poly>>
                    ::eval(N, M, ptr_x, ptr_out_2);

    double t3   = test_function<T,TS, Func_horn_cond<T, Poly>>
                    ::eval(N, M, ptr_x, ptr_out_2);

    double t4   = test_function<T,TS, Func_horn_comp_cond<T, Poly>>
                        ::eval(N, M, ptr_x, ptr_out_2);
    double d4;
    bool ok4    = test_equal_cond_comp<T, Poly>(N, ptr_x, ptr_out_1, d4);

    std::string type    = get_type_name<TS>::eval(size);

    bool ok             = ok4;
    double err          = d4;

    std::string status  = (ok == true) ? "OK" : "FAIL"; 
    dm.disp_row(type, t2, t3/t2, t4/t2, err, status);

    if (ok == false)
    {
        std::cout << d4 << "\n";
    }
};

template<class T, class TS>
void test_polynomials::make_dynamic(formatted_disp& dm, int size)
{   
    int N           = 160;

    #ifdef _DEBUG
        int M       = 1;
    #else
        int M       = 1000*1000*100/N/size;
    #endif

    using TMP       = mp_float;

    std::vector<T>  x;
    std::vector<TMP> out_1;
    std::vector<T> out_2;

    x.resize(N);
    out_1.resize(N);
    out_2.resize(N);

    T* ptr_x        = x.data();
    TMP* ptr_out_1  = out_1.data();
    T* ptr_out_2    = out_2.data();

    for (int i = 0; i < N; ++i)
        ptr_x[i]    = rand_scalar<T>() / T(2.0);

    static const int max_size   = 5000;
    using Poly  = random_polynomial_dyn<T, max_size>;
    Poly::rand(size);

    if (size > max_size)
        throw std::runtime_error("size too big");

    test_function<T,T, Func_horn_mp<Poly>>::eval(N, 1, ptr_x, ptr_out_1);

    double max_dist = 1.0;

    double t3   = test_function<T,TS, Func_horn_2<T, Poly>>
                    ::eval(N, M, ptr_x, ptr_out_2);
    double d3;
    bool ok3    = test_equal_cond<T, Poly>(N, ptr_x, ptr_out_2, ptr_out_1, max_dist, d3);

    double t6   = test_function<T,TS, Func_estrin2<T, Poly>>
                    ::eval(N, M, ptr_x, ptr_out_2);
    double d6;
    bool ok6    = test_equal_cond<T, Poly>(N, ptr_x, ptr_out_2, ptr_out_1, 10.0 * max_dist, d6);

    double d7;
    bool ok7    = test_equal_cond_comp<T, Poly>(N, ptr_x, ptr_out_1, d7);

    std::string type    = get_type_name<TS>::eval(size);

    bool ok             = ok3 && ok6 && ok7;
    double err          = std::max(d3, d6);

    std::string status  = (ok == true) ? "OK" : "FAIL"; 
    dm.disp_row(type, t3, t6/t3, err, status);

    if (ok == false)
    {
        std::cout << d3 << " " << d6 << "\n";
    }
};

template<class T, class Poly>
bool test_polynomials::test_equal_cond(int size, const T* x, const T* res1, const mp_float* res2,
                                  double max_err, double& max_error)
{
    bool ok     = true;
    max_error   = 0.0;

    for (int i = 0; i < size; ++i)
    {
        T cond1     = simd::horner_apriori_cond(x[i], Poly::get_size(), Poly::polynomial)
                    * T(2) * Poly::get_size();

        T cond2     = simd::horner_aposteriori_cond(x[i], Poly::get_size(), Poly::polynomial);
        /*
        T mu;
        T res       = simd::horner_and_error(x[i], Poly::get_size(), Poly::polynomial, mu);
        T u         = constants::eps<T>()/T(2);
        mu          = mu/u;
        T cond2     = mu / T(abs(res2[i].cast_float()));
        
        (void)res;
        */

        if (cond1 < cond2 && std::isfinite(cond2))
        {  
            std::cout << "wrong cond number " << cond1 << " " << cond2 << "\n";
        }

        T cond      = cond2;

        double err;
        bool ok_loc = test_equal_cond(res1[i], res2[i], cond, max_err, err);

        ok          = ok && ok_loc;
        max_error   = std::max(max_error, err);
    };

    return ok;
}

template<class T, class Poly>
bool test_polynomials::test_equal_cond_comp(int size, const T* x, const mp_float* res2,
                                  double& max_error)
{
    bool ok     = true;
    max_error   = 0.0;

    for (int i = 0; i < size; ++i)
    {
        T poly_err;
        bool exact;
        T res       = simd::compensated_horner_and_error(x[i], Poly::get_size(), 
                        Poly::polynomial, poly_err, exact);

        double err;
        bool ok_loc = test_equal_cond_comp(res, res2[i], poly_err, err);

        ok          = ok && ok_loc;
        max_error   = std::max(max_error, err);
    };

    return ok;
}

template<class T, class Poly>
bool test_polynomials::test_equal_ulp(int size, const T* x, const T* res1, const mp_float* res2,
                                  double max_err, double& max_error)
{
    (void)x;
    bool ok     = true;
    max_error   = 0.0;

    for (int i = 0; i < size; ++i)
    {
        double err;
        bool ok_loc = test_equal_ulp(res1[i], res2[i], max_err, err);

        ok          = ok && ok_loc;
        max_error   = std::max(max_error, err);
    };

    return ok;
}

template<class T>
bool test_polynomials::test_equal_cond(T res1, const mp_float& res2, T cond, double max_err, double& error)
{
    if (is_finite(res1) == false)
        return true;

    double dif  = abs(res1 - res2).cast_float() / abs(res2.cast_float());
    dif         = dif / cond;
    double eps  = constants::eps<T>()/2.0;
    double rel  = dif/eps;

    bool ok = (rel <= max_err);
    error   = rel;
    return ok;
};

template<class T>
struct get_precision{};

template<>
struct get_precision<float>
{
    static precision eval()
    {
        return precision::precision_float();
    };
};

template<>
struct get_precision<double>
{
    static precision eval()
    {
        return precision::precision_double();
    };
};

template<class T>
bool test_polynomials::test_equal_cond_comp(T res1, const mp_float& res2, T res_err, 
                                            double& error)
{
    if (is_finite(res1) == false)
        return true;

    mp_float dif    = mp_float(res1) - res2;
    bool ok         = abs(dif) <= mp_float(res_err);

    precision prec  = get_precision<T>::eval();
    double dist     = ulp_distance(res2, res1, prec).cast_float();
    error           = dist;

    return ok;
};

template<class T>
bool test_polynomials::test_equal_ulp(T res1, const mp_float& res2, double max_err, double& error)
{
    if (is_finite(res1) == false)
        return true;

    precision prec  = get_precision<T>::eval();
    double dist = ulp_distance(res2, res1, prec).cast_float();
    error       = dist;

    bool ok = (dist <= max_err);
    return ok;

};

}}
