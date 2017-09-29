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

#pragma warning(push)
#pragma warning(disable: 4244)  //conversion from 'const int' to 'matcl::Float', possible loss of data

#include "test_performance.h"
#include "matcl-mp/matcl_mp.h"
#include "rand_scalars.h"
#include "matcl-scalar/lib_functions/scalar_gen.h"
#include "matcl-scalar/lib_functions/utils.h"

#pragma warning(pop)

#include <iostream>
#include <set>

namespace matcl { namespace test
{

using OInteger          = dynamic::OInteger;
using OFloat            = dynamic::OFloat;
using OReal             = dynamic::OReal;
using OFloat_complex    = dynamic::OFloat_complex;
using OComplex          = dynamic::OComplex;

void test_performance()
{
    performance_tester test;
    test.make();
};

template<class T, class TR>
void performance_tester::mult(const T* A, const T* B, T* C, Integer M, Integer K, Integer N)
{
    Integer pos = 0;

    for (Integer i = 0; i < N; ++i)
    for (Integer k = 0; k < M; ++k)
        C[pos++] = T(TR());

    Integer A_ld    = M;
    Integer B_ld    = K;
    Integer C_ld    = M;
    const T* A_sav  = A;

    for (Integer i = 0; i < N; ++i)
    {
        A   = A_sav;
    
        for (Integer j = 0; j < K; ++j)
        {
            const T& B_loc  = B[j];

            for (Integer k = 0; k < M; ++k)
                C[k]    = C[k] + A[k] * B_loc;

            A   += A_ld;
        }

        C   += C_ld;
        B   += B_ld;
    };
};

template<class T> struct random;
template<class T> struct random<dynamic::object_type<T>>
{
    static dynamic::object_type<T> eval(Integer prec)
    {
        return dynamic::object_type<T>(random<T>::eval(prec));
    };
};

template<> struct random<Integer>
{
    static Integer eval(Integer)
    {
        Integer max_val = 100;
        Integer val     = matcl::irand();

        return val % max_val;
    };
};
template<> struct random<Float>
{
    static Float eval(Integer)
    {
        return frandn();
    };
};
template<> struct random<Real>
{
    static Real eval(Integer)
    {
        return randn();
    };
};
template<> struct random<Float_complex>
{
    static Float_complex eval(Integer)
    {
        return fcrandn();
    };
};
template<> struct random<Complex>
{
    static Complex eval(Integer)
    {
        return crandn();
    };
};

template<> struct random<mp_int>
{
    static mp_int eval(Integer prec)
    {
        return mp_int(random<Integer>::eval(prec));
    };
};
template<> struct random<mp_float>
{
    static mp_float eval(Integer prec)
    {
        return mp_float(random<Real>::eval(prec), precision(prec));
    };
};
template<> struct random<mp_complex>
{
    static mp_complex eval(Integer prec)
    {
        return mp_complex(random<mp_float>::eval(prec), random<mp_float>::eval(prec));
    };
};
template<> struct random<mp_rational>
{
    static mp_rational eval(Integer prec)
    {
        return mp_rational(random<mp_int>::eval(prec), random<mp_int>::eval(prec));
    };
};

template<class T>
void performance_tester :: init(T* data, Integer size, Integer prec)
{
    for (Integer i = 0; i < size; ++i)
        data[i] = random<T>::eval(prec);
}

template<class T>
double performance_tester::test_mult(T& res, Integer M, Integer K, Integer N, Integer prec)
{
    std::vector<T> data_A;
    std::vector<T> data_B;
    std::vector<T> data_C;

    data_A.resize(M*K);
    data_B.resize(K*N);
    data_C.resize(M*N);

    init(data_A.data(), M*K, prec);
    init(data_B.data(), M*K, prec);

    matcl::tic();
    mult<T>(data_A.data(), data_B.data(), data_C.data(), M, K, N);
    double t = matcl::toc();

    res = data_C[0];
    return t;
};

template<class T>
double performance_tester::test_mult_obj(T& res, Integer M, Integer K, Integer N, Integer prec)
{
    std::vector<T> data_A;
    std::vector<T> data_B;
    std::vector<T> data_C;

    data_A.resize(M*K);
    data_B.resize(K*N);
    data_C.resize(M*N);

    init(data_A.data(), M*K, prec);
    init(data_B.data(), M*K, prec);

    using object = matcl::dynamic::object;
    const object* ptr_A = reinterpret_cast<object*>(data_A.data());
    const object* ptr_B = reinterpret_cast<object*>(data_B.data());
    object* ptr_C       = reinterpret_cast<object*>(data_C.data());

    matcl::tic();
    mult<object, T>(ptr_A, ptr_B, ptr_C, M, K, N);
    double t = matcl::toc();

    res = data_C[0];
    return t;
};

void performance_tester::make()
{
    Integer M       = 100;
    Integer K       = 100;
    Integer N       = 100;
    Integer prec    = 53;

    bool int_only   = false;

    {
        Integer res;
        OInteger ores;
        double time_int     = test_mult<Integer>(res, M, K, N, prec);
        double otime_int    = test_mult<OInteger>(ores, M, K, N, prec);
        double otime_int2   = test_mult_obj<OInteger>(ores, M, K, N, prec);
        disp_res("Integer", time_int, otime_int, otime_int2);
    };

    if (int_only == true)
        return;

    {
        Float res;
        OFloat ores;
        double time_float   = test_mult<Float>(res, M, K, N, prec);
        double otime_float  = test_mult<OFloat>(ores, M, K, N, prec);
        double otime_float2 = test_mult_obj<OFloat>(ores, M, K, N, prec);
        disp_res("Float", time_float, otime_float, otime_float2);
    }

    {
        Real res;
        OReal ores;
        double time_real    = test_mult<Real>(res, M, K, N, prec);
        double otime_real   = test_mult<OReal>(ores, M, K, N, prec);
        double otime_real2  = test_mult_obj<OReal>(ores, M, K, N, prec);
        disp_res("Real", time_real, otime_real, otime_real2);
    }

    {
        Float_complex res;
        OFloat_complex ores;
        double time_fcompl  = test_mult<Float_complex>(res, M, K, N, prec);
        double otime_fcompl = test_mult<OFloat_complex>(ores, M, K, N, prec);
        double otime_fcompl2= test_mult_obj<OFloat_complex>(ores, M, K, N, prec);
        disp_res("Float_complex", time_fcompl, otime_fcompl, otime_fcompl2);
    }

    {
        Complex res;
        OComplex ores;
        double time_compl   = test_mult<Complex>(res, M, K, N, prec);
        double otime_compl  = test_mult<OComplex>(ores, M, K, N, prec);
        double otime_compl2 = test_mult_obj<OComplex>(ores, M, K, N, prec);
        disp_res("Complex", time_compl, otime_compl, otime_compl2);
    };

    {
        mp_int res;
        MP_int ores;
        double time_mpi     = test_mult<mp_int>(res, M, K, N, prec);
        double otime_mpi    = test_mult<MP_int>(ores, M, K, N, prec);
        double otime_mpi2   = test_mult_obj<MP_int>(ores, M, K, N, prec);
        disp_res("mp_int", time_mpi, otime_mpi, otime_mpi2);
    };

    {
        mp_rational res;
        MP_rational ores;
        double time_mpq     = test_mult<mp_rational>(res, M, K, N, prec);
        double otime_mpq    = test_mult<MP_rational>(ores, M, K, N, prec);
        double otime_mpq2   = test_mult_obj<MP_rational>(ores, M, K, N, prec);
        disp_res("mp_rational", time_mpq, otime_mpq, otime_mpq2);
    };

    {
        mp_float res;
        MP_float ores;
        double time_mpf     = test_mult<mp_float>(res, M, K, N, prec);
        double otime_mpf    = test_mult<MP_float>(ores, M, K, N, prec);    
        double otime_mpf2   = test_mult_obj<MP_float>(ores, M, K, N, prec);    
        disp_res("mp_float 1", time_mpf, otime_mpf, otime_mpf2);
        std::cout << "precision: " << res.get_precision() << "\n";
    };
    {
        mp_float res;
        MP_float ores;
        double time_mpf     = test_mult<mp_float>(res, M, K, N, prec*10);
        double otime_mpf    = test_mult<MP_float>(ores, M, K, N, prec*10);    
        double otime_mpf2   = test_mult_obj<MP_float>(ores, M, K, N, prec);    

        disp_res("mp_float 2", time_mpf, otime_mpf, otime_mpf2);
        std::cout << "precision: " << res.get_precision() << "\n";
    };
    {
        mp_float res;
        MP_float ores;
        double time_mpf     = test_mult<mp_float>(res, M, K, N, prec*100);
        double otime_mpf    = test_mult<MP_float>(ores, M, K, N, prec*100);    
        double otime_mpf2   = test_mult_obj<MP_float>(ores, M, K, N, prec);    

        disp_res("mp_float 3", time_mpf, otime_mpf, otime_mpf2);
        std::cout << "precision: " << res.get_precision() << "\n";
    };
    {
        mp_complex res;
        MP_complex ores;
        double time_mpc     = test_mult<mp_complex>(res, M, K, N, prec);            
        double otime_mpc    = test_mult<MP_complex>(ores, M, K, N, prec);
        double otime_mpc2   = test_mult_obj<MP_complex>(ores, M, K, N, prec);

        disp_res("mp_complex 1", time_mpc, otime_mpc, otime_mpc2);
        std::cout << "precision: " << res.get_precision() << "\n";
    }
    {
        mp_complex res;
        MP_complex ores;
        double time_mpc     = test_mult<mp_complex>(res, M, K, N, prec * 10);
        double otime_mpc    = test_mult<MP_complex>(ores, M, K, N, prec * 10);
        double otime_mpc2   = test_mult_obj<MP_complex>(ores, M, K, N, prec);

        disp_res("mp_complex 2", time_mpc, otime_mpc, otime_mpc2);
        std::cout << "precision: " << res.get_precision() << "\n";
    }
    {
        mp_complex res;
        MP_complex ores;
        double time_mpc     = test_mult<mp_complex>(res, M, K, N, prec * 100);
        double otime_mpc    = test_mult<MP_complex>(ores, M, K, N, prec * 100);
        double otime_mpc2   = test_mult_obj<MP_complex>(ores, M, K, N, prec);

        disp_res("mp_complex 3", time_mpc, otime_mpc, otime_mpc2);
        std::cout << "precision: " << res.get_precision() << "\n";
    }
};

void performance_tester::disp_res(const std::string& type, double time, double otime, double otime2)
{
    std::cout << type << ": " << "raw: " << time << ", obj: " << otime 
              << ", obj2:" << otime2 << ", rel: " << otime/time << " " << otime2/otime << " " << otime2/time << "\n";
};

}};
