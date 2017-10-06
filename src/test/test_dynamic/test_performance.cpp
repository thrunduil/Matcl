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
#include "matcl-scalar/IO/formatted_disp.h"

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
  #ifdef _DEBUG
    Integer M       = 50;
    Integer K       = 50;
    Integer N       = 50;    
  #else
    Integer M       = 100;
    Integer K       = 100;
    Integer N       = 100;    
  #endif

    bool int_only   = false;

    formatted_disp dm;

    dm.set_row_label("type", align_type::right, 14);
    dm.add_column("T", align_type::left, 7);
    dm.add_column("obj<T>", align_type::left, 7);
    dm.add_column("obj", align_type::left, 7);
    dm.add_column("obj<T>/T", align_type::left, 6);
    dm.add_column("obj/obj<T>", align_type::left, 6);
    dm.add_column("obj/T", align_type::left, 6);

    dm.disp_header();

    //TODO
    for (int i = 0; i < 200; ++i)
    {
        Real res;
        OReal ores;
        Integer prec        = 53;

        double time         = test_mult<Real>(res, M, K, N, prec);
        double otime        = test_mult<OReal>(ores, M, K, N, prec);
        double otime2       = test_mult_obj<OReal>(ores, M, K, N, prec);

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;
        dm.disp_row("Real", time, otime, otime2, rel1, rel2, rel3);
    }

    return;
    {
        Integer res;
        OInteger ores;
        Integer prec        = 53;

        double time         = test_mult<Integer>(res, M, K, N, prec);
        double otime        = test_mult<OInteger>(ores, M, K, N, prec);
        double otime2       = test_mult_obj<OInteger>(ores, M, K, N, prec);

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;
        dm.disp_row("Integer", time, otime, otime2, rel1, rel2, rel3);
    };

    if (int_only == true)
        return;

    {
        Float res;
        OFloat ores;
        Integer prec        = 53;

        double time         = test_mult<Float>(res, M, K, N, prec);
        double otime        = test_mult<OFloat>(ores, M, K, N, prec);
        double otime2       = test_mult_obj<OFloat>(ores, M, K, N, prec);

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;
        dm.disp_row("Float", time, otime, otime2, rel1, rel2, rel3);
    }

    {
        Real res;
        OReal ores;
        Integer prec        = 53;

        double time         = test_mult<Real>(res, M, K, N, prec);
        double otime        = test_mult<OReal>(ores, M, K, N, prec);
        double otime2       = test_mult_obj<OReal>(ores, M, K, N, prec);

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;
        dm.disp_row("Real", time, otime, otime2, rel1, rel2, rel3);
    }

    {
        Float_complex res;
        OFloat_complex ores;
        Integer prec        = 53;

        double time         = test_mult<Float_complex>(res, M, K, N, prec);
        double otime        = test_mult<OFloat_complex>(ores, M, K, N, prec);
        double otime2       = test_mult_obj<OFloat_complex>(ores, M, K, N, prec);

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;
        dm.disp_row("Float_complex", time, otime, otime2, rel1, rel2, rel3);
    }

    {
        Complex res;
        OComplex ores;
        Integer prec        = 53;

        double time         = test_mult<Complex>(res, M, K, N, prec);
        double otime        = test_mult<OComplex>(ores, M, K, N, prec);
        double otime2       = test_mult_obj<OComplex>(ores, M, K, N, prec);

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;
        dm.disp_row("Complex", time, otime, otime2, rel1, rel2, rel3);
    };

    {
        mp_int res;
        MP_int ores;
        Integer prec        = 53;

        double time         = test_mult<mp_int>(res, M, K, N, prec);
        double otime        = test_mult<MP_int>(ores, M, K, N, prec);
        double otime2       = test_mult_obj<MP_int>(ores, M, K, N, prec);

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;
        dm.disp_row("mp_int", time, otime, otime2, rel1, rel2, rel3);
    };

    {
        mp_rational res;
        MP_rational ores;
        Integer prec        = 53;

        double time         = test_mult<mp_rational>(res, M, K, N, prec);
        double otime        = test_mult<MP_rational>(ores, M, K, N, prec);
        double otime2       = test_mult_obj<MP_rational>(ores, M, K, N, prec);

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;
        dm.disp_row("mp_rational", time, otime, otime2, rel1, rel2, rel3);
    };

    {
        mp_float res;
        MP_float ores;
        Integer prec        = 53;

        double time         = test_mult<mp_float>(res, M, K, N, prec);
        double otime        = test_mult<MP_float>(ores, M, K, N, prec);    
        double otime2       = test_mult_obj<MP_float>(ores, M, K, N, prec);    

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;

        std::ostringstream type;
        type << "mp_float (" << res.get_precision() << ")";
        dm.disp_row(type.str(), time, otime, otime2, rel1, rel2, rel3);
    };
    {
        mp_float res;
        MP_float ores;
        Integer prec        = 530;

        double time         = test_mult<mp_float>(res, M, K, N, prec);
        double otime        = test_mult<MP_float>(ores, M, K, N, prec);    
        double otime2       = test_mult_obj<MP_float>(ores, M, K, N, prec);    

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;

        std::ostringstream type;
        type << "mp_float (" << res.get_precision() << ")";
        dm.disp_row(type.str(), time, otime, otime2, rel1, rel2, rel3);
    };
    {
        mp_float res;
        MP_float ores;
        Integer prec        = 5300;

        double time         = test_mult<mp_float>(res, M, K, N, prec);
        double otime        = test_mult<MP_float>(ores, M, K, N, prec);    
        double otime2       = test_mult_obj<MP_float>(ores, M, K, N, prec);    

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;

        std::ostringstream type;
        type << "mp_float (" << res.get_precision() << ")";
        dm.disp_row(type.str(), time, otime, otime2, rel1, rel2, rel3);
    };
    {
        mp_complex res;
        MP_complex ores;
        Integer prec        = 53;

        double time         = test_mult<mp_complex>(res, M, K, N, prec);            
        double otime        = test_mult<MP_complex>(ores, M, K, N, prec);
        double otime2       = test_mult_obj<MP_complex>(ores, M, K, N, prec);

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;

        std::ostringstream type;
        type << "mp_compl (" << res.get_precision() << ")";
        dm.disp_row(type.str(), time, otime, otime2, rel1, rel2, rel3);
    }
    {
        mp_complex res;
        MP_complex ores;
        Integer prec        = 530;

        double time         = test_mult<mp_complex>(res, M, K, N, prec);
        double otime        = test_mult<MP_complex>(ores, M, K, N, prec);
        double otime2       = test_mult_obj<MP_complex>(ores, M, K, N, prec);

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;

        std::ostringstream type;
        type << "mp_compl (" << res.get_precision() << ")";
        dm.disp_row(type.str(), time, otime, otime2, rel1, rel2, rel3);
    }
    {
        mp_complex res;
        MP_complex ores;
        Integer prec        = 5300;

        double time         = test_mult<mp_complex>(res, M, K, N, prec);
        double otime        = test_mult<MP_complex>(ores, M, K, N, prec);
        double otime2       = test_mult_obj<MP_complex>(ores, M, K, N, prec);

        double rel1         = otime/time;
        double rel2         = otime2/otime;
        double rel3         = otime2/time;

        std::ostringstream type;
        type << "mp_compl (" << res.get_precision() << ")";
        dm.disp_row(type.str(), time, otime, otime2, rel1, rel2, rel3);
    }
};

void performance_tester::disp_res(const std::string& type, double time, double otime, double otime2)
{
    std::cout << type << ": " << "raw: " << time << ", obj: " << otime 
              << ", obj2:" << otime2 << ", rel: " << otime/time << " " << otime2/otime << " " << otime2/time << "\n";
};

}};
