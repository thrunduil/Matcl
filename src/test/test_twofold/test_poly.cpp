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

#pragma once

#include "test_poly.h"
#include "matcl-simd/poly/poly_eval.h"

#include "matcl-scalar/lib_functions/utils.h"
#include "matcl-scalar/lib_functions/scalar_gen.h"

namespace matcl { namespace test
{

void test::test_poly()
{
    test_polynomials().make();
}

void test_polynomials::make()
{   
    std::cout << "start" << "\n";
    
    int N           = 200;
    int M           = 10000*1000/N;

    std::vector<double> x;
    std::vector<double> out;
    x.resize(N);
    out.resize(N);

    double* ptr_x   = x.data();
    double* ptr_out = out.data();

    for (int i = 0; i < N; ++i)
    {
        ptr_x[i]    = randn() / 3.0;
    };

    double t2   = test_horner(N, M, ptr_x, ptr_out);
    double res2 = ptr_out[0];

    double t1   = test_short_horner(N, M, ptr_x, ptr_out);    
    double res1 = ptr_out[0];

    std::cout << t1 << " " << t2 << " " << t2/t1 << "\n";
    std::cout << res2 << " " << res1 << "\n";
};

double test_polynomials::test_short_horner(int N, int M, const double* x, double* out)
{
    volatile double dum = 0.0;

    tic();

    for (int j = 0; j < M; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            double res  = simd::short_horner(x[i], 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                                                   11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0);
            out[i]  = res;
        };

        dum = dum + out[0];
    };

    double t1 = toc();
    return t1;
};

double test_polynomials::test_horner(int N, int M, const double* x, double* out)
{
    volatile double dum = 0.0;

    double poly[]   = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                       11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0};

    static const int size   = sizeof(poly) /sizeof(poly[0]);

    tic();

    for (int j = 0; j < M; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            double res  = simd::horner<size>(x[i], poly);
            out[i]  = res;
        };

        dum = dum + out[0];
    };

    double t1 = toc();
    return t1;
};

}}
