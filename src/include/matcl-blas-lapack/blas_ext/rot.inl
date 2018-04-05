/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2013-2016
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

#include "matcl-blas-lapack/blas/details/blas_simd_complex.h"
#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-lapack/blas/details/config_blas_lib.h"

#include "matcl-simd/simd.h"
#include "matcl-simd/simd_complex.h"

#include <algorithm>

namespace matcl { namespace lapack
{

template<class Val>
struct rot_blas
{};

template<>
struct rot_blas<s_type>
{
    static void eval(i_type n, s_type *x, i_type incx, s_type *y, i_type incy, s_type C, s_type S)
    {
        BLAS_NAME(srot)(&n,x,&incx,y,&incy,&C,&S);
    }
};
template<>
struct rot_blas<d_type>
{
    static void eval(i_type n, d_type *x, i_type incx, d_type *y, i_type incy, d_type C, d_type S)
    {
        BLAS_NAME(drot)(&n,x,&incx,y,&incy,&C,&S);
    }
};
template<>
struct rot_blas<z_type>
{
    static void eval(i_type n, z_type *x, i_type incx, z_type *y, i_type incy, d_type C, z_type S)
    {
        BLAS_NAME(zrot)(&n,_rc(x),&incx,_rc(y),&incy,&C,_rc(&S));
    }
};
template<>
struct rot_blas<c_type>
{
    static void eval(i_type n, c_type *x, i_type incx, c_type *y, i_type incy, s_type C, c_type S)
    {
        BLAS_NAME(crot)(&n,_rc(x),&incx,_rc(y),&incy,&C,_rc(&S));
    }
};

template<class Val>
struct rot2
{
    static void eval(i_type n, Val *x, i_type incx, Val *y, i_type incy, Val C, Val S)
    {
        if (incx == 1 && incy == 1)
            return eval_col(n, x, y, C, S);

        return rot_blas<Val>::eval(n,x,incx,y,incy,C,S);
    };
    static void eval_col(i_type N, Val* X, Val* Y, Val C, Val S)
    {
        using simd_type = typename simd::default_simd_type<Val>::type;

        static const i_type VS  = simd_type::vector_size;

        if (N == 0)
            return;

        if (N < VS)
            goto lab_tail;
        
        //TODO C might be real for complex V
        simd_type vC        = simd_type::broadcast(&C);
        simd_type vS        = simd_type::broadcast(&S);

        i_type size_1       = N / VS;

        for (i_type i = 0; i < size_1; ++i)
        {
            simd_type x     = simd_type::load(X + i * VS, std::false_type());
            simd_type y     = simd_type::load(Y + i * VS, std::false_type());

            //simd_type x2  = vC * x  + vS * y;
            //simd_type y2  = vC * y  - vS * x;
            simd_type x2    = fma_f(vC, x, vS * y);
            simd_type y2    = fms_f(vC, y, vS * x);

            x2.store(X + i * VS, std::false_type());
            y2.store(Y + i * VS, std::false_type());
        };

        i_type pos          = size_1 * VS;

        if (pos == N)
            return;

        N                   = N - pos;
        X                   = X + pos;
        Y                   = Y + pos;

      lab_tail:

        for (i_type i = 0; i < N; ++i)
        {
            Val x2          = C * X[i] + S * Y[i];
            Val y2          = C * Y[i] - S * X[i];

            X[i]            = x2;
            Y[i]            = y2;
        };

        return;
    };
};

template<class Val, bool Use_simd = true>
struct rot_2
{
    using VR = typename details::real_type<Val>::type;

    static void eval(i_type N, Val *x_11, Val *x_12_21, Val *x_22, VR C1, Val S1, VR C2, Val S2)
    {
        using simd_type = typename simd::default_simd_type<Val>::type;
        using simd_re   = typename simd::default_simd_type<VR>::type;

        static const i_type VS  = simd_type::vector_size;

        if (N == 0)
            return;

        if (N < VS)
            goto lab_tail;
    
        simd_re vC1         = simd_re::broadcast(&C1);
        simd_re vC2         = simd_re::broadcast(&C2);
        simd_type vS1       = simd_type::broadcast(&S1);
        simd_type vS2       = simd_type::broadcast(&S2);

        i_type size_1       = N / VS;

        for (i_type i = 0; i < size_1; ++i)
        {
            simd_type x1    = simd_type::load(x_11 + i * VS, std::false_type());
            simd_type x2    = simd_type::load(x_12_21 + i * VS, std::false_type());
            
            simd_type y1    = vC1 * x1  + vS1 * x2;
            y1.store(x_11 + i * VS, std::false_type());

            x1              = vC1 * x2  - vS1 * x1;
            x2              = simd_type::load(x_22 + i * VS, std::false_type());

            y1              = vC2 * x1  + vS2 * x2;
            simd_type y2    = vC2 * x2  - vS2 * x1;

            y1.store(x_12_21 + i * VS, std::false_type());
            y2.store(x_22 + i * VS, std::false_type());
        };

        i_type pos          = size_1 * VS;

        if (pos == N)
            return;

        N                   = N - pos;
        x_11                = x_11 + pos;
        x_12_21             = x_12_21 + pos;
        x_22                = x_22 + pos;

      lab_tail:

        for (i_type i = 0; i < N; ++i)
        {
            Val x1          = x_11[i];
            Val x2          = x_12_21[i];
            Val y1          = C1 * x1 + S1 * x2;
            Val y2          = C1 * x2 - S1 * x1;

            x_11[i]         = y1;

            x1              = y2;
            x2              = x_22[i];
            y1              = C2 * x1 + S2 * x2;
            y2              = C2 * x2 - S2 * x1;

            x_12_21[i]      = y1;
            x_22[i]         = y2;
        };

        return;
    };
};

template<class Val>
struct rot_2<Val, false>
{
    using VR = typename details::real_type<Val>::type;

    static void eval(i_type n, Val *x_11, Val *x_12_21, Val *x_22, VR C1, Val S1, VR C2, Val S2)
    {
        rot_blas<Val>::eval(n, x_11, 1, x_12_21, 1, C1, S1);
        rot_blas<Val>::eval(n, x_12_21, 1, x_22, 1, C2, S2);
    };
};

};};