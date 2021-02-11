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

#pragma once

#include "matcl-blas-ext/lapack_ext/blas_ext.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-lapack/blas/details/blas_simd_complex.h"
#include "matcl-blas-lapack/level1/level1.h"

#include <float.h>
#include <cstring>

#ifndef INLINE_TYPE
#define INLINE_TYPE inline
#endif

namespace matcl { namespace lapack
{

//--------------------------------------------------------------------
//                      SET
//--------------------------------------------------------------------

// Compute Y(1:n) = alpha

template<class V> BLAS_EXT_EXPORT INLINE_TYPE
typename details::enable_if_valid<void,V>::type
lapack::set(i_type n, V alpha, V* y, i_type incy)
{    
    if (n == 0)
        return;

    if (incy == 1)
    {
        if (alpha == V(0.0))
        {
            memset(y,0,n*sizeof(V));
            return;
        };

        for (i_type i = 0; i < n; ++i)
            y[i] = alpha;

        return;
    };

    for (i_type i = 0; i < n; ++i)
        y[i*incy] = alpha;
};

//--------------------------------------------------------------------
//                      YAX
//--------------------------------------------------------------------
// Compute Y = alpha * X

template<class V> BLAS_EXT_EXPORT INLINE_TYPE
typename details::enable_if_valid<void,V>::type
lapack::yax(i_type n, V alpha, const V* x, i_type incx, V* y, i_type incy)
{
    if ( n == 0)
        return;

    if ( incx == 0 )
    {   
        V ax = alpha * x[0];

        return level1::set_val<V,0>::eval(y, incy, n, ax );
    }
    
    if ( alpha == V(0.0) )
        return level1::set_val<V,0>::eval(y, incy, n, V(0.0) );

    if ( alpha == V(1.0) )
        return level1::copy<V,V,0>::eval(y, incy, x, incx, n);

    level1::ax<V,V,V,0>::eval(y,incy,x,incx,n,alpha);

    return;
};

#if 0


//=======================   YAX_CONJ    =====================================================
// Compute Y = alpha * conj(X)

template<class V> BLAS_EXT_EXPORT INLINE_TYPE
typename details::enable_if_valid<void,V>::type
yax_conj(i_type n, V alpha, const V* x,i_type incx, V* y, i_type incy)
{
    //integer     n           Last entry to update in Y
    //FLOAT       alpha       Scalar multiplier
    //FLOAT       x(*)        Vector multiplicand
    //integer     incx        Increment for X vector
    //FLOAT       y(*)        Result vector
    //integer     incy        Increment for Y vector

    V ax;
    V one = 1.;
    V zero = 0.;

    if ( n == 0)
        return;

    if ( incx == 0 )
    {   
        ax = alpha * conj(x[0]);
        lapack::set2(n,ax,y,incy);
        return;
    }
    if ( alpha == zero )
    {
        lapack::set2(n,zero,y,incy);
        return;
    };
    if ( alpha == one )
    {
        lapack::copy(n,x,incx,y,incy);
        lapack::lacgv(n,y,incy);
        return;
    };
    if ( incx == 1 && incy == 1 )
    {   
        i_type m = n%8;
        if (m != 0)
        {
            for (i_type i = 0; i < m; ++i)
                y[i] = alpha * conj(x[i]);

            y += m;
            x += m;
        };

        if (n < 8) return;

        for (i_type i = m; i < n; i += 8, x += 8, y += 8)
        {
            y[0] = alpha * conj(x[0]);
            y[1] = alpha * conj(x[1]);
            y[2] = alpha * conj(x[2]);
            y[3] = alpha * conj(x[3]);
            y[4] = alpha * conj(x[4]);
            y[5] = alpha * conj(x[5]);
            y[6] = alpha * conj(x[6]);
            y[7] = alpha * conj(x[7]);
        };
    }
    else if ( incx == 1 )
    {   
        for (i_type i = 0; i < n; ++i, ++x, y+=incy)
            *y = alpha * conj(*x);
    }
    else if ( incy == 1 )
    {   
        for (i_type i = 0; i < n; ++i, x+=incx, ++y)
            *y = alpha * conj(*x);
    }
    else
    {   
        for (i_type i = 0; i < n; ++i, x+=incx, y+=incx)
            *y = alpha * conj(*x);
    };

    return;
};

BLAS_EXT_EXPORT INLINE_TYPE
void syax_conj(const i_type* n, const s_type* alpha,const s_type* x,const i_type* incx, 
                             s_type* y, const i_type* incy )
{
    return lapack::yax_conj(*n,*alpha,x,*incx,y,*incy);
};
BLAS_EXT_EXPORT INLINE_TYPE
void dyax_conj(const i_type* n, const d_type* alpha,const d_type* x,const i_type* incx, 
                             d_type* y, const i_type* incy )
{
    return lapack::yax_conj(*n,*alpha,x,*incx,y,*incy);
};
BLAS_EXT_EXPORT INLINE_TYPE
void cyax_conj(const i_type* n, const c_type* alpha,const c_type* x,const i_type* incx, 
                             c_type* y, const i_type* incy )
{
    return lapack::yax_conj(*n,*alpha,x,*incx,y,*incy);
};
BLAS_EXT_EXPORT INLINE_TYPE
void zyax_conj(const i_type* n, const z_type* alpha,const z_type* x,const i_type* incx, 
                             z_type* y, const i_type* incy )
{
    return lapack::yax_conj(*n,*alpha,x,*incx,y,*incy);
};

#endif
};};

#undef INLINE_TYPE
