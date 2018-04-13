/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"


namespace matcl { namespace lapack
{

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
laswpc(i_type n, V *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx)
{
    i_type ix0, i1, i2, inc;

    if ( incx > 0)
    {
         ix0    = k1;
         i1     = k1;
         i2     = k2;
         inc    = 1;
    }
    else if ( incx < 0 )
    {
         ix0    = 1 + ( 1-k2 )*incx;
         i1     = k2;
         i2     = k1;
         inc    = -1;
    }
    else
    {
        return;
    };

    i_type n32  = ( n / 32 )*32;

    if (n32 != 0)
    {
        for (i_type j = 0; j < n32; j += 32)
        {
            i_type ix   = ix0;   

            if (inc > 0)
            {
                for (i_type i = i1; i <= i2; i += inc)
                {
                    i_type ip   = ipiv[ix - 1];

                    if (ip != i)
                    {
                        V* A_i      = a + (i-1) * lda;
                        V* A_ip     = a + (ip-1) * lda;

                        for (i_type k = j; k < j + 32; ++k)
                        {
                            V tmp       = A_i[k];
                            A_i[k]      = A_ip[k];
                            A_ip[k]     = tmp;
                        };
                    };

                    ix = ix + incx;
                };
            }
            else
            {
                for (i_type i = i1; i >= i2; i += inc)
                {
                    i_type ip   = ipiv[ix - 1];

                    if (ip != i)
                    {
                        V* A_i      = a + (i-1) * lda;
                        V* A_ip     = a + (ip-1) * lda;

                        for (i_type k = j; k < j + 32; ++k)
                        {
                            V tmp       = A_i[k];
                            A_i[k]      = A_ip[k];
                            A_ip[k]     = tmp;
                        };
                    };

                    ix = ix + incx;
                };
            };
        };
    };

    if ( n32 != n)
    {
         n32        = n32;
         i_type ix  = ix0;

         if (inc > 0)
         {
             for (i_type i = i1; i <= i2; i+= inc)
             {
                i_type ip   = ipiv[ix - 1];

                if ( ip != i)
                {
                    V* A_i          = a + (i-1) * lda;
                    V* A_ip         = a + (ip-1) * lda;

                    for (i_type k = n32; k < n; ++k)
                    {
                        V tmp       = A_i[k];
                        A_i[k]      = A_ip[k];
                        A_ip[k]     = tmp;
                    };
                };

                ix = ix + incx;
             };
         }
         else
         {
             for (i_type i = i1; i >= i2; i+= inc)
             {
                i_type ip   = ipiv[ix - 1];

                if ( ip != i)
                {
                    V* A_i          = a + (i-1) * lda;
                    V* A_ip         = a + (ip-1) * lda;

                    for (i_type k = n32; k < n; ++k)
                    {
                        V tmp       = A_i[k];
                        A_i[k]      = A_ip[k];
                        A_ip[k]     = tmp;
                    };
                };

                ix = ix + incx;
             };
         };
    };
};

template void BLAS_EXT_EXPORT
laswpc<s_type>(i_type n, s_type *a, i_type lda, i_type k1, i_type k2, const i_type *ipiv, i_type incx);
template void BLAS_EXT_EXPORT
laswpc<d_type>(i_type n, d_type *a, i_type lda, i_type k1, i_type k2, const i_type *ipiv, i_type incx);
template void BLAS_EXT_EXPORT
laswpc<c_type>(i_type n, c_type *a, i_type lda, i_type k1, i_type k2, const i_type *ipiv, i_type incx);
template void BLAS_EXT_EXPORT
laswpc<z_type>(i_type n, z_type *a, i_type lda, i_type k1, i_type k2, const i_type *ipiv, i_type incx);

BLAS_EXT_EXPORT void claswpc(i_type n, c_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx)
{
    laswpc<c_type>(n,a,lda,k1,k2,ipiv,incx);
};
BLAS_EXT_EXPORT void slaswpc(i_type n, s_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx)
{
    laswpc<s_type>(n,a,lda,k1,k2,ipiv,incx);
};
BLAS_EXT_EXPORT void dlaswpc(i_type n, d_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx)
{
    laswpc<d_type>(n,a,lda,k1,k2,ipiv,incx);
};
BLAS_EXT_EXPORT void zlaswpc(i_type n, z_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx)
{
    laswpc<z_type>(n,a,lda,k1,k2,ipiv,incx);
};


};};