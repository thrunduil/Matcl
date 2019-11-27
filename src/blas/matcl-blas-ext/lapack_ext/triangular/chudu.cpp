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
chudu(V* R, i_type LDR, i_type p, const typename details::real_type<V>::type& SIGMA, V* X, i_type INCX, 
    V* WORK, i_type* INFO)
{
    // Test the input parameters.
    *INFO = 0;

    if (LDR < 1 || LDR < p)
    {
        *INFO = -2;
    }
    else if (p < 1)
    {
        *INFO = -3;
    }     
    else if (INCX == 0)
    {
        *INFO = -6;
    }     

    if (*INFO != 0)
        return;

    using VR = typename details::real_type<V>::type;

    VR* work_c  = reinterpret_cast<VR*>(WORK);
    V* work_s   = WORK + p;

    // Quick return if possible
    if ( p == 0 || SIGMA == 0)
        return;

    for(i_type j = 0; j < p; ++j)
    {
        V xj            = *X;

        if (SIGMA == 1.)
        {
            for (i_type i = 0; i < j; ++i)
            {
                V t         = work_c[i] * R[i] + work_s[i] * xj;
                xj          = work_c[i] * xj   - lapack::conj(work_s[i]) * R[i];
                R[i]        = t;
            };
        }
        else if (SIGMA == -1.)
        {
            for (i_type i = 0; i < j; ++i)
            {
                V t         = work_c[i] * R[i] - work_s[i] * xj;
                xj          = work_c[i] * xj   - lapack::conj(work_s[i]) * R[i];
                R[i]        = t;
            };		
        }
        else
        {
            for (i_type i = 0; i < j; ++i)
            {
                V t         = work_c[i] * R[i] + SIGMA * work_s[i] * xj;
                xj          = work_c[i] * xj   - lapack::conj(work_s[i]) * R[i];
                R[i]        = t;
            };		
        };                    
        
        if (xj == V(0.))
        {
            work_c[j]       = 1.;
            work_s[j]       = 0.;
        }
        else
        {
            VR scale        = abs(R[j]) + sqrt(abs(SIGMA)) * abs(xj);
            VR r            = abs2(R[j]/scale) + SIGMA * abs2(xj/scale);

            if (r <= 0)
            {
                *INFO       = j + 1;
                return;
            };
            r               = scale * sqrt(r);

            work_c[j]       = real(R[j]) / r;
            work_s[j]       = xj / r;
            R[j]            = r;
        };

        R               	= R + LDR;
        X               	= X + INCX;
    };

    return;
};

template void BLAS_EXT_EXPORT
chudu<s_type>(s_type* R, i_type LDR, i_type p, const s_type& SIGMA, s_type* X, i_type INCX, s_type* WORK, i_type* INFO);

template void BLAS_EXT_EXPORT
chudu<d_type>(d_type* R, i_type LDR, i_type p, const d_type& SIGMA, d_type* X, i_type INCX, d_type* WORK, i_type* INFO);

template void BLAS_EXT_EXPORT
chudu<c_type>(c_type* R, i_type LDR, i_type p, const s_type& SIGMA, c_type* X, i_type INCX, c_type* WORK, i_type* INFO);

template void BLAS_EXT_EXPORT
chudu<z_type>(z_type* R, i_type LDR, i_type p, const d_type& SIGMA, z_type* X, i_type INCX, z_type* WORK, i_type* INFO);

BLAS_EXT_EXPORT void schudu(s_type* R, i_type LDR, i_type p, const s_type& SIGMA, s_type* X, i_type INCX, 
                        s_type* WORK, i_type* INFO)
{
    return chudu<s_type>(R,LDR,p,SIGMA, X, INCX, WORK, INFO);
};
BLAS_EXT_EXPORT void dchudu(d_type* R, i_type LDR, i_type p, const d_type& SIGMA, d_type* X, i_type INCX, 
                        d_type* WORK, i_type* INFO)
{
    return chudu<d_type>(R,LDR,p,SIGMA, X, INCX, WORK, INFO);
};
BLAS_EXT_EXPORT void cchudu(c_type* R, i_type LDR, i_type p, const s_type& SIGMA, c_type* X, i_type INCX, 
                        c_type* WORK, i_type* INFO)
{
    return chudu<c_type>(R,LDR,p,SIGMA, X, INCX, WORK, INFO);
};
BLAS_EXT_EXPORT void zchudu(z_type* R, i_type LDR, i_type p, const d_type& SIGMA, z_type* X, i_type INCX, 
                        z_type* WORK, i_type* INFO)
{
    return chudu<z_type>(R,LDR,p,SIGMA, X, INCX, WORK, INFO);
};

};};