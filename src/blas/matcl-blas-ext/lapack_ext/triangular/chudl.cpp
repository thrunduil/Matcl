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

template<class V, class VR>
static bool get_scaling(const VR& SIGMA, const VR& dj, const V& pj, VR& djbar, VR& C, V& S0, V& S1)
{
    if (pj == V(0.))
    {
        djbar   = dj;
        S0      = V(0.0);
        S1      = V(0.0);
        C       = 1.0;

        return true;
    };

	VR scale    = dj + sqrt(abs(SIGMA)) * abs(pj);
    djbar   	= abs2(dj/scale) + SIGMA * abs2(pj/scale);		
    
    if (djbar <= 0.)
        return false;
    
    djbar  		= scale * sqrt(djbar);
    S0          = pj / djbar;
    S1          = SIGMA * S0;
    C           = dj / djbar;

    return true;
};

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
chudl(V* R, i_type LDR, i_type p, const typename details::real_type<V>::type& SIGMA, V* X, i_type INCX, 
      i_type* INFO)
{
    // this is the Gill, Golub, Murray, Saunders C1 algorithm adapted from LDL to Cholesky
    // P. E. Gill, G. H. Golub, W. Murray, and M. A. Saunders, Methods for modifying matrix
    // factorizations, Math. Comp., 28 (1974), pp. 505–535.

    // let L_N - factor, D - diagonal of factor, w - vector initialized with X
    // original recurrence for Cholesky factors:
    //  w_N[p] 		<- w_O[p] - s0 * L_O[p,j]
    //  L_N[p,j] 	<- c1 * L_O[p,j] + s1 * w_O[p]
    //
    //  s0 			:= w_O[j]/D_O[j]					
    //  c1 			:= D_O[j]/D_N[j];					
    //  s1 			:= sig * (w_O[j] * A^2)/D_N[j]		
    //  
    //  D_N[j]^2 	<- D_O[j]^2 + sing * w_O[j]^2*A^2	
    //  Abar  		<- A * D_O[j] / D_N[j]				
    //
    // we do rescalling: W	= w * A, then
    //
    //  W_N[p] 		<- C * W_O[p] - S0 * L_O[p,j]
    //  L_N[p,j] 	<- C * L_O[p,j] + S1 * W_O[p]
    //  D_N[j]^2 	<- D_O[j]^2 + sing * W_O[j]^2
    //
    //  C  			:= D_O[j] / D_N[j]
    //  S0 			:= W_O[j] / D_N[j]
    //  S1 			:= sig * W_O[j]/D_N[j]

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

    // Quick return if possible
    if ( p == 0 || SIGMA == 0)
        return;
    
    for(i_type j = 0; j < p; ++j)
    {        
        V pj        = X[0];

        if (pj == V(0.))
        {
            R       = R + LDR;
            X       = X + INCX;

            continue;
        };
	
        VR dj       = real(R[j]);

        VR djbar, c;
        V s0, s1;

        bool info   = get_scaling(SIGMA, dj, pj, djbar, c, s0, s1);

        if (info == false)
        {
            *INFO   = j + 1;
            return;
        };

        R[j]        = djbar;
        V* Y        = X + INCX;
        
        for (i_type i = j + 1; i < p; ++i)
        {
            V tmp   = c * Y[0] - s0 * R[i];
            R[i]    = c * R[i] + s1 * Y[0];   
            Y[0]    = tmp;
            Y       += INCX;
        };

        R           = R + LDR;
        X           = X + INCX;
    };

    return;
};

template void BLAS_EXT_EXPORT
chudl<s_type>(s_type* R, i_type LDR, i_type p, const s_type& SIGMA, s_type* X, i_type INCX, i_type* INFO);

template void BLAS_EXT_EXPORT
chudl<d_type>(d_type* R, i_type LDR, i_type p, const d_type& SIGMA, d_type* X, i_type INCX, i_type* INFO);

template void BLAS_EXT_EXPORT
chudl<c_type>(c_type* R, i_type LDR, i_type p, const s_type& SIGMA, c_type* X, i_type INCX, i_type* INFO);

template void BLAS_EXT_EXPORT
chudl<z_type>(z_type* R, i_type LDR, i_type p, const d_type& SIGMA, z_type* X, i_type INCX, i_type* INFO);

BLAS_EXT_EXPORT void schudl(s_type* R, i_type LDR, i_type p, const s_type& SIGMA, s_type* X, i_type INCX, 
                        i_type* INFO)
{
    return chudl<s_type>(R,LDR,p,SIGMA, X, INCX, INFO);
};
BLAS_EXT_EXPORT void dchudl(d_type* R, i_type LDR, i_type p, const d_type& SIGMA, d_type* X, i_type INCX, 
                        i_type* INFO)
{
    return chudl<d_type>(R,LDR,p,SIGMA, X, INCX, INFO);
};
BLAS_EXT_EXPORT void cchudl(c_type* R, i_type LDR, i_type p, const s_type& SIGMA, c_type* X, i_type INCX, 
                        i_type* INFO)
{
    return chudl<c_type>(R,LDR,p,SIGMA, X, INCX, INFO);
};
BLAS_EXT_EXPORT void zchudl(z_type* R, i_type LDR, i_type p, const d_type& SIGMA, z_type* X, i_type INCX, 
                        i_type* INFO)
{
    return chudl<z_type>(R,LDR,p,SIGMA, X, INCX, INFO);
};

};};
