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

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-blas-ext/lapack_ext/blas_ext.h"

#include <algorithm>

namespace matcl { namespace lapack
{

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::gsg_ort(void (*CALLBACK)(void* CTX, i_type N, V* X, V* WORK), void* CTX, 
        const char* TRANS, i_type N, i_type K, const V* Q, i_type LDQ, V* X, i_type INCX, 
        V* s, i_type INCS, typename details::real_type<V>::type& x_norm, V* WORK_B, V* work, i_type& INFO)
{
    using VR = typename details::real_type<V>::type;

    INFO                    = 0;

    //test arguments
    if (N < 0)
        INFO                = -2;
    else if (K < 0 || K > N)
        INFO                = -3;
    else if (LDQ < N)
        INFO                = -5;
    else if (INCX == 0)
        INFO                = -7;
    else if (INCS == 0)
        INFO                = -9;

    if (INFO < 0)
        return;

    //fast exit
    if (N == 0)
    {
        x_norm              = VR(0.0);
        return;
    };    

    i_type max_iter         = 5;
    bool is_complex         = lapack::details::is_complex<V>::value;
    bool no_trans           = TRANS[0] == 'N' || TRANS[0] == 'n';
    i_type INCX_pos         = INCX < 0 ? -INCX : INCX;

    const char* tr_char_1;
    const char* tr_char_2;
    i_type NV, KV;

    if (no_trans)
    {
        tr_char_1           = is_complex ? "C" : "T";
        tr_char_2           = "N";
        NV                  = N;
        KV                  = K;
    }
    else
    {
        tr_char_1           = "N";
        tr_char_2           = is_complex ? "C" : "T";
        NV                  = K;
        KV                  = N;
    };

    // compute norm of starting vector with respect to inner product
    // given by B
    (*CALLBACK)(CTX, N, X, WORK_B);

    VR x_norm0              = abs( dotc(N, X, 1, WORK_B, 1) );
	x_norm0                 = sqrt(x_norm0);

    if (K == 0)
    {
        x_norm              = x_norm0;
        return;
    };    
    
    // Otherwise need to orthogonalize the starting vector against
    // the current basis using Gram-Schmidt with iter. ref
    //       s = Q^{T} * B * x;   x = x - Q * s;
    //
    // Stopping criteria used for iter. ref. is discussed in
    // Parlett's book, page 107 and in Gragg & Reichel TOMS paper    
    {
        lapack::gemv(tr_char_1, NV, KV, V(1.0), Q, LDQ, WORK_B, 1, V(0.0), s, INCS);
        lapack::gemv(tr_char_2, NV, KV, -V(1.0), Q, LDQ, s, INCS, V(1.0), X, INCX);

        // Compute the norm of the orthogonalized starting vector
        (*CALLBACK)(CTX, N, X, WORK_B);

        x_norm      = abs( dotc(N, X, 1, WORK_B, 1) );
	    x_norm      = sqrt(x_norm);

        // Check for further orthogonalization
        if (x_norm > 0.717 * x_norm0) 
            return;

        x_norm0     = x_norm;
    };

    bool succ = false;
    for (i_type iter = 1; iter < max_iter; ++iter)
    {
        lapack::gemv(tr_char_1, NV, KV, V(1.0), Q, LDQ, WORK_B, 1, V(0.0), work, 1);
        lapack::gemv(tr_char_2, NV, KV, -V(1.0), Q, LDQ, work, 1, V(1.0), X, INCX);

        lapack::axpy(K, V(1.0), work, 1, s, INCS);

        // Compute the norm of the orthogonalized starting vector
        (*CALLBACK)(CTX, N, X, WORK_B);

        x_norm      = abs( dotc(N, X, 1, WORK_B, 1) );
	    x_norm      = sqrt(x_norm);

        // Check for further orthogonalization
        if (x_norm > 0.717 * x_norm0) 
        {
            succ    = true;
            break;
        }

        x_norm0     = x_norm;
    };

    if (succ == false)
    {
        // iterative refinement step "failed"
        lapack::set(N, V(0.0), X, INCX_pos);

        x_norm          = VR(0.0);
        INFO            = 1;
    };

    return;
};

template BLAS_EXT_EXPORT
void gsg_ort<d_type>(void (*CALLBACK)(void* CTX, i_type N, d_type* X, d_type* WORK), void* CTX,
        const char* TR, i_type N, i_type K, const d_type* Q, i_type LDQ, d_type* X, 
        i_type incx, d_type* s, i_type incs, d_type& x_norm, d_type* work_b, d_type* work, i_type& info);
template BLAS_EXT_EXPORT
void gsg_ort<s_type>(void (*CALLBACK)(void* CTX, i_type N, s_type* X, s_type* WORK), void* CTX, 
        const char* TR, i_type N, i_type K, const s_type* Q, i_type LDQ, s_type* X, 
        i_type incx, s_type* s, i_type incs, s_type& x_norm, s_type* work_b, s_type* work, i_type& info);
template BLAS_EXT_EXPORT
void gsg_ort<z_type>(void (*CALLBACK)(void* CTX, i_type N, z_type* X, z_type* WORK), void* CTX, 
        const char* TR, i_type N, i_type K, const z_type* Q, i_type LDQ, z_type* X, 
        i_type incx, z_type* s, i_type incs, d_type& x_norm, z_type* work_b, z_type* work, i_type& info);
template BLAS_EXT_EXPORT
void gsg_ort<c_type>(void (*CALLBACK)(void* CTX, i_type N, c_type* X, c_type* WORK), void* CTX, 
        const char* TR, i_type N, i_type K, const c_type* Q, i_type LDQ, c_type* X, 
        i_type incx, c_type* s, i_type incs, s_type& x_norm, c_type* work_b, c_type* work, i_type& info);

}};