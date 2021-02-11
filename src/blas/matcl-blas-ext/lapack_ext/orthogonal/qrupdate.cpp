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

#include <algorithm>

namespace matcl { namespace lapack
{

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::qruphr(i_type N, i_type K, V* R, i_type LDR, const V& SIGMA, const V* U, i_type inc_u, 
        const V* W, i_type inc_w, typename details::real_type<V>::type* C, V* S, 
        i_type* IND1, i_type* IND2, i_type& SLEN, i_type& INFO)
{
    using VR    = typename details::real_type<V>::type;

    // test arguments    
    INFO        = 0;
    
    if (N < 0)
        INFO    = -1;
    else if (K < 0)
        INFO    = -2;
    else if (LDR < std::max(N, 1))
        INFO    = -4;
    else if (inc_u == 0)
        INFO    = -7;
    else if (inc_w == 0)
        INFO    = -9;

    if (INFO != 0)
        return;    

    SLEN        = 0;

    // fast exit
    if (SIGMA == V(0.0) || N == 0)
        return;

    VR  c;
    V   s, r;

    // generate plane rotations, that triangularize RR from right

    V W_mod            = conj(W[0]);

    // anihilate V and apply rotations to R and Z; result: upper hessenberg R
    i_type pos         = 0;
    for (i_type i = 0; i < N - 1 - (K + 1); ++i)
    {
        lartg<V>(conj(W[(i+1)*inc_w]), W_mod, &c, &s,  &r);
        
        if (s != V(0.0) )
        {
            C[pos]      = c;
            S[pos]      = -conj(s);
            IND1[pos]   = i + 1;
            IND2[pos]   = i + 2;
            ++pos;
        };
        W_mod           = r;
    };

    SLEN                = pos;
    K                   = std::min(K, N - 1);

    //apply rotations to R
    i_type UDIAGS       = N;
    i_type LDIAGS       = K;
    lapack::rotseq("Right", "Right", "No", SLEN, C, S, IND1, IND2, N, N, LDIAGS, UDIAGS, 
                   R, LDR, INFO);

    // make R + U*V'
    lapack::gerc(N, K + 1, SIGMA, U, inc_u, W + (N - 1 - K)*inc_w, inc_w, R + (N - 1 - K) * LDR, LDR);

    if (K < N - 1 && W_mod != V(0.0))
        lapack::geru(N, 1, SIGMA, U, inc_u, &W_mod, 1, R + (N - 1 - K - 1) * LDR, LDR);

    return;
};

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::qruphl(i_type M, i_type N, i_type K, V* R, i_type LDR, const V& SIGMA, const V* U, i_type inc_u, 
        const V* W, i_type inc_w, typename details::real_type<V>::type* C, V* S, i_type* IND1, 
        i_type* IND2, i_type& SLEN, i_type& INFO)
{
    using VR    = typename details::real_type<V>::type;

    // test arguments
    INFO        = 0;
    
    if (M < 0)
        INFO    = -1;
    else if (N < 0)
        INFO    = -2;
    else if (K < 0)
        INFO    = -3;
    else if (LDR < std::max(M, 1))
        INFO    = -5;
    else if (inc_u == 0)
        INFO    = -8;
    else if (inc_w == 0)
        INFO    = -10;

    SLEN        = 0;

    if (INFO != 0)
        return;    

    // fast exit
    if (SIGMA == V(0.0) || M == 0 || N == 0)
        return;

    VR  c;
    V   s, r;

    // generate plane rotations, that triangularize RR from right
    V U_mod             = U[(M - 1)*inc_u];
    i_type pos          = 0;

    // anihilate V and apply rotations to R and Z; result: upper hessenberg R
    for (i_type i = M - 2; i >= 1 + K; --i, ++pos)
    {
        lartg<V>(U[i*inc_u], U_mod, &c, &s,  &r);
        
        if (s != V(0.0))
        {
            C[pos]      = c;
            S[pos]      = s;
            IND1[pos]   = i + 1;
            IND2[pos]   = i + 2;            
        };

        U_mod           = r;
    };
    
    SLEN                = pos;
    K                   = std::min(K, M - 1);

    //apply rotation to R
    i_type LDIAGS  = K;
    i_type UDIAGS  = N;
    lapack::rotseq("left", "left", "No", SLEN, C, S, IND1, IND2, M, N, LDIAGS, UDIAGS,
                   R, LDR, INFO);

    // make R + U*V'
    lapack::gerc(K + 1, N, SIGMA, U, inc_u, W, inc_w, R, LDR);

    if (K < M - 1 && U_mod != V(0.0))
        lapack::gerc(1, N, SIGMA, &U_mod, inc_u, W, inc_w, R + K + 1, LDR);

    return;
};

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::huundr(i_type N, i_type D, i_type KR, V* R, i_type LDR, typename details::real_type<V>::type* C, 
        V* S, i_type* IND1, i_type* IND2, i_type& SLEN, i_type& INFO)
{
    using VR    = typename details::real_type<V>::type;

    // test arguments
    INFO        = 0;

    if (N < 1)
        INFO    = -1;
    else if (D <= 0 || D > N - 1)
        INFO    = -2;
    else if (KR < 0)
        INFO    = -3;
    else if (LDR < std::max(N+KR, 1))
        INFO    = -5;

    if (INFO != 0)
        return;

    SLEN        = 0;

    // fast exit
    if (N <= 1)
        return;

    V s, r;
    VR c;

    i_type pos          = 0;

    //anihilate subdiagonal
    for (i_type i = N - 1; i >= 1; --i)
    {
        for (i_type d = std::min(D,i); d >= 1; --d)
        {
            lartg(R[i + KR + (i-d+1)*LDR], R[i + KR + (i-d)*LDR], &c, &s, &r);

            R[i+KR+(i-d+1)*LDR] = r;
            R[i+KR+(i-d)*LDR]   = V(0.0);

            if (s == V(0.0))
                continue;

            C[pos]          = c;
            S[pos]          = s;
            IND1[pos]       = (i - d + 1) + 1;
            IND2[pos]       = (i - d) + 1;
            ++pos;

            //apply rotation
            i_type K        = i + KR;
            lapack::rot(K, R + (i-d+1)*LDR, 1, R + (i-d)*LDR, 1, c, s);
        };
    };

    SLEN                = pos;

    return;
};

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::hlundl(i_type N, i_type D, V* R, i_type LDR, typename details::real_type<V>::type* C, 
        V* S, i_type* IND1, i_type* IND2, i_type& SLEN, i_type& INFO)
{
    using VR    = typename details::real_type<V>::type;

    // test arguments
    INFO        = 0;

    if (N < 1)
        INFO    = -1;
    else if (D <= 0 || D > N - 1)
        INFO    = -2;
    else if (LDR < std::max(N, 1))
        INFO    = -4;

    if (INFO != 0)
        return;

    SLEN        = 0;

    // fast exit
    if (N <= 1)
        return;

    V s, r;
    VR c;

    i_type pos          = 0;

    //anihilate subdiagonal
    for (i_type i = N - 1; i >= 1; --i)
    {
        for (i_type d = std::min(D,i); d >= 1; --d)
        {
            lartg(R[i - d + 1 + i*LDR], R[i - d + i*LDR], &c, &s, &r);

            R[i-d+1+i*LDR]  = r;
            R[i-d+i*LDR]    = V(0.0);

            if (s == VR(0.0))
                continue;

            C[pos]          = c;
            S[pos]          = s;
            IND1[pos]       = (i - d + 1) + 1;
            IND2[pos]       = (i - d) + 1;
            ++pos;

            //apply rotation
            i_type K        = i;
            lapack::rot(K, R + i - d + 1, LDR, R + i - d, LDR, c, s);
        };
    };

    SLEN                = pos;
    return;
};

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::huundl(i_type M, i_type N, i_type D, V* R, i_type LDR, typename details::real_type<V>::type* C, 
        V* S, i_type* IND1, i_type* IND2, i_type& SLEN, i_type& INFO)
{
    using VR    = typename details::real_type<V>::type;

    // test arguments
    INFO        = 0;

    if (M < 1)
        INFO    = -1;
    else if (N < 1)
        INFO    = -2;
    else if (D <= 0 || D > M - 1)
        INFO    = -3;
    else if (LDR < std::max(M, 1))
        INFO    = -5;

    if (INFO != 0)
        return;

    SLEN        = 0;

    // fast exit
    if (M <= 1 || N == 0)
        return;

    V s, r;
    VR c;

    i_type pos          = 0;
    i_type i_last       = std::min(M - 1, N - 1);

    //anihilate subdiagonal
    for (i_type i = 0; i < i_last; ++i)
    {
        for (i_type d = std::min(M - i - 1, D); d >= 1; --d)
        {
            lartg(R[i + d - 1 + i*LDR], R[i + d + i*LDR], &c, &s, &r);

            R[i+d-1+i*LDR]  = r;
            R[i+d+i*LDR]    = V(0.0);

            if (s == V(0.0))
                continue;

            C[pos]          = c;
            S[pos]          = s;
            IND1[pos]       = (i + d - 1) + 1;
            IND2[pos]       = (i + d) + 1;
            ++pos;

            //apply rotation
            i_type K        = N - i - 1;
            lapack::rot(K, R + i + d - 1+ (i+1)*LDR, LDR, R + i + d + (i+1)*LDR, LDR, c, s);
        };
    };

    SLEN                = pos;

    return;
};

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::hlundr(i_type M, i_type N, i_type D, V* R, i_type LDR, typename details::real_type<V>::type* C, 
        V* S, i_type* IND1, i_type* IND2, i_type& SLEN, i_type& INFO)
{
    using VR    = typename details::real_type<V>::type;

    // test arguments
    INFO        = 0;

    if (M < 1)
        INFO    = -1;
    else if (N < 1)
        INFO    = -2;
    else if (D <= 0 || D > N - 1)
        INFO    = -3;
    else if (LDR < std::max(M, 1))
        INFO    = -5;

    if (INFO != 0)
        return;

    SLEN        = 0;

    // fast exit
    if (N <= 1 || M == 0)
        return;

    V s, r;
    VR c;

    i_type pos          = 0;
    i_type i_last       = std::min(M - 1, N - 1);

    //anihilate subdiagonal
    for (i_type i = 0; i < i_last; ++i)
    {
        for (i_type d = std::min(N - i - 1, D); d >= 1; --d)
        {
            lartg(R[i + (i + d - 1)*LDR], R[i + (i + d)*LDR], &c, &s, &r);

            R[i+(i+d-1)*LDR]    = r;
            R[i+(i+d)*LDR]      = V(0.0);

            if (s == V(0.0))
                continue;

            C[pos]          = c;
            S[pos]          = s;
            IND1[pos]       = (i + d - 1) + 1;
            IND2[pos]       = (i + d) + 1;
            ++pos;

            //apply rotation
            i_type K        = M - i - 1;
            lapack::rot(K, R + i + 1 + (i+d-1)*LDR, 1, R + i + 1+ (i+d)*LDR, 1, c, s);
        };
    };

    SLEN                = pos;

    return;
};

template void BLAS_EXT_EXPORT
qruphr<d_type>(i_type N, i_type K, d_type* R, i_type LDR, const d_type& SIGMA, const d_type* U, i_type inc_u,
        const d_type* W, i_type inv_w, d_type* C, d_type* S, i_type* IND1, i_type* IND2, i_type& SLEN, 
        i_type& INFO);

template void BLAS_EXT_EXPORT
qruphr<s_type>(i_type N, i_type K, s_type* R, i_type LDR, const s_type& SIGMA, const s_type* U, i_type inc_u,
        const s_type* W, i_type inv_w, s_type* C, s_type* S, i_type* IND1, i_type* IND2, i_type& SLEN, 
        i_type& INFO);

template void BLAS_EXT_EXPORT
qruphr<c_type>(i_type N, i_type K, c_type* R, i_type LDR, const c_type& SIGMA, const c_type* U, i_type inc_u,
        const c_type* W, i_type inv_w, s_type* C, c_type* S, i_type* IND1, i_type* IND2, i_type& SLEN,
        i_type& INFO);

template void BLAS_EXT_EXPORT
qruphr<z_type>(i_type N, i_type K, z_type* R, i_type LDR, const z_type& SIGMA, const z_type* U, i_type inc_u,
        const z_type* W, i_type inv_w, d_type* C, z_type* S, i_type* IND1, i_type* IND2, i_type& SLEN, 
        i_type& INFO);


template void BLAS_EXT_EXPORT
qruphl<d_type>(i_type M, i_type N, i_type K, d_type* R, i_type LDR, const d_type& SIGMA, const d_type* U, 
        i_type inc_u, const d_type* W, i_type inv_w, d_type* C, d_type* S, i_type* IND1, i_type* IND2, 
        i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
qruphl<s_type>(i_type M, i_type N, i_type K, s_type* R, i_type LDR, const s_type& SIGMA, const s_type* U,
        i_type inc_u, const s_type* W, i_type inv_w, s_type* C, s_type* S, i_type* IND1, i_type* IND2,
        i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
qruphl<c_type>(i_type M, i_type N, i_type K, c_type* R, i_type LDR, const c_type& SIGMA, const c_type* U,
        i_type inc_u, const c_type* W, i_type inv_w, s_type* C, c_type* S, i_type* IND1, i_type* IND2,
        i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
qruphl<z_type>(i_type M, i_type N, i_type K, z_type* R, i_type LDR, const z_type& SIGMA, const z_type* U,
        i_type inc_u, const z_type* W, i_type inv_w, d_type* C, z_type* S, i_type* IND1, i_type* IND2,
        i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
huundr<d_type>(i_type N, i_type D_from, i_type D_to, d_type* R, i_type LDR, d_type* C, d_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
huundr<s_type>(i_type N, i_type D_from, i_type D_to, s_type* R, i_type LDR, s_type* C, s_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
huundr<c_type>(i_type N, i_type D_from, i_type D_to, c_type* R, i_type LDR, s_type* C, c_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
huundr<z_type>(i_type N, i_type D_from, i_type D_to, z_type* R, i_type LDR, d_type* C, z_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
hlundl<d_type>(i_type N, i_type D, d_type* R, i_type LDR, d_type* C, d_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
hlundl<s_type>(i_type N, i_type D, s_type* R, i_type LDR, s_type* C, s_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
hlundl<c_type>(i_type N, i_type D, c_type* R, i_type LDR, s_type* C, c_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
hlundl<z_type>(i_type N, i_type D, z_type* R, i_type LDR, d_type* C, z_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
huundl<d_type>(i_type M, i_type N, i_type D, d_type* R, i_type LDR, d_type* C, d_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
huundl<s_type>(i_type M, i_type N, i_type D, s_type* R, i_type LDR, s_type* C, s_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
huundl<c_type>(i_type M, i_type N, i_type D, c_type* R, i_type LDR, s_type* C, c_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
huundl<z_type>(i_type M, i_type N, i_type D, z_type* R, i_type LDR, d_type* C, z_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
hlundr<d_type>(i_type M, i_type N, i_type D, d_type* R, i_type LDR, d_type* C, d_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
hlundr<s_type>(i_type M, i_type N, i_type D, s_type* R, i_type LDR, s_type* C, s_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
hlundr<c_type>(i_type M, i_type N, i_type D, c_type* R, i_type LDR, s_type* C, c_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

template void BLAS_EXT_EXPORT
hlundr<z_type>(i_type M, i_type N, i_type D, z_type* R, i_type LDR, d_type* C, z_type* S, i_type* IND1,
        i_type* IND2, i_type& SLEN, i_type& INFO);

};};