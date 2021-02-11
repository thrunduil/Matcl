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

namespace matcl { namespace lapack
{

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::ptinorm(i_type N, const typename details::real_type<V>::type* D, const V* E, 
                typename details::real_type<V>::type& NORM, typename details::real_type<V>::type* WORK, 
                i_type& INFO )
{
    using VR    = typename details::real_type<V>::type;

    INFO = 0;
    if ( N < 0 )
        INFO    = -1;

    if (INFO != 0)
        return;

    // Quick return if possible
    NORM    = VR(0);
    if (N == 0)
    {
         NORM   = VR(1.0);
         return;
    }

    // Check that D(1:N) is positive.
    for (i_type I = 0; I < N; ++I)
    {
        if (D[I] <= VR(0))
            return;
    }

    // Solve M(A) * x = e, where M(A) = (m(i,j)) is given by
    //        m(i,j) =  abs(A(i,j)), i = j,
    //        m(i,j) = -abs(A(i,j)), i .ne. j,
    //     and e = [ 1, 1, ..., 1 ]'.  Note M(A) = M(L)*D*M(L)'.
    //
    //    Solve M(L) * x = e.
    WORK[0]         = VR(1.0);

    for (i_type I = 1; I < N; ++I)
         WORK[I]    = VR(1.0) + WORK[I-1] * abs( E[I-1] );

    // Solve D * M(L)' * x = b.
    WORK[N-1]       = WORK[N-1] / D[N-1];

    for (i_type I = N - 2; I >= 0; --I)
         WORK[I]    = WORK[I] / D[I] + WORK[I+1] * abs( E[I] );

    // Compute NORM = max(x(i)), 1<=i<=n.
    i_type IX       = lapack::amax(N, WORK, 1 );
    NORM            = abs(WORK[IX-1] );
};

template BLAS_EXT_EXPORT void
lapack::ptinorm<d_type>(i_type N, const d_type* D, const d_type* E, d_type& NORM, d_type* WORK, i_type& INFO );

template BLAS_EXT_EXPORT void
lapack::ptinorm<s_type>(i_type N, const s_type* D, const s_type* E, s_type& NORM, s_type* WORK, i_type& INFO );

template BLAS_EXT_EXPORT void
lapack::ptinorm<z_type>(i_type N, const d_type* D, const z_type* E, d_type& NORM, d_type* WORK, i_type& INFO );

template BLAS_EXT_EXPORT void
lapack::ptinorm<c_type>(i_type N, const s_type* D, const c_type* E, s_type& NORM, s_type* WORK, i_type& INFO );

}};
