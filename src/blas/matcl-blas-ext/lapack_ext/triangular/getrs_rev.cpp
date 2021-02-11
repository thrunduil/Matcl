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
getrs_rev(const char *trans,i_type n,i_type nrhs,const V *a,i_type lda,const i_type *ipiv,V *b,i_type ldb,
                        i_type *info)
{
    // Test the input parameters.
    *info = 0;
    bool notran = (trans[0] == 'N' || trans[0] == 'n');
     
    if (notran == false && (trans[0] == 'T' || trans[0] == 't') == false
                        && (trans[0] == 'C' || trans[0] == 'c') == false)
    {
        *info = -1;
    }
    else if (n < 0)
    {
        *info = -2;
    }
    else if (nrhs < 0)
    {
        *info = -3;
    }
    else if (lda < lapack::maximum((matcl::lapack::i_type)1, n))
    {
        *info = -5;
    }
    else if (ldb < lapack::maximum((matcl::lapack::i_type)1, nrhs))
    {
        *info = -8;
    }
    if (*info != 0)
        return;

    // Quick return if possible
    if ( n == 0 || nrhs == 0)
        return;

    V one = 1;

    if( notran )
    {
        // Solve X * A = B.

        // Solve X * U = B, overwriting B with X.
        lapack::trsm("Right", "Upper", "No transpose", "Non-unit", nrhs, n, one, a, lda, b, ldb);

        // Solve X * L = B, overwriting B with X.
        lapack::trsm("Right", "Lower", "No trans", "Unit", nrhs, n, one, a, lda, b, ldb);

        // Apply column interchanges to X
        lapack::laswpc(nrhs, b, ldb, 1, n, ipiv, -1);
    }
    else
    {
        // Solve X * A' = B.

        // Apply column interchanges to the solution vectors.
        lapack::laswpc(nrhs, b, ldb, 1, n, ipiv, 1);

        // Solve X * L' = B, overwriting B with X.
        lapack::trsm("Right", "Lower", trans, "Unit", nrhs, n, one, a, lda, b, ldb);

        // Solve X * U' = B, overwriting B with X.
        lapack::trsm("Right", "Upper", trans, "Non-unit", nrhs, n, one, a, lda, b, ldb);
    };
};

template void BLAS_EXT_EXPORT
getrs_rev<s_type>(const char *trans, i_type n, i_type nrhs,const s_type *a,i_type lda,const i_type *ipiv,s_type *b,i_type ldb,
                        i_type *info);
template void BLAS_EXT_EXPORT
getrs_rev<d_type>(const char *trans, i_type n, i_type nrhs,const d_type *a,i_type lda,const i_type *ipiv,d_type *b,i_type ldb,
                        i_type *info);
template void BLAS_EXT_EXPORT
getrs_rev<c_type>(const char *trans, i_type n, i_type nrhs,const c_type *a,i_type lda,const i_type *ipiv,c_type *b,i_type ldb,
                        i_type *info);
template void BLAS_EXT_EXPORT
getrs_rev<z_type>(const char *trans, i_type n, i_type nrhs,const z_type *a,i_type lda,const i_type *ipiv,z_type *b,i_type ldb,
                        i_type *info);

BLAS_EXT_EXPORT void cgetrs_rev(const char *trans,i_type n,i_type nrhs,const c_type *a,i_type lda,
                            const i_type *ipiv,c_type *b,i_type ldb,i_type *info)
{
    return getrs_rev<c_type>(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
};
BLAS_EXT_EXPORT void sgetrs_rev(const char *trans,i_type n,i_type nrhs,const s_type *a,i_type lda,
                           const i_type *ipiv,s_type *b,i_type ldb,i_type *info)
{
    return getrs_rev<s_type>(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
};
BLAS_EXT_EXPORT void dgetrs_rev(const char *trans,i_type n,i_type nrhs,const d_type *a,i_type lda,
                           const i_type *ipiv,d_type *b,i_type ldb,i_type *info)
{
    return getrs_rev<d_type>(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
};
BLAS_EXT_EXPORT void zgetrs_rev(const char *trans,i_type n,i_type nrhs,const z_type *a,i_type lda,
                           const i_type *ipiv,z_type *b,i_type ldb,i_type *info)
{
    return getrs_rev<z_type>(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
};

};};