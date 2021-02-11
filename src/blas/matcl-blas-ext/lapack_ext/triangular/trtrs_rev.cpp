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
trtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const V *a,i_type lda,V *b,i_type ldb,i_type *info)
{
    *info = 0;

    V zero  = 0;
    V one   = 1;

    bool nounit = (diag[0] == 'N' || diag[0] == 'n');
    if ( (uplo[0] == 'U' || uplo[0] == 'u') == false &&  
         (uplo[0] == 'L' || uplo[0] == 'l') == false)
    {
        *info = -1;
    }
    else if ( (trans[0] == 'N' || trans[0] == 'n') == false && 
              (trans[0] == 'T' || trans[0] == 't') == false && 
              (trans[0] == 'C' || trans[0] == 'c') == false)
    {
        *info = -2;
    }
    else if (nounit == false && (diag[0] == 'U' || diag[0] == 'u') == false)
    {
        *info = -3;
    }
    else if (n < 0)
    {
        *info = -4;
    }
    else if (nrhs < 0)
    {
        *info = -5;
    }
    else if (lda < lapack::maximum((matcl::lapack::i_type)1,n))
    {
        *info = -7;
    }
    else if (ldb < lapack::maximum((matcl::lapack::i_type)1,nrhs))
    {
        *info = -9;
    }
    if (*info != 0)
        return;

    //Quick return if possible
    if (n == 0 || nrhs == 0)
        return;

    //Check for singularity.
    if (nounit == true)
    {
         for (i_type i = 0; i < n; ++i)
         {
            if( a[i + i * lda] == zero)
            {
                *info = i + 1;
                return;
            }
         };
    };

    // Solve x * A = b  or  x * A**T = b.
    lapack::trsm("Right", uplo, trans, diag, nrhs, n, one, a, lda, b, ldb);
};

template void BLAS_EXT_EXPORT
trtrs_rev<s_type>(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const s_type *a,i_type lda,s_type *b,i_type ldb,i_type *info);
template void BLAS_EXT_EXPORT
trtrs_rev<d_type>(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const d_type *a,i_type lda,d_type *b,i_type ldb,i_type *info);
template void BLAS_EXT_EXPORT
trtrs_rev<c_type>(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const c_type *a,i_type lda,c_type *b,i_type ldb,i_type *info);
template void BLAS_EXT_EXPORT
trtrs_rev<z_type>(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const z_type *a,i_type lda,z_type *b,i_type ldb,i_type *info);

BLAS_EXT_EXPORT void dtrtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const d_type *a,i_type lda,d_type *b,i_type ldb,i_type *info)
{
    trtrs_rev<d_type>(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
};
BLAS_EXT_EXPORT void ctrtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const c_type *a,i_type lda,c_type *b,i_type ldb,i_type *info)
{
    trtrs_rev<c_type>(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
};
BLAS_EXT_EXPORT void strtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const s_type *a,i_type lda,s_type *b,i_type ldb,i_type *info)
{
    trtrs_rev<s_type>(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
};
BLAS_EXT_EXPORT void ztrtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const z_type *a,i_type lda,z_type *b,i_type ldb,i_type *info)
{
    trtrs_rev<z_type>(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
};

};};