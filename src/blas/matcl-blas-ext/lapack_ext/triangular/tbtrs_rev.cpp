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
tbtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs,
                           const V *ab,i_type ldab,V *b,i_type ldb,i_type *info)
{
    *info   = 0;

    V zero  = 0;
    V one   = 1;

    bool nounit     = (diag[0] == 'N' || diag[0] == 'n');
    bool upper      = (uplo[0] == 'U' || uplo[0] == 'u');
    bool trans_n    = (trans[0] == 'N' || trans[0] == 'n');
    bool trans_t    = (trans[0] == 'T' || trans[0] == 't');
    bool trans_c    = (trans[0] == 'C' || trans[0] == 'c');

    if ( upper == false &&  (uplo[0] == 'L' || uplo[0] == 'l') == false)
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
    else if (kd < 0)
    {
        *info = -5;
    }
    else if (nrhs < 0)
    {
        *info = -6;
    }
    else if (ldab < kd + 1)
    {
        *info = -8;
    }
    else if (ldb < lapack::maximum((matcl::lapack::i_type)1,nrhs))
    {
        *info = -10;
    }

    if (*info != 0)
        return;

    //Quick return if possible
    if (n == 0 || nrhs == 0)
        return;

    // Check for singularity.

    if (nounit)
    {
        if (upper)
        {
            for (i_type i = 0; i < n; ++i)
            {
                if( ab[kd + i * ldab] == zero)
                {
                    *info = i + 1;
                    return;
                };
            };
        }
        else
        {
            for (i_type i = 0; i < n; ++i)
            {            
                if ( ab[0 + i * ldab] == zero )
                {
                    *info = i + 1;
                    return;
                };
            };
        };
    };

    // Solve X * A = B  or  X* A' = B.
    if (trans_n == true)
    {
        for (i_type j = 0; j < nrhs; ++j)
            lapack::tbsv( uplo, "T", diag, n, kd, ab, ldab, b + j, ldb );
    }
    else if (trans_t == true)
    {
        for (i_type j = 0; j < nrhs; ++j)
            lapack::tbsv( uplo, "N", diag, n, kd, ab, ldab, b + j, ldb );
    }
    else if (trans_c == true)
    {
        for (i_type j = 0; j < nrhs; ++j)
        {
            lapack::lacgv(n,b + j, ldb);
            lapack::tbsv( uplo, "N", diag, n, kd, ab, ldab, b + j, ldb );
            lapack::lacgv(n,b + j, ldb);
        };
    };
};

template void BLAS_EXT_EXPORT
tbtrs_rev<s_type>(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs,
                           const s_type *ab,i_type ldab, s_type *b,i_type ldb,i_type *info);
template void BLAS_EXT_EXPORT
tbtrs_rev<d_type>(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs,
                           const d_type *ab,i_type ldab, d_type *b,i_type ldb,i_type *info);
template void BLAS_EXT_EXPORT
tbtrs_rev<c_type>(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs,
                           const c_type *ab,i_type ldab, c_type *b,i_type ldb,i_type *info);
template void BLAS_EXT_EXPORT
tbtrs_rev<z_type>(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs,
                           const z_type *ab,i_type ldab, z_type *b,i_type ldb,i_type *info);

BLAS_EXT_EXPORT void ctbtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs, 
                           const c_type *ab,i_type ldab,c_type *b,i_type ldb,i_type *info)
{
    tbtrs_rev<c_type>(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, info);
};
BLAS_EXT_EXPORT void stbtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs,
                           const s_type *ab,i_type ldab,s_type *b,i_type ldb,i_type *info)
{
    tbtrs_rev<s_type>(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, info);
};
BLAS_EXT_EXPORT void dtbtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs,
                           const d_type *ab,i_type ldab,d_type *b,i_type ldb,i_type *info)
{
    tbtrs_rev<d_type>(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, info);
};
BLAS_EXT_EXPORT void ztbtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs, 
                           const z_type *ab,i_type ldab,z_type *b,i_type ldb,i_type *info)
{
    tbtrs_rev<z_type>(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, info);
};

};};