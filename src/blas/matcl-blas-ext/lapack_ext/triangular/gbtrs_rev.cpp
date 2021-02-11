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
gbtrs_rev(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const V *ab,
                           i_type ldab,const i_type *ipiv,V *b,i_type ldb,i_type *info)
{
    // Test the input parameters.
    *info = 0;

    bool notran = (trans[0] == 'N' || trans[0] == 'n');
    bool ctran = (trans[0] == 'C' || trans[0] == 'c');
     
    if (notran == false && (trans[0] == 'T' || trans[0] == 't') == false
                        && (trans[0] == 'C' || trans[0] == 'c') == false)
    {
        *info = -1;
    }
    else if (n < 0)
    {
        *info = -2;
    }
    else if (kl < 0)
    {
        *info = -3;
    }
    else if (ku < 0)
    {
        *info = -4;
    }
    else if (nrhs < 0)
    {
        *info = -5;
    }
    else if (ldab < 2*kl + ku + 1)
    {
        *info = -7;
    }
    else if (ldb < lapack::maximum((matcl::lapack::i_type)1, nrhs))
    {
        *info = -10;
    }
    if (*info != 0)
        return;

    // Quick return if possible
    if ( n == 0 || nrhs == 0)
        return;

    i_type KD   = ku + kl + 1;
    bool LNOTI  = (kl > 0);
    V one       = 1;

    if( notran )
    {
        // Solve  X * A = B.

        // Solve X * U = B, overwriting B with X.
        i_type info2 = 0;
        lapack::tbtrs_rev("Upper", "No transpose", "Non-unit", n, kl+ku, nrhs, ab, ldab, b, ldb, &info2);

        // Solve X * L = B, overwriting B with X.
        //
        // L is represented as a product of permutations and unit lower
        // triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
        // where each transformation L(i) is a rank-one modification of
        // the identity matrix.
        
        if (LNOTI == true)
        {            
            for (i_type j = n - 1; j >= 0; --j)
            {
                i_type LM   = lapack::minimum(kl, n - j - 1);
                i_type l    = ipiv[j] - 1;
                
                lapack::gemv("No trans", nrhs, LM, -one, b + (j + 1)*ldb, ldb, ab + KD + j*ldab, 1, 
                             one, b + j*ldb, 1);

                if (l != j)
                {
                    lapack::swap(nrhs, b + l*ldb, 1, b + j*ldb, 1);
                };
            };
        };
    }
    else
    {
        // Solve X * A' = B.

        // Solve X * L' = B, overwriting B with X.
        if( LNOTI == true )
        {
            for (i_type j = 0; j < n-1; ++j)
            {
                i_type LM   = lapack::minimum(kl, n - j - 1);
                i_type l    = ipiv[j] - 1;

                if ( l != j )
                {
                    lapack::swap(nrhs, b + l * ldb, 1, b + j*ldb, 1);
                };

                if (ctran == true)
                {
                    lapack::lacgv( LM,  const_cast<V*>(ab) + KD + j * ldab, 1);
                }

                lapack::geru(nrhs, LM, -one, b + j * ldb, 1, ab + KD + j * ldab, 1,  b + (j + 1)*ldb, ldb);

                if (ctran == true)
                {
                    lapack::lacgv( LM,  const_cast<V*>(ab) + KD + j * ldab, 1);
                }
            };
        };

        // Solve X * U' = B, overwriting B with X.
        i_type info2 = 0;
        lapack::tbtrs_rev("Upper", trans, "Non-unit", n, kl+ku, nrhs, ab, ldab, b, ldb, &info2);
    };
};

template void BLAS_EXT_EXPORT
gbtrs_rev<s_type>(const char *trans, i_type n, i_type kl, i_type ku, i_type nrhs,const s_type *ab,
                           i_type ldab, const i_type *ipiv, s_type *b, i_type ldb, i_type *info);
template void BLAS_EXT_EXPORT
gbtrs_rev<d_type>(const char *trans, i_type n, i_type kl, i_type ku, i_type nrhs,const d_type *ab,
                           i_type ldab, const i_type *ipiv, d_type *b, i_type ldb, i_type *info);
template void BLAS_EXT_EXPORT
gbtrs_rev<c_type>(const char *trans, i_type n, i_type kl, i_type ku, i_type nrhs,const c_type *ab,
                           i_type ldab, const i_type *ipiv, c_type *b, i_type ldb, i_type *info);
template void BLAS_EXT_EXPORT
gbtrs_rev<z_type>(const char *trans, i_type n, i_type kl, i_type ku, i_type nrhs,const z_type *ab,
                           i_type ldab, const i_type *ipiv, z_type *b, i_type ldb, i_type *info);

BLAS_EXT_EXPORT void cgbtrs_rev(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const c_type *ab,
                           i_type ldab,const i_type *ipiv,c_type *b,i_type ldb,i_type *info)
{
    return gbtrs_rev<c_type>(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
}
BLAS_EXT_EXPORT void sgbtrs_rev(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const s_type *ab,
                           i_type ldab,const i_type *ipiv,s_type *b,i_type ldb,i_type *info)
{
    return gbtrs_rev<s_type>(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
}
BLAS_EXT_EXPORT void dgbtrs_rev(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const d_type *ab,
                           i_type ldab,const i_type *ipiv,d_type *b,i_type ldb,i_type *info)
{
    return gbtrs_rev<d_type>(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
}
BLAS_EXT_EXPORT void zgbtrs_rev(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const z_type *ab,
                           i_type ldab,const i_type *ipiv,z_type *b,i_type ldb,i_type *info)
{
    return gbtrs_rev<z_type>(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
}

};};