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

template<class V, bool Conj>
struct make_conj
{
    static V eval(const V& x)       { return conj(x); };
};
template<class V>
struct make_conj<V,false>
{
    static V eval(const V& x)       { return x; };
};

// return static value V if V != -1 or dynamic value v otherwise
template<i_type V>
struct get_int
{
    static i_type eval(i_type v)    { return V; };
};
template<>
struct get_int<-1>
{
    static i_type eval(i_type v)    { return v; };
};

template<class V, bool Upper, i_type T_NRHS, i_type T_LDB>
struct pttrs_rev_impl
{
    using VR    = typename details::real_type<V>::type;

    static void eval(i_type N, i_type V_NRHS, const VR* D, const V* E, V* B, i_type V_LDB)
    {
        // Solve X * A = B using the factorization A = U'*D*U,
        // overwriting each right hand side vector with its solution.
                
        // U'*D*U * X^CT = B^CT
         
        // Solve U' * x' = b'.

        i_type LDB          = get_int<T_LDB>::eval(V_LDB);
        i_type NRHS         = get_int<T_NRHS>::eval(V_NRHS);

        V* TMP_BP           = B;
        V* TMP_B            = B + 1*LDB;

        for (i_type I = 1; I < N; ++I)
        {
            V TMP_E         = make_conj<V, Upper == true>::eval(E[I-1]);

            for (i_type J = 0; J < NRHS; ++J)
                TMP_B[J]    = TMP_B[J] - TMP_BP[J] * TMP_E;

            TMP_B           += LDB;
            TMP_BP          += LDB;
        };

        // Solve D * U * x' = b'.

        TMP_B               = B + (N-1)*LDB;
        VR TMP_D            = VR(1.0)/D[N-1];

        for (i_type J = 0; J < NRHS; ++J)
            TMP_B[J]        = TMP_B[J] * TMP_D;

        TMP_B               = B + (N-2) * LDB;
        TMP_BP              = B + (N-1) * LDB;

        for (i_type I = N - 2; I >= 0; --I)
        {
            V TMP_E         = make_conj<V, Upper == false>::eval(E[I]);
            TMP_D           = VR(1.0) / D[I];

            for (i_type J = 0; J < NRHS; ++J)
                TMP_B[J]    = TMP_B[J] * TMP_D - TMP_BP[J] * TMP_E;

            TMP_B           -= LDB;
            TMP_BP          -= LDB;
        };
    };
};

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
pttrs_rev(const char* UPLO, i_type N, i_type NRHS, const typename details::real_type<V>::type* D, 
      const V* E, V* B, i_type LDB, i_type& INFO )
{
    using VR    = typename details::real_type<V>::type;

    // Test the input arguments.
    INFO        = 0;
    bool UPPER  = ( UPLO[0] == 'U' || UPLO[0] == 'u' );

    if (UPPER == false && ( UPLO[0] == 'L' || UPLO[0] == 'l' ) == false )
         INFO   = -1;
    else if (N < 0)
         INFO   = -2;
    else if (NRHS < 0)
         INFO   = -3;
    else if (LDB < std::max( 1, NRHS ) )
         INFO   = -7;

    if (INFO != 0)
        return;

    // Quick return if possible
    if ( N == 0 || NRHS == 0 )
        return;

    if( N == 1 )
    {
        lapack::scal(NRHS, VR(1.0) / D[0], B, 1);
        return;
    };
    

    if (UPPER == true)
    {
        if (NRHS == 1)
        {
            if (LDB == 1)
                pttrs_rev_impl<V,true,1,1>::eval(N, NRHS, D, E, B, LDB);
            else
                pttrs_rev_impl<V,true,1,-1>::eval(N, NRHS, D, E, B, LDB);
        }
        else
        {
            pttrs_rev_impl<V,true,-1,-1>::eval(N, NRHS, D, E, B, LDB);
        }
    }
    else
    {
        if (NRHS == 1)
        {
            if (LDB == 1)
                pttrs_rev_impl<V,false,1,1>::eval(N, NRHS, D, E, B, LDB);
            else
                pttrs_rev_impl<V,false,1,-1>::eval(N, NRHS, D, E, B, LDB);
        }
        else
        {
            pttrs_rev_impl<V,false,-1,-1>::eval(N, NRHS, D, E, B, LDB);
        }
    };
};

template BLAS_EXT_EXPORT void
pttrs_rev<d_type>(const char* UPLO, i_type N, i_type NRHS, const d_type* D, 
      const d_type* E, d_type* B, i_type LDB, i_type& INFO );

template BLAS_EXT_EXPORT void
pttrs_rev<s_type>(const char* UPLO, i_type N, i_type NRHS, const s_type* D, 
      const s_type* E, s_type* B, i_type LDB, i_type& INFO );

template BLAS_EXT_EXPORT void
pttrs_rev<z_type>(const char* UPLO, i_type N, i_type NRHS, const d_type* D, 
      const z_type* E, z_type* B, i_type LDB, i_type& INFO );

template BLAS_EXT_EXPORT void
pttrs_rev<c_type>(const char* UPLO, i_type N, i_type NRHS, const s_type* D, 
      const c_type* E, c_type* B, i_type LDB, i_type& INFO );

}};
