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

template<bool Conj>
struct make_conj
{
    template<class V>
    static V eval(const V& a) { return a; };
};
template<>
struct make_conj<true>
{
    template<class V>
    static V eval(const V& a) { return conj(a); };
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

template<class V, bool Conj>
struct gttrs_rev_impl
{
    static V cn(const V& x)
    {
        return make_conj<Conj>::eval(x);
    };

    template<i_type T_NRHS, i_type T_LDB>
    static void eval(i_type N, i_type V_NRHS, const V* DL, const V* D, const V* DU, const V* DU2, 
                            const i_type* IPIV, V* B, i_type V_LDB)
    {
        // Solve X * A = B using the LU factorization of A,
        // overwriting each right hand side vector with its solution.

        i_type NRHS = get_int<T_NRHS>::eval(V_NRHS);
        i_type LDB  = get_int<T_LDB>::eval(V_LDB);

        // Solve U'*x' = b'
        for (i_type J = 0; J < NRHS; ++J)
            B[J]            = B[J] / D[0];

        if (N > 1)
        {
            V* B1           = B + LDB;
            V d1            = V(1.0) / D[1];
            V du0           = DU[0];

            for (i_type J = 0; J < NRHS; ++J)
                B1[J]       = (B1[J] - du0 * B[J]) * d1;
        };

        V* BI               = B + 2*LDB;
        V* BI1              = B + 1*LDB;
        V* BI2              = B + 0*LDB;

        for (i_type I = 2; I < N; ++I)
        {
            V du1           = DU[I-1];
            V du2           = DU2[I-2];
            V d             = V(1.0) / D[I];

            for (i_type J = 0; J < NRHS; ++J)
                BI[J]       = (BI[J] - du1 * BI1[J] - du2 * BI2[J]) * d;

            BI              += LDB;
            BI1             += LDB;
            BI2             += LDB;
        };

        // Solve L' * x' = b'.
        BI                  = B + (N-2)*LDB;
        BI1                 = B + (N-1)*LDB;

        for (i_type I = N - 2; I >= 0; --I)
        {
            if (IPIV[I] == I+1)
            {
                for (i_type J = 0; J < NRHS; ++J)
                    BI[J]   = BI[J] - DL[I] * BI1[J];
            }
            else
            {
                for (i_type J = 0; J < NRHS; ++J)
                {            
                     V TEMP = BI1[J];
                     BI1[J] = BI[J] - DL[I] * TEMP;
                     BI[J]  = TEMP;
                };
            };

            BI              -= LDB;
            BI1             -= LDB;
        };
    };

    template<i_type T_NRHS, i_type T_LDB>
    static void eval_trans(i_type N, i_type V_NRHS, const V* DL, const V* D, const V* DU, const V* DU2, 
                            const i_type* IPIV, V* B, i_type V_LDB)
    {
        i_type NRHS = get_int<T_NRHS>::eval(V_NRHS);
        i_type LDB  = get_int<T_LDB>::eval(V_LDB);

        // Solve X * A**T = B.

        //L * U * X**T = B**T

        // Solve L * x**T = b**T
        V* BI           = B;
        V* BI1          = B + LDB;

        for (i_type I = 0; I < N - 1; ++I)
        {
            if (IPIV[I] == I+1)
            {
                for (i_type J = 0; J < NRHS; ++J)
                    BI1[J]  = BI1[J] - cn(DL[I]) * BI[J];
            }
            else
            {
                for (i_type J = 0; J < NRHS; ++J)
                {
                    V TEMP  = BI[J];
                    BI[J]   = BI1[J];
                    BI1[J]  = TEMP - cn(DL[I]) * BI[J];
                };
            };

            BI              += LDB;
            BI1             += LDB;
        };

        // Solve U * x**T = b**T
        V* BN               = B + (N-1)*LDB;

        V d                 = V(1.0)/cn(D[N-1]);
        for (i_type J = 0; J < NRHS; ++J)
            BN[J]           = BN[J] *d;

        if (N > 1)
        {
            V* BN1          = B + (N-2)*LDB;
            V dnu_1         = cn(DU[N-2]);
            V dn1           = V(1.0) / cn(D[N-2]);

            for (i_type J = 0; J < NRHS; ++J)
                BN1[J]      = (BN1[J] - dnu_1 * BN[J]) * dn1;
        };

        BI                  = B + (N-3)*LDB;
        BI1                 = B + (N-2)*LDB;
        V* BI2              = B + (N-1)*LDB;

        for (i_type I = N - 3; I >= 0; --I)
        {
            V dui           = cn(DU[I]);
            V du2i          = cn(DU2[I]);
            V di            = V(1.0) / cn(D[I]);

            for (i_type J = 0; J < NRHS; ++J)
                BI[J]       = (BI[J] - dui * BI1[J] - du2i * BI2[J]) * di;

            BI              -= LDB;
            BI1             -= LDB;
            BI2             -= LDB;
        };
    };
};

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::gttrs_rev(const char* TRANS, i_type N, i_type NRHS, const V* DL, const V* D, const V* DU, const V* DU2, 
      const i_type* IPIV, V* B, i_type LDB, i_type& INFO)
{
    INFO            = 0;
    bool NOTRAN     = (TRANS[0] == 'N' || TRANS[0] == 'n');
    bool TRAN_T     = (TRANS[0] == 'T' || TRANS[0] == 't');
    bool TRAN_C     = (TRANS[0] == 'C' || TRANS[0] == 'c');

    if (!NOTRAN && !TRAN_T && !TRAN_C)
        INFO        = -1;
    else if (N < 0)
        INFO        = -2;
    else if (NRHS < 0)
        INFO        = -3;
    else if (LDB < std::max(NRHS, 1))
        INFO        = -10;

    if (INFO != 0)
        return;
          
    // Quick return if possible
    if (N == 0 || NRHS == 0)
        return;

    if (NOTRAN)
    {
        if (NRHS == 1)
        {
            if (LDB == 1)
                gttrs_rev_impl<V,false>::eval<1,1>(N, NRHS, DL, D, DU, DU2, IPIV, B, LDB);
            else
                gttrs_rev_impl<V,false>::eval<1,-1>(N, NRHS, DL, D, DU, DU2, IPIV, B, LDB);
        }
        else
        {
            gttrs_rev_impl<V,false>::eval<-1,-1>(N, NRHS, DL, D, DU, DU2, IPIV, B, LDB);
        };
    }
    else
    {
        bool conj_A   = (TRAN_T)? false : true;

        if (conj_A)
        {
            if (NRHS == 1)
            {
                if (LDB == 1)
                    gttrs_rev_impl<V,true>::eval_trans<1,1>(N, NRHS, DL, D, DU, DU2, IPIV, B, LDB);
                else
                    gttrs_rev_impl<V,true>::eval_trans<1,-1>(N, NRHS, DL, D, DU, DU2, IPIV, B, LDB);
            }
            else
            {
                gttrs_rev_impl<V,true>::eval_trans<-1,-1>(N, NRHS, DL, D, DU, DU2, IPIV, B, LDB);
            };
        }
        else
        {
            if (NRHS == 1)
            {
                if (LDB == 1)
                    gttrs_rev_impl<V,false>::eval_trans<1,1>(N, NRHS, DL, D, DU, DU2, IPIV, B, LDB);
                else
                    gttrs_rev_impl<V,false>::eval_trans<1,-1>(N, NRHS, DL, D, DU, DU2, IPIV, B, LDB);
            }
            else
            {
                gttrs_rev_impl<V,false>::eval_trans<-1,-1>(N, NRHS, DL, D, DU, DU2, IPIV, B, LDB);
            };
        };
    };
};

template BLAS_EXT_EXPORT void
gttrs_rev<d_type>(const char* TRANS, i_type N, i_type NRHS, const d_type* DL, const d_type* D, 
        const d_type* DU, const d_type* DU2, const i_type* IPIV, d_type* B, i_type LDB, i_type& INFO);

template BLAS_EXT_EXPORT void
gttrs_rev<s_type>(const char* TRANS, i_type N, i_type NRHS, const s_type* DL, const s_type* D, 
        const s_type* DU, const s_type* DU2, const i_type* IPIV, s_type* B, i_type LDB, i_type& INFO);

template BLAS_EXT_EXPORT void
gttrs_rev<z_type>(const char* TRANS, i_type N, i_type NRHS, const z_type* DL, const z_type* D, 
        const z_type* DU, const z_type* DU2, const i_type* IPIV, z_type* B, i_type LDB, i_type& INFO);

template BLAS_EXT_EXPORT void
gttrs_rev<c_type>(const char* TRANS, i_type N, i_type NRHS, const c_type* DL, const c_type* D, 
        const c_type* DU, const c_type* DU2, const i_type* IPIV, c_type* B, i_type LDB, i_type& INFO);

}};