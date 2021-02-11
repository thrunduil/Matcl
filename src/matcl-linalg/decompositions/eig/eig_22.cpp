/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-linalg/general/config_linalg.h"

#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"

#include "matcl-linalg/decompositions/givens.h"
#include "matcl-linalg/decompositions/svd.h"
#include "matcl-linalg/decompositions/schur.h"

namespace matcl { namespace details
{

template<class V, bool Is_compl = details::is_complex<V>::value>
struct schur_22_sym_impl
{
    using VR    = typename details::real_type<V>::type;

    static void eval(const V& A, const V& B, const V& C, VR& cos, V& sin, V& eig_1, V& eig_2)
    {
        return matcl::lapack::laev2(A,B,C, eig_1, eig_2,  cos, sin);
    };
};

template<class V>
struct schur_22_sym_impl<V,true>
{
    using VR    = typename details::real_type<V>::type;

    static void eval(const VR& A, const V& B, const VR& C, VR& cos, V& sin, VR& eig_1, VR& eig_2)
    {
        V W;

        if (matcl::abs(B) == VR(0.0))
            W   = V(1.0);
        else
            W   = B / abs(B);

        VR sinl;

        matcl::lapack::laev2(A, abs(B), C, eig_1, eig_2,  cos, sinl);
      
        sin = W*sinl;
    };
};

template<class V, bool Is_compl = details::is_complex<V>::value>
struct schur_22_impl
{
    using VR    = typename details::real_type<V>::type;
    using VC    = typename details::complex_type<V>::type;

    static void eval(V& A, V& B, V& C, V& D, VR& CS, V& SN, VC& eig_1, VC& eig_2)
    {
        VR RT1R, RT1I, RT2R, RT2I;
        
        lapack::lanv2(A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN);

        eig_1 = VC(RT1R, RT1I);
        eig_2 = VC(RT2R, RT2I);
    };
};

template<class V>
struct schur_22_impl<V,true>
{
    using VR    = typename details::real_type<V>::type;
    using VC    = typename details::complex_type<V>::type;
    using VL    = typename details::lapack_value_type<V>::type;

    static void eval(V& A, V& B, V& C, V& D, VR& CS, V& SN, VC& eig_1, VC& eig_2)
    {
        lapack::lanv2z<VL>(*lap(&A), *lap(&B), *lap(&C), *lap(&D), *lap(&eig_1), *lap(&eig_2), *lap(&CS), *lap(&SN));
    };
};

}};

namespace matcl
{

template<class V, class Enable, class VR>
MATCL_LINALG_EXPORT void 
matcl::schur_22_sym(const VR& A, const V& B, const VR& C, VR& cosl, V& sinl, VR& eig_1, VR& eig_2)
{
    return details::schur_22_sym_impl<V>::eval(A,B,C,cosl,sinl,eig_1, eig_2);
};

template<class V, class Enable, class VR, class VC>
MATCL_LINALG_EXPORT void 
matcl::schur_22(V& A, V& B, V& C, V& D, VR& CS, V& SN, VC& eig_1, VC& eig_2)
{
    return details::schur_22_impl<V>::eval(A,B,C,D,CS,SN,eig_1, eig_2);
};

template<class V, class Enable, class VR>
MATCL_LINALG_EXPORT void 
matcl::eig_22_sym(const VR& A, const V& B, const VR& C, VR& eig_1, VR& eig_2)
{
    using details::lap;
    return lapack::lae2(lap(A), lap(B), lap(C), *lap(&eig_1), *lap(&eig_2));
};

template<class V, class Enable, class VR, class VC>
MATCL_LINALG_EXPORT void 
matcl::eig_22(const V& A, const V& B, const V& C, const V& D, VC& eig_1, VC& eig_2)
{    
    V Al = A;
    V Bl = B;
    V Cl = C;
    V Dl = D;
    VR CS;
    V SN;

    schur_22(Al, Bl, Cl, Dl, CS, SN, eig_1, eig_2);    
};

template MATCL_LINALG_EXPORT void 
matcl::schur_22_sym<Real>(const Real& A, const Real& B, const Real& C, Real& cosl, Real& sinl, 
                        Real& eig_1, Real& eig_2);

template MATCL_LINALG_EXPORT void 
matcl::schur_22_sym<Float>(const Float& A, const Float& B, const Float& C, Float& cosl, Float& sinl, 
                        Float& eig_1, Float& eig_2);

template MATCL_LINALG_EXPORT void 
matcl::schur_22_sym<Complex>(const Real& A, const Complex& B, const Real& C, Real& cosl, 
                           Complex& sinl, Real& eig_1, Real& eig_2);

template MATCL_LINALG_EXPORT void 
matcl::schur_22_sym<Float_complex>(const Float& A, const Float_complex& B, const Float& C, Float& cosl, 
                                 Float_complex& sinl, Float& eig_1, Float& eig_2);

template MATCL_LINALG_EXPORT void 
matcl::schur_22(Real& A, Real& B, Real& C, Real& D, Real& CS, Real& SN, Complex& eig_1, Complex& eig_2);

template MATCL_LINALG_EXPORT void 
matcl::schur_22(Float& A, Float& B, Float& C, Float& D, Float& CS, Float& SN, Float_complex& eig_1, 
                Float_complex& eig_2);

template MATCL_LINALG_EXPORT void 
matcl::schur_22(Complex& A, Complex& B, Complex& C, Complex& D, Real& CS, Complex& SN, Complex& eig_1, 
                Complex& eig_2);
template MATCL_LINALG_EXPORT void 
matcl::schur_22(Float_complex& A, Float_complex& B, Float_complex& C, Float_complex& D, Float& CS, 
                Float_complex& SN, Float_complex& eig_1, Float_complex& eig_2);

template MATCL_LINALG_EXPORT void 
matcl::eig_22_sym<Real>(const Real& A, const Real& B, const Real& C, Real& eig_1, Real& eig_2);

template MATCL_LINALG_EXPORT void 
matcl::eig_22_sym<Float>(const Float& A, const Float& B, const Float& C, Float& eig_1, Float& eig_2);

template MATCL_LINALG_EXPORT void 
matcl::eig_22_sym<Complex>(const Real& A, const Complex& B, const Real& C, Real& eig_1, Real& eig_2);

template MATCL_LINALG_EXPORT void 
matcl::eig_22_sym<Float_complex>(const Float& A, const Float_complex& B, const Float& C, Float& eig_1, Float& eig_2);


template MATCL_LINALG_EXPORT void 
matcl::eig_22<Real>(const Real& A, const Real& B, const Real& C, const Real& D, Complex& eig_1, Complex& eig_2);

template MATCL_LINALG_EXPORT void 
matcl::eig_22(const Float& A, const Float& B, const Float& C, const Float& D, Float_complex& eig_1, 
                Float_complex& eig_2);
template MATCL_LINALG_EXPORT void 
matcl::eig_22(const Complex& A, const Complex& B, const Complex& C, const Complex& D, Complex& eig_1, 
                Complex& eig_2);
template MATCL_LINALG_EXPORT void 
matcl::eig_22(const Float_complex& A, const Float_complex& B, const Float_complex& C, const Float_complex& D, 
              Float_complex& eig_1, Float_complex& eig_2);

};