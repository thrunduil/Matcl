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

namespace matcl { namespace details
{

template<class V, bool Is_compl = details::is_complex<V>::value>
struct svd_22_impl
{
    using VR    = typename details::real_type<V>::type;

    static void eval(const V& A_11, const V& A_12, const V& A_21, const V& A_22, 
                     VR& cos_l,  V& sin_l, VR& cos_r, V& sin_r,V& sig_max, V& sig_min)
    {
        if (A_21 == V(0.0))
        {
            //svd of 2x2 upper triangular block
            lapack::lasv2(A_11, A_12, A_22, sig_min, sig_max, sin_r, cos_r, sin_l, cos_l);
        }
        else
        {
            //construct rotation to anihilate A_22
            V sin_g, A_11_r;
            VR cos_g;

            construct_givens2(A_11, A_12, cos_g, sin_g, A_11_r);

            V A_12_r = cos_g * A_12 + sin_g * A_22;
            V A_22_r = cos_g * A_22 - sin_g * A_12;

            V sin_ll, cos_ll;

            //svd of 2x2 upper triangular block
            lapack::lasv2(A_11_r, A_12_r, A_22_r, sig_min, sig_max, sin_r, cos_r, sin_ll, cos_ll);

            //combine left rotations
            cos_l   = cos_ll * cos_g - sin_ll * sin_g;
            sin_l   = cos_ll * sin_g + sin_ll * cos_g;
        };
    };

    static void eval_sing(const V& A_11, const V& A_12, const V& A_21, const V& A_22, 
                        V& sig_max, V& sig_min)
    {
        if (A_21 == V(0.0))
        {
            //svd of 2x2 upper triangular block
            lapack::las2(A_11, A_12, A_22, sig_min, sig_max);
        }
        else
        {
            //construct rotation to anihilate A_22
            V sin_g, A_11_r;
            VR cos_g;

            construct_givens2(A_11, A_12, cos_g, sin_g, A_11_r);

            V A_12_r = cos_g * A_12 + sin_g * A_22;
            V A_22_r = cos_g * A_22 - sin_g * A_12;

            //svd of 2x2 upper triangular block
            lapack::las2(A_11_r, A_12_r, A_22_r, sig_min, sig_max);
        };
    };

};

template<class V>
struct svd_22_impl<V,true>
{
    using VR    = typename details::real_type<V>::type;

    static void eval(const V& A_11, const V& A_12, const V& A_21, const V& A_22, 
                     VR& cos_l, V& sin_l, VR& cos_r, V& sin_r, V& sig_max, V& sig_min)
    {
        //construct rotation to anihilate A_22
        V sin_g, A_11_r, A_12_r, A_22_r;
        VR cos_g;
        bool need_g = false;

        if (A_21 != V(0.0))
        {
            need_g  = true;
            construct_givens2(A_11, A_12, cos_g, sin_g, A_11_r);

            A_12_r  = cos_g * A_12 + sin_g * A_22;
            A_22_r  = cos_g * A_22 - sin_g * A_12;
        }
        else
        {
            A_12_r  = A_12;
            A_22_r  = A_22;
        };

        // reduce A to real

        VR A_11_rr, A_12_rr, A_22_rr;
        
        V L1, L2, R2;

        if (A_11_r == V(0.0))
        {
            L1      = V(1.0);
            A_11_rr = VR(0.0);
        }
        else
        {
            A_11_rr = abs(A_11_r);
            L1      = A_11_rr / A_11_r;
            A_12_r  = L1 * A_12_r;
        }

        if (A_12_r == V(0.0))
        {
            R2      = V(1.0);
            A_12_rr = VR(0.0);
        }
        else
        {
            A_12_rr = abs(A_12_r);
            R2      = A_12_rr / A_12_r;
            A_22_r  = A_22_r * R2;
        };

        if (A_22_r == V(0.0))
        {
            L2      = V(1.0);
            A_22_rr = VR(0.0);
        }
        else
        {
            A_22_rr = abs(A_22_r);
            L2      = A_22_rr / A_22_r;
        };

        //svd of 2x2 upper triangular block
        VR sig_min_l, sig_max_l, sin_r_l, cos_r_l, sin_l_l, cos_l_l;
        lapack::lasv2<VR>(A_11_rr, A_12_rr, A_22_rr, sig_min_l, sig_max_l, sin_r_l, cos_r_l, sin_l_l, cos_l_l);
 
        //combine left rotations
        V U_11          = V(cos_l_l) * L1;
        V U_21          = V(-sin_l_l)* L1;
        V U_12          = V(sin_l_l) * L2;
        V U_22          = V(cos_l_l) * L2;

        if (need_g == true)
        {
            V U_11_c    = U_11 * cos_g - U_12 * conj(sin_g);
            V U_12_c    = U_11 * sin_g + U_12 * cos_g;
            V U_21_c    = U_21 * cos_g - U_22 * conj(sin_g);
            V U_22_c    = U_21 * sin_g + U_22 * cos_g;

            U_11        = U_11_c;
            U_21        = U_21_c;
            U_12        = U_12_c;
            U_22        = U_22_c;
        };

        //combine right rotarions
        VR V_11         = cos_r_l;
        V V_21          = V(sin_r_l) * R2;
        //V V_12        = V(-sin_r_l);
        V V_22          = V(cos_r_l) * R2;

        //convert rotation matrices to plain rotations
        VR U_11_a       = abs(U_11);
        L1              = U_11 == V(0.0) ? V(1.0) : U_11_a / U_11;
        L2              = U_22 == V(0.0) ? V(1.0) : abs(U_22) / U_22;
        VR R1           = V_11 == VR(0.0)? VR(1.0): copysign(VR(1.0), V_11);
        R2              = V_22 == V(0.0) ? V(1.0) : abs(V_22) / V_22;

        //construct final plain rotations
        cos_l           = U_11_a;
        sin_l           = U_12 * L1;

        cos_r           = V_11 * R1;
        sin_r           = conj(V_21) * R1;

        //scale singular values
        sig_max         = sig_max_l * L1 * R1;
        sig_min         = sig_min_l * L2 * R2;
    };

    static void eval_sing(const V& A_11, const V& A_12, const V& A_21, const V& A_22, 
                     VR& sig_max, VR& sig_min)
    {
        //construct rotation to anihilate A_22
        V A_11_r, A_12_r, A_22_r;      

        if (A_21 != V(0.0))
        {
            VR cos_g;
            V sin_g;

            construct_givens2(A_11, A_12, cos_g, sin_g, A_11_r);

            A_12_r  = cos_g * A_12 + sin_g * A_22;
            A_22_r  = cos_g * A_22 - sin_g * A_12;
        }
        else
        {
            A_11_r  = A_11;
            A_12_r  = A_12;
            A_22_r  = A_22;
        };

        // reduce A to real        

        VR A_11_rr, A_12_rr, A_22_rr;
        
        V L1, L2, R2;

        if (A_11_r == V(0.0))
        {
            L1      = V(1.0);
            A_11_rr = VR(0.0);
        }
        else
        {
            A_11_rr = abs(A_11_r);
            L1      = A_11_rr / A_11_r;
            A_12_r  = L1 * A_12_r;
        }

        if (A_12_r == V(0.0))
        {
            R2      = V(1.0);
            A_12_rr = VR(0.0);
        }
        else
        {
            A_12_rr = abs(A_12_r);
            R2      = A_12_rr / A_12_r;
            A_22_r  = A_22_r * R2;
        };

        if (A_22_r == V(0.0))
        {
            L2      = V(1.0);
            A_22_rr = VR(0.0);
        }
        else
        {
            A_22_rr = abs(A_22_r);
            L2      = A_22_rr / A_22_r;
        };

        //svd of 2x2 upper triangular block
        lapack::las2<VR>(A_11_rr, A_12_rr, A_22_rr, sig_min, sig_max); 
    };

};

}}

namespace matcl
{

template<class V, class Enable, class VR>
MATCL_LINALG_EXPORT void matcl::svd_22(const V& A_11, const V& A_12, const V& A_21, const V& A_22, 
                                VR& cos_l, V& sin_l, VR& cos_r, V& sin_r, V& sig_max, V& sig_min)
{
    return details::svd_22_impl<V>::eval(A_11, A_12, A_21, A_22, cos_l, sin_l, cos_r, sin_r, 
                                         sig_max, sig_min);
};

template<class V, class Enable, class VR>
MATCL_LINALG_EXPORT void matcl::svd_22(const V& A_11, const V& A_12, const V& A_21, const V& A_22, 
                                VR& sig_max, VR& sig_min)
{
    return details::svd_22_impl<V>::eval_sing(A_11, A_12, A_21, A_22, sig_max, sig_min);
};

template MATCL_LINALG_EXPORT void 
matcl::svd_22<Real>(const Real& A_11, const Real& A_12, const Real& A_21, const Real& A_22, 
                    Real& cos_l, Real& sin_l, Real& cos_r, Real& sin_r, Real& sig_max, Real& sig_min);
template MATCL_LINALG_EXPORT void 
matcl::svd_22<Float>(const Float& A_11, const Float& A_12, const Float& A_21, const Float& A_22, 
                    Float& cos_l, Float& sin_l, Float& cos_r, Float& sin_r, Float& sig_max, Float& sig_min);
template MATCL_LINALG_EXPORT void 
matcl::svd_22<Complex>(const Complex& A_11, const Complex& A_12, const Complex& A_21, const Complex& A_22, 
                    Real& cos_l, Complex& sin_l, Real& cos_r, Complex& sin_r, Complex& sig_max, Complex& sig_min);
template MATCL_LINALG_EXPORT void 
matcl::svd_22<Float_complex>(const Float_complex& A_11, const Float_complex& A_12, const Float_complex& A_21, 
                    const Float_complex& A_22, Float& cos_l, Float_complex& sin_l, Float& cos_r, Float_complex& sin_r,
                    Float_complex& sig_max, Float_complex& sig_min);

template MATCL_LINALG_EXPORT void 
matcl::svd_22<Real>(const Real& A_11, const Real& A_12, const Real& A_21, const Real& A_22, 
                    Real& sig_max, Real& sig_min);
template MATCL_LINALG_EXPORT void 
matcl::svd_22<Float>(const Float& A_11, const Float& A_12, const Float& A_21, const Float& A_22, 
                    Float& sig_max, Float& sig_min);
template MATCL_LINALG_EXPORT void 
matcl::svd_22<Complex>(const Complex& A_11, const Complex& A_12, const Complex& A_21, const Complex& A_22, 
                    Real& sig_max, Real& sig_min);
template MATCL_LINALG_EXPORT void 
matcl::svd_22<Float_complex>(const Float_complex& A_11, const Float_complex& A_12, const Float_complex& A_21, 
                    const Float_complex& A_22, Float& sig_max, Float& sig_min);

}