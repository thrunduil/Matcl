/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/decompositions/svd.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matrep/lib_functions/func_binary.h"

namespace matcl { namespace details
{

template<class V>
struct inv_22_impl
{
    using VR = typename details::real_type<V>::type;

    static void eval(V& A_11, V& A_12, V& A_21, V& A_22, VR& sig_max_r, VR& sig_min_r)
    {
        VR cos_l, cos_r;
        V sin_l, sin_r;
        V sig_max, sig_min;

        ///U * A * V' = S => A = U'*S*V, A * X = I => X = V'*inv(S)*U
        ///   [ cos_l  sin_l ] [ A_11  A_12  ] [ cos_r -sin_r ]  =  [ sig_max   0      ]
        ///   [-sin_l' cos_l ] [ A_21  A_22  ] [ sin_r' cos_r ]     [  0       sig_min ]
        matcl::svd_22(A_11, A_12, A_21, A_22, cos_l, sin_l, cos_r, sin_r, sig_max, sig_min);

        if (sig_min != V(0.0))
        {
            V sig_max_inv   = V(1.0) / sig_max;
            V sig_min_inv   = V(1.0) / sig_min;

            V sin_rc        = conj(sin_r);

            V B_11  = sig_max_inv * cos_l;
            V B_12  = sig_max_inv * sin_l;
            V B_21  = sig_min_inv * conj(-sin_l);
            V B_22  = sig_min_inv * cos_l;

            A_11    = cos_r  * B_11 - sin_r * B_21;
            A_12    = cos_r  * B_12 - sin_r * B_22;
            A_21    = sin_rc * B_11 + cos_r * B_21;
            A_22    = sin_rc * B_12 + cos_r * B_22;
        };

        sig_max_r   = abs(sig_max);
        sig_min_r   = abs(sig_min);
    };
};

};};

namespace matcl
{

template<class V, class Enable, class VR>
MATCL_LINALG_EXPORT void
matcl::inv_22(V& A_11, V& A_12, V& A_21, V& A_22, VR& sig_max, VR& sig_min)
{
    return details::inv_22_impl<V>::eval(A_11, A_12, A_21, A_22, sig_max, sig_min);
};

template MATCL_LINALG_EXPORT void 
matcl::inv_22<Real>(Real& A_11, Real& A_12, Real& A_21, Real& A_22, Real& sig_max, Real& sig_min);

template MATCL_LINALG_EXPORT void 
matcl::inv_22<Float>(Float& A_11, Float& A_12, Float& A_21, Float& A_22, Float& sig_max, Float& sig_min);

template MATCL_LINALG_EXPORT void 
matcl::inv_22<Complex>(Complex& A_11, Complex& A_12, Complex& A_21, Complex& A_22, Real& sig_max, Real& sig_min);

template MATCL_LINALG_EXPORT void 
matcl::inv_22<Float_complex>(Float_complex& A_11, Float_complex& A_12, Float_complex& A_21, Float_complex& A_22, 
                             Float& sig_max, Float& sig_min);

};