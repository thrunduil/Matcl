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

#include "matcl-specfunc/objects/object_func.h"
#include "matcl-dynamic/register_function.h"
#include "matcl-dynamic/object_type.h"
#include "matcl-dynamic/details/register_function_macro.h"

#include "matcl-specfunc/matcl_specfunc.h"

namespace matcl 
{

using Type = Real;

namespace mdd = matcl::dynamic::details;

#define MATCL_REGISTER_SCALAR_FUNC_SIMP(fn)                                 \
MATCL_REGISTER_SCALAR_FUNC(f_##fn, matcl::fn, Type, function_name::fn)   \

MATCL_REGISTER_SCALAR_FUNC_SIMP(erf)
MATCL_REGISTER_SCALAR_FUNC_SIMP(erfc)
MATCL_REGISTER_SCALAR_FUNC_SIMP(erf_inv)
MATCL_REGISTER_SCALAR_FUNC_SIMP(erfc_inv)

MATCL_REGISTER_SCALAR_FUNC_SIMP(gamma)
MATCL_REGISTER_SCALAR_FUNC_SIMP(gammaln)
MATCL_REGISTER_SCALAR_FUNC_SIMP(gamma1pm1)
MATCL_REGISTER_SCALAR_FUNC_SIMP(digamma)
MATCL_REGISTER_SCALAR_FUNC_SIMP(trigamma)
MATCL_REGISTER_BIN_FUNC(f_gamma_ratio, matcl::gamma_ratio, Type, Type, function_name::gamma_ratio)
MATCL_REGISTER_BIN_FUNC(f_gamma_delta_ratio, matcl::gamma_delta_ratio, Type, Type, function_name::gamma_delta_ratio)
MATCL_REGISTER_BIN_FUNC(f_igamma_lower, matcl::igamma_lower, Type, Type, function_name::igamma_lower)
MATCL_REGISTER_BIN_FUNC(f_igamma_upper, matcl::igamma_upper, Type, Type, function_name::igamma_upper)
MATCL_REGISTER_BIN_FUNC(f_igamma_lower_norm, matcl::igamma_lower_norm, Type, Type, function_name::igamma_lower_norm)
MATCL_REGISTER_BIN_FUNC(f_igamma_upper_norm, matcl::igamma_upper_norm, Type, Type, function_name::igamma_upper_norm)
MATCL_REGISTER_BIN_FUNC(f_igamma_lower_inv, matcl::igamma_lower_inv, Type, Type, function_name::igamma_lower_inv)
MATCL_REGISTER_BIN_FUNC(f_igamma_upper_inv, matcl::igamma_upper_inv, Type, Type, function_name::igamma_upper_inv)
MATCL_REGISTER_BIN_FUNC(f_igamma_lower_inva, matcl::igamma_lower_inva, Type, Type, function_name::igamma_lower_inva)
MATCL_REGISTER_BIN_FUNC(f_igamma_upper_inva, matcl::igamma_upper_inva, Type, Type, function_name::igamma_upper_inva)
MATCL_REGISTER_BIN_FUNC(f_igamma_lower_dif, matcl::igamma_lower_dif, Type, Type, function_name::igamma_lower_dif)
MATCL_REGISTER_BIN_FUNC(f_polygamma2, matcl::polygamma, Type, Integer, function_name::polygamma)

MATCL_REGISTER_SCALAR_FUNC_SIMP(zeta)

MATCL_REGISTER_BIN_FUNC(f_beta, matcl::beta, Type, Type, function_name::beta)
MATCL_REGISTER_3_FUNC(f_ibeta, matcl::ibeta, Type, Type, Type, function_name::ibeta)
MATCL_REGISTER_3_FUNC(f_ibetac, matcl::ibetac, Type, Type, Type, function_name::ibetac)
MATCL_REGISTER_3_FUNC(f_ibeta_norm, matcl::ibeta_norm, Type, Type, Type, function_name::ibeta_norm)
MATCL_REGISTER_3_FUNC(f_ibetac_norm, matcl::ibetac_norm, Type, Type, Type, function_name::ibetac_norm)
MATCL_REGISTER_3_FUNC(f_ibeta_inv, matcl::ibeta_inv, Type, Type, Type, function_name::ibeta_inv)
MATCL_REGISTER_3_FUNC(f_ibetac_inv, matcl::ibetac_inv, Type, Type, Type, function_name::ibetac_inv)
MATCL_REGISTER_3_FUNC(f_ibeta_inva, matcl::ibeta_inva, Type, Type, Type, function_name::ibeta_inva)
MATCL_REGISTER_3_FUNC(f_ibetac_inva, matcl::ibetac_inva, Type, Type, Type, function_name::ibetac_inva)
MATCL_REGISTER_3_FUNC(f_ibeta_invb, matcl::ibeta_invb, Type, Type, Type, function_name::ibeta_invb)
MATCL_REGISTER_3_FUNC(f_ibetac_invb, matcl::ibetac_invb, Type, Type, Type, function_name::ibetac_invb)
MATCL_REGISTER_3_FUNC(f_ibeta_dif, matcl::ibeta_dif, Type, Type, Type, function_name::ibeta_dif)

MATCL_REGISTER_SCALAR_FUNC_SIMP(airy_ai)
MATCL_REGISTER_SCALAR_FUNC_SIMP(airy_bi)
MATCL_REGISTER_SCALAR_FUNC_SIMP(airy_ai_dif)
MATCL_REGISTER_SCALAR_FUNC_SIMP(airy_bi_dif)

MATCL_REGISTER_BIN_FUNC(f_cyl_bessel_j, matcl::cyl_bessel_j, Type, Type, function_name::cyl_bessel_j)
MATCL_REGISTER_BIN_FUNC(f_cyl_neumann, matcl::cyl_neumann, Type, Type, function_name::cyl_neumann)
MATCL_REGISTER_BIN_FUNC(f_cyl_bessel_j_dif, matcl::cyl_bessel_j_dif, Type, Type, function_name::cyl_bessel_j_dif)
MATCL_REGISTER_BIN_FUNC(f_cyl_neumann_dif, matcl::cyl_neumann_dif, Type, Type, function_name::cyl_neumann_dif)
MATCL_REGISTER_BIN_FUNC(f_cyl_bessel_i, matcl::cyl_bessel_i, Type, Type, function_name::cyl_bessel_i)
MATCL_REGISTER_BIN_FUNC(f_cyl_bessel_k, matcl::cyl_bessel_k, Type, Type, function_name::cyl_bessel_k)
MATCL_REGISTER_BIN_FUNC(f_cyl_bessel_i_dif, matcl::cyl_bessel_i_dif, Type, Type, function_name::cyl_bessel_i_dif)
MATCL_REGISTER_BIN_FUNC(f_cyl_bessel_k_dif, matcl::cyl_bessel_k_dif, Type, Type, function_name::cyl_bessel_k_dif)
MATCL_REGISTER_BIN_FUNC(f_cyl_bessel_j_zero, matcl::cyl_bessel_j_zero, Type, Integer, function_name::cyl_bessel_j_zero)
MATCL_REGISTER_BIN_FUNC(f_cyl_neumann_zero, matcl::cyl_neumann_zero, Type, Integer, function_name::cyl_neumann_zero)
MATCL_REGISTER_BIN_FUNC(f_sph_bessel, matcl::sph_bessel, Integer, Type, function_name::sph_bessel)
MATCL_REGISTER_BIN_FUNC(f_sph_neumann, matcl::sph_neumann, Integer, Type, function_name::sph_neumann)
MATCL_REGISTER_BIN_FUNC(f_sph_bessel_dif, matcl::sph_bessel_dif, Integer, Type, function_name::sph_bessel_dif)
MATCL_REGISTER_BIN_FUNC(f_sph_neumann_dif, matcl::sph_neumann_dif, Integer, Type, function_name::sph_neumann_dif)

MATCL_REGISTER_BIN_FUNC(f_cyl_hankel_1, matcl::cyl_hankel_1, Type, Type, function_name::cyl_hankel_1)
MATCL_REGISTER_BIN_FUNC(f_cyl_hankel_2, matcl::cyl_hankel_2, Type, Type, function_name::cyl_hankel_2)
MATCL_REGISTER_BIN_FUNC(f_sph_hankel_1, matcl::sph_hankel_1, Type, Type, function_name::sph_hankel_1)
MATCL_REGISTER_BIN_FUNC(f_sph_hankel_2, matcl::sph_hankel_2, Type, Type, function_name::sph_hankel_2)

MATCL_REGISTER_SCALAR_FUNC_SIMP(ellint_1)
MATCL_REGISTER_SCALAR_FUNC_SIMP(ellint_2)
MATCL_REGISTER_SCALAR_FUNC_SIMP(ellint_d)
MATCL_REGISTER_BIN_FUNC(f_ellint_rc, matcl::ellint_rc, Type, Type, function_name::ellint_rc)
MATCL_REGISTER_BIN_FUNC(f_ellint_1_2, matcl::ellint_1, Type, Type, function_name::ellint_1)
MATCL_REGISTER_BIN_FUNC(f_ellint_2_2, matcl::ellint_2, Type, Type, function_name::ellint_2)
MATCL_REGISTER_BIN_FUNC(f_ellint_3, matcl::ellint_3, Type, Type, function_name::ellint_3)
MATCL_REGISTER_BIN_FUNC(f_ellint_d_2, matcl::ellint_d, Type, Type, function_name::ellint_d)
MATCL_REGISTER_BIN_FUNC(f_jacobi_sn, matcl::jacobi_sn, Type, Type, function_name::jacobi_sn)
MATCL_REGISTER_BIN_FUNC(f_jacobi_cn, matcl::jacobi_cn, Type, Type, function_name::jacobi_cn)
MATCL_REGISTER_BIN_FUNC(f_jacobi_dn, matcl::jacobi_dn, Type, Type, function_name::jacobi_dn)
MATCL_REGISTER_3_FUNC(f_ellint_rf, matcl::ellint_rf, Type, Type, Type, function_name::ellint_rf)
MATCL_REGISTER_3_FUNC(f_ellint_rd, matcl::ellint_rd, Type, Type, Type, function_name::ellint_rd)
MATCL_REGISTER_3_FUNC(f_ellint_rg, matcl::ellint_rg, Type, Type, Type, function_name::ellint_rg)
MATCL_REGISTER_3_FUNC(f_ellint_3_2, matcl::ellint_3, Type, Type, Type, function_name::ellint_3)
MATCL_REGISTER_4_FUNC(f_ellint_rj, matcl::ellint_rj, Type, Type, Type, Type, function_name::ellint_rj)

MATCL_REGISTER_SCALAR_FUNC_SIMP(expint)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sinc)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sinhc)
MATCL_REGISTER_BIN_FUNC(f_expint2, matcl::expint, Type, Integer, function_name::expint)
MATCL_REGISTER_BIN_FUNC(f_owens_t, matcl::owens_t, Type, Type, function_name::owens_t)
};