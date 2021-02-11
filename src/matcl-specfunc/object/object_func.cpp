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
#include "matcl-dynamic/predefined_functions.h"
#include "matcl-dynamic/function.h"
#include "matcl-dynamic/object_type.h"
#include "matcl-scalar/object.h"

namespace matcl { namespace function_name
{

dynamic::function_name function_name::beta::eval()
{
    static dynamic::function_name f("beta");
    return f;
};
dynamic::function_name function_name::ibeta::eval()
{
    static dynamic::function_name f("ibeta");
    return f;
};
dynamic::function_name function_name::ibetac::eval()
{
    static dynamic::function_name f("ibetac");
    return f;
};
dynamic::function_name function_name::ibeta_norm::eval()
{
    static dynamic::function_name f("ibeta_norm");
    return f;
};
dynamic::function_name function_name::ibetac_norm::eval()
{
    static dynamic::function_name f("ibetac_norm");
    return f;
};
dynamic::function_name function_name::ibeta_inv::eval()
{
    static dynamic::function_name f("ibeta_inv");
    return f;
};
dynamic::function_name function_name::ibetac_inv::eval()
{
    static dynamic::function_name f("ibetac_inv");
    return f;
};
dynamic::function_name function_name::ibeta_inva::eval()
{
    static dynamic::function_name f("ibeta_inva");
    return f;
};
dynamic::function_name function_name::ibetac_inva::eval()
{
    static dynamic::function_name f("ibetac_inva");
    return f;
};
dynamic::function_name function_name::ibeta_invb::eval()
{
    static dynamic::function_name f("ibeta_invb");
    return f;
};
dynamic::function_name function_name::ibetac_invb::eval()
{
    static dynamic::function_name f("ibetac_invb");
    return f;
};
dynamic::function_name function_name::ibeta_dif::eval()
{
    static dynamic::function_name f("ibeta_dif");
    return f;
};

dynamic::function_name function_name::erf::eval()
{
    static dynamic::function_name f("erf");
    return f;
};
dynamic::function_name function_name::erfc::eval()
{
    static dynamic::function_name f("erfc");
    return f;
};
dynamic::function_name function_name::erf_inv::eval()
{
    static dynamic::function_name f("erf_inv");
    return f;
};
dynamic::function_name function_name::erfc_inv::eval()
{
    static dynamic::function_name f("erfc_inv");
    return f;
};

dynamic::function_name function_name::gamma::eval()
{
    static dynamic::function_name f("gamma");
    return f;
};
dynamic::function_name function_name::gammaln::eval()
{
    static dynamic::function_name f("gammaln");
    return f;
};

dynamic::function_name function_name::gamma1pm1::eval()
{
    static dynamic::function_name f("gamma1pm1");
    return f;
};
dynamic::function_name function_name::digamma::eval()
{
    static dynamic::function_name f("digamma");
    return f;
};
dynamic::function_name function_name::trigamma::eval()
{
    static dynamic::function_name f("trigamma");
    return f;
};
dynamic::function_name function_name::polygamma::eval()
{
    static dynamic::function_name f("polygamma");
    return f;
};
dynamic::function_name function_name::gamma_ratio::eval()
{
    static dynamic::function_name f("gamma_ratio");
    return f;
};
dynamic::function_name function_name::gamma_delta_ratio::eval()
{
    static dynamic::function_name f("gamma_delta_ratio");
    return f;
};
dynamic::function_name function_name::igamma_lower::eval()
{
    static dynamic::function_name f("igamma_lower");
    return f;
};
dynamic::function_name function_name::igamma_upper::eval()
{
    static dynamic::function_name f("igamma_upper");
    return f;
};
dynamic::function_name function_name::igamma_lower_norm::eval()
{
    static dynamic::function_name f("igamma_lower_norm");
    return f;
};
dynamic::function_name function_name::igamma_upper_norm::eval()
{
    static dynamic::function_name f("igamma_upper_norm");
    return f;
};
dynamic::function_name function_name::igamma_lower_inv::eval()
{
    static dynamic::function_name f("igamma_lower_inv");
    return f;
};
dynamic::function_name function_name::igamma_upper_inv::eval()
{
    static dynamic::function_name f("igamma_upper_inv");
    return f;
};
dynamic::function_name function_name::igamma_lower_inva::eval()
{
    static dynamic::function_name f("igamma_lower_inva");
    return f;
};
dynamic::function_name function_name::igamma_upper_inva::eval()
{
    static dynamic::function_name f("igamma_upper_inva");
    return f;
};
dynamic::function_name function_name::igamma_lower_dif::eval()
{
    static dynamic::function_name f("igamma_lower_dif");
    return f;
};

dynamic::function_name function_name::zeta::eval()
{
    static dynamic::function_name f("zeta");
    return f;
};

dynamic::function_name function_name::airy_ai::eval()
{
    static dynamic::function_name f("airy_ai");
    return f;
};
dynamic::function_name function_name::airy_bi::eval()
{
    static dynamic::function_name f("airy_bi");
    return f;
};
dynamic::function_name function_name::airy_ai_dif::eval()
{
    static dynamic::function_name f("airy_ai_dif");
    return f;
};
dynamic::function_name function_name::airy_bi_dif::eval()
{
    static dynamic::function_name f("airy_bi_dif");
    return f;
};

dynamic::function_name function_name::cyl_bessel_j::eval()
{
    static dynamic::function_name f("cyl_bessel_j");
    return f;
};
dynamic::function_name function_name::cyl_neumann::eval()
{
    static dynamic::function_name f("cyl_neumann");
    return f;
};
dynamic::function_name function_name::cyl_bessel_j_zero::eval()
{
    static dynamic::function_name f("cyl_bessel_j_zero");
    return f;
};
dynamic::function_name function_name::cyl_neumann_zero::eval()
{
    static dynamic::function_name f("cyl_neumann_zero");
    return f;
};
dynamic::function_name function_name::cyl_bessel_j_dif::eval()
{
    static dynamic::function_name f("cyl_bessel_j_dif");
    return f;
};
dynamic::function_name function_name::cyl_neumann_dif::eval()
{
    static dynamic::function_name f("cyl_neumann_dif");
    return f;
};
dynamic::function_name function_name::cyl_bessel_i::eval()
{
    static dynamic::function_name f("cyl_bessel_i");
    return f;
};
dynamic::function_name function_name::cyl_bessel_k::eval()
{
    static dynamic::function_name f("cyl_bessel_k");
    return f;
};
dynamic::function_name function_name::cyl_bessel_i_dif::eval()
{
    static dynamic::function_name f("cyl_bessel_i_dif");
    return f;
};
dynamic::function_name function_name::cyl_bessel_k_dif::eval()
{
    static dynamic::function_name f("cyl_bessel_k_dif");
    return f;
};
dynamic::function_name function_name::sph_bessel::eval()
{
    static dynamic::function_name f("sph_bessel");
    return f;
};
dynamic::function_name function_name::sph_neumann::eval()
{
    static dynamic::function_name f("sph_neumann");
    return f;
};
dynamic::function_name function_name::sph_bessel_dif::eval()
{
    static dynamic::function_name f("sph_bessel_dif");
    return f;
};
dynamic::function_name function_name::sph_neumann_dif::eval()
{
    static dynamic::function_name f("sph_neumann_dif");
    return f;
};

dynamic::function_name function_name::cyl_hankel_1::eval()
{
    static dynamic::function_name f("cyl_hankel_1");
    return f;
};
dynamic::function_name function_name::cyl_hankel_2::eval()
{
    static dynamic::function_name f("cyl_hankel_2");
    return f;
};
dynamic::function_name function_name::sph_hankel_1::eval()
{
    static dynamic::function_name f("sph_hankel_1");
    return f;
};
dynamic::function_name function_name::sph_hankel_2::eval()
{
    static dynamic::function_name f("sph_hankel_2");
    return f;
};

dynamic::function_name function_name::ellint_rf::eval()
{
    static dynamic::function_name f("ellint_rf");
    return f;
};
dynamic::function_name function_name::ellint_rd::eval()
{
    static dynamic::function_name f("ellint_rd");
    return f;
};
dynamic::function_name function_name::ellint_rj::eval()
{
    static dynamic::function_name f("ellint_rj");
    return f;
};
dynamic::function_name function_name::ellint_rc::eval()
{
    static dynamic::function_name f("ellint_rc");
    return f;
};
dynamic::function_name function_name::ellint_rg::eval()
{
    static dynamic::function_name f("ellint_rg");
    return f;
};
dynamic::function_name function_name::ellint_1::eval()
{
    static dynamic::function_name f("ellint_1");
    return f;
};
dynamic::function_name function_name::ellint_2::eval()
{
    static dynamic::function_name f("ellint_2");
    return f;
};
dynamic::function_name function_name::ellint_3::eval()
{
    static dynamic::function_name f("ellint_3");
    return f;
};
dynamic::function_name function_name::ellint_d::eval()
{
    static dynamic::function_name f("ellint_d");
    return f;
};
dynamic::function_name function_name::jacobi_sn::eval()
{
    static dynamic::function_name f("jacobi_sn");
    return f;
};
dynamic::function_name function_name::jacobi_cn::eval()
{
    static dynamic::function_name f("jacobi_cn");
    return f;
};
dynamic::function_name function_name::jacobi_dn::eval()
{
    static dynamic::function_name f("jacobi_dn");
    return f;
};
dynamic::function_name function_name::jacobi_elliptic::eval()
{
    static dynamic::function_name f("jacobi_elliptic");
    return f;
};
dynamic::function_name function_name::expint::eval()
{
    static dynamic::function_name f("expint");
    return f;
};
dynamic::function_name function_name::owens_t::eval()
{
    static dynamic::function_name f("owens_t");
    return f;
};
dynamic::function_name function_name::sinc::eval()
{
    static dynamic::function_name f("sinc");
    return f;
};
dynamic::function_name function_name::sinhc::eval()
{
    static dynamic::function_name f("sinhc");
    return f;
};

};};

namespace matcl
{

Object object_func::beta(const Object& arg1, const Object& arg2)
{
    dynamic::function_name fn = function_name::beta::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, arg1, arg2);
    return ret;
};
Object object_func::ibeta(const Object& A, const Object& B, const Object& C)
{
    dynamic::function_name fn = function_name::ibeta::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B, C);
    return ret;
};
Object object_func::ibetac(const Object& A, const Object& B, const Object& C)
{
    dynamic::function_name fn = function_name::ibetac::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B, C);
    return ret;
};
Object object_func::ibeta_norm(const Object& A, const Object& B, const Object& C)
{
    dynamic::function_name fn = function_name::ibeta_norm::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B, C);
    return ret;
};
Object object_func::ibetac_norm(const Object& A, const Object& B, const Object& C)
{
    dynamic::function_name fn = function_name::ibetac_norm::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B, C);
    return ret;
};
Object object_func::ibeta_inv(const Object& A, const Object& B, const Object& C)
{
    dynamic::function_name fn = function_name::ibeta_inv::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B, C);
    return ret;
};
Object object_func::ibetac_inv(const Object& A, const Object& B, const Object& C)
{
    dynamic::function_name fn = function_name::ibetac_inv::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B, C);
    return ret;
};
Object object_func::ibeta_inva(const Object& A, const Object& B, const Object& C)
{
    dynamic::function_name fn = function_name::ibeta_inva::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B, C);
    return ret;
};
Object object_func::ibetac_inva(const Object& A, const Object& B, const Object& C)
{
    dynamic::function_name fn = function_name::ibetac_inva::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B, C);
    return ret;
};
Object object_func::ibeta_invb(const Object& A, const Object& B, const Object& C)
{
    dynamic::function_name fn = function_name::ibeta_invb::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B, C);
    return ret;
};
Object object_func::ibetac_invb(const Object& A, const Object& B, const Object& C)
{
    dynamic::function_name fn = function_name::ibetac_invb::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B, C);
    return ret;
};
Object object_func::ibeta_dif(const Object& A, const Object& B, const Object& C)
{
    dynamic::function_name fn = function_name::ibeta_dif::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B, C);
    return ret;
};

Object object_func::erf(const Object& arg1)
{
    dynamic::function_name fn = function_name::erf::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, arg1);
    return ret;
};
Object object_func::erfc(const Object& arg1)
{
    dynamic::function_name fn = function_name::erfc::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, arg1);
    return ret;
};
Object object_func::erf_inv(const Object& arg1)
{
    dynamic::function_name fn = function_name::erf_inv::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, arg1);
    return ret;
};
Object object_func::erfc_inv(const Object& arg1)
{
    dynamic::function_name fn = function_name::erfc_inv::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, arg1);
    return ret;
};

Object object_func::gamma(const Object& arg1)
{
    dynamic::function_name fn = function_name::gamma::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, arg1);
    return ret;
};
Object object_func::gammaln(const Object& arg1)
{
    dynamic::function_name fn = function_name::gammaln::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, arg1);
    return ret;
};
Object object_func::gamma1pm1(const Object& arg1)
{
    dynamic::function_name fn = function_name::gamma1pm1::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, arg1);
    return ret;
};
Object object_func::digamma(const Object& arg1)
{
    dynamic::function_name fn = function_name::digamma::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, arg1);
    return ret;
};
Object object_func::trigamma(const Object& arg1)
{
    dynamic::function_name fn = function_name::trigamma::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, arg1);
    return ret;
};
    
Object object_func::gamma_ratio(const Object& A, const Object& B)
{
    dynamic::function_name fn = function_name::gamma_ratio::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B);
    return ret;
};

Object object_func::gamma_delta_ratio(const Object& A, const Object& B)
{
    dynamic::function_name fn = function_name::gamma_delta_ratio::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B);
    return ret;
};

Object object_func::igamma_lower(const Object& A, const Object& B)
{
    dynamic::function_name fn = function_name::igamma_lower::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B);
    return ret;
};

Object object_func::igamma_upper(const Object& A, const Object& B)
{
    dynamic::function_name fn = function_name::igamma_upper::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B);
    return ret;
};

Object object_func::igamma_lower_norm(const Object& A, const Object& B)
{
    dynamic::function_name fn = function_name::igamma_lower_norm::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B);
    return ret;
};

Object object_func::igamma_upper_norm(const Object& A, const Object& B)
{
    dynamic::function_name fn = function_name::igamma_upper_norm::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B);
    return ret;
};

Object object_func::igamma_lower_inv(const Object& A, const Object& B)
{
    dynamic::function_name fn = function_name::igamma_lower_inv::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B);
    return ret;
};

Object object_func::igamma_upper_inv(const Object& A, const Object& B)
{
    dynamic::function_name fn = function_name::igamma_upper_inv::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B);
    return ret;
};

Object object_func::igamma_lower_inva(const Object& A, const Object& B)
{
    dynamic::function_name fn = function_name::igamma_lower_inva::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B);
    return ret;
};

Object object_func::igamma_upper_inva(const Object& A, const Object& B)
{
    dynamic::function_name fn = function_name::igamma_upper_inva::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B);
    return ret;
};

Object object_func::igamma_lower_dif(const Object& A, const Object& B)
{
    dynamic::function_name fn = function_name::igamma_lower_dif::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B);
    return ret;
};

Object object_func::polygamma(const Object& A, Integer n)
{
    dynamic::function_name fn = function_name::polygamma::eval();
    OInteger B(n);
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A, B);
    return ret;
};

Object object_func::zeta(const Object& A)
{
    dynamic::function_name fn = function_name::zeta::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A);
    return ret;
};

Object object_func::airy_ai(const Object& A)
{
    dynamic::function_name fn = function_name::airy_ai::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A);
    return ret;
};

Object object_func::airy_bi(const Object& A)
{
    dynamic::function_name fn = function_name::airy_bi::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A);
    return ret;
};

Object object_func::airy_ai_dif(const Object& A)
{
    dynamic::function_name fn = function_name::airy_ai_dif::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A);
    return ret;
};

Object object_func::airy_bi_dif(const Object& A)
{
    dynamic::function_name fn = function_name::airy_bi_dif::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, A);
    return ret;
};

Object object_func::cyl_bessel_j(const Object& v, const Object& x)
{
    dynamic::function_name fn = function_name::cyl_bessel_j::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, x);
    return ret;
};
Object object_func::cyl_neumann(const Object& v, const Object& x)
{
    dynamic::function_name fn = function_name::cyl_neumann::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, x);
    return ret;
};
Object object_func::cyl_bessel_j_zero(const Object& v, Integer m)
{
    dynamic::function_name fn = function_name::cyl_bessel_j_zero::eval();
    OInteger B(m);
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, B);
    return ret;
};
Object object_func::cyl_neumann_zero(const Object& v, Integer m)
{
    dynamic::function_name fn = function_name::cyl_neumann_zero::eval();
    OInteger B(m);
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, B);
    return ret;
};
Object object_func::cyl_bessel_j_dif(const Object& v, const Object& x)
{
    dynamic::function_name fn = function_name::cyl_bessel_j_dif::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, x);
    return ret;
};
Object object_func::cyl_neumann_dif(const Object& v, const Object& x)
{
    dynamic::function_name fn = function_name::cyl_neumann_dif::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, x);
    return ret;
};
Object object_func::cyl_bessel_i(const Object& v, const Object& x)
{
    dynamic::function_name fn = function_name::cyl_bessel_i::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, x);
    return ret;
};
Object object_func::cyl_bessel_k(const Object& v, const Object& x)
{
    dynamic::function_name fn = function_name::cyl_bessel_k::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, x);
    return ret;
};
Object object_func::cyl_bessel_i_dif(const Object& v, const Object& x)
{
    dynamic::function_name fn = function_name::cyl_bessel_i_dif::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, x);
    return ret;
};
Object object_func::cyl_bessel_k_dif(const Object& v, const Object& x)
{
    dynamic::function_name fn = function_name::cyl_bessel_k_dif::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, x);
    return ret;
};
Object object_func::sph_bessel(Integer n, const Object& x)
{
    dynamic::function_name fn = function_name::sph_bessel::eval();
    OInteger N(n);
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, N, x);
    return ret;
};
Object object_func::sph_neumann(Integer n, const Object& x)
{
    dynamic::function_name fn = function_name::sph_neumann::eval();
    OInteger N(n);
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, N, x);
    return ret;
};
Object object_func::sph_bessel_dif(Integer n, const Object& x)
{
    dynamic::function_name fn = function_name::sph_bessel_dif::eval();
    OInteger N(n);
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, N, x);
    return ret;
};
Object object_func::sph_neumann_dif(Integer n, const Object& x)
{
    dynamic::function_name fn = function_name::sph_neumann_dif::eval();
    OInteger N(n);
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, N, x);
    return ret;
};

Object object_func::cyl_hankel_1(const Object& v, const Object& x)
{
    dynamic::function_name fn = function_name::cyl_hankel_1::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, x);
    return ret;
};
Object object_func::cyl_hankel_2(const Object& v, const Object& x)
{
    dynamic::function_name fn = function_name::cyl_hankel_2::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, x);
    return ret;
};
Object object_func::sph_hankel_1(const Object& v, const Object& x)
{
    dynamic::function_name fn = function_name::sph_hankel_1::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, x);
    return ret;
};
Object object_func::sph_hankel_2(const Object& v, const Object& x)
{
    dynamic::function_name fn = function_name::sph_hankel_2::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, v, x);
    return ret;
};

Object object_func::ellint_rf(const Object& x, const Object& y, const Object& z)
{
    dynamic::function_name fn = function_name::ellint_rf::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, x, y, z);
    return ret;
};
Object object_func::ellint_rd(const Object& x, const Object& y, const Object& z)
{
    dynamic::function_name fn = function_name::ellint_rd::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, x, y, z);
    return ret;
};
Object object_func::ellint_rj(const Object& x, const Object& y, const Object& z, const Object& p)
{
    dynamic::function_name fn = function_name::ellint_rj::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, x, y, z, p);
    return ret;
};
Object object_func::ellint_rc(const Object& x, const Object& y)
{
    dynamic::function_name fn = function_name::ellint_rc::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, x, y);
    return ret;
};
Object object_func::ellint_rg(const Object& x, const Object& y, const Object& z)
{
    dynamic::function_name fn = function_name::ellint_rg::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, x, y, z);
    return ret;
};
Object object_func::ellint_1(const Object& k, const Object& phi)
{
    dynamic::function_name fn = function_name::ellint_1::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, k, phi);
    return ret;
};
Object object_func::ellint_1(const Object& k)
{
    dynamic::function_name fn = function_name::ellint_1::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, k);
    return ret;
};
Object object_func::ellint_2(const Object& k, const Object& phi)
{
    dynamic::function_name fn = function_name::ellint_2::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, k, phi);
    return ret;
};
Object object_func::ellint_2(const Object& k)
{
    dynamic::function_name fn = function_name::ellint_2::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, k);
    return ret;
};
Object object_func::ellint_3(const Object& k, const Object& n, const Object& phi)
{
    dynamic::function_name fn = function_name::ellint_3::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, k, n, phi);
    return ret;
};

Object object_func::ellint_3(const Object& k, const Object& n)
{
    dynamic::function_name fn = function_name::ellint_3::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, k, n);
    return ret;
};
Object object_func::ellint_d(const Object& k, const Object& phi)
{
    dynamic::function_name fn = function_name::ellint_d::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, k, phi);
    return ret;
};
Object object_func::ellint_d(const Object& k)
{
    dynamic::function_name fn = function_name::ellint_d::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, k);
    return ret;
};
Object object_func::jacobi_sn(const Object& k, const Object& u)
{
    dynamic::function_name fn = function_name::jacobi_sn::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, k, u);
    return ret;
};
Object object_func::jacobi_cn(const Object& k, const Object& u)
{
    dynamic::function_name fn = function_name::jacobi_cn::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, k, u);
    return ret;
};
Object object_func::jacobi_dn(const Object& k, const Object& u)
{
    dynamic::function_name fn = function_name::jacobi_dn::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, k, u);
    return ret;
};

Object object_func::expint(const Object& k, Integer n)
{
    dynamic::function_name fn = function_name::expint::eval();
    OInteger no(n);
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, k, no);
    return ret;
};
Object object_func::expint(const Object& k)
{
    dynamic::function_name fn = function_name::expint::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, k);
    return ret;
};

Object object_func::owens_t(const Object& h, const Object& a)
{
    dynamic::function_name fn = function_name::owens_t::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, h, a);
    return ret;
};

Object object_func::sinc(const Object& x)
{
    dynamic::function_name fn = function_name::sinc::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, x);
    return ret;
};
Object object_func::sinhc(const Object& x)
{
    dynamic::function_name fn = function_name::sinhc::eval();
    dynamic::object ret;
    dynamic::eval_function::eval(fn, ret, x);
    return ret;
};

};
