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

#pragma once

#include "matcl-specfunc/general/config.h"
#include "matcl-scalar/object.h"

namespace matcl 
{

namespace function_name
{    
    // predefined function names

    struct MATCL_SF_EXPORT erf                  { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT erfc                 { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT erf_inv              { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT erfc_inv             { static dynamic::function_name eval(); };

    struct MATCL_SF_EXPORT gamma                { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT gammaln              { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT gamma1pm1            { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT digamma              { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT trigamma             { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT polygamma            { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT gamma_ratio          { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT gamma_delta_ratio    { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT igamma_lower         { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT igamma_upper         { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT igamma_lower_norm    { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT igamma_upper_norm    { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT igamma_lower_inv     { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT igamma_upper_inv     { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT igamma_lower_inva    { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT igamma_upper_inva    { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT igamma_lower_dif     { static dynamic::function_name eval(); };

    struct MATCL_SF_EXPORT zeta                 { static dynamic::function_name eval(); };    

    struct MATCL_SF_EXPORT beta                 { static dynamic::function_name eval(); };    
    struct MATCL_SF_EXPORT ibeta                { static dynamic::function_name eval(); };    
    struct MATCL_SF_EXPORT ibetac               { static dynamic::function_name eval(); };    
    struct MATCL_SF_EXPORT ibeta_norm           { static dynamic::function_name eval(); };    
    struct MATCL_SF_EXPORT ibetac_norm          { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ibeta_inv            { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ibetac_inv           { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ibeta_inva           { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ibetac_inva          { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ibeta_invb           { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ibetac_invb          { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ibeta_dif            { static dynamic::function_name eval(); };

    struct MATCL_SF_EXPORT airy_ai              { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT airy_bi              { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT airy_ai_dif          { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT airy_bi_dif          { static dynamic::function_name eval(); };

    struct MATCL_SF_EXPORT cyl_bessel_j         { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT cyl_neumann          { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT cyl_bessel_j_zero    { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT cyl_neumann_zero     { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT cyl_bessel_j_dif     { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT cyl_neumann_dif      { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT cyl_bessel_i         { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT cyl_bessel_k         { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT cyl_bessel_i_dif     { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT cyl_bessel_k_dif     { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT sph_bessel           { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT sph_neumann          { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT sph_bessel_dif       { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT sph_neumann_dif      { static dynamic::function_name eval(); };

    struct MATCL_SF_EXPORT cyl_hankel_1         { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT cyl_hankel_2         { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT sph_hankel_1         { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT sph_hankel_2         { static dynamic::function_name eval(); };

    struct MATCL_SF_EXPORT ellint_rf            { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ellint_rd            { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ellint_rj            { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ellint_rc            { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ellint_rg            { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ellint_1             { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ellint_2             { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ellint_3             { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT ellint_d             { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT jacobi_sn            { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT jacobi_cn            { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT jacobi_dn            { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT jacobi_elliptic      { static dynamic::function_name eval(); };

    struct MATCL_SF_EXPORT expint               { static dynamic::function_name eval(); };

    struct MATCL_SF_EXPORT owens_t              { static dynamic::function_name eval(); };

    struct MATCL_SF_EXPORT sinc                 { static dynamic::function_name eval(); };
    struct MATCL_SF_EXPORT sinhc                { static dynamic::function_name eval(); };
};

namespace object_func
{            
    MATCL_SF_EXPORT Object  erf(const Object& A);
    MATCL_SF_EXPORT Object  erfc(const Object& A);
    MATCL_SF_EXPORT Object  erf_inv(const Object& A);
    MATCL_SF_EXPORT Object  erfc_inv(const Object& A);

    MATCL_SF_EXPORT Object  gamma(const Object& A);
    MATCL_SF_EXPORT Object  gammaln(const Object& A);    
    MATCL_SF_EXPORT Object  gamma1pm1(const Object& A);
    MATCL_SF_EXPORT Object  digamma(const Object& A);
    MATCL_SF_EXPORT Object  trigamma(const Object& A);
    MATCL_SF_EXPORT Object  polygamma(const Object& A, Integer n);
    MATCL_SF_EXPORT Object  gamma_ratio(const Object& A, const Object& B);    
    MATCL_SF_EXPORT Object  gamma_delta_ratio(const Object& A, const Object& B);    
    MATCL_SF_EXPORT Object  igamma_lower(const Object& A, const Object& B);
    MATCL_SF_EXPORT Object  igamma_upper(const Object& A, const Object& B);
    MATCL_SF_EXPORT Object  igamma_lower_norm(const Object& A, const Object& B);
    MATCL_SF_EXPORT Object  igamma_upper_norm(const Object& A, const Object& B);
    MATCL_SF_EXPORT Object  igamma_lower_inv(const Object& A, const Object& B);
    MATCL_SF_EXPORT Object  igamma_upper_inv(const Object& A, const Object& B);
    MATCL_SF_EXPORT Object  igamma_lower_inva(const Object& A, const Object& B);
    MATCL_SF_EXPORT Object  igamma_upper_inva(const Object& A, const Object& B);
    MATCL_SF_EXPORT Object  igamma_lower_dif(const Object& A, const Object& B);
    
    MATCL_SF_EXPORT Object  zeta(const Object& A);

    MATCL_SF_EXPORT Object  beta(const Object& A, const Object& B);
    MATCL_SF_EXPORT Object  ibeta(const Object& A, const Object& B, const Object& C);
    MATCL_SF_EXPORT Object  ibetac(const Object& A, const Object& B, const Object& C);
    MATCL_SF_EXPORT Object  ibeta_norm(const Object& A, const Object& B, const Object& C);
    MATCL_SF_EXPORT Object  ibetac_norm(const Object& A, const Object& B, const Object& C);
    MATCL_SF_EXPORT Object  ibeta_inv(const Object& A, const Object& B, const Object& C);
    MATCL_SF_EXPORT Object  ibetac_inv(const Object& A, const Object& B, const Object& C);
    MATCL_SF_EXPORT Object  ibeta_inva(const Object& A, const Object& B, const Object& C);
    MATCL_SF_EXPORT Object  ibetac_inva(const Object& A, const Object& B, const Object& C);
    MATCL_SF_EXPORT Object  ibeta_invb(const Object& A, const Object& B, const Object& C);
    MATCL_SF_EXPORT Object  ibetac_invb(const Object& A, const Object& B, const Object& C);
    MATCL_SF_EXPORT Object  ibeta_dif(const Object& A, const Object& B, const Object& C);

    MATCL_SF_EXPORT Object  airy_ai(const Object& A);
    MATCL_SF_EXPORT Object  airy_bi(const Object& A);
    MATCL_SF_EXPORT Object  airy_ai_dif(const Object& A);
    MATCL_SF_EXPORT Object  airy_bi_dif(const Object& A);

    MATCL_SF_EXPORT Object  cyl_bessel_j(const Object& v, const Object& x);
    MATCL_SF_EXPORT Object  cyl_neumann(const Object& v, const Object& x);
    MATCL_SF_EXPORT Object  cyl_bessel_j_zero(const Object& v, Integer m);
    MATCL_SF_EXPORT Object  cyl_neumann_zero(const Object& v, Integer m);
    MATCL_SF_EXPORT Object  cyl_bessel_j_dif(const Object& v, const Object& x);
    MATCL_SF_EXPORT Object  cyl_neumann_dif(const Object& v, const Object& x);
    MATCL_SF_EXPORT Object  cyl_bessel_i(const Object& v, const Object& x);
    MATCL_SF_EXPORT Object  cyl_bessel_k(const Object& v, const Object& x);
    MATCL_SF_EXPORT Object  cyl_bessel_i_dif(const Object& v, const Object& x);
    MATCL_SF_EXPORT Object  cyl_bessel_k_dif(const Object& v, const Object& x);
    MATCL_SF_EXPORT Object  sph_bessel(Integer n, const Object& x);
    MATCL_SF_EXPORT Object  sph_neumann(Integer n, const Object& x);
    MATCL_SF_EXPORT Object  sph_bessel_dif(Integer n, const Object& x);
    MATCL_SF_EXPORT Object  sph_neumann_dif(Integer n, const Object& x);

    MATCL_SF_EXPORT Object  cyl_hankel_1(const Object& v, const Object& x);
    MATCL_SF_EXPORT Object  cyl_hankel_2(const Object& v, const Object& x);
    MATCL_SF_EXPORT Object  sph_hankel_1(const Object& v, const Object& x);
    MATCL_SF_EXPORT Object  sph_hankel_2(const Object& v, const Object& x);

    MATCL_SF_EXPORT Object  ellint_rf(const Object& x, const Object& y, const Object& z);
    MATCL_SF_EXPORT Object  ellint_rd(const Object& x, const Object& y, const Object& z);
    MATCL_SF_EXPORT Object  ellint_rj(const Object& x, const Object& y, const Object& z, const Object& p);
    MATCL_SF_EXPORT Object  ellint_rc(const Object& x, const Object& y);
    MATCL_SF_EXPORT Object  ellint_rg(const Object& x, const Object& y, const Object& z);
    MATCL_SF_EXPORT Object  ellint_1(const Object& k, const Object& phi);
    MATCL_SF_EXPORT Object  ellint_1(const Object& k);
    MATCL_SF_EXPORT Object  ellint_2(const Object& k, const Object& phi);
    MATCL_SF_EXPORT Object  ellint_2(const Object& k);
    MATCL_SF_EXPORT Object  ellint_3(const Object& k, const Object& n, const Object& phi);
    MATCL_SF_EXPORT Object  ellint_3(const Object& k, const Object& n);
    MATCL_SF_EXPORT Object  ellint_d(const Object& k, const Object& phi);
    MATCL_SF_EXPORT Object  ellint_d(const Object& k);
    MATCL_SF_EXPORT Object  jacobi_sn(const Object& k, const Object& u);
    MATCL_SF_EXPORT Object  jacobi_cn(const Object& k, const Object& u);
    MATCL_SF_EXPORT Object  jacobi_dn(const Object& k, const Object& u);

    MATCL_SF_EXPORT Object  expint(const Object& k, Integer n);
    MATCL_SF_EXPORT Object  expint(const Object& k);

    MATCL_SF_EXPORT Object  owens_t(const Object& h, const Object& a);

    MATCL_SF_EXPORT Object  sinc(const Object& x);
    MATCL_SF_EXPORT Object  sinhc(const Object& x);

    // TODO
    //MATCL_SF_EXPORT Object  gammaln(const Object& A, int& sign);
    //MATCL_SF_EXPORT void    jacobi_elliptic(const Object& k, const Object& u, 
    //                                        Object& sn, Object& cn, Object& dn);

};}
