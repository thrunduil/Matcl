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

#include "matcl-specfunc/lib_functions/elliptic.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-specfunc/objects/object_func.h"

#pragma warning(push)
#pragma warning(disable:4127)   //conditional expression is constant
#include <boost/math/special_functions.hpp>
#pragma warning(pop)

namespace matcl { namespace details
{

namespace md = matcl::details;

namespace ignore_errors 
{
    using boost__ignore_policy
        = ::boost::math::policies::policy 
          <
            ::boost::math::policies::domain_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::pole_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::overflow_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::underflow_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::rounding_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::denorm_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::indeterminate_result_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::evaluation_error<::boost::math::policies::ignore_error>
         >;
    BOOST_MATH_DECLARE_SPECIAL_FUNCTIONS(boost__ignore_policy)
}

//-------------------------------------------------------------------------------
//                                  REAL
//-------------------------------------------------------------------------------
template<>
typename elliptic_helper<Real>::return_type 
elliptic_helper<Real>::eval_ellint_rf(const Real& x, const Real& y, const Real& z)
{
    return ignore_errors::ellint_rf(x,y,z);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_ellint_rd(const Real& x, const Real& y, const Real& z)
{
    return ignore_errors::ellint_rd(x,y,z);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_ellint_rj(const Real& x, const Real& y, const Real& z,
                                      const Real& p)
{
    return ignore_errors::ellint_rj(x,y,z,p);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_ellint_rc(const Real& x, const Real& y)
{
    return ignore_errors::ellint_rc(x,y);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_ellint_rg(const Real& x, const Real& y, const Real& z)
{
    return ignore_errors::ellint_rg(x,y,z);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_ellint_1(const Real& k, const Real& phi)
{
    return ignore_errors::ellint_1(k,phi);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_ellint_1(const Real& k)
{
    return ignore_errors::ellint_1(k);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_ellint_2(const Real& k, const Real& phi)
{
    return ignore_errors::ellint_2(k,phi);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_ellint_2(const Real& k)
{
    return ignore_errors::ellint_2(k);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_ellint_3(const Real& k, const Real& n, const Real& phi)
{
    return ignore_errors::ellint_3(k,n,phi);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_ellint_3(const Real& k, const Real& n)
{
    return ignore_errors::ellint_3(k,n);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_ellint_d(const Real& k, const Real& phi)
{
    return ignore_errors::ellint_d(k,phi);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_ellint_d(const Real& k)
{
    return ignore_errors::ellint_d(k);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_jacobi_sn(const Real& k, const Real& u)
{
    return ignore_errors::jacobi_sn(k,u);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_jacobi_cn(const Real& k, const Real& u)
{
    return ignore_errors::jacobi_cn(k,u);
};

template<>
typename elliptic_helper<Real>::return_type
elliptic_helper<Real>::eval_jacobi_dn(const Real& k, const Real& u)
{
    return ignore_errors::jacobi_dn(k,u);
};

template<>
void elliptic_helper<Real>::eval_jacobi_elliptic(const Real& k, const Real& u, 
                                            Real& sn, Real& cn, Real& dn)
{
    sn = ignore_errors::jacobi_elliptic(k,u,&cn,&dn);
};

//-------------------------------------------------------------------------------
//                              FLOAT
//-------------------------------------------------------------------------------
template<>
typename elliptic_helper<Float>::return_type 
elliptic_helper<Float>::eval_ellint_rf(const Float& x, const Float& y, const Float& z)
{
    return ignore_errors::ellint_rf(x,y,z);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_ellint_rd(const Float& x, const Float& y, const Float& z)
{
    return ignore_errors::ellint_rd(x,y,z);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_ellint_rj(const Float& x, const Float& y, const Float& z,
                                       const Float& p)
{
    return ignore_errors::ellint_rj(x,y,z,p);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_ellint_rc(const Float& x, const Float& y)
{
    return ignore_errors::ellint_rc(x,y);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_ellint_rg(const Float& x, const Float& y, const Float& z)
{
    return ignore_errors::ellint_rg(x,y,z);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_ellint_1(const Float& k, const Float& phi)
{
    return ignore_errors::ellint_1(k,phi);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_ellint_1(const Float& k)
{
    return ignore_errors::ellint_1(k);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_ellint_2(const Float& k, const Float& phi)
{
    return ignore_errors::ellint_2(k,phi);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_ellint_2(const Float& k)
{
    return ignore_errors::ellint_2(k);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_ellint_3(const Float& k, const Float& n, const Float& phi)
{
    return ignore_errors::ellint_3(k,n,phi);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_ellint_3(const Float& k, const Float& n)
{
    return ignore_errors::ellint_3(k,n);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_ellint_d(const Float& k, const Float& phi)
{
    return ignore_errors::ellint_d(k,phi);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_ellint_d(const Float& k)
{
    return ignore_errors::ellint_d(k);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_jacobi_sn(const Float& k, const Float& u)
{
    return ignore_errors::jacobi_sn(k,u);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_jacobi_cn(const Float& k, const Float& u)
{
    return ignore_errors::jacobi_cn(k,u);
};

template<>
typename elliptic_helper<Float>::return_type
elliptic_helper<Float>::eval_jacobi_dn(const Float& k, const Float& u)
{
    return ignore_errors::jacobi_dn(k,u);
};

template<>
void elliptic_helper<Float>::eval_jacobi_elliptic(const Float& k, const Float& u, 
                                             Float& sn, Float& cn, Float& dn)
{
    sn = ignore_errors::jacobi_elliptic(k,u,&cn,&dn);
};

//-------------------------------------------------------------------------------
//                              COMPLEX
//-------------------------------------------------------------------------------
template<>
typename elliptic_helper<Complex>::return_type 
elliptic_helper<Complex>::eval_ellint_rf(const Complex& x, const Complex& y, const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return elliptic_helper<Real>::eval_ellint_rf(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_rf"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_ellint_rd(const Complex& x, const Complex& y, const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return elliptic_helper<Real>::eval_ellint_rd(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_rd"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_ellint_rj(const Complex& x, const Complex& y, 
                                         const Complex& z, const Complex& p)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0 && imag(p) == 0.0)
        return elliptic_helper<Real>::eval_ellint_rj(real(x),real(y),real(z),real(p));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_rj"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_ellint_rc(const Complex& x, const Complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return elliptic_helper<Real>::eval_ellint_rc(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_rc"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_ellint_rg(const Complex& x, const Complex& y, const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return elliptic_helper<Real>::eval_ellint_rg(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_rg"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_ellint_1(const Complex& k, const Complex& phi)
{
    if (imag(k) == 0.0 && imag(phi))
        return elliptic_helper<Real>::eval_ellint_1(real(k),real(phi));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_1"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_ellint_1(const Complex& k)
{
    if (imag(k) == 0.0)
        return elliptic_helper<Real>::eval_ellint_1(real(k));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_1"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_ellint_2(const Complex& k, const Complex& phi)
{
    if (imag(k) == 0.0 && imag(phi))
        return elliptic_helper<Real>::eval_ellint_2(real(k),real(phi));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_2"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_ellint_2(const Complex& k)
{
    if (imag(k) == 0.0)
        return elliptic_helper<Real>::eval_ellint_2(real(k));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_2"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_ellint_3(const Complex& k, const Complex& n,
                                        const Complex& phi)
{
    if (imag(k) == 0.0 && imag(phi) && imag(n))
        return elliptic_helper<Real>::eval_ellint_3(real(k),real(n),real(phi));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_3"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_ellint_3(const Complex& k, const Complex& n)
{
    if (imag(k) == 0.0 && imag(n))
        return elliptic_helper<Real>::eval_ellint_3(real(k),real(n));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_3"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_ellint_d(const Complex& k, const Complex& phi)
{
    if (imag(k) == 0.0 && imag(phi))
        return elliptic_helper<Real>::eval_ellint_d(real(k),real(phi));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_d"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_ellint_d(const Complex& k)
{
    if (imag(k) == 0.0)
        return elliptic_helper<Real>::eval_ellint_d(real(k));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_d"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_jacobi_sn(const Complex& k, const Complex& u)
{
    if (imag(k) == 0.0 && imag(u))
        return elliptic_helper<Real>::eval_jacobi_sn(real(k),real(u));
    else
        throw matcl::error::function_not_defined_for_complex("jacobi_sn"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_jacobi_cn(const Complex& k, const Complex& u)
{
    if (imag(k) == 0.0 && imag(u))
        return elliptic_helper<Real>::eval_jacobi_cn(real(k),real(u));
    else
        throw matcl::error::function_not_defined_for_complex("jacobi_cn"); 
};

template<>
typename elliptic_helper<Complex>::return_type
elliptic_helper<Complex>::eval_jacobi_dn(const Complex& k, const Complex& u)
{
    if (imag(k) == 0.0 && imag(u))
        return elliptic_helper<Real>::eval_jacobi_dn(real(k),real(u));
    else
        throw matcl::error::function_not_defined_for_complex("jacobi_dn"); 
};

template<>
void elliptic_helper<Complex>::eval_jacobi_elliptic(const Complex& k, const Complex& u, 
                                               Complex& sn, Complex& cn, Complex& dn)
{
    if (imag(k) == 0.0 && imag(u))
    {
        Real ret_sn, ret_cn, ret_dn;
        elliptic_helper<Real>::eval_jacobi_elliptic(real(k),real(u),ret_sn,ret_cn,ret_dn);
        sn = ret_sn;
        cn = ret_cn;
        dn = ret_dn;
        return;
    }
    else
    {
        throw matcl::error::function_not_defined_for_complex("jacobi_elliptic"); 
    };
};

//-------------------------------------------------------------------------------
//                              FLOAT_COMPLEX
//-------------------------------------------------------------------------------
template<>
typename elliptic_helper<Float_complex>::return_type 
elliptic_helper<Float_complex>::eval_ellint_rf(const Float_complex& x, const Float_complex& y,
                                               const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return elliptic_helper<Float>::eval_ellint_rf(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_rf"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_ellint_rd(const Float_complex& x, const Float_complex& y,
                                               const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return elliptic_helper<Float>::eval_ellint_rd(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_rd"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_ellint_rj(const Float_complex& x, const Float_complex& y, 
                                         const Float_complex& z, const Float_complex& p)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0 && imag(p) == 0.0)
        return elliptic_helper<Float>::eval_ellint_rj(real(x),real(y),real(z),real(p));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_rj"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_ellint_rc(const Float_complex& x, const Float_complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return elliptic_helper<Float>::eval_ellint_rc(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_rc"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_ellint_rg(const Float_complex& x, const Float_complex& y,
                                               const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return elliptic_helper<Float>::eval_ellint_rg(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_rg"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_ellint_1(const Float_complex& k, const Float_complex& phi)
{
    if (imag(k) == 0.0 && imag(phi))
        return elliptic_helper<Float>::eval_ellint_1(real(k),real(phi));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_1"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_ellint_1(const Float_complex& k)
{
    if (imag(k) == 0.0)
        return elliptic_helper<Float>::eval_ellint_1(real(k));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_1"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_ellint_2(const Float_complex& k, const Float_complex& phi)
{
    if (imag(k) == 0.0 && imag(phi))
        return elliptic_helper<Float>::eval_ellint_2(real(k),real(phi));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_2"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_ellint_2(const Float_complex& k)
{
    if (imag(k) == 0.0)
        return elliptic_helper<Float>::eval_ellint_2(real(k));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_2"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_ellint_3(const Float_complex& k, const Float_complex& n,
                                        const Float_complex& phi)
{
    if (imag(k) == 0.0 && imag(phi) && imag(n))
        return elliptic_helper<Float>::eval_ellint_3(real(k),real(n),real(phi));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_3"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_ellint_3(const Float_complex& k, const Float_complex& n)
{
    if (imag(k) == 0.0 && imag(n))
        return elliptic_helper<Float>::eval_ellint_3(real(k),real(n));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_3"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_ellint_d(const Float_complex& k, const Float_complex& phi)
{
    if (imag(k) == 0.0 && imag(phi))
        return elliptic_helper<Float>::eval_ellint_d(real(k),real(phi));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_d"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_ellint_d(const Float_complex& k)
{
    if (imag(k) == 0.0)
        return elliptic_helper<Float>::eval_ellint_d(real(k));
    else
        throw matcl::error::function_not_defined_for_complex("ellint_d"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_jacobi_sn(const Float_complex& k, const Float_complex& u)
{
    if (imag(k) == 0.0 && imag(u))
        return elliptic_helper<Float>::eval_jacobi_sn(real(k),real(u));
    else
        throw matcl::error::function_not_defined_for_complex("jacobi_sn"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_jacobi_cn(const Float_complex& k, const Float_complex& u)
{
    if (imag(k) == 0.0 && imag(u))
        return elliptic_helper<Float>::eval_jacobi_cn(real(k),real(u));
    else
        throw matcl::error::function_not_defined_for_complex("jacobi_cn"); 
};

template<>
typename elliptic_helper<Float_complex>::return_type
elliptic_helper<Float_complex>::eval_jacobi_dn(const Float_complex& k, const Float_complex& u)
{
    if (imag(k) == 0.0 && imag(u))
        return elliptic_helper<Float>::eval_jacobi_dn(real(k),real(u));
    else
        throw matcl::error::function_not_defined_for_complex("jacobi_dn"); 
};

template<>
void elliptic_helper<Float_complex>::eval_jacobi_elliptic(const Float_complex& k, const Float_complex& u, 
                                               Float_complex& sn, Float_complex& cn, Float_complex& dn)
{
    if (imag(k) == 0.0 && imag(u))
    {
        Float ret_sn, ret_cn, ret_dn;
        elliptic_helper<Float>::eval_jacobi_elliptic(real(k),real(u),ret_sn,ret_cn,ret_dn);
        sn = ret_sn;
        cn = ret_cn;
        dn = ret_dn;
        return;
    }
    else
    {
        throw matcl::error::function_not_defined_for_complex("jacobi_elliptic"); 
    };
};

//-------------------------------------------------------------------------------
//                              OBJECT
//-------------------------------------------------------------------------------
template<>
typename elliptic_helper<Object>::return_type 
elliptic_helper<Object>::eval_ellint_rf(const Object& x, const Object& y, const Object& z)
{
    return object_func::ellint_rf(x,y,z);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_ellint_rd(const Object& x, const Object& y, const Object& z)
{
    return object_func::ellint_rd(x,y,z);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_ellint_rj(const Object& x, const Object& y, const Object& z,
                                       const Object& p)
{
    return object_func::ellint_rj(x,y,z,p);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_ellint_rc(const Object& x, const Object& y)
{
    return object_func::ellint_rc(x,y);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_ellint_rg(const Object& x, const Object& y, const Object& z)
{
    return object_func::ellint_rg(x,y,z);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_ellint_1(const Object& k, const Object& phi)
{
    return object_func::ellint_1(k,phi);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_ellint_1(const Object& k)
{
    return object_func::ellint_1(k);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_ellint_2(const Object& k, const Object& phi)
{
    return object_func::ellint_2(k,phi);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_ellint_2(const Object& k)
{
    return object_func::ellint_2(k);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_ellint_3(const Object& k, const Object& n, const Object& phi)
{
    return object_func::ellint_3(k,n,phi);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_ellint_3(const Object& k, const Object& n)
{
    return object_func::ellint_3(k,n);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_ellint_d(const Object& k, const Object& phi)
{
    return object_func::ellint_d(k,phi);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_ellint_d(const Object& k)
{
    return object_func::ellint_d(k);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_jacobi_sn(const Object& k, const Object& u)
{
    return object_func::jacobi_sn(k,u);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_jacobi_cn(const Object& k, const Object& u)
{
    return object_func::jacobi_cn(k,u);
};

template<>
typename elliptic_helper<Object>::return_type
elliptic_helper<Object>::eval_jacobi_dn(const Object& k, const Object& u)
{
    return object_func::jacobi_dn(k,u);
};

template<>
void elliptic_helper<Object>::eval_jacobi_elliptic(const Object& k, const Object& u, 
                                             Object& sn, Object& cn, Object& dn)
{
    (void)k;
    (void)u;
    (void)sn;
    (void)cn;
    (void)dn;
    //TODO: impl for object
    throw matcl::error::object_value_type_not_allowed("jacobi_elliptic");
};

//template struct elliptic_helper<Integer>;
template struct elliptic_helper<Real>;
template struct elliptic_helper<Float>;
template struct elliptic_helper<Complex>;
template struct elliptic_helper<Float_complex>;
template struct elliptic_helper<Object>;

};};

namespace matcl
{

void test_ellipt_inst()
{
    Real x = 0.0;
    Real y = 0.0;
    Real z = 0.0;
    Real p = 0.0;
    Real k = 0.0;
    Real phi = 0.0;
    Real n = 0.0;
    Real u = 0.0;

    ellint_rf(x, y, z);
    ellint_rd(x, y, z);
    ellint_rj(x, y, z, p);
    ellint_rc(x, y);
    ellint_rg(x, y, z);
    ellint_1(k, phi);
    ellint_1(k);
    ellint_2(k, phi);
    ellint_2(k);
    ellint_3(k, n, phi);
    ellint_3(k, n);
    ellint_d(k, phi);
    ellint_d(k);
    jacobi_sn(k, u);
    jacobi_cn(k, u);
    jacobi_dn(k, u);

    Real sn, cn, dn;
    jacobi_elliptic(k, u, sn, cn, dn);
};

};
