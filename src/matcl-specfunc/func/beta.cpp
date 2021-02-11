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

#include "matcl-specfunc/lib_functions/beta.h"
#include "matcl-specfunc/lib_functions/gamma.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-specfunc/func/raw/complex_func.h"
#include "matcl-specfunc/objects/object_func.h"

#pragma warning(push)
#pragma warning(disable:4127)   //conditional expression is constant
#include <boost/math/special_functions.hpp>
#pragma warning(pop)

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace details
{

namespace md = matcl::details;
namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;

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
typename beta_helper<Real>::return_type beta_helper<Real>::eval_beta(const Real& x, const Real& y)
{
    return ignore_errors::beta(x,y);
};
template<>
typename beta_helper<Real>::return_type
beta_helper<Real>::eval_ibeta(const Real& x,const Real& y,const Real& z)
{
    return ignore_errors::beta(x,y,z);
};
template<>
typename beta_helper<Real>::return_type
beta_helper<Real>::eval_ibetac(const Real& x,const Real& y,const Real& z)
{
    return ignore_errors::betac(x,y,z);
};
template<>
typename beta_helper<Real>::return_type
beta_helper<Real>::eval_ibeta_norm(const Real& x,const Real& y,const Real& z)
{
    return ignore_errors::ibeta(x,y,z);
};
template<>
typename beta_helper<Real>::return_type
beta_helper<Real>::eval_ibetac_norm(const Real& x,const Real& y,const Real& z)
{
    return ignore_errors::ibetac(x,y,z);
};
template<>
typename beta_helper<Real>::return_type
beta_helper<Real>::eval_ibeta_inv(const Real& x,const Real& y,const Real& z)
{
    return ignore_errors::ibeta_inv(x,y,z);
};
template<>
typename beta_helper<Real>::return_type
beta_helper<Real>::eval_ibetac_inv(const Real& x,const Real& y,const Real& z)
{
    return ignore_errors::ibetac_inv(x,y,z);
};
template<>
typename beta_helper<Real>::return_type
beta_helper<Real>::eval_ibeta_inva(const Real& x,const Real& y,const Real& z)
{
    return ignore_errors::ibeta_inva(x,y,z);
};
template<>
typename beta_helper<Real>::return_type
beta_helper<Real>::eval_ibetac_inva(const Real& x,const Real& y,const Real& z)
{
    return ignore_errors::ibetac_inva(x,y,z);
};
template<>
typename beta_helper<Real>::return_type
beta_helper<Real>::eval_ibeta_invb(const Real& x,const Real& y,const Real& z)
{
    return ignore_errors::ibeta_invb(x,y,z);
};
template<>
typename beta_helper<Real>::return_type
beta_helper<Real>::eval_ibetac_invb(const Real& x,const Real& y,const Real& z)
{
    return ignore_errors::ibetac_invb(x,y,z);
};
template<>
typename beta_helper<Real>::return_type
beta_helper<Real>::eval_ibeta_dif(const Real& x,const Real& y,const Real& z)
{
    return ignore_errors::ibeta_derivative(x,y,z);
};

//-------------------------------------------------------------------------------
//                                  FLOAT
//-------------------------------------------------------------------------------
template<>
typename beta_helper<Float>::return_type beta_helper<Float>::eval_beta(const Float& x, const Float& y)
{
    return ignore_errors::beta(x,y);
};
template<>
typename beta_helper<Float>::return_type
beta_helper<Float>::eval_ibeta(const Float& x,const Float& y,const Float& z)
{
    return ignore_errors::beta(x,y,z);
};
template<>
typename beta_helper<Float>::return_type
beta_helper<Float>::eval_ibetac(const Float& x,const Float& y,const Float& z)
{
    return ignore_errors::betac(x,y,z);
};
template<>
typename beta_helper<Float>::return_type
beta_helper<Float>::eval_ibeta_norm(const Float& x,const Float& y,const Float& z)
{
    return ignore_errors::ibeta(x,y,z);
};
template<>
typename beta_helper<Float>::return_type
beta_helper<Float>::eval_ibetac_norm(const Float& x,const Float& y,const Float& z)
{
    return ignore_errors::ibetac(x,y,z);
};
template<>
typename beta_helper<Float>::return_type
beta_helper<Float>::eval_ibeta_inv(const Float& x,const Float& y,const Float& z)
{
    return ignore_errors::ibeta_inv(x,y,z);
};
template<>
typename beta_helper<Float>::return_type
beta_helper<Float>::eval_ibetac_inv(const Float& x,const Float& y,const Float& z)
{
    return ignore_errors::ibetac_inv(x,y,z);
};
template<>
typename beta_helper<Float>::return_type
beta_helper<Float>::eval_ibeta_inva(const Float& x,const Float& y,const Float& z)
{
    return ignore_errors::ibeta_inva(x,y,z);
};
template<>
typename beta_helper<Float>::return_type
beta_helper<Float>::eval_ibetac_inva(const Float& x,const Float& y,const Float& z)
{
    return ignore_errors::ibetac_inva(x,y,z);
};
template<>
typename beta_helper<Float>::return_type
beta_helper<Float>::eval_ibeta_invb(const Float& x,const Float& y,const Float& z)
{
    return ignore_errors::ibeta_invb(x,y,z);
};
template<>
typename beta_helper<Float>::return_type
beta_helper<Float>::eval_ibetac_invb(const Float& x,const Float& y,const Float& z)
{
    return ignore_errors::ibetac_invb(x,y,z);
};
template<>
typename beta_helper<Float>::return_type
beta_helper<Float>::eval_ibeta_dif(const Float& x,const Float& y,const Float& z)
{
    return ignore_errors::ibeta_derivative(x,y,z);
};

//-------------------------------------------------------------------------------
//                                  COMPLEX
//-------------------------------------------------------------------------------
template<>
typename beta_helper<Complex>::return_type 
beta_helper<Complex>::eval_beta(const Complex& x, const Complex& y)
{
    return exp( gammaln(x) + gammaln(y) - gammaln(x+y));
};
template<>
typename beta_helper<Complex>::return_type
beta_helper<Complex>::eval_ibeta(const Complex& x,const Complex& y,const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Real>::eval_ibeta(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibeta"); 
};
template<>
typename beta_helper<Complex>::return_type
beta_helper<Complex>::eval_ibetac(const Complex& x,const Complex& y,const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Real>::eval_ibetac(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibetac"); 
};
template<>
typename beta_helper<Complex>::return_type
beta_helper<Complex>::eval_ibeta_norm(const Complex& x,const Complex& y,const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Real>::eval_ibeta_norm(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibeta_norm"); 
};
template<>
typename beta_helper<Complex>::return_type
beta_helper<Complex>::eval_ibetac_norm(const Complex& x,const Complex& y,const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Real>::eval_ibetac_norm(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibetac_norm"); 
};
template<>
typename beta_helper<Complex>::return_type
beta_helper<Complex>::eval_ibeta_inv(const Complex& x,const Complex& y,const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Real>::eval_ibeta_inv(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibeta_inv"); 
};
template<>
typename beta_helper<Complex>::return_type
beta_helper<Complex>::eval_ibetac_inv(const Complex& x,const Complex& y,const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Real>::eval_ibetac_inv(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibetac_inv"); 
};
template<>
typename beta_helper<Complex>::return_type
beta_helper<Complex>::eval_ibeta_inva(const Complex& x,const Complex& y,const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Real>::eval_ibeta_inva(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibeta_inva"); 
};
template<>
typename beta_helper<Complex>::return_type
beta_helper<Complex>::eval_ibetac_inva(const Complex& x,const Complex& y,const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Real>::eval_ibetac_inva(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibetac_inva"); 
};
template<>
typename beta_helper<Complex>::return_type
beta_helper<Complex>::eval_ibeta_invb(const Complex& x,const Complex& y,const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Real>::eval_ibeta_invb(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibeta_invb"); 
};
template<>
typename beta_helper<Complex>::return_type
beta_helper<Complex>::eval_ibetac_invb(const Complex& x,const Complex& y,const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Real>::eval_ibetac_invb(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibetac_invb"); 
};
template<>
typename beta_helper<Complex>::return_type
beta_helper<Complex>::eval_ibeta_dif(const Complex& x,const Complex& y,const Complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Real>::eval_ibeta_dif(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibeta_dif"); 
};

//-------------------------------------------------------------------------------
//                                  FLOAT_COMPLEX
//-------------------------------------------------------------------------------
template<>
typename beta_helper<Float_complex>::return_type 
beta_helper<Float_complex>::eval_beta(const Float_complex& x, const Float_complex& y)
{
    return exp( gammaln(x) + gammaln(y) - gammaln(x+y));
};

template<>
typename beta_helper<Float_complex>::return_type
beta_helper<Float_complex>::eval_ibeta(const Float_complex& x,const Float_complex& y,
                                       const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Float>::eval_ibeta(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibeta"); 
};
template<>
typename beta_helper<Float_complex>::return_type
beta_helper<Float_complex>::eval_ibetac(const Float_complex& x,const Float_complex& y,
                                        const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Float>::eval_ibetac(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibetac"); 
};
template<>
typename beta_helper<Float_complex>::return_type
beta_helper<Float_complex>::eval_ibeta_norm(const Float_complex& x,const Float_complex& y,
                                            const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Float>::eval_ibeta_norm(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibeta_norm"); 
};
template<>
typename beta_helper<Float_complex>::return_type
beta_helper<Float_complex>::eval_ibetac_norm(const Float_complex& x,const Float_complex& y,
                                             const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Float>::eval_ibetac_norm(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibetac_norm"); 
};
template<>
typename beta_helper<Float_complex>::return_type
beta_helper<Float_complex>::eval_ibeta_inv(const Float_complex& x,const Float_complex& y,
                                           const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Float>::eval_ibeta_inv(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibeta_inv"); 
};
template<>
typename beta_helper<Float_complex>::return_type
beta_helper<Float_complex>::eval_ibetac_inv(const Float_complex& x,const Float_complex& y,
                                            const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Float>::eval_ibetac_inv(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibetac_inv"); 
};
template<>
typename beta_helper<Float_complex>::return_type
beta_helper<Float_complex>::eval_ibeta_inva(const Float_complex& x,const Float_complex& y,
                                            const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Float>::eval_ibeta_inva(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibeta_inva"); 
};
template<>
typename beta_helper<Float_complex>::return_type
beta_helper<Float_complex>::eval_ibetac_inva(const Float_complex& x,const Float_complex& y,
                                             const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Float>::eval_ibetac_inva(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibetac_inva"); 
};
template<>
typename beta_helper<Float_complex>::return_type
beta_helper<Float_complex>::eval_ibeta_invb(const Float_complex& x,const Float_complex& y,
                                            const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Float>::eval_ibeta_invb(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibeta_invb"); 
};
template<>
typename beta_helper<Float_complex>::return_type
beta_helper<Float_complex>::eval_ibetac_invb(const Float_complex& x,const Float_complex& y,
                                             const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Float>::eval_ibetac_invb(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibetac_invb"); 
};
template<>
typename beta_helper<Float_complex>::return_type
beta_helper<Float_complex>::eval_ibeta_dif(const Float_complex& x,const Float_complex& y,
                                           const Float_complex& z)
{
    if (imag(x) == 0.0 && imag(y) == 0.0 && imag(z) == 0.0)
        return beta_helper<Float>::eval_ibeta_dif(real(x),real(y),real(z));
    else
        throw matcl::error::function_not_defined_for_complex("ibeta_dif"); 
};

//-------------------------------------------------------------------------------
//                                  INTEGER
//-------------------------------------------------------------------------------
template<>
typename beta_helper<Integer>::return_type
beta_helper<Integer>::eval_beta(const Integer& x,const Integer& y)
{
    return beta_helper<Real>::eval_beta((Real)x, (Real)y);
};
template<>
typename beta_helper<Integer>::return_type
beta_helper<Integer>::eval_ibeta(const Integer& x,const Integer& y,const Integer& z)
{
    return beta_helper<Real>::eval_ibeta((Real)x, (Real)y, (Real)z);
};
template<>
typename beta_helper<Integer>::return_type
beta_helper<Integer>::eval_ibetac(const Integer& x,const Integer& y,const Integer& z)
{
    return beta_helper<Real>::eval_ibetac((Real)x, (Real)y, (Real)z);
};
template<>
typename beta_helper<Integer>::return_type
beta_helper<Integer>::eval_ibeta_norm(const Integer& x,const Integer& y,const Integer& z)
{
    return beta_helper<Real>::eval_ibeta_norm((Real)x, (Real)y, (Real)z);
};
template<>
typename beta_helper<Integer>::return_type
beta_helper<Integer>::eval_ibetac_norm(const Integer& x,const Integer& y,const Integer& z)
{
    return beta_helper<Real>::eval_ibetac_norm((Real)x, (Real)y, (Real)z);
};
template<>
typename beta_helper<Integer>::return_type
beta_helper<Integer>::eval_ibeta_inv(const Integer& x,const Integer& y,const Integer& z)
{
    return beta_helper<Real>::eval_ibeta_inv((Real)x, (Real)y, (Real)z);
};
template<>
typename beta_helper<Integer>::return_type
beta_helper<Integer>::eval_ibetac_inv(const Integer& x,const Integer& y,const Integer& z)
{
    return beta_helper<Real>::eval_ibetac_inv((Real)x, (Real)y, (Real)z);
};
template<>
typename beta_helper<Integer>::return_type
beta_helper<Integer>::eval_ibeta_inva(const Integer& x,const Integer& y,const Integer& z)
{
    return beta_helper<Real>::eval_ibeta_inva((Real)x, (Real)y, (Real)z);
};
template<>
typename beta_helper<Integer>::return_type
beta_helper<Integer>::eval_ibetac_inva(const Integer& x,const Integer& y,const Integer& z)
{
    return beta_helper<Real>::eval_ibetac_inva((Real)x, (Real)y, (Real)z);
};
template<>
typename beta_helper<Integer>::return_type
beta_helper<Integer>::eval_ibeta_invb(const Integer& x,const Integer& y,const Integer& z)
{
    return beta_helper<Real>::eval_ibeta_invb((Real)x, (Real)y, (Real)z);
};
template<>
typename beta_helper<Integer>::return_type
beta_helper<Integer>::eval_ibetac_invb(const Integer& x,const Integer& y,const Integer& z)
{
    return beta_helper<Real>::eval_ibetac_invb((Real)x, (Real)y, (Real)z);
};
template<>
typename beta_helper<Integer>::return_type
beta_helper<Integer>::eval_ibeta_dif(const Integer& x,const Integer& y,const Integer& z)
{
    return beta_helper<Real>::eval_ibeta_dif((Real)x, (Real)y, (Real)z);
};

//-------------------------------------------------------------------------------
//                                  OBJECT
//-------------------------------------------------------------------------------
template<>
typename beta_helper<Object>::return_type 
beta_helper<Object>::eval_beta(const Object& x, const Object& y)
{
    return object_func::beta(x,y);
};
template<>
typename beta_helper<Object>::return_type
beta_helper<Object>::eval_ibeta(const Object& x,const Object& y,const Object& z)
{
    return object_func::ibeta(x,y,z);
};
template<>
typename beta_helper<Object>::return_type
beta_helper<Object>::eval_ibetac(const Object& x,const Object& y,const Object& z)
{
    return object_func::ibetac(x,y,z);
};
template<>
typename beta_helper<Object>::return_type
beta_helper<Object>::eval_ibeta_norm(const Object& x,const Object& y,const Object& z)
{
    return object_func::ibeta_norm(x,y,z);
};
template<>
typename beta_helper<Object>::return_type
beta_helper<Object>::eval_ibetac_norm(const Object& x,const Object& y,const Object& z)
{
    return object_func::ibetac_norm(x,y,z);
};
template<>
typename beta_helper<Object>::return_type
beta_helper<Object>::eval_ibeta_inv(const Object& x,const Object& y,const Object& z)
{
    return object_func::ibeta_inv(x,y,z);
};
template<>
typename beta_helper<Object>::return_type
beta_helper<Object>::eval_ibetac_inv(const Object& x,const Object& y,const Object& z)
{
    return object_func::ibetac_inv(x,y,z);
};
template<>
typename beta_helper<Object>::return_type
beta_helper<Object>::eval_ibeta_inva(const Object& x,const Object& y,const Object& z)
{
    return object_func::ibeta_inva(x,y,z);
};
template<>
typename beta_helper<Object>::return_type
beta_helper<Object>::eval_ibetac_inva(const Object& x,const Object& y,const Object& z)
{
    return object_func::ibetac_inva(x,y,z);
};
template<>
typename beta_helper<Object>::return_type
beta_helper<Object>::eval_ibeta_invb(const Object& x,const Object& y,const Object& z)
{
    return object_func::ibeta_invb(x,y,z);
};
template<>
typename beta_helper<Object>::return_type
beta_helper<Object>::eval_ibetac_invb(const Object& x,const Object& y,const Object& z)
{
    return object_func::ibetac_invb(x,y,z);
};
template<>
typename beta_helper<Object>::return_type
beta_helper<Object>::eval_ibeta_dif(const Object& x,const Object& y,const Object& z)
{
    return object_func::ibeta_dif(x,y,z);
};

template struct beta_helper<Integer>;
template struct beta_helper<Real>;
template struct beta_helper<Float>;
template struct beta_helper<Complex>;
template struct beta_helper<Float_complex>;
//template struct beta_helper<Object>;

//-------------------------------------------------------------------------------
//                                  EVAL
//-------------------------------------------------------------------------------
struct eval_beta
{
    template<class T1, class T2>
    static auto eval(const T1& A, const T2& B)
                    -> typename md::unify_types2<T1,T2,Float>::type
    {
        using ret = typename md::unify_types2<T1,T2,Float>::type;
        return md::beta_helper<ret>::eval_beta(ret(A),ret(B));
    };
};

}}

matcl::Matrix matcl::beta(const Matrix &A,const Matrix &B)
{
    return eval_binary_func(A,B,details::eval_beta());
};

matcl::Matrix matcl::beta(Matrix&& A, const Matrix& B)
{
    return eval_binary_func(std::move(A),B,details::eval_beta());
};
matcl::Matrix matcl::beta(Matrix&& A, Matrix&& B)
{
    return eval_binary_func(std::move(A),std::move(B),details::eval_beta());
};
matcl::Matrix matcl::beta(const Matrix& A, Matrix&& B)
{
    return eval_binary_func(A,std::move(B),details::eval_beta());
};


#pragma warning( pop )