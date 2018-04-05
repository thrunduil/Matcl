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

#include "matcl-specfunc/lib_functions/sin_cardinal.h"
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
//                          REAL
//-------------------------------------------------------------------------------
template<>
typename sin_cardinal_helper<Real>::return_type 
sin_cardinal_helper<Real>::eval_sinc(const Real& x)
{
    return return_type(ignore_errors::sinc_pi(x));
};

template<>
typename sin_cardinal_helper<Real>::return_type 
sin_cardinal_helper<Real>::eval_sinhc(const Real& x)
{
    return return_type(ignore_errors::sinhc_pi(x));
}

//-------------------------------------------------------------------------------
//                          FLOAT
//-------------------------------------------------------------------------------
template<>
typename sin_cardinal_helper<Float>::return_type 
sin_cardinal_helper<Float>::eval_sinc(const Float& x)
{
    return return_type(ignore_errors::sinc_pi(x));
};

template<>
typename sin_cardinal_helper<Float>::return_type 
sin_cardinal_helper<Float>::eval_sinhc(const Float& x)
{
    return return_type(ignore_errors::sinhc_pi(x));
}

//-------------------------------------------------------------------------------
//                          COMPLEX
//-------------------------------------------------------------------------------
template<>
typename sin_cardinal_helper<Complex>::return_type 
sin_cardinal_helper<Complex>::eval_sinc(const Complex& x)
{
    return return_type(ignore_errors::sinc_pi(x.value));
};

template<>
typename sin_cardinal_helper<Complex>::return_type 
sin_cardinal_helper<Complex>::eval_sinhc(const Complex& x)
{
    return return_type(ignore_errors::sinhc_pi(x.value));
}

//-------------------------------------------------------------------------------
//                          FLOAT_COMPLEX
//-------------------------------------------------------------------------------
template<>
typename sin_cardinal_helper<Float_complex>::return_type 
sin_cardinal_helper<Float_complex>::eval_sinc(const Float_complex& x)
{
    return return_type(ignore_errors::sinc_pi(x.value));
};

template<>
typename sin_cardinal_helper<Float_complex>::return_type 
sin_cardinal_helper<Float_complex>::eval_sinhc(const Float_complex& x)
{
    return return_type(ignore_errors::sinhc_pi(x.value));
}

//-------------------------------------------------------------------------------
//                          INTEGER
//-------------------------------------------------------------------------------
template<>
typename sin_cardinal_helper<Integer>::return_type 
sin_cardinal_helper<Integer>::eval_sinc(const Integer& x)
{
    return sin_cardinal_helper<Real>::eval_sinc(Real(x));
};

template<>
typename sin_cardinal_helper<Integer>::return_type 
sin_cardinal_helper<Integer>::eval_sinhc(const Integer& x)
{
    return sin_cardinal_helper<Real>::eval_sinhc(Real(x));
}

//-------------------------------------------------------------------------------
//                          OBJECT
//-------------------------------------------------------------------------------
template<>
typename sin_cardinal_helper<Object>::return_type 
sin_cardinal_helper<Object>::eval_sinc(const Object& x)
{
    return object_func::sinc(x);
};

template<>
typename sin_cardinal_helper<Object>::return_type 
sin_cardinal_helper<Object>::eval_sinhc(const Object& x)
{
    return object_func::sinhc(x);
}

template struct sin_cardinal_helper<Real>;
template struct sin_cardinal_helper<Float>;
template struct sin_cardinal_helper<Complex>;
template struct sin_cardinal_helper<Float_complex>;
template struct sin_cardinal_helper<Object>;
template struct sin_cardinal_helper<Integer>;

struct func_sinc
{
    template<class T>
    auto eval(const T& v) const -> decltype(sin_cardinal_helper<T>::eval_sinc(v))
    {
        return sin_cardinal_helper<T>::eval_sinc(v);
    }
};
struct func_sinhc
{
    template<class T>
    auto eval(const T& v) const -> decltype(sin_cardinal_helper<T>::eval_sinhc(v))
    {
        return sin_cardinal_helper<T>::eval_sinhc(v);
    }
};

};};

namespace matcl
{

void test_sin_cardinal_inst()
{
    Real x = 0.0;

    sinc(x);
    sinhc(x);
};

Matrix matcl::sinc(const Matrix &A)
{
    return eval_scalar_func(A, details::func_sinc());
};
Matrix matcl::sinhc(const Matrix &A)
{
    return eval_scalar_func(A, details::func_sinhc());
};

};
