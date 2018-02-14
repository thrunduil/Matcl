/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017-2018
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

#include "matcl-core/config.h"
#include "matcl-scalar/lib_functions/sequence/objects.h"

#include <functional>

#pragma warning(push)
#pragma warning(disable: 4251) //needs to have dll-interface to be used by clients

namespace matcl { namespace details
{

template<class Float>
class limit_estimator_impl;

}}

namespace matcl
{

// limit type
enum class limit_type
{
    left,   // compute left-handed limit
    right,  // compute rigth-handed limit
    both    // assume, that left-handed and right-handed limit exists
};

// estimate limit of a function f(x) for x -> x0; x0 can a finite number
// or +- Inf.  An estimate of the absolute error is also given. It is 
// assumed that the limit exists.The algorithm builds a sequence of points
// dx_i converging to zero:
//      dx_{i+1} = scale_i * dx_{i}, given dx_0 
// a sequence of function values y_i
//      y_i = f(x0 + mult_i * dx_i)
// where mult_i = +1 if limit_type = right, mult_i = -1 if limit_type = left,
// and mult_i = (-1)^i otherwise. The limit is approximated using a sequence
// accelerator. When x0 = +- Inf, then the limit f(1/x) for x -> 0 is 
// computed
template<class Float>
class MATCL_SCALAR_EXPORT limit_estimator
{
    private:
        using impl_type         = std::shared_ptr<details::limit_estimator_impl<Float>>;

    public:
        // type of function evaluating f(x)
        using function_type     = std::function<Float (const Float&)>;

        // type of an extrapolator
        using extrapolator_type = extrapolator_ptr<Float>;

    private:
        impl_type   m_impl;

    public:
        // initialize the algorithm; precision is working precision; when
        // precision = 0, then default precision is used; nonzero value has 
        // effect only when Float is a multiprecision type
        limit_estimator(int precision = 0);

        // destructor
        ~limit_estimator();

        // approximate the limit for a function f at point x0 of type lim_type
        // estimated absolute error of computed limit is returned in abs_err.
        // It is not guaranteed, that | res - lim| <= abs_err, where res is 
        // estimated limit, and lim is the true limit; estimated error abs_err
        // is the value returned by the extrapolator used 
        Float       eval(const function_type& f, const Float& x0, 
                        limit_type lim_type, Float& abs_err);

        // initialize the algorithm; precision is working precision; when
        // precision = 0, then default precision is used; nonzero value has effect
        // only when Float is a multiprecision type
        void        clear(int precision = 0);

        // return number of iterations performed
        int         number_iterations() const;

        // set extrapolator; on default the wynn_epsilon algorithm
        // is used
        void        set_extrapolator(const extrapolator_type& extr);

        // return current extrapolator used by the algorithm
        extrapolator_type
                    get_extrapolator() const;

        // set maximum number of iterations allowed; when number of 
        // iterations reach this value, then algorithm is stopped and
        // current best limit estimation is returned; if max_iter <= 0
        // then default value (100) is used
        void        set_max_iterations(int max_iter);

        // get current maximum iterations limit
        int         get_max_iterations() const;

        // set initial scaling factor scale_0; if scaling <= 0, then 
        // default value (0.5) is used
        void        set_initial_scaling(const Float& scaling);

        // get current scaling factor scale_0
        Float       get_initial_scaling() const;

        // initial step is determined as dx_0 = |x0| * scaling if x0 != 0,
        // and dx_0 = scaling if x0 = 0. If scaling <= 0, then default 
        // value (1.0) is used
        void        set_initial_point_scaling(const Float& scaling);

        // return current point scaling
        Float       get_initial_point_scaling() const;
};

}

#pragma warning(pop)