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
#include "matcl-scalar/details/sequence/seq_transform_base.h"

namespace matcl
{

template<class Float>
class extrapolator
{
    public:
        virtual ~extrapolator(){};

        // estimate the limit of a given sequence;
        // x        - a new point different than all previous points
        // y        - value of a function at point x
        // lim      - estimation of the limit
        // abs_err  - estimation of the absolute error
        virtual void    eval(const Float& x, const Float& val, Float& lim, 
                            Float& err) = 0;

        // estimate the limit of a given sequence
        // x        - a new point different than all previous points
        // y        - value of a function at point x
        // y_err    - estimation of absolute error of y, if y_err <= 0, then 
        //            y_err is set to |y| * epsilon
        // lim      - estimation of the limit
        // abs_err  - estimation of the absolute error
        virtual void    eval(const Float& x, const Float& val, const Float& val_abs_err, 
                            Float& lim, Float& err) = 0;

        // return estimation of limit computed in last call of the eval function
        virtual const Float&
                        last_result() const = 0;

        // return estimation of round-off errors of computed limit in the last
        // call of the eval function
        virtual const Float&    
                        numerical_error() const = 0;

        // return estimation of |result - lim| ignoring round-off errors, where
        // result is computed in the last call of the eval function
        virtual const Float&    
                        limit_error() const = 0;

        // return estimation of absolute error of computed result in the last call
        // of the eval function; returned value is equal to numerical_error() + 
        // limit_error()
        virtual Float   total_error() const = 0;

        // initialize the algorithm; precision is working precision; when
        // precision = 0, then default precision is used; nonzero value has effect
        // only when Float is a multiprecision type
        virtual void    clear(int precision = 0) = 0;
};

// shared pointer storing an extrapolator object
template<class Float>
using extrapolator_ptr      = std::shared_ptr<extrapolator<Float>>;

// make new extrapolator based on wynn_epsilon algorithm;
// see extrapolator::clear for description of the precision argument
template<class Float>
extrapolator_ptr<Float> make_extrapolator_wynn_epsilon(int precision = 0);

// make new extrapolator based on aitken_delta algorithm;
// see extrapolator::clear for description of the precision argument
template<class Float>
extrapolator_ptr<Float> make_extrapolator_aitken_delta(int precision = 0);

// make new extrapolator based on richardson algorithm;
// see extrapolator::clear for description of the precision argument
template<class Float>
extrapolator_ptr<Float> make_extrapolator_richardson(int precision = 0);

};
