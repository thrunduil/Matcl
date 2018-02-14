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

#include "matcl-scalar/lib_functions/sequence/objects.h"
#include "matcl-scalar/lib_functions/sequence/wynn_epsilon.h"
#include "matcl-scalar/lib_functions/sequence/aitken_delta.h"
#include "matcl-scalar/lib_functions/sequence/richardson.h"

namespace matcl
{

// impementation of extrapolator based on wynn_epsilon algorithm;
// see also wynn_epsilon class
template<class Float>
class extrapolator_wynn_epsilon : public extrapolator<Float>
{
    private:
        wynn_epsilon<Float> m_impl;

    public:
        // initialize the algorithm; precision is working precision; when
        // precision = 0, then default precision is used; nonzero value has effect
        // only when Float is a multiprecision type
        extrapolator_wynn_epsilon(int precision = 0);

        // destructor
        ~extrapolator_wynn_epsilon() override{};

        // estimate the limit of a given sequence;
        // x        - a new point different than all previous points
        // y        - value of a function at point x
        // lim      - estimation of the limit
        // abs_err  - estimation of the absolute error
        void        eval(const Float& x, const Float& val, Float& lim, 
                            Float& err) override;

        // estimate the limit of a given sequence
        // x        - a new point different than all previous points
        // y        - value of a function at point x
        // y_err    - estimation of absolute error of y, if y_err <= 0, then 
        //            y_err is set to |y| * epsilon
        // lim      - estimation of the limit
        // abs_err  - estimation of the absolute error
        void        eval(const Float& x, const Float& val, const Float& val_abs_err, 
                            Float& lim, Float& err) override;

        // return estimation of limit computed in last call of the eval function
        const Float& last_result() const override;

        // return estimation of round-off errors of computed limit in the last
        // call of the eval function
        const Float& numerical_error() const override;

        // return estimation of |result - lim| ignoring round-off errors, where
        // result is computed in the last call of the eval function
        const Float& limit_error() const override;

        // return estimation of absolute error of computed result in the last call
        // of the eval function; returned value is equal to numerical_error() + 
        // limit_error()
        Float       total_error() const override;

        // initialize the algorithm; precision is working precision; when
        // precision = 0, then default precision is used; nonzero value has effect
        // only when Float is a multiprecision type
        void        clear(int precision = 0) override;
};

// impementation of extrapolator based on aitken_delta algorithm;
// see also aitken_delta class
template<class Float>
class extrapolator_aitken_delta : public extrapolator<Float>
{
    private:
        aitken_delta<Float> m_impl;

    public:
        // initialize the algorithm; precision is working precision; when
        // precision = 0, then default precision is used; nonzero value has effect
        // only when Float is a multiprecision type
        extrapolator_aitken_delta(int precision = 0);

        // destructor
        ~extrapolator_aitken_delta() override{};

        // estimate the limit of a given sequence;
        // x        - a new point different than all previous points
        // y        - value of a function at point x
        // lim      - estimation of the limit
        // abs_err  - estimation of the absolute error
        void        eval(const Float& x, const Float& val, Float& lim, 
                            Float& err) override;

        // estimate the limit of a given sequence
        // x        - a new point different than all previous points
        // y        - value of a function at point x
        // y_err    - estimation of absolute error of y, if y_err <= 0, then 
        //            y_err is set to |y| * epsilon
        // lim      - estimation of the limit
        // abs_err  - estimation of the absolute error
        void        eval(const Float& x, const Float& val, const Float& val_abs_err, 
                            Float& lim, Float& err) override;

        // return estimation of limit computed in last call of the eval function
        const Float& last_result() const override;

        // return estimation of round-off errors of computed limit in the last
        // call of the eval function
        const Float& numerical_error() const override;

        // return estimation of |result - lim| ignoring round-off errors, where
        // result is computed in the last call of the eval function
        const Float& limit_error() const override;

        // return estimation of absolute error of computed result in the last call
        // of the eval function; returned value is equal to numerical_error() + 
        // limit_error()
        Float       total_error() const override;

        // initialize the algorithm; precision is working precision; when
        // precision = 0, then default precision is used; nonzero value has effect
        // only when Float is a multiprecision type
        void        clear(int precision = 0) override;
};

// impementation of extrapolator based on richardson algorithm;
// see also richardson class
template<class Float>
class extrapolator_richardson : public extrapolator<Float>
{
    private:
        richardson<Float>   m_impl;

    public:
        // initialize the algorithm; precision is working precision; when
        // precision = 0, then default precision is used; nonzero value has effect
        // only when Float is a multiprecision type
        extrapolator_richardson(int precision = 0);

        // destructor
        ~extrapolator_richardson() override{};

        // estimate the limit of a given sequence;
        // x        - a new point different than all previous points
        // y        - value of a function at point x
        // lim      - estimation of the limit
        // abs_err  - estimation of the absolute error
        void        eval(const Float& x, const Float& val, Float& lim, 
                            Float& err) override;

        // estimate the limit of a given sequence
        // x        - a new point different than all previous points
        // y        - value of a function at point x
        // y_err    - estimation of absolute error of y, if y_err <= 0, then 
        //            y_err is set to |y| * epsilon
        // lim      - estimation of the limit
        // abs_err  - estimation of the absolute error
        void        eval(const Float& x, const Float& val, const Float& val_abs_err, 
                            Float& lim, Float& err) override;

        // return estimation of limit computed in last call of the eval function
        const Float& last_result() const override;

        // return estimation of round-off errors of computed limit in the last
        // call of the eval function
        const Float& numerical_error() const override;

        // return estimation of |result - lim| ignoring round-off errors, where
        // result is computed in the last call of the eval function
        const Float& limit_error() const override;

        // return estimation of absolute error of computed result in the last call
        // of the eval function; returned value is equal to numerical_error() + 
        // limit_error()
        Float       total_error() const override;

        // initialize the algorithm; precision is working precision; when
        // precision = 0, then default precision is used; nonzero value has effect
        // only when Float is a multiprecision type
        void        clear(int precision = 0) override;
};

};

#include "matcl-scalar/details/sequence/objects_impl.inl"
