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

#include "matcl-core/config.h"
#include "matcl-scalar/details/sequence/seq_transform_base.h"

namespace matcl
{
    
// The routine determines the limit of a given sequence of approximations
// at points x0, x1, ..., xn, xi-> 0, for i -> INF by means of the  
// Richardson algorithm (Neville’s scheme). An estimate of the absolute
// error is also given.
// The algorithm is exact for a sequence of the form:
//      s_n = s + Sum_{j=0}^k c_j x_n^{j+1}
// provided that |x_{i} / x_{i+1}| >= a for some a > 1.
//
// It is not guaranteed, that | res - lim| <= abs_err, where res is estimated
// limit, lim is the true limit, and abs_err is the estimated error; abs_err 
// is a sum of estimated round-off errors due to the algorithm (first order 
// approximation) and a rough estimation of deviation of res from the true 
// limit based on history of computed approximations
//
// References:
//  [1] Weniger, E. J. (1989). Nonlinear sequence transformations for the 
//      acceleration of convergence and the summation of divergent series. 
//      Computer Physics Reports, 10(5-6), 189-371.
template<class Float>
class richardson : protected details::seq_transform_base<Float>
{
    private:        
        vector_points   m_points;

    public:
        // initialize the algorithm; precision is working precision; when
        // precision = 0, then default precision is used; nonzero value has effect
        // only when Float is a multiprecision type
        richardson(int precision = 0);

        // estimate the limit of a given sequence;
        // x        - a new point different than all previous points
        // y        - value of a function at point x
        // res      - estimation of the limit
        // abs_err  - estimation of the absolute error
        void            eval(const Float& x, const Float& y, Float& res, Float& abs_err);

        // estimate the limit of a given sequence
        // x        - a new point different than all previous points
        // y        - value of a function at point x
        // y_err    - estimation of absolute error of y, if y_err <= 0, then 
        //            y_err is set to |y| * epsilon
        // res      - estimation of the limit
        // abs_err  - estimation of the absolute error
        void            eval(const Float& x, const Float& y, const Float& y_err, 
                            Float& res, Float& abs_err);

        // return estimation of limit computed in last call of the eval function
        const Float&    last_result() const;

        // return estimation of round-off errors of computed limit in the last
        // call of the eval function
        const Float&    numerical_error() const;

        // return estimation of |result - lim| ignoring round-off errors, where
        // result is computed in the last call of the eval function
        const Float&    limit_error() const;

        // return estimation of absolute error of computed result in the last call
        // of the eval function; returned value is equal to numerical_error() + 
        // limit_error()
        Float           total_error() const;

        // initialize the algorithm; precision is working precision; when
        // precision = 0, then default precision is used; nonzero value has effect
        // only when Float is a multiprecision type
        void            clear(int precision = 0);

    private:
        void            make();
};

}

#include "matcl-scalar/details/sequence/richardson.inl"
