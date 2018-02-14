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
    
// The routine determines the limit of a given sequence of approximations,
// by means of the epsilon algorithm of P. Wynn. An estimate of the absolute
// error is also given.
// The epsilon algorithm accelerates convergence of sequence of the form:
//      s_n ~ s + Sum_{j = 0}^inf c_j l_j^n; 1 > l_0 > l_1 > ... > 0, n -> inf
//      s_n ~ s + (-1)^n Sum_{j=0}^inf c_j / (n+b)^{j+1}, b > 0, n -> inf
// The epsilon algorithm is not able to accelerate converge of logarithmically
// convergent sequences, for example
//      s_n ~ s + Sum_{j=0}^inf c_j/(n+b)^{j+1}, b > 0, n ->inf
// The epsilon algorithm is able to find the antilimit of a divergent sequence
//
// It is not guaranteed, that | res - lim| <= abs_err, where res is estimated
// limit, lim is the true limit (or antilimit), and abs_err is the estimated
// error; abs_err is a sum of estimated round-off errors due to the algorithm
// (first order approximation) and a rough estimation of deviation of res from
// the true limit based on history of computed approximations
//
// References:
//  [1] Weniger, E. J. (1989). Nonlinear sequence transformations for the 
//      acceleration of convergence and the summation of divergent series. 
//      Computer Physics Reports, 10(5-6), 189-371.
template<class Float>
class wynn_epsilon : protected details::seq_transform_base<Float>
{
    public:
        // initialize the algorithm; precision is working precision; when
        // precision = 0, then default precision is used; nonzero value has 
        // effect only when Float is a multiprecision type
        wynn_epsilon(int precision = 0);

        // estimate the limit of a given sequence;
        // elem     - new element of the sequence
        // res      - estimation of the limit
        // abs_err  - estimation of the absolute error
        void            eval(const Float& elem, Float& res, Float& abs_err);

        // estimate the limit of a given sequence;
        // elem     - new element of the sequence
        // elem_err - estimation of absolute error of elem, if elem_err <= 0, 
        //            then elem_err is set to |elem| * epsilon
        // res      - estimation of the limit
        // abs_err  - estimation of the absolute error
        void            eval(const Float& elem, const Float& elem_err, Float& res, 
                            Float& abs_err);

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

};

#include "matcl-scalar/details/sequence/wynn_epsilon.inl"