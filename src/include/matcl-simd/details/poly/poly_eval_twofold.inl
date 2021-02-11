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

#include "matcl-simd/poly/poly_eval_twofold.h"
#include "matcl-simd/details/poly/utils.h"

namespace matcl { namespace simd { namespace details
{

//-----------------------------------------------------------------------
//                      HORNER COMPENSATED
//-----------------------------------------------------------------------
template<class Arg_type, class Coef_type>
struct eval_horner2_compensated
{
    using twofold_type          = twofold<Arg_type>;

    force_inline
    static twofold_type eval(const Arg_type& x, int size, const Coef_type* poly)
    {
        Arg_type comp   = Arg_type(0.0);
        Arg_type res    = broadcast<Arg_type>::eval(poly[size - 1]);

        // algorithm from "Faithful Polynomial Evaluation with Compensated
        // Horner Algorithm", P. Langlois, N. Louvet
        for (int i = size - 2; i >= 0; --i)
        {
            Arg_type c          = broadcast<Arg_type>::eval(poly[i]);
            twofold_type prod   = twofold_mult(res, x);            
            twofold_type sum    = twofold_sum(prod.value, c);

            Arg_type err        = prod.error + sum.error;
            comp                = eval_fma<Arg_type, Arg_type>::eval(comp, x, err);
            res                 = sum.value;
        }

        // instead of returning res + comp we return the twofold number
        return twofold_sum(res, comp);
    };

    force_inline
    static twofold_type eval(const Arg_type& x, int size, const twofold<Coef_type>* poly)
    {
        twofold_type res    = poly[size - 1];

        // algorithm from "Faithful Polynomial
        // Evaluation with Compensated Horner Algorithm", P. Langlois, N. Louvet

        for (int i = size - 2; i >= 0; --i)
        {
            twofold_type c      = twofold_type(broadcast<Arg_type>::eval(poly[i].value), 
                                               broadcast<Arg_type>::eval(poly[i].error));

            // 3 u^2 relative error (2 if FMA is available)
            twofold_type prod   = res * x;

            // 3 u^2 relative error in quad precision
            res                 = prod + c;
        }

        // repeating error counting for double precision we obtain
        // |p(x) - p_app(x)|/|p(x)| <= c * cond(p,x) * u^2
        // m - accuracy of multiplication (in u^2 units)
        // a - accuracy of addition (in u^2 units)
        // c = gam((N-1)(m+a))
        // gam(k) = k*u^2/(1-k*u^2)
        // thus m + a = 6, c = 6*(N-1)/(1-(N-1)*6*u^2)
        return res;
    };

    force_inline
    static twofold_type eval(const twofold_type& x, int size, const twofold<Coef_type>* poly)
    {
        twofold_type res    = poly[size - 1];

        // algorithm from "Faithful Polynomial
        // Evaluation with Compensated Horner Algorithm", P. Langlois, N. Louvet

        for (int i = size - 2; i >= 0; --i)
        {
            twofold_type c      = twofold_type(broadcast<Arg_type>::eval(poly[i].value), 
                                               broadcast<Arg_type>::eval(poly[i].error));

            // 7 u^2 relative error
            twofold_type prod   = res * x;

            // 3 u^2 relative error in quad precision
            res                 = prod + c;
        }

        // repeating error counting for double precision we obtain
        // |p(x) - p_app(x)|/|p(x)| <= c * cond(p,x) * u^2
        // m - accuracy of multiplication (in u^2 units)
        // a - accuracy of addition (in u^2 units)
        // c = gam((N-1)(m+a))
        // gam(k) = k*u^2/(1-k*u^2)
        // thus m + a = 10, c = 10*(N-1)/(1-(N-1)*10*u^2)
        return res;
    };
};

//-----------------------------------------------------------------------
//                      HORNER COMPENSATED AND ERROR
//-----------------------------------------------------------------------
template<class Arg_type, class Coef_type>
struct eval_horner2_compensated_error
{
    force_inline
    static Arg_type eval(const Arg_type& x, const Arg_type& res0, int size, 
                         const Coef_type* poly, Arg_type& res_err, 
                         bool& is_faithfully_rounded)
    {
        using twofold_type = twofold<Arg_type>;

        // modified version of the algorithm from "Faithful Polynomial
        // Evaluation with Compensated Horner Algorithm", P. Langlois, N. Louvet

        Arg_type abs_x      = eval_abs<Arg_type>::eval(x);
        Arg_type comp_b     = Arg_type(0.0);
        Arg_type comp_c     = Arg_type(0.0);
        Arg_type res        = res0;

        for (int i = size - 1; i >= 0; --i)
        {
            // evaluate one term using compensated horner:
            // [res, err] = res * x + a
            Coef_type c         = poly[i];
            twofold_type prod   = twofold_mult(res, x);            
            twofold_type sum    = twofold_sum(prod.value, Arg_type(c));
            res                 = sum.value;
            
            // new coefficient of residual polynomial horner
            Arg_type err        = prod.error + sum.error;

            // new coefficient of |residual polynomial horner|
            Arg_type abs_err    = eval_abs<Arg_type>::eval(prod.error)
                                + eval_abs<Arg_type>::eval(sum.error);

            // evaluate one term of the residual polynomial using horner:
            comp_c              = eval_fma<Arg_type, Arg_type>
                                    ::eval(comp_c, x, err);

            // evaluate one term of the |residual polynomial| using horner:
            comp_b              = eval_fma<Arg_type, Arg_type>
                                    ::eval(comp_b, abs_x, abs_err);
        }

        twofold_type sum        = twofold_sum(res, comp_c);

        res         = sum.value;
        
        // N is the order of polynomial (the highest power), i.e. N = size - 1
        // but we have additional term (res)
        int N           = size;

        const Arg_type eps    = eval_eps<Arg_type>::eval();
        const Arg_type half   = Arg_type(0.5);
        const Arg_type one    = Arg_type(1);
        const Arg_type u      = eps * half;
        const Arg_type u_half = u * half;        

        Arg_type g      = eval_gamma(2 * N - 1, u);        
        Arg_type N1     = Arg_type(float(N + 1));

        Arg_type a      = (g * comp_b) / (one - N1 * eps);

        res_err         = a + eval_abs<Arg_type>::eval(sum.error);
        res_err         = res_err / (one - eps);

        Arg_type max_a  = u_half * eval_abs<Arg_type>::eval(res);
        
        if (eval_lt<Arg_type>::eval(a, max_a) == true)
            is_faithfully_rounded = true;
        else
            is_faithfully_rounded = false;

        return res;
    };

    force_inline
    static Arg_type eval_gamma(int k, const Arg_type& u)
    {
        Arg_type ku = Arg_type(float(k)) * u;
        Arg_type r  = ku / (Arg_type(1) - ku);
        return r;
    };
};

template<class Arg_type, class Coef_type>
Arg_type horner_abs(const Arg_type& x, int poly_size, const Coef_type* poly)
{
    Arg_type res = broadcast<Arg_type>::eval(poly[poly_size - 1]);
    res          = details::eval_abs<Arg_type>::eval(res);

    return details::eval_horner2<details::trans_abs, Arg_type, Coef_type>
                ::eval(x, res, poly_size - 1, poly);
};

}}}

namespace matcl
{

template<class Arg_type, class Coef_type>
twofold<Arg_type>
simd::compensated_horner(const Arg_type& x, int poly_size, const Coef_type* poly)
{
    return details::eval_horner2_compensated<Arg_type, Coef_type>
                ::eval(x, poly_size, poly);
};

template<class Arg_type, class Coef_type>
twofold<Arg_type>
simd::compensated_horner(const Arg_type& x, int poly_size, const twofold<Coef_type>* poly)
{
    return details::eval_horner2_compensated<Arg_type, Coef_type>
                ::eval(x, poly_size, poly);
};

template<class Arg_type, class Coef_type>
twofold<Arg_type>
simd::compensated_horner(const twofold<Arg_type>& x, int poly_size, 
                         const twofold<Coef_type>* poly)
{
    return details::eval_horner2_compensated<Arg_type, Coef_type>
                ::eval(x, poly_size, poly);
};

template<class Arg_type, class Coef_type>
Arg_type simd::compensated_horner_and_error(const Arg_type& x, int N, const Coef_type* poly,
                Arg_type& error, bool& is_faithfully_rounded)
{
    Arg_type res = details::broadcast<Arg_type>::eval(poly[N - 1]);
    res          = details::eval_horner2_compensated_error<Arg_type, Coef_type>
                    ::eval(x, res, N - 1, poly, error, is_faithfully_rounded);
    return res;
};

}
