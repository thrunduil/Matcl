/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-simd/poly/poly_eval.h"
#include "matcl-simd/details/poly/utils.h"

namespace matcl { namespace simd { namespace details
{

//-----------------------------------------------------------------------
//                   SHORT HORNER
//-----------------------------------------------------------------------
template<class Trans, class Arg_type, class ... Args>
struct eval_small_horner
{};

template<class Trans, class Arg_type, class Arg1, class ... Args>
struct eval_small_horner<Trans, Arg_type, Arg1, Args...>
{
    force_inline
    static Arg_type eval(const Arg_type& x, const Arg1& coef1, 
                         const Args& ... coef)
    {
        Arg_type res    = eval_small_horner<Trans, Arg_type, Args...>
                            ::eval(x, coef...);
        Arg1 coef_tr    = Trans::eval(coef1);
        res             = eval_fma<Arg_type, Arg1>::eval(res, x, coef_tr);
        return res;
    }
};

template<class Trans, class Arg_type, class Arg1>
struct eval_small_horner<Trans, Arg_type, Arg1>
{
    force_inline
    static Arg_type eval(const Arg_type& x, const Arg1& coef1)
    {
        (void)x;
        Arg1 coef_tr    = Trans::eval(coef1);
        return Arg_type(coef_tr);
    }
};

//-----------------------------------------------------------------------
//                   HORNER
//-----------------------------------------------------------------------
template<class Trans, int Poly_size, class Arg_type, class Coef_type>
struct eval_horner
{
    force_inline
    static Arg_type eval(const Arg_type& x, const Arg_type& res0, 
                         const Coef_type* poly)
    {
        static const int size_1 = Poly_size/2;
        static const int size_2 = Poly_size - size_1;

        Arg_type res = res0;
        res = eval_horner<Trans, size_1, Arg_type, Coef_type>
                ::eval(x, res, poly + size_2);
        res = eval_horner<Trans, size_2, Arg_type, Coef_type>
                ::eval(x, res, poly);
        return res;
    }

    force_inline
    static Arg_type eval2(const Arg_type& x, const Coef_type* poly)
    {
        Arg_type res    = Arg_type(Trans::eval(poly[Poly_size - 1]));

        return eval_horner<Trans, Poly_size - 1, Arg_type, Coef_type>
                ::eval(x, res, poly);
    }
};

template<class Trans, class Arg_type, class Coef_type>
struct eval_horner<Trans, 4, Arg_type, Coef_type>
{
    force_inline
    static Arg_type eval(const Arg_type& x, const Arg_type& res, 
                         const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        Arg_type c1   = Arg_type(Trans::eval(poly[1]));
        Arg_type c2   = Arg_type(Trans::eval(poly[2]));
        Arg_type c3   = Arg_type(Trans::eval(poly[3]));

        return small_horner(x, c0, c1, c2, c3, res);
    }

    force_inline
    static Arg_type eval2(const Arg_type& x, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        Arg_type c1   = Arg_type(Trans::eval(poly[1]));
        Arg_type c2   = Arg_type(Trans::eval(poly[2]));
        Arg_type c3   = Arg_type(Trans::eval(poly[3]));

        return small_horner(x, c0, c1, c2, c3);
    }
};

template<class Trans, class Arg_type, class Coef_type>
struct eval_horner<Trans, 3, Arg_type, Coef_type>
{
    force_inline
    static Arg_type eval(const Arg_type& x, const Arg_type& res, 
                         const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        Arg_type c1   = Arg_type(Trans::eval(poly[1]));
        Arg_type c2   = Arg_type(Trans::eval(poly[2]));

        return small_horner(x, c0, c1, c2, res);
    }

    force_inline
    static Arg_type eval2(const Arg_type& x, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        Arg_type c1   = Arg_type(Trans::eval(poly[1]));
        Arg_type c2   = Arg_type(Trans::eval(poly[2]));

        return small_horner(x, c0, c1, c2);
    }
};

template<class Trans, class Arg_type, class Coef_type>
struct eval_horner<Trans, 2, Arg_type, Coef_type>
{
    force_inline
    static Arg_type eval(const Arg_type& x, const Arg_type& res, 
                         const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        Arg_type c1   = Arg_type(Trans::eval(poly[1]));

        return small_horner(x, c0, c1, res);
    }

    force_inline
    static Arg_type eval2(const Arg_type& x, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        Arg_type c1   = Arg_type(Trans::eval(poly[1]));

        return small_horner(x, c0, c1);
    }
};

template<class Trans, class Arg_type, class Coef_type>
struct eval_horner<Trans, 1, Arg_type, Coef_type>
{
    force_inline
    static Arg_type eval(const Arg_type& x, const Arg_type& res, 
                         const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        return small_horner(x, c0, res);
    }

    force_inline
    static Arg_type eval2(const Arg_type& x, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));

        return small_horner(x, c0);
    }
};

template<class Trans, class Arg_type, class Coef_type>
struct eval_horner<Trans, 0, Arg_type, Coef_type>
{
    force_inline
    static Arg_type eval(const Arg_type& x, const Arg_type& res, 
                         const Coef_type* poly)
    {
        (void)x;
        (void)poly;
        return res;
    }
};

//-----------------------------------------------------------------------
//                   HORNER + ERROR
//-----------------------------------------------------------------------
template<int Poly_size, class Arg_type, class Coef_type>
struct eval_horner_error
{
    force_inline
    static Arg_type eval(const Arg_type& x, const Arg_type& abs_x, 
                         const Arg_type& res0, const Coef_type* poly, Arg_type& mu)
    {
        Arg_type res= res0;

        Arg_type cl = Arg_type(poly[Poly_size - 1]);
        res         = small_horner(x, cl, res);

        Arg_type ar = details::eval_abs<Arg_type>::eval(res);
        mu          = small_horner(abs_x, ar, mu);

        return eval_horner_error<Poly_size - 1, Arg_type, Coef_type>
                    ::eval(x, abs_x, res, poly, mu);
    }
};

template<class Arg_type, class Coef_type>
struct eval_horner_error<0, Arg_type, Coef_type>
{
    force_inline
    static Arg_type eval(const Arg_type& x, const Arg_type& abs_x, const Arg_type& res, 
                         const Coef_type* poly, Arg_type& mu)
    {
        (void)x;
        (void)abs_x;
        (void)poly;
        (void)mu;
        return res;
    }
};

//-----------------------------------------------------------------------
//                      HORNER DYNAMIC
//-----------------------------------------------------------------------
template<class Trans, class Arg_type, class Coef_type>
struct eval_horner2
{
    static const int unroll_size    = 8;

    force_inline
    static Arg_type eval(const Arg_type& x, const Arg_type& res0, int size, 
                         const Coef_type* poly)
    {
        int n1      = size / unroll_size;
        int n2      = size % unroll_size;

        const Coef_type* poly_tail   = poly + size;
        Arg_type res    = res0;

        for (int j = 0; j < n1; ++j)
        {
            poly_tail   = poly_tail - unroll_size;
            res         = eval_horner<Trans, unroll_size, Arg_type, Coef_type>
                                ::eval(x, res, poly_tail);            
        };

        for (int i = n2 - 1; i >= 0; --i)
        {
            Coef_type c = Trans::eval(poly[i]);
            res = eval_fma<Arg_type, Coef_type>::eval(res, x, c);
        }

        return res;
    }
};

//-----------------------------------------------------------------------
//                      HORNER APOSTERIORI COND
//-----------------------------------------------------------------------
template<class Arg_type, class Coef_type>
struct eval_horner_post_cond
{
    static const int unroll_size    = 8;

    using Trans         = trans_id;

    force_inline
    static Arg_type eval(const Arg_type& x, int size, const Coef_type* poly,
                         Arg_type& mu)
    {
        Arg_type res;
        eval_base(x, size, poly, res, mu);

        Arg_type u      = Arg_type(0.5) * details::eval_eps<Arg_type>::eval();
        mu              = mu * u;

        return res;
    };

    force_inline
    static Arg_type eval_cond(const Arg_type& x, int size, const Coef_type* poly,
                              Arg_type& p_app, Arg_type& mu)
    {
        Arg_type err;

        eval_base(x, size, poly, p_app, err);

        Arg_type u      = Arg_type(0.5) * details::eval_eps<Arg_type>::eval();
        mu              = err * u;
        Arg_type p_abs  = eval_abs<Arg_type>::eval(p_app);

        // take lower bound for |p(x)| based on |p_ap(x)| and mu
        Arg_type cond   = err/(p_abs - mu);

        return cond;
    };

    force_inline
    static void eval_base(const Arg_type& x, int size, const Coef_type* poly,
                         Arg_type& res, Arg_type& mu)
    {
        Arg_type abs_x  = details::eval_abs<Arg_type>::eval(x);
        res             = Arg_type(poly[size-1]);
        mu              = details::eval_abs<Arg_type>::eval(res)/Arg_type(2);
        size            -= 1;

        int n1          = size / unroll_size;
        int n2          = size % unroll_size;

        const Coef_type* poly_tail   = poly + size;

        for (int j = 0; j < n1; ++j)
        {
            poly_tail   = poly_tail - unroll_size;
            res         = eval_horner_error<unroll_size, Arg_type, Coef_type>
                                ::eval(x, abs_x, res, poly_tail, mu);            
        };

        for (int i = n2 - 1; i >= 0; --i)
        {
            res         = eval_horner_error<1, Arg_type, Coef_type>
                                ::eval(x, abs_x, res, poly + i, mu);            
        }

        Arg_type abs_r  = details::eval_abs<Arg_type>::eval(res);
        mu              = fms_f(Arg_type(2.0), mu, abs_r);
    };
};

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
        Arg_type res    = Arg_type(poly[size - 1]);

        // algorithm from "Faithful Polynomial Evaluation with Compensated
        // Horner Algorithm", P. Langlois, N. Louvet
        for (int i = size - 2; i >= 0; --i)
        {
            Arg_type c          = Arg_type(poly[i]);
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
            twofold_type c      = twofold_type(Arg_type(poly[i].value), 
                                    Arg_type(poly[i].error));

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
            twofold_type c      = twofold_type(Arg_type(poly[i].value), 
                                    Arg_type(poly[i].error));

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
    Arg_type res = Arg_type(poly[poly_size - 1]);
    res          = details::eval_abs<Arg_type>::eval(res);

    return details::eval_horner2<details::trans_abs, Arg_type, Coef_type>
                ::eval(x, res, poly_size - 1, poly);
};

}}}

namespace matcl
{

template<class Arg_type, class ... Args>
force_inline
Arg_type simd::small_horner(const Arg_type& x, const Args& ... coef)
{
    return details::eval_small_horner<details::trans_id, Arg_type, Args ...>
                        ::eval(x, coef...);
};

template<int Poly_size, class Arg_type, class Coef_type>
force_inline
Arg_type simd::horner(const Arg_type& x, const Coef_type* poly)
{
    Arg_type res    = Arg_type(poly[Poly_size - 1]);

    return details::eval_horner<details::trans_id, Poly_size - 1, Arg_type, Coef_type>
                        ::eval(x, res, poly);
}

template<class Arg_type, class Coef_type>
Arg_type simd::horner(const Arg_type& x, int poly_size, const Coef_type* poly)
{
    Arg_type res = Arg_type(poly[poly_size - 1]);

    return details::eval_horner2<details::trans_id, Arg_type, Coef_type>
                ::eval(x, res, poly_size - 1, poly);
};

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
    Arg_type res = Arg_type(poly[N - 1]);
    res          = details::eval_horner2_compensated_error<Arg_type, Coef_type>
                    ::eval(x, res, N - 1, poly, error, is_faithfully_rounded);
    return res;
};

template<class Arg_type, class Coef_type>
Arg_type simd::horner_apriori_cond(const Arg_type& x, int N, const Coef_type* poly)
{
    Arg_type val;
    Arg_type abs_val;

    return horner_apriori_cond(x, N, poly, val, abs_val);
};

template<class Arg_type, class Coef_type>
Arg_type simd::horner_apriori_cond(const Arg_type& x, int N, const Coef_type* poly,
                Arg_type& val, Arg_type& abs_val)
{
    Arg_type abs_x  = details::eval_abs<Arg_type>::eval(x);
    abs_val         = details::horner_abs(abs_x, N, poly);
    val             = horner(x, N, poly);    
    Arg_type aval   = details::eval_abs<Arg_type>::eval(val);

    return abs_val / aval;
};

template<class Arg_type, class Coef_type>
Arg_type simd::horner_and_error(const Arg_type& x, int poly_size, const Coef_type* poly,
                                       Arg_type& error)
{
    Arg_type res = details::eval_horner_post_cond<Arg_type, Coef_type>
                    ::eval(x, poly_size, poly, error);
    return res;
};

template<class Arg_type, class Coef_type>
Arg_type simd::horner_aposteriori_cond(const Arg_type& x, int N, const Coef_type* poly)
{
    Arg_type val, err;
    return horner_aposteriori_cond(x, N, poly, val, err);
};

template<class Arg_type, class Coef_type>
Arg_type simd::horner_aposteriori_cond(const Arg_type& x, int N, const Coef_type* poly,
                Arg_type& val, Arg_type& err)
{
    Arg_type res = details::eval_horner_post_cond<Arg_type, Coef_type>
                    ::eval_cond(x, N, poly, val, err);
    return res;
};

}
