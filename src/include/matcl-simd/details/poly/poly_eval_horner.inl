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
    static Arg_type eval(Arg_type x, Arg1 coef1, Args ... coef)
    {
        Arg_type res  = eval_small_horner<Trans, Arg_type, Args...>
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
    static Arg_type eval(Arg_type x, Arg1 coef1)
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
    static Arg_type eval(Arg_type x, Arg_type res, const Coef_type* poly)
    {
        static const int size_1 = Poly_size/2;
        static const int size_2 = Poly_size - size_1;

        res = eval_horner<Trans, size_1, Arg_type, Coef_type>
                ::eval(x, res, poly + size_2);
        res = eval_horner<Trans, size_2, Arg_type, Coef_type>
                ::eval(x, res, poly);
        return res;
    }

    force_inline
    static Arg_type eval2(Arg_type x, const Coef_type* poly)
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
    static Arg_type eval(Arg_type x, Arg_type res, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        Arg_type c1   = Arg_type(Trans::eval(poly[1]));
        Arg_type c2   = Arg_type(Trans::eval(poly[2]));
        Arg_type c3   = Arg_type(Trans::eval(poly[3]));

        return small_horner(x, c0, c1, c2, c3, res);
    }

    force_inline
    static Arg_type eval2(Arg_type x, const Coef_type* poly)
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
    static Arg_type eval(Arg_type x, Arg_type res, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        Arg_type c1   = Arg_type(Trans::eval(poly[1]));
        Arg_type c2   = Arg_type(Trans::eval(poly[2]));

        return small_horner(x, c0, c1, c2, res);
    }

    force_inline
    static Arg_type eval2(Arg_type x, const Coef_type* poly)
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
    static Arg_type eval(Arg_type x, Arg_type res, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        Arg_type c1   = Arg_type(Trans::eval(poly[1]));

        return small_horner(x, c0, c1, res);
    }

    force_inline
    static Arg_type eval2(Arg_type x, const Coef_type* poly)
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
    static Arg_type eval(Arg_type x, Arg_type res, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        return small_horner(x, c0, res);
    }

    force_inline
    static Arg_type eval2(Arg_type x, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));

        return small_horner(x, c0);
    }
};

template<class Trans, class Arg_type, class Coef_type>
struct eval_horner<Trans, 0, Arg_type, Coef_type>
{
    force_inline
    static Arg_type eval(Arg_type x, Arg_type res, const Coef_type* poly)
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
    static Arg_type eval(Arg_type x, Arg_type abs_x, Arg_type res, 
                         const Coef_type* poly, Arg_type& mu)
    {
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
    static Arg_type eval(Arg_type x, Arg_type abs_x, Arg_type res, 
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
    static Arg_type eval(Arg_type x, Arg_type res, int size, const Coef_type* poly)
    {
        int n1      = size / unroll_size;
        int n2      = size % unroll_size;

        const Coef_type* poly_tail   = poly + size;

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
    static Arg_type eval(Arg_type x, int size, const Coef_type* poly,
                         Arg_type& mu)
    {
        Arg_type abs_x  = details::eval_abs<Arg_type>::eval(x);
        Arg_type res    = Arg_type(poly[size-1]);
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
        Arg_type u      = details::eval_eps<Arg_type>::eval();
        mu              = (mu - Arg_type(0.5) * abs_r) * u;

        return res;
    };
};

//-----------------------------------------------------------------------
//                      HORNER COMPENSATED
//-----------------------------------------------------------------------
template<class Trans, class Arg_type, class Coef_type>
struct eval_horner2_compensated
{
    force_inline
    static Arg_type eval(Arg_type x, Arg_type res, int size, const Coef_type* poly)
    {
        using twofold_type = twofold<Arg_type>;

        Arg_type comp   = Arg_type(0.0);

        for (int i = size - 1; i >= 0; --i)
        {
            twofold_type prod   = twofold_mult(res, x);
            Coef_type c         = Trans::eval(poly[i]);
            twofold_type sum    = twofold_sum(prod.value, Arg_type(c));

            Arg_type err        = prod.error + sum.error;
            comp                = eval_fma<Arg_type, Arg_type>::eval(comp, x, err);
            res                 = sum.value;
        }

        return res + comp;
    };
};

//-----------------------------------------------------------------------
//                      HORNER COMPENSATED AND ERROR
//-----------------------------------------------------------------------
template<class Arg_type, class Coef_type>
struct eval_horner2_compensated_error
{
    force_inline
    static Arg_type eval(Arg_type x, Arg_type res, int size, const Coef_type* poly,
                         Arg_type& res_err, bool& is_exactly_rounded)
    {
        using twofold_type = twofold<Arg_type>;

        // modified version of the algorithm from "Faithful Polynomial
        // Evaluation with Compensated Horner Algorithm", P. Langlois, N. Louvet

        Arg_type abs_x      = eval_abs<Arg_type>::eval(x);
        Arg_type comp_b     = Arg_type(0.0);
        Arg_type comp_c     = Arg_type(0.0);

        for (int i = size - 1; i >= 0; --i)
        {
            // evaluate one term using compensated horner:
            // [res, err] = res * x + a
            twofold_type prod   = twofold_mult(res, x);
            Coef_type c         = poly[i];
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
        int N       = size;

        Arg_type eps    = eval_eps<Arg_type>::eval();
        Arg_type u      = eps / Arg_type(2);
        Arg_type g      = eval_gamma(2 * N - 1, u);
        Arg_type one    = Arg_type(1);
        Arg_type N1     = Arg_type(float(N + 1));

        Arg_type a      = (g * comp_b) / (one - N1 * eps);

        res_err         = a + eval_abs<Arg_type>::eval(sum.error);
        res_err         = res_err / (one - eps);

        Arg_type max_a  = (u / Arg_type(2)) * eval_abs<Arg_type>::eval(res);
        
        if (eval_lt<Arg_type>::eval(a, max_a) == true)
            is_exactly_rounded = true;
        else
            is_exactly_rounded = false;

        return res;
    };

    force_inline
    static Arg_type eval_gamma(int k, Arg_type u)
    {
        Arg_type ku = Arg_type(float(k)) * u;
        Arg_type r  = ku / (Arg_type(1) - ku);
        return r;
    };
};

template<class Arg_type, class Coef_type>
Arg_type horner_abs(Arg_type x, int poly_size, const Coef_type* poly)
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
Arg_type simd::small_horner(Arg_type x, Args ... coef)
{
    return details::eval_small_horner<details::trans_id, Arg_type, Args ...>
                        ::eval(x, coef...);
};

template<int Poly_size, class Arg_type, class Coef_type>
force_inline
Arg_type simd::horner(Arg_type x, const Coef_type* poly)
{
    Arg_type res    = Arg_type(poly[Poly_size - 1]);

    return details::eval_horner<details::trans_id, Poly_size - 1, Arg_type, Coef_type>
                        ::eval(x, res, poly);
}

template<class Arg_type, class Coef_type>
Arg_type simd::horner(Arg_type x, int poly_size, const Coef_type* poly)
{
    Arg_type res = Arg_type(poly[poly_size - 1]);

    return details::eval_horner2<details::trans_id, Arg_type, Coef_type>
                ::eval(x, res, poly_size - 1, poly);
};

template<class Arg_type, class Coef_type>
Arg_type simd::compensated_horner(Arg_type x, int poly_size, const Coef_type* poly)
{
    Arg_type res = Arg_type(poly[poly_size - 1]);

    return details::eval_horner2_compensated<details::trans_id, Arg_type, Coef_type>
                ::eval(x, res, poly_size - 1, poly);
};

template<class Arg_type, class Coef_type>
Arg_type simd::compensated_horner_and_error(Arg_type x, int N, const Coef_type* poly,
                Arg_type& error, bool& is_exactly_rounded)
{
    Arg_type res = Arg_type(poly[N - 1]);
    res          = details::eval_horner2_compensated_error<Arg_type, Coef_type>
                    ::eval(x, res, N - 1, poly, error, is_exactly_rounded);
    return res;
};

template<class Arg_type, class Coef_type>
Arg_type simd::horner_apriori_cond(Arg_type x, int poly_size, const Coef_type* poly)
{
    Arg_type val;
    Arg_type val_abs;

    return horner_apriori_cond(x, poly_size, poly, val, val_abs);
};

template<class Arg_type, class Coef_type>
Arg_type simd::horner_apriori_cond(Arg_type x, int N, const Coef_type* poly,
                                Arg_type& val, Arg_type& val_abs)
{
    Arg_type abs_x  = details::eval_abs<Arg_type>::eval(x);
    val_abs         = details::horner_abs(abs_x, N, poly);
    val             = horner(x, N, poly);    
    Arg_type aval   = details::eval_abs<Arg_type>::eval(val);

    return val_abs / aval;
};

template<class Arg_type, class Coef_type>
Arg_type simd::horner_and_error(Arg_type x, int poly_size, const Coef_type* poly,
                                       Arg_type& error)
{
    Arg_type res = details::eval_horner_post_cond<Arg_type, Coef_type>
                    ::eval(x, poly_size, poly, error);
    return res;
};

}
