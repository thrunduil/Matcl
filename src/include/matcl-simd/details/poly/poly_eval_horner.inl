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
        Arg1 coef_tr    = Arg1(Trans::eval(coef1));
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
        Arg1 coef_tr    = Arg1(Trans::eval(coef1));
        return coef_tr;
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
        Arg_type res    = broadcast<Arg1>::eval(Trans::eval(poly[Poly_size - 1]));

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
        Arg_type c0   = broadcast<Arg_type>::eval(Trans::eval(poly[0]));
        Arg_type c1   = broadcast<Arg_type>::eval(Trans::eval(poly[1]));
        Arg_type c2   = broadcast<Arg_type>::eval(Trans::eval(poly[2]));
        Arg_type c3   = broadcast<Arg_type>::eval(Trans::eval(poly[3]));

        return small_horner(x, c0, c1, c2, c3, res);
    }

    force_inline
    static Arg_type eval2(const Arg_type& x, const Coef_type* poly)
    {
        Arg_type c0   = broadcast<Arg_type>::eval(Trans::eval(poly[0]));
        Arg_type c1   = broadcast<Arg_type>::eval(Trans::eval(poly[1]));
        Arg_type c2   = broadcast<Arg_type>::eval(Trans::eval(poly[2]));
        Arg_type c3   = broadcast<Arg_type>::eval(Trans::eval(poly[3]));

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
        Arg_type c0   = broadcast<Arg_type>::eval(Trans::eval(poly[0]));
        Arg_type c1   = broadcast<Arg_type>::eval(Trans::eval(poly[1]));
        Arg_type c2   = broadcast<Arg_type>::eval(Trans::eval(poly[2]));
                        
        return small_horner(x, c0, c1, c2, res);
    }

    force_inline
    static Arg_type eval2(const Arg_type& x, const Coef_type* poly)
    {
        Arg_type c0   = broadcast<Arg_type>::eval(Trans::eval(poly[0]));
        Arg_type c1   = broadcast<Arg_type>::eval(Trans::eval(poly[1]));
        Arg_type c2   = broadcast<Arg_type>::eval(Trans::eval(poly[2]));

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
        Arg_type c0   = broadcast<Arg_type>::eval(Trans::eval(poly[0]));
        Arg_type c1   = broadcast<Arg_type>::eval(Trans::eval(poly[1]));

        return small_horner(x, c0, c1, res);
    }

    force_inline
    static Arg_type eval2(const Arg_type& x, const Coef_type* poly)
    {
        Arg_type c0   = broadcast<Arg_type>::eval(Trans::eval(poly[0]));
        Arg_type c1   = broadcast<Arg_type>::eval(Trans::eval(poly[1]));

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
        Arg_type c0   = broadcast<Arg_type>::eval(Trans::eval(poly[0]));
        return small_horner(x, c0, res);
    }

    force_inline
    static Arg_type eval2(const Arg_type& x, const Coef_type* poly)
    {
        Arg_type c0   = broadcast<Arg_type>::eval(Trans::eval(poly[0]));

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

        Arg_type cl = broadcast<Arg_type>::eval(poly[Poly_size - 1]);
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
            Arg_type c  = broadcast<Arg_type>::eval(Trans::eval(poly[i]));
            res = eval_fma<Arg_type, Arg_type>::eval(res, x, c);
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
        res             = broadcast<Arg_type>::eval((poly[size-1]));
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
    Arg_type res    = details::broadcast<Arg_type>::eval(poly[Poly_size - 1]);

    return details::eval_horner<details::trans_id, Poly_size - 1, Arg_type, Coef_type>
                        ::eval(x, res, poly);
}

template<class Arg_type, class Coef_type>
Arg_type simd::horner(const Arg_type& x, int poly_size, const Coef_type* poly)
{
    Arg_type res = details::broadcast<Arg_type>::eval(poly[poly_size - 1]);

    return details::eval_horner2<details::trans_id, Arg_type, Coef_type>
                ::eval(x, res, poly_size - 1, poly);
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
