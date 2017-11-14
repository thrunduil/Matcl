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
template<class Trans, class Float_type, class ... Args>
struct eval_short_horner
{};

template<class Trans, class Float_type, class Arg1, class ... Args>
struct eval_short_horner<Trans, Float_type, Arg1, Args...>
{
    force_inline
    static Float_type eval(Float_type x, Arg1 coef1, Args ... coef)
    {
        Float_type res  = eval_short_horner<Trans, Float_type, Args...>
                            ::eval(x, coef...);
        Arg1 coef_tr    = Trans::eval(coef1);
        res             = eval_fma<Float_type, Arg1>::eval(res, x, coef_tr);
        return res;
    }
};

template<class Trans, class Float_type, class Arg1>
struct eval_short_horner<Trans, Float_type, Arg1>
{
    force_inline
    static Float_type eval(Float_type x, Arg1 coef1)
    {
        (void)x;
        Arg1 coef_tr    = Trans::eval(coef1);
        return Float_type(coef_tr);
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

        return short_horner(x, c0, c1, c2, c3, res);
    }

    force_inline
    static Arg_type eval2(Arg_type x, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        Arg_type c1   = Arg_type(Trans::eval(poly[1]));
        Arg_type c2   = Arg_type(Trans::eval(poly[2]));
        Arg_type c3   = Arg_type(Trans::eval(poly[3]));

        return short_horner(x, c0, c1, c2, c3);
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

        return short_horner(x, c0, c1, c2, res);
    }

    force_inline
    static Arg_type eval2(Arg_type x, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        Arg_type c1   = Arg_type(Trans::eval(poly[1]));
        Arg_type c2   = Arg_type(Trans::eval(poly[2]));

        return short_horner(x, c0, c1, c2);
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

        return short_horner(x, c0, c1, res);
    }

    force_inline
    static Arg_type eval2(Arg_type x, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));
        Arg_type c1   = Arg_type(Trans::eval(poly[1]));

        return short_horner(x, c0, c1);
    }
};

template<class Trans, class Arg_type, class Coef_type>
struct eval_horner<Trans, 1, Arg_type, Coef_type>
{
    force_inline
    static Arg_type eval(Arg_type x, Arg_type res, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));

        return short_horner(x, c0, res);
    }

    force_inline
    static Arg_type eval2(Arg_type x, const Coef_type* poly)
    {
        Arg_type c0   = Arg_type(Trans::eval(poly[0]));

        return short_horner(x, c0);
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
//                      HORNER TWOFOLD
//-----------------------------------------------------------------------
template<class Trans, class Arg_type, class Coef_type>
struct eval_horner2_twofold
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
            twofold_type sum    = twofold_plus(prod.value, Arg_type(c));

            Arg_type err        = prod.error + sum.error;
            comp                = eval_fma<Arg_type, Arg_type>::eval(comp, x, err);
            res                 = sum.value;
        }

        return res + comp;
    };
};

}}}

namespace matcl
{

template<class Float_type, class ... Args>
force_inline
Float_type simd::short_horner(Float_type x, Args ... coef)
{
    return details::eval_short_horner<details::trans_id, Float_type, Args ...>
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
force_inline
Arg_type simd::horner(Arg_type x, int poly_size, const Coef_type* poly)
{
    Arg_type res = Arg_type(poly[poly_size - 1]);

    return details::eval_horner2<details::trans_id, Arg_type, Coef_type>
                ::eval(x, res, poly_size - 1, poly);
};

template<class Arg_type, class Coef_type>
force_inline
Arg_type simd::horner_abs(Arg_type x, int poly_size, const Coef_type* poly)
{
    Arg_type res = Arg_type(poly[poly_size - 1]);
    res          = details::eval_abs<Arg_type>::eval(res);

    return details::eval_horner2<details::trans_abs, Arg_type, Coef_type>
                ::eval(x, res, poly_size - 1, poly);
};

template<class Arg_type, class Coef_type>
force_inline
Arg_type simd::twofold_horner(Arg_type x, int poly_size, const Coef_type* poly)
{
    Arg_type res = Arg_type(poly[poly_size - 1]);

    return details::eval_horner2_twofold<details::trans_id, Arg_type, Coef_type>
                ::eval(x, res, poly_size - 1, poly);
};

template<class Arg_type, class Coef_type>
Arg_type simd::horner_cond(Arg_type x, int poly_size, const Coef_type* poly)
{
    Arg_type abs_x  = details::eval_abs<Arg_type>::eval(x);
    Arg_type num    = horner_abs(abs_x, poly_size, poly);
    Arg_type den    = horner(x, poly_size, poly);    
    Arg_type aden   = details::eval_abs<Arg_type>::eval(den);

    return num / aden;
};

}
