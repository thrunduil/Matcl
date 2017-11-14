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

template<int Poly_size, int Max_pow>
struct get_estrin_max_size
{
    static const int exp    = eval_pow2(Max_pow);
    static const int size   = (exp == Poly_size)? exp/2 : exp;
    static const int max_pow= (exp == Poly_size)? Max_pow : Max_pow + 1;
};

//-----------------------------------------------------------------------
//                   ESTRIN
//-----------------------------------------------------------------------

template<class Trans, int Poly_size, class Arg, class Coef>
struct eval_estrin
{
    static const int log2       = eval_log2(Poly_size);
    static const int max_pow    = get_estrin_max_size<Poly_size, log2>::max_pow;
    static const int size1      = get_estrin_max_size<Poly_size, log2>::size;
    static const int size2      = Poly_size - size1;

    static_assert(size2 > 0, "invalid size");

    force_inline
    static Arg eval(Arg x, const Coef* poly)
    {
        Arg xpow[max_pow];
        xpow[0] = x;

        for (int i = 1; i < max_pow; ++i)
            xpow[i] = xpow[i-1] * xpow[i-1];

        return fma(eval_estrin<Trans, size2, Arg, Coef>::eval_rec(xpow, poly + size1), 
                       xpow[max_pow-1], 
                       eval_estrin<Trans, size1, Arg, Coef>::eval_rec(xpow, poly));
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        return fma(eval_estrin<Trans, size2, Arg, Coef>::eval_rec(xpow, poly + size1), 
                       xpow[max_pow-1], 
                       eval_estrin<Trans, size1, Arg, Coef>::eval_rec(xpow, poly));
    }
};

template<class Trans, class Arg, class Coef>
struct eval_estrin<Trans, 8, Arg, Coef>
{
    force_inline
    static Arg eval(Arg x, const Coef* poly)
    {
        Arg x2      = x * x;
        Arg x4      = x2 * x2;
        Arg xpow[3] = {x, x2, x4};
        
        return eval_rec(xpow, poly);
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg c0   = Arg(Trans::eval(poly[0]));
        Arg c1   = Arg(Trans::eval(poly[1]));
        Arg c2   = Arg(Trans::eval(poly[2]));
        Arg c3   = Arg(Trans::eval(poly[3]));
        Arg c4   = Arg(Trans::eval(poly[4]));
        Arg c5   = Arg(Trans::eval(poly[5]));
        Arg c6   = Arg(Trans::eval(poly[6]));
        Arg c7   = Arg(Trans::eval(poly[7]));

        Arg x    = xpow[0];
        Arg x2   = xpow[1];
        Arg x4   = xpow[2];

        return fma(fma(fma(c7, x, c6), x2, fma(c5, x, c4)), 
                   x4, 
                   fma(fma(c3, x, c2), x2, fma(c1, x, c0)));
    }
};

template<class Trans, class Arg, class Coef>
struct eval_estrin<Trans, 7, Arg, Coef>
{
    force_inline
    static Arg eval(Arg x, const Coef* poly)
    {
        Arg x2      = x * x;
        Arg x4      = x2 * x2;
        Arg xpow[3] = {x, x2, x4};
        
        return eval_rec(xpow, poly);
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg c0   = Arg(Trans::eval(poly[0]));
        Arg c1   = Arg(Trans::eval(poly[1]));
        Arg c2   = Arg(Trans::eval(poly[2]));
        Arg c3   = Arg(Trans::eval(poly[3]));
        Arg c4   = Arg(Trans::eval(poly[4]));
        Arg c5   = Arg(Trans::eval(poly[5]));
        Arg c6   = Arg(Trans::eval(poly[6]));

        Arg x    = xpow[0];
        Arg x2   = xpow[1];
        Arg x4   = xpow[2];

        return fma(fma(c6, x2, fma(c5, x, c4)), 
                   x4, 
                   fma(fma(c3, x, c2), x2, fma(c1, x, c0)));
    }
};

template<class Trans, class Arg, class Coef>
struct eval_estrin<Trans, 6, Arg, Coef>
{
    force_inline
    static Arg eval(Arg x, const Coef* poly)
    {
        Arg x2      = x * x;
        Arg x4      = x2 * x2;
        Arg xpow[3] = {x, x2, x4};
        
        return eval_rec(xpow, poly);
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg c0   = Arg(Trans::eval(poly[0]));
        Arg c1   = Arg(Trans::eval(poly[1]));
        Arg c2   = Arg(Trans::eval(poly[2]));
        Arg c3   = Arg(Trans::eval(poly[3]));
        Arg c4   = Arg(Trans::eval(poly[4]));
        Arg c5   = Arg(Trans::eval(poly[5]));

        Arg x    = xpow[0];
        Arg x2   = xpow[1];
        Arg x4   = xpow[2];

        return fma(fma(c5, x, c4), 
                   x4, 
                   fma(fma(c3, x, c2), x2, fma(c1, x, c0)));
    }
};

template<class Trans, class Arg, class Coef>
struct eval_estrin<Trans, 5, Arg, Coef>
{
    force_inline
    static Arg eval(Arg x, const Coef* poly)
    {        
        Arg x2      = x * x;
        Arg x4      = x2 * x2;
        Arg xpow[3] = {x, x2, x4};
        Arg res     = eval_rec(xpow, poly);

        return res;
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg c0   = Arg(Trans::eval(poly[0]));
        Arg c1   = Arg(Trans::eval(poly[1]));
        Arg c2   = Arg(Trans::eval(poly[2]));
        Arg c3   = Arg(Trans::eval(poly[3]));
        Arg c4   = Arg(Trans::eval(poly[4]));

        Arg x    = xpow[0];
        Arg x2   = xpow[1];
        Arg x4   = xpow[2];

        return fma(c4, x4, fma(fma(c3, x, c2), x2, fma(c1, x, c0)));
    }
};

template<class Trans, class Arg, class Coef>
struct eval_estrin<Trans, 4, Arg, Coef>
{
    force_inline
    static Arg eval(Arg x, const Coef* poly)
    {
        Arg c0   = Arg(Trans::eval(poly[0]));
        Arg c1   = Arg(Trans::eval(poly[1]));
        Arg c2   = Arg(Trans::eval(poly[2]));
        Arg c3   = Arg(Trans::eval(poly[3]));

        Arg x2   = x * x;

        return fma(fma(c3, x, c2), x2, fma(c1, x, c0));
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg c0   = Arg(Trans::eval(poly[0]));
        Arg c1   = Arg(Trans::eval(poly[1]));
        Arg c2   = Arg(Trans::eval(poly[2]));
        Arg c3   = Arg(Trans::eval(poly[3]));

        Arg x    = xpow[0];
        Arg x2   = xpow[1];

        return fma(fma(c3, x, c2), x2, fma(c1, x, c0));
    }
};

template<class Trans, class Arg, class Coef>
struct eval_estrin<Trans, 3, Arg, Coef>
{
    force_inline
    static Arg eval(Arg x, const Coef* poly)
    {
        Arg c0   = Arg(Trans::eval(poly[0]));
        Arg c1   = Arg(Trans::eval(poly[1]));
        Arg c2   = Arg(Trans::eval(poly[2]));

        Arg x2   = x * x;

        return fma(c2, x2, fma(c1, x, c0));
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg c0   = Arg(Trans::eval(poly[0]));
        Arg c1   = Arg(Trans::eval(poly[1]));
        Arg c2   = Arg(Trans::eval(poly[2]));

        Arg x    = xpow[0];
        Arg x2   = xpow[1];

        return fma(c2, x2, fma(c1, x, c0));
    }
};

template<class Trans, class Arg, class Coef>
struct eval_estrin<Trans, 2, Arg, Coef>
{
    force_inline
    static Arg eval(Arg x, const Coef* poly)
    {
        Arg c0   = Arg(Trans::eval(poly[0]));
        Arg c1   = Arg(Trans::eval(poly[1]));

        return fma(c1, x, c0);
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg x    = xpow[0];
        Arg c0   = Arg(Trans::eval(poly[0]));
        Arg c1   = Arg(Trans::eval(poly[1]));

        return fma(c1, x, c0);
    }
};

template<class Trans, class Arg, class Coef>
struct eval_estrin<Trans, 1, Arg, Coef>
{
    force_inline
    static Arg eval(Arg x, const Coef* poly)
    {
        (void)x;
        Arg c0   = Arg(Trans::eval(poly[0]));
        return c0;
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        (void)xpow;
        Arg c0   = Arg(Trans::eval(poly[0]));
        return c0;
    }
};

//-----------------------------------------------------------------------
//                      ESTRIN DYNAMIC
//-----------------------------------------------------------------------
template<class Trans, class Arg, class Coef>
struct eval_estrin2
{
    //TODO
    static const int max_pow    = 10;

    force_inline
    static Arg eval(Arg x, int size, const Coef* poly)
    {
        Arg xpow[max_pow];
        xpow[0] = x;

        switch (size)
        {
            case 1:
                return Arg(Trans::eval(poly[0]));
            case 2:
                return eval_estrin<Trans, 2, Arg, Coef>::eval_rec(xpow, poly);
        };        
        
        int pow, level;
        get_pow_level(size, pow, level);

        for (int i = 1; i <= level; ++i)
            xpow[i] = xpow[i-1] * xpow[i-1];

        int size1   = pow;
        int size2   = size - size1;

        if (size1 == size2)
        {
            return fma(eval_estrin2<Trans, Arg, Coef>::eval_rec2(xpow, size1, level - 1, poly + size1), 
                       xpow[level], 
                       eval_estrin2<Trans, Arg, Coef>::eval_rec2(xpow, size1, level - 1, poly));
        }
        else
        {
            return fma(eval_estrin2<Trans, Arg, Coef>::eval_rec(xpow, size2, poly + size1), 
                       xpow[level], 
                       eval_estrin2<Trans, Arg, Coef>::eval_rec2(xpow, size1, level - 1, poly));
        };
    }

    static Arg eval_rec(const Arg* xpow, int size, const Coef* poly)
    {        
        switch (size)
        {
            case 1:
                return Arg(Trans::eval(poly[0]));
            case 2:
                return eval_estrin<Trans, 2, Arg, Coef>::eval_rec(xpow, poly);
        };        
        
        int pow, level;
        get_pow_level(size, pow, level);

        int size1   = pow;
        int size2   = size - size1;

        if (size1 == size2)
        {
            return fma(eval_estrin2<Trans, Arg, Coef>::eval_rec2(xpow, size1, level - 1, poly + size1), 
                       xpow[level], 
                       eval_estrin2<Trans, Arg, Coef>::eval_rec2(xpow, size1, level - 1, poly));
        }
        else
        {
            return fma(eval_estrin2<Trans, Arg, Coef>::eval_rec(xpow, size2, poly + size1), 
                       xpow[level], 
                       eval_estrin2<Trans, Arg, Coef>::eval_rec2(xpow, size1, level - 1, poly));
        };
    }

    static Arg eval_rec2(const Arg* xpow, int size, int level, const Coef* poly)
    {
        //this function is recursive and cannot be inlined

        switch(level)
        {
            case 0:
                return eval_estrin<Trans, 2, Arg, Coef>::eval_rec(xpow, poly);
            case 1:
                return eval_estrin<Trans, 4, Arg, Coef>::eval_rec(xpow, poly);
            case 2:
                return eval_estrin<Trans, 8, Arg, Coef>::eval_rec(xpow, poly);
            case 3:
                return eval_estrin<Trans, 16, Arg, Coef>::eval_rec(xpow, poly);
            case 4:
                return eval_estrin<Trans, 32, Arg, Coef>::eval_rec(xpow, poly);

            case 5:
            {
                int size_part   = size/2;
                Arg r[2];
                
                const Coef* p   = poly;

                for (int i = 1; i >= 0; --i)
                {
                    r[i]    = eval_estrin<Trans, 32, Arg, Coef>::eval_rec(xpow, p);
                    p       = p + size_part;
                }

                return fma(r[0], xpow[level], r[1]);
            }

            case 6:
            {
                int size_part   = size/4;
                Arg r[4];
                
                const Coef* p   = poly;

                for (int i = 3; i >= 0; --i)
                {
                    r[i]    = eval_estrin<Trans, 32, Arg, Coef>::eval_rec(xpow, p);
                    p       = p + size_part;
                }

                return fma(fma(r[0], xpow[level-1], r[1]), 
                           xpow[level], 
                           fma(r[2], xpow[level-1], r[3]));
            }

            case 7:
            {
                int size_part   = size/8;
                Arg r[8];
                
                const Coef* p   = poly;

                for (int i = 7; i >= 0; --i)
                {
                    r[i]    = eval_estrin<Trans, 32, Arg, Coef>::eval_rec(xpow, p);
                    p       = p + size_part;
                }

                Arg r1      = fma(fma(r[0], xpow[level-2], r[1]), 
                                xpow[level-1], 
                                fma(r[2], xpow[level-2], r[3]));

                Arg r2      = fma(fma(r[4], xpow[level-2], r[5]), 
                                xpow[level-1], 
                                fma(r[6], xpow[level-2], r[7]));

                return fma(r1, xpow[level], r2);
            }
        };

        int size_half   = size/2;

        Arg r1  = eval_estrin2<Trans, Arg, Coef>
                    ::eval_rec2(xpow, size_half, level - 1, poly + size_half);

        Arg r2  = eval_estrin2<Trans, Arg, Coef>
                    ::eval_rec2(xpow, size_half, level - 1, poly);

        return fma(r1, xpow[level], r2);
    }

    static void get_pow_level(int size, int& pow, int& level)
    {
        //TODO

        pow     = 2;
        level   = 1;

        while (pow < size)
        {
            level   += 1;
            pow     = pow * 2;
        };
        
        if (pow >= size)
        {
            pow     = pow / 2;
            level   = level - 1;
        };
    };
};

}}}

namespace matcl
{

template<int Poly_size, class Arg, class Coef>
force_inline
Arg simd::estrin(Arg x, const Coef* poly)
{
    return details::eval_estrin<details::trans_id, Poly_size, Arg, Coef>
                        ::eval(x, poly);
}

template<class Arg, class Coef>
Arg simd::estrin(Arg x, int poly_size, const Coef* poly)
{
    return details::eval_estrin2<details::trans_id, Arg, Coef>
                ::eval(x, poly_size, poly);
}

}
