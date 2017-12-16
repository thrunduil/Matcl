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

template<int Poly_size, class Arg, class Coef>
struct eval_estrin
{
    static const int log2       = eval_log2(Poly_size);
    static const int max_pow    = get_estrin_max_size<Poly_size, log2>::max_pow;
    static const int size1      = get_estrin_max_size<Poly_size, log2>::size;
    static const int size2      = Poly_size - size1;

    static_assert(size2 > 0, "invalid size");

    force_inline
    static Arg eval(const Arg& x, const Coef* poly)
    {
        Arg xpow[max_pow];
        xpow[0] = x;

        for (int i = 1; i < max_pow; ++i)
            xpow[i] = xpow[i-1] * xpow[i-1];

        return fma(eval_estrin<size2, Arg, Coef>::eval_rec(xpow, poly + size1), 
                       xpow[max_pow-1], 
                       eval_estrin<size1, Arg, Coef>::eval_rec(xpow, poly));
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        return fma(eval_estrin<size2, Arg, Coef>::eval_rec(xpow, poly + size1), 
                       xpow[max_pow-1], 
                       eval_estrin<size1, Arg, Coef>::eval_rec(xpow, poly));
    }
};

template<class Arg, class Coef>
struct eval_estrin<9, Arg, Coef>
{
    force_inline
    static Arg eval(const Arg& x, const Coef* poly)
    {
        Arg x2  = x * x;
        Arg x4  = x2 * x2;

        Arg p11 = broadcast<Arg>::eval(poly + 8);
        Arg p12 = broadcast<Arg>::eval(poly + 5);
        Arg p21 = broadcast<Arg>::eval(poly + 3);
        Arg p22 = broadcast<Arg>::eval(poly + 1);
        
        p11     = fma(p11, x, broadcast<Arg>::eval(poly + 7));
        p12     = fma(p12, x, broadcast<Arg>::eval(poly + 4));
        p21     = fma(p21, x, broadcast<Arg>::eval(poly + 2));
        p22     = fma(p22, x, broadcast<Arg>::eval(poly + 0));

        p11     = fma(p11, x, broadcast<Arg>::eval(poly + 6));

        Arg p2  = fma(p21, x2, p22);
        Arg p1  = fma(p11, x2, p12);
        
        Arg p   = fma(p1, x4, p2);

        return p;
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg x   = xpow[0];
        Arg x2  = xpow[1];
        Arg x4  = xpow[2];

        Arg p11 = broadcast<Arg>::eval(poly + 8);
        Arg p12 = broadcast<Arg>::eval(poly + 5);
        Arg p21 = broadcast<Arg>::eval(poly + 3);
        Arg p22 = broadcast<Arg>::eval(poly + 1);
        
        p11     = fma(p11, x, broadcast<Arg>::eval(poly + 7));
        p12     = fma(p12, x, broadcast<Arg>::eval(poly + 4));
        p21     = fma(p21, x, broadcast<Arg>::eval(poly + 2));
        p22     = fma(p22, x, broadcast<Arg>::eval(poly + 0));

        p11     = fma(p11, x, broadcast<Arg>::eval(poly + 6));

        Arg p2  = fma(p21, x2, p22);
        Arg p1  = fma(p11, x2, p12);
        
        Arg p   = fma(p1, x4, p2);
        return p;
    }
};

template<class Arg, class Coef>
struct eval_estrin<8, Arg, Coef>
{
    force_inline
    static Arg eval(const Arg& x, const Coef* poly)
    {
        Arg x2  = x * x;
        Arg x4  = x2 * x2;

        Arg p11 = broadcast<Arg>::eval(poly + 7);
        Arg p12 = broadcast<Arg>::eval(poly + 5);
        Arg p21 = broadcast<Arg>::eval(poly + 3);
        Arg p22 = broadcast<Arg>::eval(poly + 1);
        
        p11     = fma(p11, x, broadcast<Arg>::eval(poly + 6));
        p12     = fma(p12, x, broadcast<Arg>::eval(poly + 4));
        p21     = fma(p21, x, broadcast<Arg>::eval(poly + 2));
        p22     = fma(p22, x, broadcast<Arg>::eval(poly + 0));

        Arg p1  = fma(p11, x2, p12);
        Arg p2  = fma(p21, x2, p22);

        Arg p   = fma(p1, x4, p2);

        return p;
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg x   = xpow[0];
        Arg x2  = xpow[1];
        Arg x4  = xpow[2];

        Arg p11 = broadcast<Arg>::eval(poly + 7);
        Arg p12 = broadcast<Arg>::eval(poly + 5);
        Arg p21 = broadcast<Arg>::eval(poly + 3);
        Arg p22 = broadcast<Arg>::eval(poly + 1);
        
        p11     = fma(p11, x, broadcast<Arg>::eval(poly + 6));
        p12     = fma(p12, x, broadcast<Arg>::eval(poly + 4));
        p21     = fma(p21, x, broadcast<Arg>::eval(poly + 2));
        p22     = fma(p22, x, broadcast<Arg>::eval(poly + 0));

        Arg p1  = fma(p11, x2, p12);
        Arg p2  = fma(p21, x2, p22);

        Arg p   = fma(p1, x4, p2);
        return p;
    }
};

template<class Arg, class Coef>
struct eval_estrin<7, Arg, Coef>
{
    force_inline
    static Arg eval(const Arg& x, const Coef* poly)
    {
        Arg x2  = x * x;
        Arg x4  = x2 * x2;

        Arg p1  = broadcast<Arg>::eval(poly + 6);
        Arg p21 = broadcast<Arg>::eval(poly + 3);
        Arg p22 = broadcast<Arg>::eval(poly + 1);
        
        p1      = fma(p1, x, broadcast<Arg>::eval(poly + 5));
        p21     = fma(p21, x, broadcast<Arg>::eval(poly + 2));
        p22     = fma(p22, x, broadcast<Arg>::eval(poly + 0));

        p1      = fma(p1, x, broadcast<Arg>::eval(poly + 4));

        Arg p2  = fma(p21, x2, p22);
        Arg p   = fma(p1, x4, p2);

        return p;
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg x   = xpow[0];
        Arg x2  = xpow[1];
        Arg x4  = xpow[2];

        Arg p1  = broadcast<Arg>::eval(poly + 6);
        Arg p21 = broadcast<Arg>::eval(poly + 3);
        Arg p22 = broadcast<Arg>::eval(poly + 1);
        
        p1      = fma(p1, x, broadcast<Arg>::eval(poly + 5));
        p21     = fma(p21, x, broadcast<Arg>::eval(poly + 2));
        p22     = fma(p22, x, broadcast<Arg>::eval(poly + 0));

        p1      = fma(p1, x, broadcast<Arg>::eval(poly + 4));

        Arg p2  = fma(p21, x2, p22);
        Arg p   = fma(p1, x4, p2);
        return p;
    }
};

template<class Arg, class Coef>
struct eval_estrin<6, Arg, Coef>
{
    force_inline
    static Arg eval(const Arg& x, const Coef* poly)
    {
        Arg x2  = x * x;
        Arg x4  = x2 * x2;

        Arg p1  = broadcast<Arg>::eval(poly + 5);
        Arg p21 = broadcast<Arg>::eval(poly + 3);
        Arg p22 = broadcast<Arg>::eval(poly + 1);
        
        p1      = fma(p1,  x, broadcast<Arg>::eval(poly + 4));
        p21     = fma(p21, x, broadcast<Arg>::eval(poly + 2));
        p22     = fma(p22, x, broadcast<Arg>::eval(poly + 0));

        Arg p2  = fma(p21, x2, p22);
        Arg p   = fma(p1, x4, p2);

        return p;
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg x   = xpow[0];
        Arg x2  = xpow[1];
        Arg x4  = xpow[2];

        Arg p1  = broadcast<Arg>::eval(poly + 5);
        Arg p21 = broadcast<Arg>::eval(poly + 3);
        Arg p22 = broadcast<Arg>::eval(poly + 1);
        
        p1      = fma(p1, x, broadcast<Arg>::eval(poly + 4));
        p21     = fma(p21, x, broadcast<Arg>::eval(poly + 2));
        p22     = fma(p22, x, broadcast<Arg>::eval(poly + 0));

        Arg p2  = fma(p21, x2, p22);
        Arg p   = fma(p1, x4, p2);

        return p;
    }
};

template<class Arg, class Coef>
struct eval_estrin<5, Arg, Coef>
{
    force_inline
    static Arg eval(const Arg& x, const Coef* poly)
    {        
        Arg x2  = x * x;

        Arg p1  = broadcast<Arg>::eval(poly + 4);
        Arg p2  = broadcast<Arg>::eval(poly + 1);
        
        p1      = fma(p1, x, broadcast<Arg>::eval(poly + 3));
        p2      = fma(p2, x, broadcast<Arg>::eval(poly + 0));

        p1      = fma(p1, x, broadcast<Arg>::eval(poly + 2));

        Arg p   = fma(p1, x2, p2);
        return p;
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg x   = xpow[0];
        Arg x2  = xpow[1];

        Arg p1  = broadcast<Arg>::eval(poly + 4);
        Arg p2  = broadcast<Arg>::eval(poly + 1);
        
        p1      = fma(p1, x, broadcast<Arg>::eval(poly + 3));
        p2      = fma(p2, x, broadcast<Arg>::eval(poly + 0));

        p1      = fma(p1, x, broadcast<Arg>::eval(poly + 2));

        Arg p   = fma(p1, x2, p2);
        return p;
    }
};

template<class Arg, class Coef>
struct eval_estrin<4, Arg, Coef>
{
    force_inline
    static Arg eval(const Arg& x, const Coef* poly)
    {
        Arg x2  = x * x;

        Arg p1  = broadcast<Arg>::eval(poly + 3);
        Arg p2  = broadcast<Arg>::eval(poly + 1);
        
        p1      = fma(p1, x, broadcast<Arg>::eval(poly + 2));        
        p2      = fma(p2, x, broadcast<Arg>::eval(poly + 0));

        Arg p   = fma(p1, x2, p2);
        return p;
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg x   = xpow[0];
        Arg x2  = xpow[1];

        Arg p1  = broadcast<Arg>::eval(poly + 3);
        Arg p2  = broadcast<Arg>::eval(poly + 1);
        
        p1      = fma(p1, x, broadcast<Arg>::eval(poly + 2));        
        p2      = fma(p2, x, broadcast<Arg>::eval(poly + 0));

        Arg p   = fma(p1, x2, p2);
        return p;
    }
};

template<class Arg, class Coef>
struct eval_estrin<3, Arg, Coef>
{
    force_inline
    static Arg eval(const Arg& x, const Coef* poly)
    {
        Arg p   = broadcast<Arg>::eval(poly + 2);
        p       = fma(p, x, broadcast<Arg>::eval(poly + 1));
        p       = fma(p, x, broadcast<Arg>::eval(poly + 0));
        return p;
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg x    = xpow[0];

        Arg p   = broadcast<Arg>::eval(poly + 2);
        p       = fma(p, x, broadcast<Arg>::eval(poly + 1));
        p       = fma(p, x, broadcast<Arg>::eval(poly + 0));
        return p;

    }
};

template<class Arg, class Coef>
struct eval_estrin<2, Arg, Coef>
{
    force_inline
    static Arg eval(const Arg& x, const Coef* poly)
    {
        Arg p   = broadcast<Arg>::eval(poly + 1);
        p       = fma(p, x, broadcast<Arg>::eval(poly + 0));
        return p;
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        Arg x    = xpow[0];

        Arg p   = broadcast<Arg>::eval(poly + 1);
        p       = fma(p, x, broadcast<Arg>::eval(poly + 0));
        return p;
    }
};

template<class Arg, class Coef>
struct eval_estrin<1, Arg, Coef>
{
    force_inline
    static Arg eval(const Arg& x, const Coef* poly)
    {
        (void)x;
        return broadcast<Arg>::eval(poly + 0);
    }

    force_inline
    static Arg eval_rec(const Arg* xpow, const Coef* poly)
    {
        (void)xpow;
        return broadcast<Arg>::eval(poly + 0);
    }
};

//-----------------------------------------------------------------------
//                      ESTRIN DYNAMIC
//-----------------------------------------------------------------------
template<class Arg, class Coef>
struct eval_estrin2
{
    static const int max_pow            = 64;
    static const int horner_threshold   = 16;

    force_inline
    static Arg eval(const Arg& x, int size, const Coef* poly)
    {
        if (size < horner_threshold)
            return horner(x, size, poly);

        Arg xpow[max_pow];
        xpow[0] = x;

        int pow, level;
        get_pow_level(size, pow, level);

        for (int i = 1; i <= level; ++i)
            xpow[i] = xpow[i-1] * xpow[i-1];

        int size1   = pow;
        int size2   = size - size1;

        if (size1 == size2)
        {
            return fma(eval_estrin2<Arg, Coef>::eval_rec2(xpow, size1, level - 1, poly + size1), 
                       xpow[level], 
                       eval_estrin2<Arg, Coef>::eval_rec2(xpow, size1, level - 1, poly));
        }
        else
        {
            return fma(eval_estrin2<Arg, Coef>::eval_rec(xpow, size2, poly + size1), 
                       xpow[level], 
                       eval_estrin2<Arg, Coef>::eval_rec2(xpow, size1, level - 1, poly));
        };
    }

    static Arg eval_rec(const Arg* xpow, int size, const Coef* poly)
    {        
        if (size < horner_threshold)
            return horner(xpow[0], size, poly);
        
        int pow, level;
        get_pow_level(size, pow, level);

        int size1   = pow;
        int size2   = size - size1;

        if (size1 == size2)
        {
            return fma(eval_estrin2<Arg, Coef>::eval_rec2(xpow, size1, level - 1, poly + size1), 
                       xpow[level], 
                       eval_estrin2<Arg, Coef>::eval_rec2(xpow, size1, level - 1, poly));
        }
        else
        {
            return fma(eval_estrin2<Arg, Coef>::eval_rec(xpow, size2, poly + size1), 
                       xpow[level], 
                       eval_estrin2<Arg, Coef>::eval_rec2(xpow, size1, level - 1, poly));
        };
    }

    static Arg eval_rec2(const Arg* xpow, int size, int level, const Coef* poly)
    {
        switch(level)
        {
            case 0:
                return eval_estrin<2, Arg, Coef>::eval_rec(xpow, poly);
            case 1:
                return eval_estrin<4, Arg, Coef>::eval_rec(xpow, poly);
            case 2:
                return eval_estrin<8, Arg, Coef>::eval_rec(xpow, poly);
            case 3:
                return eval_estrin<16, Arg, Coef>::eval_rec(xpow, poly);
            case 4:
                return eval_estrin<32, Arg, Coef>::eval_rec(xpow, poly);

            case 5:
            {
                int size_part   = size/2;
                Arg r[2];
                
                const Coef* p   = poly;

                for (int i = 1; i >= 0; --i)
                {
                    r[i]    = eval_estrin<32, Arg, Coef>::eval_rec(xpow, p);
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
                    r[i]    = eval_estrin<32, Arg, Coef>::eval_rec(xpow, p);
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
                    r[i]    = eval_estrin<32, Arg, Coef>::eval_rec(xpow, p);
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

        Arg r1  = eval_estrin2<Arg, Coef>
                    ::eval_rec2(xpow, size_half, level - 1, poly + size_half);

        Arg r2  = eval_estrin2<Arg, Coef>
                    ::eval_rec2(xpow, size_half, level - 1, poly);

        return fma(r1, xpow[level], r2);
    }

    static void get_pow_level(int size, int& pow, int& level)
    {
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

template<int Poly_size, class Arg_type, class Coef_type>
force_inline
Arg_type simd::estrin(const Arg_type& x, const Coef_type* poly)
{
    return details::eval_estrin<Poly_size, Arg_type, Coef_type>
                        ::eval(x, poly);
}

template<class Arg_type, class Coef_type>
Arg_type simd::estrin(const Arg_type& x, int poly_size, const Coef_type* poly)
{
    return details::eval_estrin2<Arg_type, Coef_type>
                ::eval(x, poly_size, poly);
}

}
