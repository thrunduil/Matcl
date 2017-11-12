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

namespace matcl { namespace simd { namespace details
{

template<class Float_type, class ... Args>
struct eval_short_horner
{};

template<class Float_type, class Arg1, class ... Args>
struct eval_short_horner<Float_type, Arg1, Args...>
{
    force_inline
    static Float_type eval(Float_type x, Arg1 coef1, Args ... coef)
    {
        Float_type res  = eval_short_horner<Float_type, Args...>::eval(x, coef...);
        res             = eval_fma<Float_type>::eval(res, x, coef1);
        return res;
    }
};

template<class Float_type, class Arg1>
struct eval_short_horner<Float_type, Arg1>
{
    force_inline
    static Float_type eval(Float_type x, Arg1 coef1)
    {
        (void)x;
        return coef1;
    }
};

template<class Float_type>
struct eval_fma
{
    force_inline
    static Float_type eval(Float_type x, Float_type y, Float_type z)
    {
        return x * y + z;
        //return fma_f(x, y, z);
    }
};

template<int Poly_size, class Float_type>
struct eval_horner
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        static const int size_1 = Poly_size/2;
        static const int size_2 = Poly_size - size_1;

        res = eval_horner<size_1, Float_type>::eval(x, res, poly + size_2);
        res = eval_horner<size_2, Float_type>::eval(x, res, poly);
        return res;
    }
};

template<class Float_type>
struct eval_horner<15, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], poly[2], poly[3], poly[4], 
                            poly[5], poly[6], poly[8], poly[9], poly[10], 
                            poly[11], poly[12], poly[14], poly[15], res);
    }
};

template<class Float_type>
struct eval_horner<14, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], poly[2], poly[3], poly[4], 
                            poly[5], poly[6], poly[8], poly[9], poly[10], 
                            poly[11], poly[12], poly[14], res);
    }
};

template<class Float_type>
struct eval_horner<13, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], poly[2], poly[3], poly[4], 
                            poly[5], poly[6], poly[8], poly[9], poly[10], 
                            poly[11], poly[12], poly[13], res);
    }
};

template<class Float_type>
struct eval_horner<12, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], poly[2], poly[3], poly[4], 
                            poly[5], poly[6], poly[8], poly[9], poly[10], 
                            poly[11], poly[12], res);
    }
};

template<class Float_type>
struct eval_horner<11, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], poly[2], poly[3], poly[4], 
                            poly[5], poly[6], poly[8], poly[9], poly[10], 
                            poly[11], res);
    }
};

template<class Float_type>
struct eval_horner<10, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], poly[2], poly[3], poly[4], 
                            poly[5], poly[6], poly[8], poly[9], poly[10], res);
    }
};

template<class Float_type>
struct eval_horner<9, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], poly[2], poly[3], poly[4], 
                            poly[5], poly[6], poly[8], poly[9], res);
    }
};

template<class Float_type>
struct eval_horner<8, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], poly[2], poly[3], poly[4], 
                            poly[5], poly[6], poly[8], res);
    }
};

template<class Float_type>
struct eval_horner<7, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], poly[2], poly[3], poly[4], 
                            poly[5], poly[6], res);
    }
};

template<class Float_type>
struct eval_horner<6, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], poly[2], poly[3], poly[4], 
                            poly[5], res);
    }
};

template<class Float_type>
struct eval_horner<5, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], poly[2], poly[3], poly[4], res);
    }
};

template<class Float_type>
struct eval_horner<4, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], poly[2], poly[3], res);
    }
};

template<class Float_type>
struct eval_horner<3, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], poly[2], res);
    }
};

template<class Float_type>
struct eval_horner<2, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], poly[1], res);
    }
};

template<class Float_type>
struct eval_horner<1, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        return short_horner(x, poly[0], res);
    }
};

template<class Float_type>
struct eval_horner<0, Float_type>
{
    static Float_type eval(Float_type x, Float_type res, const Float_type* poly)
    {
        (void)x;
        (void)poly;
        return res;
    }
};

/*
//TODO
template<class Float_type>
struct eval_horner<2, Float_type>
{
    static Float_type eval(Float_type x, const Float_type* poly)
    {
        using FT = Float_type;
        return eval_short_horner<FT, FT,FT>::eval(x, poly[0], poly[1]);
    }
};

template<class Float_type>
struct eval_horner<3, Float_type>
{
    static Float_type eval(Float_type x, const Float_type* poly)
    {
        using FT = Float_type;
        return eval_short_horner<FT, FT,FT,FT>::eval(x, poly[0], poly[1], poly[2]);
    }
};

template<class Float_type>
struct eval_horner<4, Float_type>
{
    static Float_type eval(Float_type x, const Float_type* poly)
    {
        using FT = Float_type;
        return eval_short_horner<FT, FT,FT,FT,FT>::eval(x, poly[0], poly[1], poly[2], poly[3]);
    }
};
*/

}}}

namespace matcl
{

template<class Float_type, class ... Args>
Float_type simd::short_horner(Float_type x, Args ... coef)
{
    return details::eval_short_horner<Float_type, Args ...>::eval(x, coef...);
};

template<int Poly_size, class Float_type>
Float_type simd::horner(Float_type x, const Float_type* poly)
{
    Float_type res = poly[Poly_size - 1];
    return details::eval_horner<Poly_size - 1, Float_type>::eval(x, res, poly);
}

#if 0
template<int Poly_size, class Float_type>
Float_type simd::horner(const Float_type* poly, Float_type x)
{
    Float_type res  = poly[Poly_size - 1];

    for (int i = Poly_size - 2; i >= 0; --i)
        res = fma_f(res, x, poly[i]);

    return res;
};

template<int Poly_size, class Float_type>
Float_type simd::horner_abs(const Float_type* poly, Float_type x)
{
    Float_type res  = std::abs(poly[Poly_size - 1]);

    for (int i = Poly_size - 2; i >= 0; --i)
        res = fma_f(res, x, std::abs(poly[i]));

    return res;
};

template<int Poly_size, class Float_type>
Float_type simd::twofold_horner(const Float_type* poly, Float_type x,
                twofold<Float_type>* err)
{
    using twofold_type = twofold<Float_type>;

    Float_type res  = poly[Poly_size - 1];

    for (int i = Poly_size - 2; i >= 0; --i)
    {
        twofold_type prod   = twofold_mult(res, x);
        pi[i]               = prod.error;
        twofold_type sum    = twofold_sum(prod.value, poly[i]);

        sig[i]              = sum.error;
        res                 = sum.value;
    }

    return res;
};

template<int Poly_size, class Float_type>
Float_type simd::poly_cond(const Float_type* poly, Float_type x)
{
    Float_type num  = horner_abs<Poly_size>(poly, std::abs(x));
    Float_type den  = horner<Poly_size>(poly, x);    

    return num / std::abs(den);
};

#endif
}

#include "matcl-simd/details/poly/poly_eval.inl"
