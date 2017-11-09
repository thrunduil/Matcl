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

#include "matcl-simd/arch/simd_impl.h"
#include "matcl-simd/func/simd_func_def.h"
#include "matcl-simd/details/scalar_func.h"

namespace matcl { namespace simd
{

template<class T>
struct simd_reverse<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(reverse(x.data[1]), reverse(x.data[0]));
    };
};

template<class T>
struct simd_mult<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(x.data[0] * y.data[0], x.data[1] * y.data[1]);
    };
};

template<class T>
struct simd_div<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(x.data[0] / y.data[0], x.data[1] / y.data[1]);
    };
};

template<class T>
struct simd_plus<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(x.data[0] + y.data[0], x.data[1] + y.data[1]);
    };
};

template<class T>
struct simd_minus<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(x.data[0] - y.data[0], x.data[1] - y.data[1]);
    };
};

template<class T>
struct simd_uminus<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(-x.data[0], -x.data[1]);
    };
};

template<class T>
struct simd_sum_all<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static T eval(const simd_type& x)
    {
        return sum_all(x.data[0]) + sum_all(x.data[1]);
    };
};

template<class T>
struct simd_sub_add<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(sub_add(x.data[0], y.data[0]), sub_add(x.data[1], y.data[1]));
    };
};
template<class T>
struct simd_fma_f<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(fma_f(x.data[0], y.data[0], z.data[0]), 
                         fma_f(x.data[1], y.data[1], z.data[1]));
    };
};

template<class T>
struct simd_fms_f<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(fms_f(x.data[0], y.data[0], z.data[0]), 
                         fms_f(x.data[1], y.data[1], z.data[1]));
    };
};

template<class T>
struct simd_abs<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(abs(x.data[0]), abs(x.data[1]));
    };
};

template<class T>
struct simd_max<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(max(x.data[0], y.data[0]), max(x.data[1], y.data[1]));
    };
};

template<class T>
struct simd_min<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(min(x.data[0], y.data[0]), min(x.data[1], y.data[1]));
    };
};

template<class T>
struct simd_sqrt<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(sqrt(x.data[0]), sqrt(x.data[1]));
    };
};

template<class T>
struct simd_round<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(round(x.data[0]), round(x.data[1]));
    };
};

template<class T>
struct simd_floor<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(floor(x.data[0]), floor(x.data[1]));
    };
};

template<class T>
struct simd_ceil<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(ceil(x.data[0]), ceil(x.data[1]));
    };
};

template<class T>
struct simd_trunc<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(trunc(x.data[0]), trunc(x.data[1]));
    };
};

template<class T>
struct simd_eeq<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(eeq(x.data[0], y.data[0]), eeq(x.data[1], y.data[1]));
    };
};

template<class T>
struct simd_neq<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(neq(x.data[0], y.data[0]), neq(x.data[1], y.data[1]));
    };
};

template<class T>
struct simd_lt<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(lt(x.data[0], y.data[0]), lt(x.data[1], y.data[1]));
    };
};

template<class T>
struct simd_gt<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(gt(x.data[0], y.data[0]), gt(x.data[1], y.data[1]));
    };
};

template<class T>
struct simd_leq<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(leq(x.data[0], y.data[0]), leq(x.data[1], y.data[1]));
    };
};

template<class T>
struct simd_geq<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(geq(x.data[0], y.data[0]), geq(x.data[1], y.data[1]));
    };
};

template<class T>
struct simd_any_nan<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        bool b1 = any_nan(x.extract_low());
        bool b2 = any_nan(x.extract_high());

        return b1 || b2;
    };
};

template<class T>
struct simd_any<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        bool b1 = any(x.extract_low());
        bool b2 = any(x.extract_high());

        return b1 || b2;
    };
};

template<class T>
struct simd_all<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        bool b1 = all(x.extract_low());
        bool b2 = all(x.extract_high());

        return b1 && b2;
    };
};

}}
