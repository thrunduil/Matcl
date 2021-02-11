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

#include "matcl-simd/details/arch/simd_impl.h"
#include "matcl-simd/details/complex/simd_complex_impl.h"
#include "matcl-simd/details/func/simd_func_complex_def.h"
#include "matcl-simd/details/complex/recover_nan.h"

namespace matcl { namespace simd { namespace details
{

template<class T>
struct simd_compl_reverse<T, 256, sse_tag>
{
    using simd_type = simd_compl<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(reverse(x.extract_high()), reverse(x.extract_low()));
    };
};

template<class T>
struct simd_compl_conj<T, 256, sse_tag>
{
    using simd_type = simd_compl<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(conj(x.extract_low()), conj(x.extract_high()));
    };
};

template<class T>
struct simd_compl_mult<T, 256, sse_tag>
{
    using simd_type         = simd_compl<T, 256, sse_tag>;
    using simd_real_type    = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {        
        return simd_type(x.extract_low() * y.extract_low(), x.extract_high() * y.extract_high());
    };

    force_inline
    static simd_type eval_rc(const simd_real_type& x, const simd_type& y)
    {
        return simd_type(x.extract_low() * y.extract_low(), x.extract_high() * y.extract_high());
    }

    force_inline
    static simd_type eval_cr(const simd_type& x, const simd_real_type& y)
    {
        return simd_type(x.extract_low() * y.extract_low(), x.extract_high() * y.extract_high());
    }
};

template<class T>
struct simd_compl_div<T, 256, sse_tag>
{
    using simd_type = simd_compl<T, 256, sse_tag>;
    using simd_real = simd<T, 256, sse_tag>;

    // (a.re * b.re + a.im * b.im, b.re * a.im - a.re * b.im) / (b.re * b.re + b.im * b.im)
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {        
        return simd_type(x.extract_low() / y.extract_low(), x.extract_high() / y.extract_high());
    };

    // (a.re * b.re + a.im * b.im, b.re * a.im - a.re * b.im) / (b.re * b.re + b.im * b.im)
    force_inline
    static simd_type eval(const simd_real& x, const simd_type& y)
    {        
        return simd_type(x.extract_low() / y.extract_low(), x.extract_high() / y.extract_high());
    };

    force_inline
    static simd_type eval_rc(const simd_real& x, const simd_type& y)
    {        
        return simd_type(x.extract_low() / y.extract_low(), x.extract_high() / y.extract_high());
    };

    force_inline
    static simd_type eval_cr(const simd_type& x, const simd_real& y)
    {
        return simd_type(x.extract_low() / y.extract_low(), x.extract_high() / y.extract_high());
    }
};

template<class T>
struct simd_compl_plus<T, 256, sse_tag>
{
    using simd_type = simd_compl<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(x.data.data[0] + y.data.data[0], x.data.data[1] + y.data.data[1]);
    };
};

template<class T>
struct simd_compl_minus<T, 256, sse_tag>
{
    using simd_type = simd_compl<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(x.data.data[0] - y.data.data[0], x.data.data[1] - y.data.data[1]);
    };
};

template<class T>
struct simd_compl_uminus<T, 256, sse_tag>
{
    using simd_type = simd_compl<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(-x.data.data[0], -x.data.data[1]);
    };
};

template<class T>
struct simd_compl_horizontal_sum<T, 256, sse_tag>
{
    using simd_type = simd_compl<T, 256, sse_tag>;

    force_inline
    static simd_double_complex eval(const simd_type& x)
    {
        return horizontal_sum(x.extract_low() + x.extract_high());
    };
};

}}}
