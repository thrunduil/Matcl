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

#include "matcl-simd/details/arch/simd_impl.h"
#include "matcl-simd/details/func/simd_func_def.h"

namespace matcl { namespace simd { namespace details
{

template<int Bits, class Tag>
struct simd_sub_add<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function sub_add not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_sub_add<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function sub_add not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fma_f<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fma_f not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fma_f<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fma_f not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fms_f<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fms_f not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fms_f<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fms_f not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fnma_f<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fnma_f not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fnma_f<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fnma_f not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fnms_f<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fnms_f not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fnms_f<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fnms_f not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fma_a<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fma_a not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fma_a<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fma_a not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fms_a<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fms_a not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fms_a<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fms_a not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fnma_a<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fnma_a not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fnma_a<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fnma_a not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fnms_a<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fnms_a not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_fnms_a<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function fnms_a not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_sqrt<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function sqrt not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_sqrt<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function sqrt not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_any_nan<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function any_nan not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_any_nan<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function any_nan not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_is_nan<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function is_nan not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_is_nan<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function is_nan not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_div<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function div not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_div<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function div not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_pow2k<int32_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function pow2k not defined for integer arguments");
};

template<int Bits, class Tag>
struct simd_pow2k<int64_t, Bits, Tag>
{
    static_assert(md::dependent_false<Tag>::value, 
                "function pow2k not defined for integer arguments");
};

}}}
