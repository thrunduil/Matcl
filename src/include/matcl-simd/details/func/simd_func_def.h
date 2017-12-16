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

#include "matcl-simd/simd_general.h"

namespace matcl { namespace simd { namespace details
{

namespace ms = matcl::simd;

template<class Val, int Bits, class Simd_tag>
struct simd_reverse
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function reverse not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_signbit_base
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function signbit_base not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_bitwise_or
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function bitwise_or not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_bitwise_and
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function bitwise_and not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_bitwise_xor
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function bitwise_xor not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_bitwise_andnot
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function bitwise_andnot not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_bitwise_not
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function bitwise_not not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_shift_left
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function shift_left not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_shift_right
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function shift_right not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_shift_right_arithmetic
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function shift_right_arithmetic not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_mult
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function mult not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_div
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function div not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_plus
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function plus not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_minus
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function minus not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_uminus
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function uminus not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_abs
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function abs not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_sub_add
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function sub_add not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_fma_f
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function fma_f not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_fms_f
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function fms_f not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_fma_a
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function fma_a not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_fms_a
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function fms_a not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_fnma_f
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function fnma_f not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_fnms_f
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function fnms_f not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_fnma_a
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function fnma_a not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_fnms_a
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function fnms_a not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_sum_all
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function sum_all not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_max
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function max not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_min
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function min not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_eeq
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function eeq not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_neq
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function neq not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_lt
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function lt not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_gt
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function gt not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_leq
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function leq not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_geq
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function geq not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_is_nan
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function is_nan not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_sqrt
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function sqrt not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_round
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function round not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_floor
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function floor not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_ceil
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function ceil not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_trunc
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function trunc not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_any_nan
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function any_nan not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_any
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function any not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_all
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function all not defined for for given arguments");
};

template<class Val, int Bits, class Simd_tag>
struct simd_if_then_else
{
    static_assert(md::dependent_false<Simd_tag>::value, 
                "function if_then_else not defined for for given arguments");
};

}}}
