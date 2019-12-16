/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Paweł Kowal 2019
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

#include "matcl-core/details/mpl.h"
#include "mkgen/matrix/matrix.h"

namespace matcl { namespace mkgen { namespace colon_func
{

//------------------------------------------------------------------------------
//                      check_colon
//------------------------------------------------------------------------------
// check if colon is in range of a matrix with M elements
template<class Colon, Integer M>
struct check_colon
{
    static_assert(md::dependent_false<Colon>::value, "unknown colon type");
};

template<Integer M>
struct check_colon<colon_all, M>
{
    static_assert(M >= 0, "invalid number of elements");
    using type  = void;
};

template<Integer M, Integer M1>
struct check_colon<colon<M1>, M>
{
    static_assert(M >= 0, "invalid number of elements");
    static_assert(M1 >= 1 && M1 <= M, "colon out of range");
    
    using type  = void;
};

template<Integer M, Integer Start, Integer End>
struct check_colon<colon2<Start, End>, M>
{
    static_assert(M >= 0, "invalid number of elements");

    static_assert(Start >= 1 && Start <= M, "colon out of range");
    static_assert(End >= 1 && End <= M, "colon out of range");

    using type = void;
};

template<Integer M, Integer Start, Integer Step, Integer End>
struct check_colon<colon3<Start, Step, End>, M>
{
    static_assert(M >= 0, "invalid number of elements");

    static_assert(Start >= 1 && Start <= M, "colon out of range");
    static_assert(End >= 1 && End <= M, "colon out of range");
    static_assert(Step != 0, "invalid step");

    using type = void;
};

//------------------------------------------------------------------------------
//                      size
//------------------------------------------------------------------------------

// number of elements selected by the colon, when matrix has M elements
template<class Colon, Integer M>
struct size
{
    static_assert(md::dependent_false<Colon>::value, "unknown colon type");
};

template<Integer M>
struct size<colon_all, M>
{
    using dum   = typename check_colon<colon_all, M>::type;
    static const Integer value  = M;
};

template<Integer M, Integer M1>
struct size<colon<M1>, M>
{
    using dum   = typename check_colon<colon<M1>, M>::type;

    static const Integer value  = 1;
};

template<Integer M, Integer Start, Integer End>
struct size<colon2<Start, End>, M>
{
    using dum   = typename check_colon<colon2<Start, End>, M>::type;

    static const Integer value0 = End - Start + 1;
    static const Integer value  = (value0 < 0 ? 0 : value0);
};

template<Integer M, Integer Start, Integer Step, Integer End>
struct size<colon3<Start, Step, End>, M>
{
    using dum   = typename check_colon<colon3<Start, Step, End>, M>::type;

    static const Integer value0 = (End - Start) / Step + 1;
    static const Integer value  = (value0 < 0 ? 0 : value0);
};

//------------------------------------------------------------------------------
//                      first
//------------------------------------------------------------------------------
// index of the first element selected by the colon, when matrix has M elements
template<class Colon, Integer M>
struct first
{
    static_assert(md::dependent_false<Colon>::value, "unknown colon type");
};

template<Integer M>
struct first<colon_all, M>
{
    using dum   = typename check_colon<colon_all, M>::type;
    static const Integer value  = 1;
};

template<Integer M, Integer M1>
struct first<colon<M1>, M>
{
    using dum   = typename check_colon<colon<M1>, M>::type;
    static const Integer value  = M1;
};

template<Integer M, Integer Start, Integer End>
struct first<colon2<Start,End>, M>
{
    using dum   = typename check_colon<colon2<Start,End>, M>::type;
    static const Integer value  = Start;
};

template<Integer M, Integer Start, Integer Step, Integer End>
struct first<colon3<Start, Step, End>, M>
{
    using dum   = typename check_colon<colon3<Start, Step, End>, M>::type;
    static const Integer value  = Start;
};

//------------------------------------------------------------------------------
//                      last
//------------------------------------------------------------------------------
// index of the last element selected by the colon, when matrix has M elements
template<class Colon, Integer M>
struct last
{
    static_assert(md::dependent_false<Colon>::value, "unknown colon type");
};

template<Integer M>
struct last<colon_all, M>
{
    using dum   = typename check_colon<colon_all, M>::type;
    static const Integer value  = M;
};

template<Integer M, Integer M1>
struct last<colon<M1>, M>
{
    using dum   = typename check_colon<colon<M1>, M>::type;
    static const Integer value  = M1;
};

template<Integer M, Integer Start, Integer End>
struct last<colon2<Start,End>, M>
{
    using dum   = typename check_colon<colon2<Start,End>, M>::type;
    static const Integer value  = End;
};

template<Integer M, Integer Start, Integer Step, Integer End>
struct last<colon3<Start, Step, End>, M>
{
    using dum   = typename check_colon<colon3<Start, Step, End>, M>::type;

    static const Integer size   = colon_func::size<colon3<Start, Step, End>,M>::value;
    static const Integer value  = Start + (size-1) * Step;
};

//------------------------------------------------------------------------------
//                      offset
//------------------------------------------------------------------------------
// return offset of the first element from beginning of an array (0-based)
template<class Colon>
struct offset
{
    static_assert(md::dependent_false<Colon>::value, "unknown colon type");
};

template<>
struct offset<colon_all>
{
    static const Integer value  = 0;
};

template<Integer M1>
struct offset<colon<M1>>
{
    static const Integer value  = M1 - 1;
};

template<Integer Start, Integer End>
struct offset<colon2<Start,End>>
{
    static const Integer value  = Start - 1;
};

template<Integer Start, Integer Step, Integer End>
struct offset<colon3<Start, Step, End>>
{
    static const Integer value  = Start - 1;
};

//------------------------------------------------------------------------------
//                      step
//------------------------------------------------------------------------------
// return distance between two consecutive elements pointed by the colon
template<class Colon>
struct step
{
    static_assert(md::dependent_false<Colon>::value, "unknown colon type");
};

template<>
struct step<colon_all>
{
    static const Integer value  = 1;
};

template<Integer M1>
struct step<colon<M1>>
{
    static const Integer value  = 1;
};

template<Integer Start, Integer End>
struct step<colon2<Start,End>>
{
    static const Integer value  = 1;
};

template<Integer Start, Integer Step, Integer End>
struct step<colon3<Start, Step, End>>
{
    static const Integer value  = Step;
};

//------------------------------------------------------------------------------
//                      index
//------------------------------------------------------------------------------
// gen index of K-th element selected by the colon (index and K is 1-based)
template<Integer K, class Colon>
struct index
{
    static_assert(md::dependent_false<Colon>::value,
                  "this type should not be instantiated");
};

template<Integer K>
struct index<K,colon_all>
{
    static_assert(K >= 1, "invalid element");
    static const Integer value = K;
};

template<Integer K, Integer Start>
struct index<K, colon<Start>>
{
    static_assert(K == 1, "invalid element");
    static const Integer value = Start;
};

template<Integer K, Integer Start, Integer End>
struct index<K, colon2<Start, End>>
{
    static_assert(K >= 1, "invalid element");    
    static const Integer value = Start + (K - 1);

    static_assert(value <= End, "invalid element");
};

template<Integer K, Integer Start, Integer Step, Integer End>
struct index<K, colon3<Start, Step, End>>
{
    static_assert(K >= 1, "invalid element");    
    static const Integer value = Start + Step * (K - 1);

    static_assert((Step > 0 && value <= End) || (Step < 0 && value >= End), 
                  "invalid element");
};

}}}
