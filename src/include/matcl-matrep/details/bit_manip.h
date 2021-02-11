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

#include "matcl-matrep/details/mpl.h"

namespace matcl 
{

// different routines for maniputation of bits stored in size_t type
template<class Sizet_type, int Bits = sizeof(Sizet_type)>
struct bit_manip
{
    static_assert(details::dependent_false<Sizet_type>::value, "not implemented");
};

template<class Sizet_type>
struct bit_manip<Sizet_type, 4>
{
    // clear all bits on positions pos and higher
    static Sizet_type   bits_before_pos(Sizet_type bits, Sizet_type pos);

    // return number of bits set
    static Sizet_type   count_bits(Sizet_type bits);

    // position of the most significant bit (starting from zero), at least one
    // bit must be set
    static Sizet_type   most_significant_bit_pos(Sizet_type bits);

    // position of the least significant bit (starting from zero), at least one
    // bit must be set
    static Sizet_type   least_significant_bit_pos(Sizet_type bits);

    // clear all bits except the most significant one
    static Sizet_type   most_significant_bit(Sizet_type bits);

    // clear all bits except the least significant one
    static Sizet_type   least_significant_bit(Sizet_type bits);
};

template<class Sizet_type>
struct bit_manip<Sizet_type, 8>
{
    // clear all bits on positions pos and higher
    static Sizet_type   bits_before_pos(Sizet_type bits, Sizet_type pos);

    // return number of bits set
    static Sizet_type   count_bits(Sizet_type bits);

    // position of the most significant bit (starting from zero), at least one
    // bit must be set
    static Sizet_type   most_significant_bit_pos(Sizet_type bits);

    // position of the least significant bit (starting from zero), at least one
    // bit must be set
    static Sizet_type   least_significant_bit_pos(Sizet_type bits);

    // clear all bits except the most significant one
    static Sizet_type   most_significant_bit(Sizet_type bits);

    // clear all bits except the least significant one
    static Sizet_type   least_significant_bit(Sizet_type bits);
};

};

#include "matcl-matrep/details/bit_manip.inl"
