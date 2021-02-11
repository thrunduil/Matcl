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

#include "matcl-matrep/details/bit_manip.h"

namespace matcl 
{

//-----------------------------------------------------------------
//                      32-bit version
//-----------------------------------------------------------------
template<class Sizet_type>
inline Sizet_type 
bit_manip<Sizet_type, 4>::bits_before_pos(Sizet_type bits, Sizet_type pos)
{
    size_t mask     = ~(0xffffffff << pos);
    size_t sel      = bits & mask;
    return sel;
};

template<class Sizet_type>
inline Sizet_type 
bit_manip<Sizet_type, 4>::count_bits(Sizet_type bits)
{
    //taken from The Aggregate Magic Algorithms
    bits = bits - ((bits >> 1U) & 0x55555555);
    bits = (bits & 0x33333333) + ((bits >> 2) & 0x33333333);
    bits = (((bits + (bits >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;

    return bits;
}

template<class Sizet_type>
inline Sizet_type 
bit_manip<Sizet_type, 4>::most_significant_bit(Sizet_type x)
{
    //taken from The Aggregate Magic Algorithms
    x       |= (x >> 1);
    x       |= (x >> 2);
    x       |= (x >> 4);
    x       |= (x >> 8);
    x       |= (x >> 16);

    return(x & ~(x >> 1));
}

#pragma warning(push)
#pragma warning(disable: 4146)  // unary minus operator applied to unsigned type,
                                // result still unsigned

template<class Sizet_type>
inline Sizet_type 
bit_manip<Sizet_type, 4>::least_significant_bit(Sizet_type bits)
{
    //taken from The Aggregate Magic Algorithms
    return bits & -bits;
};

#pragma warning(pop)

template<class Sizet_type>
inline Sizet_type
bit_manip<Sizet_type, 4>::most_significant_bit_pos(Sizet_type x)
{
    //taken from The Aggregate Magic Algorithms
    x       |= (x >> 1);
    x       |= (x >> 2);
    x       |= (x >> 4);
    x       |= (x >> 8);
    x       |= (x >> 16);

    return count_bits(x) - 1;
};

template<class Sizet_type>
inline Sizet_type 
bit_manip<Sizet_type, 4>::least_significant_bit_pos(Sizet_type x)
{
    //taken from The Aggregate Magic Algorithms
    x       |= (x << 1);
    x       |= (x << 2);
    x       |= (x << 4);
    x       |= (x << 8);
    x       |= (x << 16);

    return 32 - count_bits(x);
};

//-----------------------------------------------------------------
//                      64-bit version
//-----------------------------------------------------------------
template<class Sizet_type>
inline Sizet_type 
bit_manip<Sizet_type, 8>::bits_before_pos(Sizet_type bits, Sizet_type pos)
{
	size_t mask = ~(0xffffffffffffffff << pos);
	size_t sel = bits & mask;
	return sel;
};

template<class Sizet_type>
inline Sizet_type 
bit_manip<Sizet_type, 8>::count_bits(Sizet_type x)
{
	const uint64_t m1 = 0x5555555555555555; //binary: 0101...
	const uint64_t m2 = 0x3333333333333333; //binary: 00110011..
	const uint64_t m4 = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
	const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

	x   -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
	x   = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
	x   = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
	return (x * h01) >> 56;           //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
}

#pragma warning(push)
#pragma warning(disable: 4146)  // unary minus operator applied to unsigned type, 
                                // result still unsigned

//suppress this warning on 32-bit builds
#pragma warning(disable: 4293) // shift count negative or too big, undefined behavior

template<class Sizet_type>
inline Sizet_type 
bit_manip<Sizet_type, 8>::most_significant_bit(Sizet_type x)
{
	//taken from The Aggregate Magic Algorithms
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	x |= (x >> 8);
	x |= (x >> 16);
	x |= (x >> 32);

	return(x & ~(x >> 1));
}

template<class Sizet_type>
inline Sizet_type
bit_manip<Sizet_type, 8>::least_significant_bit(Sizet_type bits)
{
	//taken from The Aggregate Magic Algorithms
	return bits & -bits;
};

template<class Sizet_type>
inline Sizet_type
bit_manip<Sizet_type, 8>::most_significant_bit_pos(Sizet_type x)
{
	//taken from The Aggregate Magic Algorithms
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	x |= (x >> 8);
	x |= (x >> 16);
	x |= (x >> 32);

	return count_bits(x) - 1;
};

template<class Sizet_type>
inline Sizet_type
bit_manip<Sizet_type, 8>::least_significant_bit_pos(Sizet_type x)
{
	//taken from The Aggregate Magic Algorithms
	x |= (x << 1);
	x |= (x << 2);
	x |= (x << 4);
	x |= (x << 8);
	x |= (x << 16);
	x |= (x << 32);

	return 64 - count_bits(x);
};

#pragma warning(pop)

};
