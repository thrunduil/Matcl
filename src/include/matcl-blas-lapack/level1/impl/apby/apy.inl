/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2021
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

#include "matcl-blas-lapack/level1/level1_apby.h"

namespace matcl { namespace level1
{

template<class TY, Integer Rows>
struct apy<TY, Rows, false, details::true_t>
{
    force_inline
    static void eval(TY* Y1, Integer rows, const TY& a)
    {
        Integer size =  get_int<Rows>::eval(rows);

        TY* MATCL_RESTRICTED Y  = Y1;

		for (Integer i = 0; i < size; ++i)
			Y[i] = a + Y[i];
    };

    force_inline
    static void eval(TY* Y1, Integer Y_step, Integer rows, const TY& a)
    {
        if (Y_step == 1)
            return eval(Y1,rows,a);

        Integer size =  get_int<Rows>::eval(rows);

        TY* MATCL_RESTRICTED Y  = Y1;

		for (Integer i = 0; i < size; ++i)
			Y[i*Y_step] = a + Y[i*Y_step];
    };
};

template<class TY, Integer Rows>
struct apy<TY, Rows, true, details::true_t>
{
    force_inline
    static void eval(TY* Y1, Integer rows, const TY& a)
    {
        using simd_type = typename simd::default_simd_type<TY>::type;
        static const int vec_size   = simd_type::vector_size;

        Integer size   =  get_int<Rows>::eval(rows);
        Integer size_1 =  get_int<Rows>::eval(rows) / vec_size;

        TY* MATCL_RESTRICTED Y = Y1;

        simd_type sa(a);

		for (Integer i = 0; i < size_1; ++i)
        {
            simd_type sY_1  = simd_type::load(Y+i*vec_size, std::false_type());
			simd_type res_1 = sa + sY_1;

            res_1.store(Y+i*vec_size, std::false_type());
        };

        for (Integer i = size_1*vec_size; i < size; ++i)
            Y[i]        = a  + Y[i];
    };

    force_inline
    static void eval(TY* Y1, Integer Y_step, Integer rows, const TY& a)
    {
        if (Y_step == 1)
            return eval(Y1,rows,a);

        return apy<TY,Rows,false,details::true_t>::eval(Y1,Y_step,rows,a);
    };
};

}}