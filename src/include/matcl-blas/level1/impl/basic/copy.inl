/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-blas/level1/level1_basic.h"

#pragma warning(push)
#pragma warning(disable: 4127)  //conditional expression is constant

namespace matcl { namespace level1
{

template<class T1, class T2, Integer Rows>
struct copy<T1, T2, Rows, false, details::true_t>
{
    force_inline
    static void eval(T1* dest0, const T2* source0, Integer rows)
    {
        Integer size =  get_int<Rows>::eval(rows);

        T1* MATCL_RESTRICTED dest         = dest0;
        const T2* MATCL_RESTRICTED source = source0;

        for (Integer pos = 0; pos < size; ++pos)
            dest[pos] = T1(source[pos]);
    };

    force_inline
    static void eval_step2(T1* dest0, const T2* source0, Integer source_step, Integer rows)
    {
        Integer size =  get_int<Rows>::eval(rows);

        T1* MATCL_RESTRICTED dest         = dest0;
        const T2* MATCL_RESTRICTED source = source0;

        for (Integer pos = 0; pos < size; ++pos)
            dest[pos] = T1(source[pos*source_step]);
    };

    force_inline
    static void eval_step1(T1* dest0, Integer dest_step, const T2* source0, Integer rows)
    {
        Integer size =  get_int<Rows>::eval(rows);

        T1* MATCL_RESTRICTED dest         = dest0;
        const T2* MATCL_RESTRICTED source = source0;

        for (Integer pos = 0; pos < size; ++pos)
            dest[pos*dest_step] = T1(source[pos]);
    };

    force_inline
    static void eval(T1* dest0, Integer dest_step, const T2* source0, Integer source_step, Integer rows)
    {
        if (dest_step == 1)
        {
            if (source_step == 1)
                return eval(dest0, source0, rows);
            else
                return eval_step2(dest0, source0, source_step, rows);
        }
        else if (source_step == 1)                
        {
            return eval_step1(dest0, dest_step, source0, rows);
        };

        Integer size =  get_int<Rows>::eval(rows);

        T1* MATCL_RESTRICTED dest         = dest0;
        const T2* MATCL_RESTRICTED source = source0;

        for (Integer pos = 0; pos < size; ++pos)
            dest[pos*dest_step] = T1(source[pos*source_step]);
    };
};

template<class T, Integer Rows>
struct copy<T, T, Rows, true, details::true_t>
{
    force_inline
    static void eval(T* dest0, const T* source0, Integer rows)
    {
        using simd_type = typename simd::default_simd_type<T>::type;
        static const int vec_size   = simd_type::vector_size;

        Integer size   =  get_int<Rows>::eval(rows);
        Integer size_1 =  get_int<Rows>::eval(rows) / vec_size;

        T* MATCL_RESTRICTED Y        = dest0;
        const T* MATCL_RESTRICTED X  = source0;

		for (Integer i = 0; i < size_1; ++i)
        {
            simd_type sX_1  = simd_type::load(X+i*vec_size, std::false_type());
            sX_1.store(Y+i*vec_size, std::false_type());
        };

        for (Integer i = size_1*vec_size; i < size; ++i)
            Y[i]        = X[i];
    };

    force_inline
    static void eval(T* dest0, Integer dest_step, const T* source0, Integer source_step, Integer rows)
    {
        if (dest_step == 1 && source_step == 1)
            return eval(dest0, source0, rows);

        return copy<T,T, Rows, false, details::true_t>::eval(dest0, dest_step, source0, source_step, rows);
    };
};

//-------------------------------------------------------------
//              remove Continuous == 0
//-------------------------------------------------------------
template<bool Select, class T1, class T2, Integer Rows, Integer Cols, bool Use_simd>
struct copy_mat<Select, T1, T2, Rows, Cols, 0, Use_simd, details::true_t>
{
    static void eval(T1* dest, Integer dest_ld, const T2* source, 
                     Integer source_ld, Integer rows, Integer cols)
    {
        if (Cols == 1)
        {
            return copy_mat<Select, T1, T2, Rows, 1, 1, Use_simd, details::true_t>
                    ::eval(dest, dest_ld, source, source_ld, rows, cols);
        };

        if (Select == false)
        {
            return copy_mat<false, T1, T2, Rows, Cols, 2, Use_simd, details::true_t>
                ::eval(dest, dest_ld, source, source_ld, rows, cols);
        }

        if (dest_ld == get_int<Rows>::eval(rows) && source_ld == get_int<Rows>::eval(rows))
        {
            return copy_mat<true, T1, T2, Rows, Cols, 1, Use_simd, details::true_t>
                    ::eval(dest, dest_ld, source, source_ld, rows, cols);
        }
        else
        {
            return copy_mat<true, T1, T2, Rows, Cols, 2, Use_simd, details::true_t>
                    ::eval(dest, dest_ld, source, source_ld, rows, cols);
        }    
    };
};

//-------------------------------------------------------------
//              Continuous == 1
//-------------------------------------------------------------
template<bool Select, class T1, class T2, Integer Rows, bool Use_simd, Integer Cols>
struct copy_mat<Select, T1, T2, Rows, Cols, 1, Use_simd, details::true_t>
{
    static void eval(T1* dest, Integer dest_ld, const T2* source,
                     Integer source_ld, Integer rows, Integer cols)
    {
        (void)source_ld;
        (void)dest_ld;
        Integer size =  get_int<Rows>::eval(rows) * get_int<Cols>::eval(cols);
        copy<T1,T2,Rows*Cols>::eval(dest,source,size);
    };
};

//-------------------------------------------------------------
//              Continuous == 2
//-------------------------------------------------------------
template<class T1, class T2, Integer Cols, bool Use_simd>
struct copy_mat<true, T1, T2, 0, Cols, 2, Use_simd, details::true_t>
{
    static void eval(T1* dest, Integer dest_ld, const T2* source,
                     Integer source_ld, Integer rows, Integer cols)
    {
        switch(rows)
        {
            case 1:
                return copy_mat<true, T1, T2, 1, Cols, 2, Use_simd, details::true_t>
                    ::eval(dest, dest_ld, source, source_ld, rows, cols);
            case 2:
                return copy_mat<true, T1, T2, 2, Cols, 2, Use_simd, details::true_t>
                    ::eval(dest, dest_ld, source, source_ld, rows, cols);
            case 3:
                return copy_mat<true, T1, T2,  3, Cols, 2, Use_simd, details::true_t>
                    ::eval(dest, dest_ld, source, source_ld, rows, cols);
            case 4:
                return copy_mat<true, T1, T2,  4, Cols, 2, Use_simd, details::true_t>
                    ::eval(dest, dest_ld, source, source_ld, rows, cols);
            default:
                return copy_mat<false, T1, T2,  0, Cols, 2, Use_simd, details::true_t>
                    ::eval(dest, dest_ld, source, source_ld, rows, cols);
        };
    };
};

template<bool Select, class T1, class T2, Integer Rows, Integer Cols, bool Use_simd>
struct copy_mat<Select, T1, T2, Rows, Cols, 2, Use_simd, details::true_t>
{
    static void eval(T1* dest, Integer dest_ld, const T2* source,
                     Integer source_ld, Integer rows, Integer cols)
    {
        for (Integer i = 0 ; i < get_int<Cols>::eval(cols); ++i)
        {
            copy<T1, T2, Rows>::eval(dest, source, rows);

            dest        += dest_ld;
            source      += source_ld;
        }
    };
};

}}

#pragma warning(pop)