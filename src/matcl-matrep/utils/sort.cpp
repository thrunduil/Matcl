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

#include "matcl-internals/base/sort.h"
#include "matcl-scalar/object.h"
#include "matcl-matrep/objects/details/type_info_object.h"
#include "matcl-core/details/integer.h"

#include <algorithm>
#include <stack>

namespace matcl { namespace utils
{

static const Integer SORT_MAX  = 50;

template <class Type> 
inline Type * partition(Type *pF, Type *pL, const Type pivot)
{
    for (pF--; ;) 
    {
        while (*(++pF) < pivot);
        while (*(--pL) > pivot);
        if (pF>=pL) { return pF; }
        std::iter_swap(pF, pL);
    }
}

template <typename Type> 
inline const Type& get_median(const Type& x, const Type& y, const Type& z)
{
    if (x < y)	{ return (y < z ? y : x < z ? z : x); }
    else		{ return (x < z ? x : y < z ? z : y); }
}

template <typename Type> 
inline const Type& get_pivot(const Type *pF, Integer n)
{	
    return get_median(pF[0], pF[n/2], pF[n-1]);	// using first, mid and last
}

template <typename Type> 
inline void sort_q(Type *pF, Type *pL)
{
    Integer n = matcl::cast_int64(pL-pF);
    if (n <= 1) 
        return;

    if (n <= SORT_MAX) 
    {
        sort_ins(pF, n);
    } 
    else 
    {
        Type *pM = partition(pF, pL, get_pivot(pF, n));
        sort_q(pF, pM);
        sort_q(pM, pL);
    };
};

template<class T>
inline void swap(T& arg1, T& arg2)
{
    T tmp = arg1;
    arg1 = arg2;
    arg2 = tmp;
};

template<>
inline void swap(Object& arg1, Object& arg2)
{
    swap(arg1,arg2);
};

template <class Ttype>
void sort_q(Ttype *ptr, Integer n)
{
    if (n <= 1)
        return;

    if (n <= SORT_MAX) 
    {
        sort_ins(ptr, n);
        return;
    };

    return sort_q(ptr,ptr+n);
}

template <class Ttype>
void sort_q_rev(Ttype *ptr, Integer n)
{
    if (n <= 1)
        return;

    if (n <= SORT_MAX) 
    {
        sort_ins_rev(ptr, n);
        return;
    };

    sort_q(ptr,ptr+n);
    std::reverse(ptr,ptr+n);
}


template <class Ttype1, class Ttype2>
void sort_q(Ttype1 *ptr1, Ttype2 *ptr2, Integer n)
{
    if (n <= 1) 
        return;

    if (n <= SORT_MAX) 
    {
        sort_ins(ptr1,ptr2, n);
        return;
    };

    //TODO: this may break when pivot is bad
    Ttype1 x1 = get_pivot(ptr1, n);
    Integer i,j;

    ti::ti_type<Ttype1> ti_1 = ti::get_ti(ptr1[0]);
    ti::ti_type<Ttype2> ti_2 = ti::get_ti(ptr2[0]);
    
    for (i = 0, j = n - 1; i <= j;)
    {
        while (ptr1[i] < x1) ++i;
        while (ptr1[j] > x1) --j;
        if (i <= j)
        {
            swap(ptr1[i],ptr1[j]);
            swap(ptr2[i],ptr2[j]);
            ++i;
            --j;
        }
    }
    if (j > 0)		sort_q(ptr1, ptr2, j + 1);
    if (i < n - 1)	sort_q(ptr1 + i, ptr2 + i, n - i);
};

template <class Ttype1, class Ttype2>
void sort_q_rev(Ttype1 *ptr1, Ttype2 *ptr2, Integer n)
{
    if (n <= 1) 
        return;

    if (n <= SORT_MAX) 
    {
        sort_ins_rev(ptr1,ptr2, n);
        return;
    };

    sort_q(ptr1, ptr2, n);
    std::reverse(ptr1,ptr1+n);
    std::reverse(ptr2,ptr2+n);
};

template <class Ttype1, class Ttype2, class Ttype3>
void sort_q(Ttype1 *ptr1, Ttype2 *ptr2, Ttype3 *ptr3, Integer n)
{
    if (n <= 1) 
        return;

    if (n <= SORT_MAX) 
    {
        sort_ins(ptr1,ptr2, ptr3, n);
        return;
    };

    Ttype1 x1 = get_pivot(ptr1, n);
    Integer i,j;

    for (i = 0, j = n - 1; i <= j;)
    {
        while (ptr1[i] < x1) ++i;
        while (ptr1[j] > x1) --j;
        if (i <= j)
        {
            swap(ptr1[i],ptr1[j]);
            swap(ptr2[i],ptr2[j]);
            swap(ptr3[i],ptr3[j]);

            ++i;
            --j;
        }
    }
    if (j > 0)		sort_q(ptr1, ptr2,ptr3, j + 1);
    if (i < n - 1)	sort_q(ptr1 + i, ptr2 + i, ptr3 + i, n - i);
}

};};


template MATCL_MATREP_EXPORT void matcl::utils::sort_ins(Integer* const, const Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_ins(Real* const, const Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_ins(Float* const, const Integer);

template MATCL_MATREP_EXPORT void matcl::utils::sort_q(Integer*, Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_q(Float*, Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_q(Real*, Integer);

template MATCL_MATREP_EXPORT void matcl::utils::sort_q_rev(Integer*, Integer);

template MATCL_MATREP_EXPORT void matcl::utils::sort_ins(Integer*, Integer*, Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_ins(Float*, Integer*, Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_ins(Real*, Integer*, Integer);

template MATCL_MATREP_EXPORT void matcl::utils::sort_ins(Integer*, Float*, Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_ins(Integer*, Real*, Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_ins(Integer*, Float_complex*, Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_ins(Integer*, Complex*, Integer);

template MATCL_MATREP_EXPORT void matcl::utils::sort_q(Integer*, Integer*, Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_q(Float*, Integer*, Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_q(Real*, Integer*, Integer);

template MATCL_MATREP_EXPORT void matcl::utils::sort_q(Integer*, Float*, Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_q(Integer*, Real*, Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_q(Integer*, Float_complex*, Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_q(Integer*, Complex*, Integer);
template MATCL_MATREP_EXPORT void matcl::utils::sort_q(Integer*, Object*, Integer);

template MATCL_MATREP_EXPORT void matcl::utils::sort_q(Integer*, Integer*, Integer*, Integer);

template MATCL_MATREP_EXPORT void matcl::utils::sort_q_rev(Integer*, Integer*, Integer);
