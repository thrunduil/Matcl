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

//TODO: compare with sym_arrow
namespace matcl { namespace utils
{

template <class T>
void sort_ins(T * const ptr, const Integer n)
{
    for(Integer j = 1; j < n; ++j)
    {
        T x(std::move(ptr[j]));
        Integer i;

        for (i = j - 1; ((i >= 0) && (ptr[i] > x)); --i)
            ptr[i + 1] = std::move(ptr[i]);

        ptr[i + 1] = std::move(x);
    };
};

template <class T>
void sort_ins_rev(T * const ptr, const Integer n)
{
    for(Integer j = 1; j < n; ++j)
    {
        T x(std::move(ptr[j]));
        Integer i;

        for (i = j - 1; ((i >= 0) && (ptr[i] < x)); --i)
            ptr[i + 1] = std::move(ptr[i]);

        ptr[i + 1] = std::move(x);
    };
};

template <class T1, class T2>
void sort_ins(T1 * const ptr1, T2 * const ptr2,const Integer n)
{
    for(Integer j = 1; j < n; ++j)
    {
        T1 x1(std::move(ptr1[j]));
        T2 x2(std::move(ptr2[j]));
        Integer i;

        for (i = j - 1; ((i >= 0) && (ptr1[i] > x1)); --i)
        {
            ptr1[i + 1] = std::move(ptr1[i]);
            ptr2[i + 1] = std::move(ptr2[i]);
        }

        ptr1[i + 1] = std::move(x1);
        ptr2[i + 1] = std::move(x2);
    };
};

template <class T1, class T2>
void sort_ins_rev(T1 * const ptr1, T2 * const ptr2, const Integer n)
{
    for(Integer j = 1; j < n; ++j)
    {
        T1 x1(std::move(ptr1[j]));
        T2 x2(std::move(ptr2[j]));

        Integer i;
        for (i = j - 1; ((i >= 0) && (ptr1[i] < x1)); --i)
        {
            ptr1[i + 1] = std::move(ptr1[i]);
            ptr2[i + 1] = std::move(ptr2[i]);
        };

        ptr1[i + 1] = std::move(x1);
        ptr2[i + 1] = std::move(x2);
    };
};

template <class T1, class T2, class T3>
void sort_ins(T1 * const ptr1, T2 * const ptr2, T3 * const ptr3,const Integer n)
{
    for(Integer j = 1; j < n; ++j)
    {
        T1 x1(std::move(ptr1[j]));
        T2 x2(std::move(ptr2[j]));
        T3 x3(std::move(ptr3[j]));
        Integer i;

        for (i = j - 1; ((i >= 0) && (ptr1[i] > x1)); --i)
        {
            ptr1[i + 1] = std::move(ptr1[i]);
            ptr2[i + 1] = std::move(ptr2[i]);
            ptr3[i + 1] = std::move(ptr3[i]);
        }

        ptr1[i + 1] = std::move(x1);
        ptr2[i + 1] = std::move(x2);
        ptr3[i + 1] = std::move(x3);
    };
};

template <class T1, class T2, class T3>
void sort_ins_rev(T1 * const ptr1, T2 * const ptr2, T3 * const ptr3,const Integer n)
{
    for(Integer j = 1; j < n; ++j)
    {
        T1 x1(std::move(ptr1[j]));
        T2 x2(std::move(ptr2[j]));
        T3 x3(std::move(ptr3[j]));
        Integer i;

        for (i = j - 1; ((i >= 0) && (ptr1[i] < x1)); --i)
        {
            ptr1[i + 1] = std::move(ptr1[i]);
            ptr2[i + 1] = std::move(ptr2[i]);
            ptr3[i + 1] = std::move(ptr3[i]);
        }

        ptr1[i + 1] = std::move(x1);
        ptr2[i + 1] = std::move(x2);
        ptr3[i + 1] = std::move(x3);
    };
};

};};