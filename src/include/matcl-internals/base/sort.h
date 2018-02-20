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

#include "matcl-core/matrix/scalar_types.h"
#include "matcl-matrep/general/config.h"

namespace matcl { namespace utils
{

template<class T> 
void sort_ins(T* ptr, Integer size);

template<class T> 
void sort_ins_rev(T* ptr, Integer size);

template<class T1,class T2> 
void sort_ins(T1* ptr1, T2* ptr2, Integer size);

template<class T1,class T2> 
void sort_ins_rev(T1* ptr1, T2* ptr2, Integer size);

template<class T1,class T2,class T3> 
void sort_ins(T1* ptr1, T2* ptr2, T3* ptr3, Integer size);

template<class T1,class T2,class T3> 
void sort_ins_rev(T1* ptr1, T2* ptr2, T3* ptr3, Integer size);

template <class T> 
void MATCL_MATREP_EXPORT sort_q(T* ptr, Integer size);

template<class T> 
void MATCL_MATREP_EXPORT sort_q_rev(T* ptr, Integer size);

template<class T1, class T2> 
void MATCL_MATREP_EXPORT sort_q(T1* ptr1, T2* ptr2, Integer size);

template<class T1, class T2>
void MATCL_MATREP_EXPORT sort_q_rev(T1* ptr1, T2* ptr2, Integer size);

template<class T1, class T2, class T3> 
void MATCL_MATREP_EXPORT sort_q(T1* ptr1, T2* ptr2, T3* ptr3, Integer size);

};};

#include "sort.inl"