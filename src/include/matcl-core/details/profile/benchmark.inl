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

#include "matcl-core/profile/benchmark.h"

namespace matcl
{

inline
benchmark::benchmark(const function_ptr& f)
    :m_func(f)
{};

template<class T>
inline 
static void benchmark::use_value(const T& v)
{
    size_t addr = *reinterpret_cast<const size_t*>(&v);
    m_dummy     += addr;
};

inline 
void benchmark::report_time(double t)
{
    m_stats.report(t);
};

};