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

#include "matcl-core/profile/time_stats.h"

#include <algorithm>

namespace matcl
{

inline
time_stats::time_stats()
{
    clear();
}

inline
void time_stats::clear()
{
    m_min   = std::numeric_limits<double>::infinity();
    m_max   = 0.0;
    m_sum   = 0.0;
    m_sum2  = 0.0;
    m_cases = 0.0;
};

inline
void time_stats::report(double t)
{
    m_min   = std::min(m_min, t);
    m_max   = std::max(m_max, t);
    m_sum   += t;
    m_sum2  += t*t;
    m_cases += 1.0;
}

inline
double time_stats::mean() const
{
    return m_sum / m_cases;
}

inline
double time_stats::std() const
{
    return std::sqrt(this->variance());
}

inline
double time_stats::variance() const
{
    return m_sum2 / m_cases - std::pow(m_sum / m_cases, 2);
}

inline
double time_stats::minimum() const
{
    return m_min;
}

inline
double time_stats::maximum() const
{
    return m_max;
}

};