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

#include "matcl-core/config.h"
#include <iostream>

namespace matcl
{

// calculate statistics of a set of values of elapsed time measured
// by the benchmark function
class time_stats
{
    private:
        double          m_min;
        double          m_max;
        double          m_sum;
        double          m_sum2;
        double          m_cases;

    public:
        // initialize with empty set of measured time
        time_stats();

        // clear the set
        void    clear();

        // add new time measurement to the set and update statistics
        void    report(double t);

        // calculate mean of the set
        double  mean() const;

        // calculate standard deviation of the set
        double  std() const;

        // calculate variance of the set
        double  variance() const;

        // return minimum of elements in the se
        double  minimum() const;

        // return maximum of elements in the se
        double  maximum() const;

        // print computed statistics
        MATCL_CORE_EXPORT
        friend  std::ostream& operator<<(std::ostream& os, const time_stats& stats);

        // approximate statistics of a distribution, which is the ration of a
        // distribution represented by the first time_stats r and the second 
        // time_stats s; resulting time_stats should not be appended with new
        // observations unless new observations represent ratios
        MATCL_CORE_EXPORT
        friend  time_stats operator/(const time_stats& r, const time_stats& s);

        // approximate statistics of a distribution, which is the difference of a
        // distribution represented by the first time_stats r and the second 
        // time_stats s; resulting time_stats should not be appended with new
        // observations unless new observations represent differences
        MATCL_CORE_EXPORT
        friend  time_stats operator-(const time_stats& r, const time_stats& s);
};

};

#include "matcl-core/details/profile/time_stats.inl"