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

#include "matcl-core/profile/time_stats.h"
#include "matcl-core/memory/global_objects.h"

#include <iomanip>

namespace matcl
{

MATCL_CORE_EXPORT
std::ostream& matcl::operator<<(std::ostream& os, const time_stats& stats)
{
    os << std::setprecision(3) << std::scientific;

    os <<  "min: " << stats.minimum() << ", "
        << "mean: "<< stats.mean() << ", "
        << "std: " << stats.std() << ", "
        << "max: " << stats.maximum();

    return os;
};

MATCL_CORE_EXPORT
time_stats matcl::operator/(const time_stats& r, const time_stats& s)
{
    // second order approximation of moments of r/s
    // taken from http://www.stat.cmu.edu/~hseltman/files/ratio.pdf

    time_stats res;

    res.m_min       = r.minimum() / s.maximum();
    res.m_max       = r.maximum() / s.minimum();

    double r_mean   = r.mean();
    double s_mean   = s.mean();
    double r_var    = r.variance();
    double s_var    = s.variance();
    double r_mean2  = r_mean * r_mean;
    double s_mean2  = s_mean * s_mean;
    res.m_cases     = (r.m_cases + s.m_cases) / 2.0;

    double res_mean = r_mean / s_mean + s_var * r_mean / std::pow(s_mean, 3);
    double res_var  = r_mean2 / s_mean2 * (r_var / r_mean2 + s_var/s_mean2);

    res.m_sum       = res_mean * res.m_cases;
    res.m_sum2      = (res_var + std::pow(res_mean, 2)) * res.m_cases;

    return res;
};

MATCL_CORE_EXPORT
time_stats matcl::operator-(const time_stats& r, const time_stats& s)
{
    time_stats res;

    res.m_min       = r.minimum() - s.maximum();
    res.m_max       = r.maximum() - s.minimum();

    double r_mean   = r.mean();
    double s_mean   = s.mean();
    double r_var    = r.variance();
    double s_var    = s.variance();
    res.m_cases     = (r.m_cases + s.m_cases) / 2.0;

    double res_mean = r_mean - s_mean;
    double res_var  = r_var + s_var;

    res.m_sum       = res_mean * res.m_cases;
    res.m_sum2      = (res_var + std::pow(res_mean, 2)) * res.m_cases;

    return res;
};

};