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
#include "matcl-core/profile/timer.h"

#include <memory>

#pragma warning(push)
#pragma warning(disable : 4251) // needs to have dll-interface to be used by clients

namespace matcl
{

// base class for functions, that can be tested by the benchmark class
struct benchmark_function
{
    virtual ~benchmark_function(){};

    // derived class must implement this function, which executes code to be 
    // profiled; evaluation of this function must require at least a few
    // microsecond (resolution of timer is ~ 1 microsecond), otherwise 
    // results will be measingless
    virtual void    eval() = 0;
};

// test execution time of a function
class MATCL_CORE_EXPORT benchmark
{
    public:
        // pointer to an object implementing 
        using function_ptr  = std::shared_ptr<benchmark_function>;

    private:
        timer           m_time;
        function_ptr    m_func;
        time_stats      m_stats;

        static volatile 
        size_t          m_dummy;

    public:
        // create benchmark object for a function f
        benchmark(const function_ptr& f);

        // run tests, function will be executed n_rep times
        time_stats  make(int n_rep);

        // function implementing a test should call this function in order
        // to suppress optimizations by the compiler (i.e. unsued value 
        // optimizations)
        template<class T>
        static void use_value(const T& v);

    private:
        void        report_time(double t);
};

};

#pragma warning(pop)

#include "matcl-core/details/profile/benchmark.inl"