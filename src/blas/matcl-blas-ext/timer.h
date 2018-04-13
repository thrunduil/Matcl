/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2011-2015
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

#include <string>

#ifdef __unix__
    #include <sys/timeb.h>
#endif

namespace matcl { namespace blas
{

namespace details
{

struct timer_base
{
    bool            tic_started;
    int64_t         tic_int64;

    #ifdef __unix__
      struct timeb  tic_timeb;
    #endif

    void            init();
    void            tic();
    double          toc();
};

};

/// calculate elapsed time
class timer : private details::timer_base
{
    public:
        /// create timer object; in order to start counting time call tic()
        timer()        { init();};

        /// restart timer
        using details::timer_base::tic;
        
        /// finish timer; return elapsed time is seconds from last call to tic()
        using details::timer_base::toc;
};

/// do nothing
class null_timer
{
    public:
        /// create timer object; in order to start counting time call tic()
        null_timer()            {};

        /// restart timer
        void            tic()   {};    
        
        /// finish timer; return elapsed time is seconds from last call to tic()
        double          toc()   { return 0.0; };
};

}};
