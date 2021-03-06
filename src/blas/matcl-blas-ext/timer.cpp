/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe� Kowal 2011
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
#include "timer.h"

#include <iostream>
#include <sstream>

#ifndef __unix__
    #include <windows.h>
#endif

//FIXME: implement *nix ver
#ifndef _MSC_VER

    #include <cstdarg>
    namespace
    {
        int sprintf_s(char * _DstBuf, size_t _SizeInBytes, const char * _Format, ...)
        {
            va_list args;
            va_start (args, _Format);
            int n = vsnprintf(_DstBuf, _SizeInBytes, _Format, args);
            va_end (args);

            return n;
        }
    }
#endif

#ifdef __unix__
    static inline matcl::double MicrosecsToSecs( int64_t microseconds )
    {
        return 1e-06*microseconds;
    }
#endif

namespace matcl { namespace blas
{

void details::timer_base::init()
{
    tic_started = false;
};

void details::timer_base::tic(void)
{
    tic_started = true;

    #ifndef __unix__
        QueryPerformanceCounter((LARGE_INTEGER*) &tic_int64);
    #else
        ftime(&tic_timeb);
    #endif
}

double details::timer_base::toc(void)
{
    #ifndef __unix__
        double t;
        int64_t toc_int64, fr_int64;

        if (tic_started)
        {
            tic_started = false;
            QueryPerformanceCounter((LARGE_INTEGER*) &toc_int64);
            QueryPerformanceFrequency((LARGE_INTEGER*) &fr_int64);
            t = (double) (toc_int64 - tic_int64) / (double) fr_int64;
            return t;
        }
        else
        {
            return 0.;
        };
    #else
        struct timeb toc_timeb;
        ftime(&toc_timeb);

        long dSec    = (long) toc_timeb.time - (long)tic_timeb.time;
        long dUSec   = 1000 * toc_timeb.millitm - 1000 * tic_timeb.millitm;

        if (dUSec < 0)
        {
            dUSec += 1000000;
            --dSec;
        }

        return MicrosecsToSecs( 1000000 * dSec + dUSec );
    #endif
}

}};

