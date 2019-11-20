/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*! @file
 * \brief Returns the time in seconds used by the process
 *
 * <pre>
 * Purpose
 * ======= 
 *	Returns the time in seconds used by the process.
 *
 * Note: the timer function call is machine dependent. Use conditional
 *       compilation to choose the appropriate function.
 * </pre>
 */

#include "superlu_defs.h"

#ifdef SUN 
/*
 * 	It uses the system call gethrtime(3C), which is accurate to 
 *	nanoseconds. 
*/
#include <sys/time.h>
 
double SuperLU_timer_() {
    return ( (double)gethrtime() / 1e9 );
}

#elif defined ( UNIX_TIMER )

#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <sys/time.h>

double SuperLU_timer_()
{
    struct tms use;
    double tmp;
    int clocks_per_sec = sysconf(_SC_CLK_TCK);

    times(&use);
    tmp = use.tms_utime;
    tmp += use.tms_stime;
    return (double)(tmp) / clocks_per_sec;
}

#elif _WIN32

#include <time.h>

double SuperLU_timer_()
{
    clock_t t;
    t=clock();

    return ((double)t)/CLOCKS_PER_SEC;
}

#else

#include "mpi/mpiuni.h"

double SuperLU_timer_()
{
    return MPI_Wtime();
}

#endif

