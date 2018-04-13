/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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
#include "crash_behaviour.h"

#ifdef __MINGW32__
#define __STDC_LIMIT_MACROS
#endif

#if !defined (_WIN32) || defined (__MINGW32__)

namespace matcl {

void disable_crash_dialog()
{
    ;
}

void report_crash_to_stderr()
{
    ;//llvm::sys::PrintStackTraceOnErrorSignal();
}

}

#else // #ifdef _WIN32

#include <windows.h>
#include <crtdbg.h>
#include <stdlib.h>

namespace matcl {


void disable_crash_dialog()
{
    SetErrorMode(SetErrorMode(0)
        | SEM_FAILCRITICALERRORS
        // | SEM_NOALIGNMENTFAULTEXCEPT
        | SEM_NOGPFAULTERRORBOX
        | SEM_NOOPENFILEERRORBOX
    );

    _set_abort_behavior(0, 0
        | _WRITE_ABORT_MSG
        | _CALL_REPORTFAULT
    );
}

void report_crash_to_stderr()
{
    _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDERR);
    _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDERR);
    _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);

}

}

#endif // #ifdef _WIN32
