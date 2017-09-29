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

#include "matcl-core/config.h"
#include "matcl-mp/config.h"

namespace matcl { namespace error_flags
{
    //-----------------------------------------------------------------
    // functions to manipulate floating point error flags set by functions
    // in mfpr and matcl-mp library; these flags are thread local
    //-----------------------------------------------------------------

    // clear the underflow, overflow, divide-by-zero, invalid, inexact
    // and erange flags.
    void    clear_underflow();
    void    clear_overflow();
    void    clear_divby0();
    void    clear_nan();
    void    clear_inexact();
    void    clear_erange();

    // set the underflow, overflow, divide-by-zero, invalid, inexact and
    // erange flags
    void    set_underflow();
    void    set_overflow();
    void    set_divby0();
    void    set_nan();
    void    set_inexact();
    void    set_erange();

    // return the corresponding (underflow, overflow, divide-by-zero, 
    // invalid, inexact, erange) flag, which is non-zero iff the flag is set.
    bool    get_underflow();
    bool    get_overflow();
    bool    get_divby0();
    bool    get_nan();
    bool    get_inexact();
    bool    get_erange();

    //-----------------------------------------------------------------
    // functions to manipulate integer error flags set by functions
    // in matcl-mp library; these flags are thread local
    //-----------------------------------------------------------------
    void    clear_integer_overflow();
    void    clear_integer_erange();

    void    set_integer_overflow();
    void    set_integer_erange();

    bool    get_integer_overflow();
    bool    get_integer_erange();

    // clear all floating point and integer flags
    void    clear_flags();

};}
