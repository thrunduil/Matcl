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

#include "matcl-mp/error_flags.h"
#include "utils/impl_types.h"
#include "matcl-core/general/thread.h"

namespace matcl
{

void error_flags::clear_underflow()
{
    mpfr_clear_underflow();
};
void error_flags::clear_overflow()
{
    mpfr_clear_overflow();
};
void error_flags::clear_divby0()
{
    mpfr_clear_divby0();
};
void error_flags::clear_nan()
{
    mpfr_clear_nanflag();
};
void error_flags::clear_inexact()
{
    mpfr_clear_inexflag();
};
void error_flags::clear_erange()
{
    mpfr_clear_erangeflag();
};

void error_flags::set_underflow()
{
    mpfr_set_underflow();
};
void error_flags::set_overflow()
{
    mpfr_set_overflow();
};
void error_flags::set_divby0()
{
    mpfr_set_divby0();
};
void error_flags::set_nan()
{
    mpfr_set_nanflag();
};
void error_flags::set_inexact()
{
    mpfr_set_inexflag();
};
void error_flags::set_erange()
{
    mpfr_set_erangeflag();
};

bool error_flags::get_underflow()
{
    return mpfr_underflow_p() != 0;
};
bool error_flags::get_overflow()
{
    return mpfr_overflow_p() != 0;
};
bool error_flags::get_divby0()
{
    return mpfr_divby0_p() != 0;
};
bool error_flags::get_nan()
{
    return mpfr_nanflag_p() != 0;
};
bool error_flags::get_inexact()
{
    return mpfr_inexflag_p() != 0;
};
bool error_flags::get_erange()
{
    return mpfr_erangeflag_p() != 0;
};

MATCL_THREAD_LOCAL static int g_integer_overflow = 0;
MATCL_THREAD_LOCAL static int g_integer_erange = 0;

void error_flags::clear_integer_overflow()
{
    g_integer_overflow = 0;
};
void error_flags::clear_integer_erange()
{
    g_integer_erange = 0;
};

void error_flags::set_integer_overflow()
{
    g_integer_overflow = 1;
};
void error_flags::set_integer_erange()
{
    g_integer_erange = 1;
};

bool error_flags::get_integer_overflow()
{
    return g_integer_overflow != 0;
};
bool error_flags::get_integer_erange()
{
    return g_integer_erange != 0;
}

void error_flags::clear_flags()
{
    mpfr_clear_flags();
    clear_integer_overflow();
    clear_integer_erange();
};

};