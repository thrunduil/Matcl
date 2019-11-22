/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "extern/clapack/blaswrap.h"
#include "extern/clapack/f2c.h"
#include "extern/clapack/clapack.h"

// include the file with the interface such plugin needs to implement
#include "matcl-blas-lapack/blas_loader/blas_plugin.h"

#define CALL_SYNTAX(x) f2c::f2c_##x
#define CALL_SYNTAX2(x) f2c_##x

//missing functions
namespace f2c
{
    static s_type_wr 
    CALL_SYNTAX2(sdsdot)(i_type_wr *n, s_type_wr *sb, s_type_wr *sx, i_type_wr *incx, 
                         s_type_wr *sy, i_type_wr *incy)
    {
        return f2c::sdsdot_(n, sb, sx, incx, sy, incy);
    };

    static d_type_wr 
    CALL_SYNTAX2(dsdot)(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
    {
        return f2c::dsdot_(n, sx, incx, sy, incy);
    };
    
    static s_type_wr CALL_SYNTAX2(scabs1)(c_type_wr *z__)
    {
        return f2c::scabs1_((f2c::complex*)z__);
    };

    static s_type_wr CALL_SYNTAX2(sasum)(i_type_wr *n, s_type_wr *sx, i_type_wr *incx)
    {
        return f2c::sasum_(n, sx, incx);
    };

    static s_type_wr CALL_SYNTAX2(scasum)(i_type_wr *n, c_type_wr *cx, i_type_wr *incx)
    {
        return f2c::scasum_(n, (f2c::complex*)cx, incx);
    };

    static s_type_wr CALL_SYNTAX2(scnrm2)(i_type_wr *n, c_type_wr *cx, i_type_wr *incx)
    {
        return f2c::scnrm2_(n, (f2c::complex*)cx, incx);
    };

    static s_type_wr 
    CALL_SYNTAX2(sdot)(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
    {
        return f2c::sdot_(n, sx, incx, sy, incy);
    };

    static s_type_wr CALL_SYNTAX2(snrm2)(i_type_wr *n, s_type_wr *sx, i_type_wr *incx)
    {
        return f2c::snrm2_(n, sx, incx);
    };

};

static const ::blas_plugin g_clapack_plugin;

extern "C"
{
    BLAS_PLUGIN_EXPORT 
    const ::blas_plugin* get_blas_plugin()
    {
        return &g_clapack_plugin;
    }
}

i_type_wr get_num_threads()
{
    return 1;
};

i_type_wr get_default_threads()
{
    return 1;
};

void set_num_threads(i_type_wr*)
{
    return;
};	

bool are_user_threads_allowed()
{
    return true;
};	

const char* get_name()
{
    return "CLAPACK";
};

void initialize()
{};

void force_initialization()
{};

// generic stuff - the constructor which initializes all pointers to blas functions
#include "matcl-blas-lapack/blas_loader/blas_plugin_common.h"
