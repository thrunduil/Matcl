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

#include "mkl_blas.h"

#include <atomic>
#include <thread>
#include <algorithm>
#include <iostream>

// include the file with the interface such plugin needs to implement
#include "matcl-blas-lapack/blas_loader/blas_plugin.h"
#include "omp.h"
#include "mkl.h"

#define CALL_SYNTAX(x) x

#define FUNCTION_NAME_cdotc wrap_cdotc
#define FUNCTION_NAME_cdotu wrap_cdotu
#define FUNCTION_NAME_zdotc wrap_zdotc
#define FUNCTION_NAME_zdotu wrap_zdotu

// WARNING: OMP spport must be enabled

class global_data
{
    private:
        std::atomic<long>   m_threads;

    public:
        ::blas_plugin       m_plugin;

    public:
        global_data();

        int     get_num_threads();
        void    set_num_threads(int n);
        bool    are_user_threads_allowed();
};

global_data::global_data()
    :m_threads(0)
{
    int max_threads = std::thread::hardware_concurrency();
    m_threads       = max_threads;

    omp_set_num_threads(max_threads);
};

int global_data::get_num_threads()
{
    return m_threads;
};

void global_data::set_num_threads(int n)
{
    int max_threads = std::thread::hardware_concurrency();
    m_threads = std::min(std::max(n,1),max_threads);
    omp_set_num_threads(m_threads);
};

bool global_data::are_user_threads_allowed()
{
    return true;
};

static global_data m_data;

extern "C"
{
    BLAS_PLUGIN_EXPORT 
    const ::blas_plugin* get_blas_plugin()
    {
        return &m_data.m_plugin;
    }
}
	
i_type_wr get_num_threads()
{
    return m_data.get_num_threads();
};

i_type_wr get_default_threads()
{
    return std::thread::hardware_concurrency();
};

void set_num_threads(i_type_wr* n)
{
    m_data.set_num_threads(*n);
};	

bool are_user_threads_allowed()
{
    return m_data.are_user_threads_allowed();
};

void initialize()
{};

void force_initialization()
{};

const char* get_name()
{
    return "MKL";
};

extern "C"
{
    static c_type_wr wrap_cdotc(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, 
                     i_type_wr *incy)
    {
        c_type_wr res;

        cdotc((MKL_Complex8*)&res, n, (MKL_Complex8*)cx, incx, (MKL_Complex8*)cy, incy);
        return res;        
    };

    static c_type_wr wrap_cdotu(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, 
                 i_type_wr *incy)
    {
        c_type_wr res;

        cdotu((MKL_Complex8*)&res, n, (MKL_Complex8*)cx, incx, (MKL_Complex8*)cy, incy);
        return res;        
    };

    static z_type_wr wrap_zdotc(i_type_wr *n, z_type_wr *cx, i_type_wr *incx, z_type_wr *cy, 
                     i_type_wr *incy)
    {
        z_type_wr res;

        zdotc((MKL_Complex16*)&res, n, (MKL_Complex16*)cx, incx, (MKL_Complex16*)cy, incy);
        return res;        
    };

    static z_type_wr wrap_zdotu(i_type_wr *n, z_type_wr *cx, i_type_wr *incx, z_type_wr *cy, 
                     i_type_wr *incy)
    {
        z_type_wr res;

        zdotu((MKL_Complex16*)&res, n, (MKL_Complex16*)cx, incx, (MKL_Complex16*)cy, incy);
        return res;        
    };

};

// generic stuff - the constructor which initializes all pointers to blas functions
#include "matcl-blas-lapack/blas_loader/blas_plugin_common.h"
