/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2011
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

#include <atomic>
#include <thread>
#include <algorithm>
#include <iostream>

// include the file with the interface such plugin needs to implement
#include "blas_wrapper.h"
#include <Windows.h>
#undef min
#undef max

class global_data
{
    public:
        ::blas_wrapper          m_wrapper;
        const ::blas_wrapper*   m_plugin_cpu;
        const ::blas_wrapper*   m_plugin_gpu;

    public:
        global_data();

        const ::blas_wrapper*
                    get_wrapper(const ::blas_wrapper* plugin_cpu, const ::blas_wrapper* plugin_gpu);
        void        init();
        int         get_num_threads();
        void        set_num_threads(int n);
        bool        are_user_threads_allowed();
};

global_data::global_data()
    :m_plugin_cpu(nullptr), m_plugin_gpu(nullptr)
{};
void global_data::init()
{}
int global_data::get_num_threads()
{
    return m_plugin_cpu->get_num_threads_fptr();
};
void global_data::set_num_threads(int n)
{
    m_plugin_cpu->set_num_threads_fptr(&n);
};
bool global_data::are_user_threads_allowed()
{
    return m_plugin_cpu->are_user_threads_allowed_fptr();
};

const ::blas_wrapper*  global_data::get_wrapper(const ::blas_wrapper* plugin_cpu, const ::blas_wrapper* plugin_gpu)
{
    this->m_plugin_cpu = plugin_cpu;
    this->m_plugin_gpu = plugin_gpu;
    this->init();

    return &this->m_wrapper;
};

static global_data m_data;

BOOL APIENTRY DllMain(HMODULE hModule,
                      DWORD  ul_reason_for_call,
                      LPVOID lpReserved
                      )
{
    switch (ul_reason_for_call)
    {
        case DLL_PROCESS_ATTACH:
        case DLL_THREAD_ATTACH:
            break;
        case DLL_THREAD_DETACH:
            break;
        case DLL_PROCESS_DETACH:
            break;
    }
    return TRUE;
}

extern "C"
{
    BLAS_PLUGIN_EXPORT 
    const ::blas_wrapper* get_blas_wrapper_gpu_cpu(const ::blas_wrapper* plugin_cpu, const ::blas_wrapper* plugin_gpu)
    {
        return m_data.get_wrapper(plugin_cpu, plugin_gpu);
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

//--------------------------------------------------------------------------
//                      BLAS I, II
//--------------------------------------------------------------------------
i_type_wr caxpy_gpu_wrap_impl(i_type_wr *n, c_type_wr *alpha, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->caxpy_fptr(n, alpha, cx, incx, cy, incy);
};
i_type_wr daxpy_gpu_wrap_impl(i_type_wr *n, d_type_wr *da, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->daxpy_fptr(n, da, dx, incx, dy, incy);
};
i_type_wr saxpy_gpu_wrap_impl(i_type_wr *n, s_type_wr *sa, s_type_wr *sx, i_type_wr *incx,
                              s_type_wr *sy, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->saxpy_fptr(n, sa, sx, incx, sy, incy);
};
i_type_wr zaxpy_gpu_wrap_impl(i_type_wr *n, z_type_wr *za, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->zaxpy_fptr(n, za, zx, incx, zy, incy);
};

i_type_wr ccopy_gpu_wrap_impl(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, i_type_wr *incy)
{
    return  m_data.m_plugin_cpu->ccopy_fptr(n, cx, incx, cy, incy);
};
i_type_wr cdotc_gpu_wrap_impl(c_type_wr *ret_val, i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->cdotc_fptr(ret_val, n, cx, incx, cy, incy);
};
i_type_wr cdotu_gpu_wrap_impl(c_type_wr *ret_val, i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->cdotu_fptr(ret_val, n, cx, incx, cy, incy);
};
i_type_wr cgbmv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku,
                              c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx,
                              c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->cgbmv_fptr(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr cgemv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *a,
                              i_type_wr *lda, c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->cgemv_fptr(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr cgerc_gpu_wrap_impl(i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                              c_type_wr *y, i_type_wr *incy, c_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->cgerc_fptr(m, n, alpha, x, incx, y, incy, a, lda);
};
i_type_wr cgeru_gpu_wrap_impl(i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                              c_type_wr *y, i_type_wr *incy, c_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->cgeru_fptr(m, n, alpha, x, incx, y, incy, a, lda);
};
i_type_wr chbmv_gpu_wrap_impl(char *uplo, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, c_type_wr *a,
                              i_type_wr *lda, c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->chbmv_fptr(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr chemv_gpu_wrap_impl(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda,
                              c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->chemv_fptr(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr cher_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                             c_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->cher_fptr(uplo, n, alpha, x, incx, a, lda);
};
i_type_wr cher2_gpu_wrap_impl(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                              c_type_wr *y, i_type_wr *incy, c_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->cher2_fptr(uplo, n, alpha, x, incx, y, incy, a, lda);
};
i_type_wr chpmv_gpu_wrap_impl(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *ap, c_type_wr *x,
                              i_type_wr *incx, c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->chpmv_fptr(uplo, n, alpha, ap, x, incx, beta, y, incy);
};
i_type_wr chpr_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                             c_type_wr *ap)
{
    return m_data.m_plugin_cpu->chpr_fptr(uplo, n, alpha, x, incx, ap);
};
i_type_wr chpr2_gpu_wrap_impl(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                              c_type_wr *y, i_type_wr *incy, c_type_wr *ap)
{
    return m_data.m_plugin_cpu->chpr2_fptr(uplo, n, alpha, x, incx, y, incy, ap);
};
i_type_wr crotg_gpu_wrap_impl(c_type_wr *ca, c_type_wr *cb, s_type_wr *c__, c_type_wr *s)
{
    return m_data.m_plugin_cpu->crotg_fptr(ca, cb, c__, s);
};
i_type_wr cscal_gpu_wrap_impl(i_type_wr *n, c_type_wr *ca, c_type_wr *cx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->cscal_fptr(n, ca, cx, incx);
};
i_type_wr csrot_gpu_wrap_impl(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, i_type_wr *incy,
                              s_type_wr *c__, s_type_wr *s)
{
    return m_data.m_plugin_cpu->csrot_fptr(n, cx, incx, cy, incy, c__, s);
};
i_type_wr csscal_gpu_wrap_impl(i_type_wr *n, s_type_wr *sa, c_type_wr *cx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->csscal_fptr(n, sa, cx, incx);
};
i_type_wr cswap_gpu_wrap_impl(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->cswap_fptr(n, cx, incx, cy, incy);
};
i_type_wr ctbmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, c_type_wr *a,
                              i_type_wr *lda, c_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->ctbmv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
};
i_type_wr ctbsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, c_type_wr *a,
                              i_type_wr *lda, c_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->ctbsv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
};
i_type_wr ctpmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *ap, c_type_wr *x,
                              i_type_wr *incx)
{
    return m_data.m_plugin_cpu->ctpmv_fptr(uplo, trans, diag, n, ap, x, incx);
};
i_type_wr ctpsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *ap, c_type_wr *x,
                              i_type_wr *incx)
{
    return m_data.m_plugin_cpu->ctpsv_fptr(uplo, trans, diag, n, ap, x, incx);
};
i_type_wr ctrmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *a, i_type_wr *lda,
                              c_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->ctrmv_fptr(uplo, trans, diag, n, a, lda, x, incx);
};
i_type_wr ctrsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *a, i_type_wr *lda,
                              c_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->ctrsv_fptr(uplo, trans, diag, n, a, lda, x, incx);
};
d_type_wr dasum_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->dasum_fptr(n, dx, incx);
};
d_type_wr dcabs1_gpu_wrap_impl(z_type_wr *z__)
{
    return m_data.m_plugin_cpu->dcabs1_fptr(z__);
};
i_type_wr dcopy_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->dcopy_fptr(n, dx, incx, dy, incy);
};
d_type_wr ddot_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->ddot_fptr(n, dx, incx, dy, incy);
};
i_type_wr dgbmv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku,
                              d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *x, i_type_wr *incx,
                              d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->dgbmv_fptr(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr dgemv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *a,
                              i_type_wr *lda, d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->dgemv_fptr(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr dger_gpu_wrap_impl(i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx,
                             d_type_wr *y, i_type_wr *incy, d_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->dger_fptr(m, n, alpha, x, incx, y, incy, a, lda);
};
d_type_wr dnrm2_gpu_wrap_impl(i_type_wr *n, d_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->dnrm2_fptr(n, x, incx);
};
i_type_wr drot_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy,
                             d_type_wr *c__, d_type_wr *s)
{
    return m_data.m_plugin_cpu->drot_fptr(n, dx, incx, dy, incy, c__, s);
};
i_type_wr drotg_gpu_wrap_impl(d_type_wr *da, d_type_wr *db, d_type_wr *c__, d_type_wr *s)
{
    return m_data.m_plugin_cpu->drotg_fptr(da, db, c__, s);
};
i_type_wr drotm_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy,
                              d_type_wr *dparam)
{
    return m_data.m_plugin_cpu->drotm_fptr(n, dx, incx, dy, incy, dparam);
};
i_type_wr drotmg_gpu_wrap_impl(d_type_wr *dd1, d_type_wr *dd2, d_type_wr *dx1, d_type_wr *dy1, d_type_wr *dparam)
{
    return m_data.m_plugin_cpu->drotmg_fptr(dd1, dd2, dx1, dy1, dparam);
};
i_type_wr dsbmv_gpu_wrap_impl(char *uplo, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, d_type_wr *a,
                              i_type_wr *lda, d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->dsbmv_fptr(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr dscal_gpu_wrap_impl(i_type_wr *n, d_type_wr *da, d_type_wr *dx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->dscal_fptr(n, da, dx, incx);
};
/*
d_type_wr dsdot_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->dsdot_fptr(n, sx, incx, sy, incy);
};
*/
i_type_wr dspmv_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *ap, d_type_wr *x,
                              i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->dspmv_fptr(uplo, n, alpha, ap, x, incx, beta, y, incy);
};
i_type_wr dspr_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx,
                             d_type_wr *ap)
{
    return m_data.m_plugin_cpu->dspr_fptr(uplo, n, alpha, x, incx, ap);
};
i_type_wr dspr2_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx,
                              d_type_wr *y, i_type_wr *incy, d_type_wr *ap)
{
    return m_data.m_plugin_cpu->dspr2_fptr(uplo, n, alpha, x, incx, y, incy, ap);
};
i_type_wr dswap_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->dswap_fptr(n, dx, incx, dy, incy);
};
i_type_wr dsymv_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda,
                              d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->dsymv_fptr(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr dsyr_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx,
                             d_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->dsyr_fptr(uplo, n, alpha, x, incx, a, lda);
};
i_type_wr dsyr2_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx,
                              d_type_wr *y, i_type_wr *incy, d_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->dsyr2_fptr(uplo, n, alpha, x, incx, y, incy, a, lda);
};
i_type_wr dtbmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, d_type_wr *a,
                              i_type_wr *lda, d_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->dtbmv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
};
i_type_wr dtbsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, d_type_wr *a,
                              i_type_wr *lda, d_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->dtbsv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
};
i_type_wr dtpmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *ap,
                              d_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->dtpmv_fptr(uplo, trans, diag, n, ap, x, incx);
};
i_type_wr dtpsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *ap,
                              d_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->dtpsv_fptr(uplo, trans, diag, n, ap, x, incx);
};
i_type_wr dtrmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *a, i_type_wr *lda,
                              d_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->dtrmv_fptr(uplo, trans, diag, n, a, lda, x, incx);
};
i_type_wr dtrsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *a, i_type_wr *lda,
                              d_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->dtrsv_fptr(uplo, trans, diag, n, a, lda, x, incx);
};
d_type_wr dzasum_gpu_wrap_impl(i_type_wr *n, z_type_wr *zx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->dzasum_fptr(n, zx, incx);
};
d_type_wr dznrm2_gpu_wrap_impl(i_type_wr *n, z_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->dznrm2_fptr(n, x, incx);
};
i_type_wr icamax_gpu_wrap_impl(i_type_wr *n, c_type_wr *cx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->icamax_fptr(n, cx, incx);
};
i_type_wr idamax_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->idamax_fptr(n, dx, incx);
};
i_type_wr isamax_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->isamax_fptr(n, sx, incx);
};
i_type_wr izamax_gpu_wrap_impl(i_type_wr *n, z_type_wr *zx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->izamax_fptr(n, zx, incx);
};
/*
i_type_wr lsame_gpu_wrap_impl(char *ca, char *cb)
{
    return m_data.m_plugin_cpu->lsame_fptr(ca, cb);
};
*/
d_type_wr sasum_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->sasum_fptr(n, sx, incx);
};
/*
d_type_wr scabs1_gpu_wrap_impl(c_type_wr *z__)
{
    return m_data.m_plugin_cpu->scabs1_fptr(z__);
};
*/
d_type_wr scasum_gpu_wrap_impl(i_type_wr *n, c_type_wr *cx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->scasum_fptr(n, cx, incx);
};
d_type_wr scnrm2_gpu_wrap_impl(i_type_wr *n, c_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->scnrm2_fptr(n, x, incx);
};
i_type_wr scopy_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->scopy_fptr(n, sx, incx, sy, incy);
};
s_type_wr sdot_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->sdot_fptr(n, sx, incx, sy, incy);
};
/*
d_type_wr sdsdot_gpu_wrap_impl(i_type_wr *n, s_type_wr *sb, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy,
                               i_type_wr *incy)
{
    return m_data.m_plugin_cpu->sdsdot_fptr(n, sb, sx, incx, sy, incy);
};
*/
i_type_wr sgbmv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku,
                              s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *x, i_type_wr *incx,
                              s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->sgbmv_fptr(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr sgemv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *a,
                              i_type_wr *lda, s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->sgemv_fptr(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr sger_gpu_wrap_impl(i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx,
                             s_type_wr *y, i_type_wr *incy, s_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->sger_fptr(m, n, alpha, x, incx, y, incy, a, lda);
};
d_type_wr snrm2_gpu_wrap_impl(i_type_wr *n, s_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->snrm2_fptr(n, x, incx);
};
i_type_wr srot_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy,
                             s_type_wr *c__, s_type_wr *s)
{
    return m_data.m_plugin_cpu->srot_fptr(n, sx, incx, sy, incy, c__, s);
};
i_type_wr srotg_gpu_wrap_impl(s_type_wr *sa, s_type_wr *sb, s_type_wr *c__, s_type_wr *s)
{
    return m_data.m_plugin_cpu->srotg_fptr(sa, sb, c__, s);
};
i_type_wr srotm_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy,
                              s_type_wr *sparam)
{
    return m_data.m_plugin_cpu->srotm_fptr(n, sx, incx, sy, incy, sparam);
};
i_type_wr srotmg_gpu_wrap_impl(s_type_wr *sd1, s_type_wr *sd2, s_type_wr *sx1, s_type_wr *sy1, s_type_wr *sparam)
{
    return m_data.m_plugin_cpu->srotmg_fptr(sd1, sd2, sx1, sy1, sparam);
};
i_type_wr ssbmv_gpu_wrap_impl(char *uplo, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, s_type_wr *a,
                              i_type_wr *lda, s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->ssbmv_fptr(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr sscal_gpu_wrap_impl(i_type_wr *n, s_type_wr *sa, s_type_wr *sx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->sscal_fptr(n, sa, sx, incx);
};
i_type_wr sspmv_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *ap, s_type_wr *x,
                              i_type_wr *incx, s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->sspmv_fptr(uplo, n, alpha, ap, x, incx, beta, y, incy);
};
i_type_wr sspr_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx,
                             s_type_wr *ap)
{
    return m_data.m_plugin_cpu->sspr_fptr(uplo, n, alpha, x, incx, ap);
};
i_type_wr sspr2_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx,
                              s_type_wr *y, i_type_wr *incy, s_type_wr *ap)
{
    return m_data.m_plugin_cpu->sspr2_fptr(uplo, n, alpha, x, incx, y, incy, ap);
};
i_type_wr sswap_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->sswap_fptr(n, sx, incx, sy, incy);
};
i_type_wr ssymv_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda,
                              s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->ssymv_fptr(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr ssyr_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx,
                             s_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->ssyr_fptr(uplo, n, alpha, x, incx, a, lda);
};
i_type_wr ssyr2_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx,
                              s_type_wr *y, i_type_wr *incy, s_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->ssyr2_fptr(uplo, n, alpha, x, incx, y, incy, a, lda);
};
i_type_wr stbmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, s_type_wr *a,
                              i_type_wr *lda, s_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->stbmv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
};
i_type_wr stbsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, s_type_wr *a,
                              i_type_wr *lda, s_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->stbsv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
};
i_type_wr stpmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *ap, s_type_wr *x,
                              i_type_wr *incx)
{
    return m_data.m_plugin_cpu->stpmv_fptr(uplo, trans, diag, n, ap, x, incx);
};
i_type_wr stpsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *ap, s_type_wr *x,
                              i_type_wr *incx)
{
    return m_data.m_plugin_cpu->stpsv_fptr(uplo, trans, diag, n, ap, x, incx);
};
i_type_wr strmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *a, i_type_wr *lda,
                              s_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->strmv_fptr(uplo, trans, diag, n, a, lda, x, incx);
};
i_type_wr strsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *a, i_type_wr *lda,
                              s_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->strsv_fptr(uplo, trans, diag, n, a, lda, x, incx);
};
/*
i_type_wr xerbla_gpu_wrap_impl(char *srname, i_type_wr *info)
{
    return m_data.m_plugin_cpu->xerbla_fptr(srname, info);
};
*/
i_type_wr zcopy_gpu_wrap_impl(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->zcopy_fptr(n, zx, incx, zy, incy);
};
i_type_wr zdotc_gpu_wrap_impl(z_type_wr *res, i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->zdotc_fptr(res, n, zx, incx, zy, incy);
};
i_type_wr zdotu_gpu_wrap_impl(z_type_wr *res, i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->zdotu_fptr(res, n, zx, incx, zy, incy);
};
i_type_wr zdrot_gpu_wrap_impl(i_type_wr *n, z_type_wr *cx, i_type_wr *incx, z_type_wr *cy, i_type_wr *incy,
                              d_type_wr *c__, d_type_wr *s)
{
    return m_data.m_plugin_cpu->zdrot_fptr(n, cx, incx, cy, incy, c__, s);
};
i_type_wr zdscal_gpu_wrap_impl(i_type_wr *n, d_type_wr *da, z_type_wr *zx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->zdscal_fptr(n, da, zx, incx);
};
i_type_wr zgbmv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku,
                              z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx,
                              z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->zgbmv_fptr(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr zgemv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *a,
                              i_type_wr *lda, z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->zgemv_fptr(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr zgerc_gpu_wrap_impl(i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                              z_type_wr *y, i_type_wr *incy, z_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->zgerc_fptr(m, n, alpha, x, incx, y, incy, a, lda);
};
i_type_wr zgeru_gpu_wrap_impl(i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                              z_type_wr *y, i_type_wr *incy, z_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->zgeru_fptr(m, n, alpha, x, incx, y, incy, a, lda);
};
i_type_wr zhbmv_gpu_wrap_impl(char *uplo, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, z_type_wr *a,
                              i_type_wr *lda, z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y,
                              i_type_wr *incy)
{
    return m_data.m_plugin_cpu->zhbmv_fptr(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr zhemv_gpu_wrap_impl(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda,
                              z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->zhemv_fptr(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
};
i_type_wr zher_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                             z_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->zher_fptr(uplo, n, alpha, x, incx, a, lda);
};
i_type_wr zher2_gpu_wrap_impl(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                              z_type_wr *y, i_type_wr *incy, z_type_wr *a, i_type_wr *lda)
{
    return m_data.m_plugin_cpu->zher2_fptr(uplo, n, alpha, x, incx, y, incy, a, lda);
};
i_type_wr zhpmv_gpu_wrap_impl(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *ap, z_type_wr *x,
                              i_type_wr *incx, z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->zhpmv_fptr(uplo, n, alpha, ap, x, incx, beta, y, incy);
};
i_type_wr zhpr_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                             z_type_wr *ap)
{
    return m_data.m_plugin_cpu->zhpr_fptr(uplo, n, alpha, x, incx, ap);
};
i_type_wr zhpr2_gpu_wrap_impl(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                              z_type_wr *y, i_type_wr *incy, z_type_wr *ap)
{
    return m_data.m_plugin_cpu->zhpr2_fptr(uplo, n, alpha, x, incx, y, incy, ap);
};
i_type_wr zrotg_gpu_wrap_impl(z_type_wr *ca, z_type_wr *cb, d_type_wr *c__, z_type_wr *s)
{
    return m_data.m_plugin_cpu->zrotg_fptr(ca, cb, c__, s);
};
i_type_wr zscal_gpu_wrap_impl(i_type_wr *n, z_type_wr *za, z_type_wr *zx, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->zscal_fptr(n, za, zx, incx);
};
i_type_wr zswap_gpu_wrap_impl(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, i_type_wr *incy)
{
    return m_data.m_plugin_cpu->zswap_fptr(n, zx, incx, zy, incy);
};
i_type_wr ztbmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, z_type_wr *a,
                              i_type_wr *lda, z_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->ztbmv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
};
i_type_wr ztbsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, z_type_wr *a,
                              i_type_wr *lda, z_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->ztbsv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
};
i_type_wr ztpmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *ap, z_type_wr *x,
                              i_type_wr *incx)
{
    return m_data.m_plugin_cpu->ztpmv_fptr(uplo, trans, diag, n, ap, x, incx);
};
i_type_wr ztpsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *ap, z_type_wr *x,
                              i_type_wr *incx)
{
    return m_data.m_plugin_cpu->ztpsv_fptr(uplo, trans, diag, n, ap, x, incx);
};
i_type_wr ztrmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *a, i_type_wr *lda,
                              z_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->ztrmv_fptr(uplo, trans, diag, n, a, lda, x, incx);
};
i_type_wr ztrsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *a, i_type_wr *lda,
                              z_type_wr *x, i_type_wr *incx)
{
    return m_data.m_plugin_cpu->ztrsv_fptr(uplo, trans, diag, n, a, lda, x, incx);
};

//--------------------------------------------------------------------------
//                      BLAS III
//--------------------------------------------------------------------------
i_type_wr dgemm_gpu_wrap_impl(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                              d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb,
                              d_type_wr *beta, d_type_wr *c__, i_type_wr *ldc)
{
    if ((transa[0] == 'n' || transa[0] == 'N') && (transb[0] == 'n' || transb[0] == 'N'))
        return m_data.m_plugin_gpu->dgemm_fptr(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    else
        return m_data.m_plugin_cpu->dgemm_fptr(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr cgemm_gpu_wrap_impl(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                              c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb,
                              c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_gpu->cgemm_fptr(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr sgemm_gpu_wrap_impl(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                              s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb,
                              s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_gpu->sgemm_fptr(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr zgemm_gpu_wrap_impl(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                              z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb,
                              z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_gpu->zgemm_fptr(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};

i_type_wr csyrk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha,
                              c_type_wr *a, i_type_wr *lda, c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->csyrk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
};
i_type_wr dsyrk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, d_type_wr *a,
                              i_type_wr *lda, d_type_wr *beta, d_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->dsyrk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
};
i_type_wr ssyrk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha,
                              s_type_wr *a, i_type_wr *lda, s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->ssyrk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
};
i_type_wr zsyrk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha,
                              z_type_wr *a, i_type_wr *lda, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->zsyrk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
};

i_type_wr cherk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha,
                              c_type_wr *a, i_type_wr *lda, s_type_wr *beta, c_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->cherk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
};
i_type_wr zherk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha,
                              z_type_wr *a, i_type_wr *lda, d_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->zherk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
};

i_type_wr csyr2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha,
                               c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta,
                               c_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->csyr2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr dsyr2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha,
                               d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta,
                               d_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->dsyr2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr ssyr2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha,
                               s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *beta,
                               s_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->ssyr2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr zsyr2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha,
                               z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta,
                               z_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->zsyr2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};

i_type_wr ctrsm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb)
{
    return m_data.m_plugin_cpu->ctrsm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};
i_type_wr dtrsm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb)
{
    return m_data.m_plugin_cpu->dtrsm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};
i_type_wr strsm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb)
{
    return m_data.m_plugin_cpu->strsm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};
i_type_wr ztrsm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb)
{
    return m_data.m_plugin_cpu->ztrsm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};

i_type_wr csymm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, c_type_wr *alpha,
                              c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta,
                              c_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->csymm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr dsymm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, d_type_wr *alpha,
                              d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta,
                              d_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->dsymm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr ssymm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, s_type_wr *alpha,
                              s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *beta,
                              s_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->ssymm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr zsymm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, z_type_wr *alpha,
                              z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta,
                              z_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->zsymm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
};

i_type_wr chemm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, c_type_wr *alpha,
                              c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta,
                              c_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->chemm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr zhemm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, z_type_wr *alpha,
                              z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta,
                              z_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->zhemm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr cher2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha,
                               c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, s_type_wr *beta,
                               c_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->cher2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr zher2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha,
                               z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, d_type_wr *beta,
                               z_type_wr *c__, i_type_wr *ldc)
{
    return m_data.m_plugin_cpu->zher2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};

i_type_wr ctrmm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb)
{
    return m_data.m_plugin_cpu->ctrmm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};
i_type_wr dtrmm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb)
{
    return m_data.m_plugin_cpu->dtrmm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};
i_type_wr strmm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb)
{
    return m_data.m_plugin_cpu->strmm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};
i_type_wr ztrmm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb)
{
    return m_data.m_plugin_cpu->ztrmm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};

// generic stuff - the constructor which initializes all pointers to blas functions
#define CALL_SYNTAX(x) x##_gpu_wrap_impl
#include "blas_plugin_common.h"
