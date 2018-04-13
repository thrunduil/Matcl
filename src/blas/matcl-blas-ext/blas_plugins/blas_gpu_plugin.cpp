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
#include "cublas.h"
#include "cublasxt.h"
#include <Windows.h>

#undef min
#undef max

_declspec(thread) cublasXtHandle_t  cublas_handle = nullptr;
_declspec(thread) bool              cublas_handle_initialized = false;

class global_data
{
    private:
        std::atomic<long>   m_threads;

    public:
        ::blas_wrapper      m_wrapper;

    public:
        global_data();

        int         get_num_threads();
        void        set_num_threads(int n);
        bool        are_user_threads_allowed();

        static void exit_thread();
};

global_data::global_data()
    :m_threads(0)
{
    int max_threads = std::thread::hardware_concurrency();
    m_threads       = max_threads;
};
int global_data::get_num_threads()
{
    return m_threads;
};
void global_data::set_num_threads(int n)
{
    int max_threads = std::thread::hardware_concurrency();
    m_threads = std::min(std::max(n,1),max_threads);
};
bool global_data::are_user_threads_allowed()
{
    return true;
};
void global_data::exit_thread()
{
    if (cublas_handle != nullptr)
        cublasXtDestroy(cublas_handle);
}

static cublasXtHandle_t get_cublas_handle()
{
    if (cublas_handle_initialized == true)
        return cublas_handle;

    cublasStatus_t stat = cublasXtCreate(&cublas_handle);

    if (stat != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "cublas initialization error" << "\n";
        return nullptr;
    }

    int devices[1] = { 0 };
    if (cublasXtDeviceSelect(cublas_handle, 1, devices) != CUBLAS_STATUS_SUCCESS)
    { 
        std::cerr << "set devices fail" << "\n";
        return nullptr;
    }

    cublas_handle_initialized = true;

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);

    return cublas_handle;
};

static void release_cublas_handle()
{
    cublasXtDestroy(cublas_handle);
    cublas_handle = nullptr;
    cublas_handle_initialized = false;
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
            //global_data::init_thread();
            break;
        case DLL_THREAD_DETACH:
            global_data::exit_thread();
            break;
        case DLL_PROCESS_DETACH:
            break;
    }
    return TRUE;
}

extern "C"
{
    BLAS_PLUGIN_EXPORT 
    const ::blas_wrapper* get_blas_wrapper()
    {
        return &m_data.m_wrapper;
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

cublasOperation_t op(const char* type)
{
    if (type[0] == 'n' || type[0] == 'N')
        return CUBLAS_OP_N;
    if (type[0] == 't' || type[0] == 'T')
        return CUBLAS_OP_T;
    if (type[0] == 'c' || type[0] == 'C')
        return CUBLAS_OP_C;

    return CUBLAS_OP_N;
}
cublasFillMode_t op_lo(const char* type)
{
    if (type[0] == 'l' || type[0] == 'L')
        return CUBLAS_FILL_MODE_LOWER;
    if (type[0] == 'u' || type[0] == 'U')
        return CUBLAS_FILL_MODE_UPPER;

    return CUBLAS_FILL_MODE_LOWER;
}
cublasSideMode_t op_side(const char* type)
{
    if (type[0] == 'l' || type[0] == 'L')
        return CUBLAS_SIDE_LEFT;
    if (type[0] == 'r' || type[0] == 'R')
        return CUBLAS_SIDE_RIGHT;

    return CUBLAS_SIDE_LEFT;
}
cublasDiagType_t op_diag(const char* type)
{
    if (type[0] == 'u' || type[0] == 'U')
        return CUBLAS_DIAG_UNIT;
    if (type[0] == 'n' || type[0] == 'N')
        return CUBLAS_DIAG_NON_UNIT;

    return CUBLAS_DIAG_NON_UNIT;
}
//--------------------------------------------------------------------------
//                      NOT IMPLEMENTED
//--------------------------------------------------------------------------
i_type_wr not_implemented()
{
    std::cerr << "not implemented" << "\n";
    return 1;
}
i_type_wr caxpy_gpu_wrap_impl(i_type_wr *n, c_type_wr *alpha, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr ccopy_gpu_wrap_impl(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr cdotc_gpu_wrap_impl(c_type_wr *ret_val, i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr cdotu_gpu_wrap_impl(c_type_wr *ret_val, i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr cgbmv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku,
                              c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx,
                              c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr cgemv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *a,
                              i_type_wr *lda, c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr cgerc_gpu_wrap_impl(i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                              c_type_wr *y, i_type_wr *incy, c_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
i_type_wr cgeru_gpu_wrap_impl(i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                              c_type_wr *y, i_type_wr *incy, c_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
i_type_wr chbmv_gpu_wrap_impl(char *uplo, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, c_type_wr *a,
                              i_type_wr *lda, c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr chemv_gpu_wrap_impl(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda,
                              c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr cher_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                             c_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
i_type_wr cher2_gpu_wrap_impl(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                              c_type_wr *y, i_type_wr *incy, c_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
i_type_wr chpmv_gpu_wrap_impl(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *ap, c_type_wr *x,
                              i_type_wr *incx, c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr chpr_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                             c_type_wr *ap)
{
    return not_implemented();
};
i_type_wr chpr2_gpu_wrap_impl(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                              c_type_wr *y, i_type_wr *incy, c_type_wr *ap)
{
    return not_implemented();
};
i_type_wr crotg_gpu_wrap_impl(c_type_wr *ca, c_type_wr *cb, s_type_wr *c__, c_type_wr *s)
{
    return not_implemented();
};
i_type_wr cscal_gpu_wrap_impl(i_type_wr *n, c_type_wr *ca, c_type_wr *cx, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr csrot_gpu_wrap_impl(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, i_type_wr *incy,
                              s_type_wr *c__, s_type_wr *s)
{
    return not_implemented();
};
i_type_wr csscal_gpu_wrap_impl(i_type_wr *n, s_type_wr *sa, c_type_wr *cx, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr cswap_gpu_wrap_impl(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr ctbmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, c_type_wr *a,
                              i_type_wr *lda, c_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr ctbsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, c_type_wr *a,
                              i_type_wr *lda, c_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr ctpmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *ap, c_type_wr *x,
                              i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr ctpsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *ap, c_type_wr *x,
                              i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr ctrmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *a, i_type_wr *lda,
                              c_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr ctrsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *a, i_type_wr *lda,
                              c_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
d_type_wr dasum_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr daxpy_gpu_wrap_impl(i_type_wr *n, d_type_wr *da, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy,
                              i_type_wr *incy)
{
    return not_implemented();
};
d_type_wr dcabs1_gpu_wrap_impl(z_type_wr *z__)
{
    return not_implemented();
};
i_type_wr dcopy_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
{
    return not_implemented();
};
d_type_wr ddot_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr dgbmv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku,
                              d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *x, i_type_wr *incx,
                              d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr dgemv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *a,
                              i_type_wr *lda, d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr dger_gpu_wrap_impl(i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx,
                             d_type_wr *y, i_type_wr *incy, d_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
d_type_wr dnrm2_gpu_wrap_impl(i_type_wr *n, d_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr drot_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy,
                             d_type_wr *c__, d_type_wr *s)
{
    return not_implemented();
};
i_type_wr drotg_gpu_wrap_impl(d_type_wr *da, d_type_wr *db, d_type_wr *c__, d_type_wr *s)
{
    return not_implemented();
};
i_type_wr drotm_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy,
                              d_type_wr *dparam)
{
    return not_implemented();
};
i_type_wr drotmg_gpu_wrap_impl(d_type_wr *dd1, d_type_wr *dd2, d_type_wr *dx1, d_type_wr *dy1, d_type_wr *dparam)
{
    return not_implemented();
};
i_type_wr dsbmv_gpu_wrap_impl(char *uplo, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, d_type_wr *a,
                              i_type_wr *lda, d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr dscal_gpu_wrap_impl(i_type_wr *n, d_type_wr *da, d_type_wr *dx, i_type_wr *incx)
{
    return not_implemented();
};
d_type_wr dsdot_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr dspmv_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *ap, d_type_wr *x,
                              i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr dspr_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx,
                             d_type_wr *ap)
{
    return not_implemented();
};
i_type_wr dspr2_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx,
                              d_type_wr *y, i_type_wr *incy, d_type_wr *ap)
{
    return not_implemented();
};
i_type_wr dswap_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr dsymv_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda,
                              d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr dsyr_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx,
                             d_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
i_type_wr dsyr2_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx,
                              d_type_wr *y, i_type_wr *incy, d_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
i_type_wr dtbmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, d_type_wr *a,
                              i_type_wr *lda, d_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr dtbsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, d_type_wr *a,
                              i_type_wr *lda, d_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr dtpmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *ap,
                              d_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr dtpsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *ap,
                              d_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr dtrmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *a, i_type_wr *lda,
                              d_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr dtrsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *a, i_type_wr *lda,
                              d_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
d_type_wr dzasum_gpu_wrap_impl(i_type_wr *n, z_type_wr *zx, i_type_wr *incx)
{
    return not_implemented();
};
d_type_wr dznrm2_gpu_wrap_impl(i_type_wr *n, z_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr icamax_gpu_wrap_impl(i_type_wr *n, c_type_wr *cx, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr idamax_gpu_wrap_impl(i_type_wr *n, d_type_wr *dx, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr isamax_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr izamax_gpu_wrap_impl(i_type_wr *n, z_type_wr *zx, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr lsame_gpu_wrap_impl(char *ca, char *cb)
{
    return not_implemented();
};
d_type_wr sasum_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr saxpy_gpu_wrap_impl(i_type_wr *n, s_type_wr *sa, s_type_wr *sx, i_type_wr *incx,
                              s_type_wr *sy, i_type_wr *incy)
{
    return not_implemented();
};
d_type_wr scabs1_gpu_wrap_impl(c_type_wr *z__)
{
    return not_implemented();
};
d_type_wr scasum_gpu_wrap_impl(i_type_wr *n, c_type_wr *cx, i_type_wr *incx)
{
    return not_implemented();
};
d_type_wr scnrm2_gpu_wrap_impl(i_type_wr *n, c_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr scopy_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
{
    return not_implemented();
};
s_type_wr sdot_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
{
    return (s_type_wr)not_implemented();
};
d_type_wr sdsdot_gpu_wrap_impl(i_type_wr *n, s_type_wr *sb, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy,
                               i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr sgbmv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku,
                              s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *x, i_type_wr *incx,
                              s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr sgemv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *a,
                              i_type_wr *lda, s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr sger_gpu_wrap_impl(i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx,
                             s_type_wr *y, i_type_wr *incy, s_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
d_type_wr snrm2_gpu_wrap_impl(i_type_wr *n, s_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr srot_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy,
                             s_type_wr *c__, s_type_wr *s)
{
    return not_implemented();
};
i_type_wr srotg_gpu_wrap_impl(s_type_wr *sa, s_type_wr *sb, s_type_wr *c__, s_type_wr *s)
{
    return not_implemented();
};
i_type_wr srotm_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy,
                              s_type_wr *sparam)
{
    return not_implemented();
};
i_type_wr srotmg_gpu_wrap_impl(s_type_wr *sd1, s_type_wr *sd2, s_type_wr *sx1, s_type_wr *sy1, s_type_wr *sparam)
{
    return not_implemented();
};
i_type_wr ssbmv_gpu_wrap_impl(char *uplo, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, s_type_wr *a,
                              i_type_wr *lda, s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr sscal_gpu_wrap_impl(i_type_wr *n, s_type_wr *sa, s_type_wr *sx, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr sspmv_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *ap, s_type_wr *x,
                              i_type_wr *incx, s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr sspr_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx,
                             s_type_wr *ap)
{
    return not_implemented();
};
i_type_wr sspr2_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx,
                              s_type_wr *y, i_type_wr *incy, s_type_wr *ap)
{
    return not_implemented();
};
i_type_wr sswap_gpu_wrap_impl(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr ssymv_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda,
                              s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr ssyr_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx,
                             s_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
i_type_wr ssyr2_gpu_wrap_impl(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx,
                              s_type_wr *y, i_type_wr *incy, s_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
i_type_wr stbmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, s_type_wr *a,
                              i_type_wr *lda, s_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr stbsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, s_type_wr *a,
                              i_type_wr *lda, s_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr stpmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *ap, s_type_wr *x,
                              i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr stpsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *ap, s_type_wr *x,
                              i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr strmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *a, i_type_wr *lda,
                              s_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr strsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *a, i_type_wr *lda,
                              s_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr xerbla_gpu_wrap_impl(char *srname, i_type_wr *info)
{
    return not_implemented();
};
i_type_wr zaxpy_gpu_wrap_impl(i_type_wr *n, z_type_wr *za, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr zcopy_gpu_wrap_impl(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr zdotc_gpu_wrap_impl(z_type_wr *res, i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr zdotu_gpu_wrap_impl(z_type_wr *res, i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr zdrot_gpu_wrap_impl(i_type_wr *n, z_type_wr *cx, i_type_wr *incx, z_type_wr *cy, i_type_wr *incy,
                              d_type_wr *c__, d_type_wr *s)
{
    return not_implemented();
};
i_type_wr zdscal_gpu_wrap_impl(i_type_wr *n, d_type_wr *da, z_type_wr *zx, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr zgbmv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku,
                              z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx,
                              z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr zgemv_gpu_wrap_impl(char *trans, i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *a,
                              i_type_wr *lda, z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr zgerc_gpu_wrap_impl(i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                              z_type_wr *y, i_type_wr *incy, z_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
i_type_wr zgeru_gpu_wrap_impl(i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                              z_type_wr *y, i_type_wr *incy, z_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
i_type_wr zhbmv_gpu_wrap_impl(char *uplo, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, z_type_wr *a,
                              i_type_wr *lda, z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y,
                              i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr zhemv_gpu_wrap_impl(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda,
                              z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr zher_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                             z_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
i_type_wr zher2_gpu_wrap_impl(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                              z_type_wr *y, i_type_wr *incy, z_type_wr *a, i_type_wr *lda)
{
    return not_implemented();
};
i_type_wr zhpmv_gpu_wrap_impl(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *ap, z_type_wr *x,
                              i_type_wr *incx, z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr zhpr_gpu_wrap_impl(char *uplo, i_type_wr *n, d_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                             z_type_wr *ap)
{
    return not_implemented();
};
i_type_wr zhpr2_gpu_wrap_impl(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                              z_type_wr *y, i_type_wr *incy, z_type_wr *ap)
{
    return not_implemented();
};
i_type_wr zrotg_gpu_wrap_impl(z_type_wr *ca, z_type_wr *cb, d_type_wr *c__, z_type_wr *s)
{
    return not_implemented();
};
i_type_wr zscal_gpu_wrap_impl(i_type_wr *n, z_type_wr *za, z_type_wr *zx, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr zswap_gpu_wrap_impl(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, i_type_wr *incy)
{
    return not_implemented();
};
i_type_wr ztbmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, z_type_wr *a,
                              i_type_wr *lda, z_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr ztbsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, z_type_wr *a,
                              i_type_wr *lda, z_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr ztpmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *ap, z_type_wr *x,
                              i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr ztpsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *ap, z_type_wr *x,
                              i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr ztrmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *a, i_type_wr *lda,
                              z_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr ztrsv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *a, i_type_wr *lda,
                              z_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};

//--------------------------------------------------------------------------
//                      HELPERS
//--------------------------------------------------------------------------
static bool is_zero(s_type_wr val)
{
    return val == 0;
};
static bool is_zero(d_type_wr val)
{
    return val == 0;
};
static bool is_zero(c_type_wr val)
{
    return val.r == 0 && val.i == 0;
};
static bool is_zero(z_type_wr val)
{
    return val.r == 0 && val.i == 0;
};

static float gpu_c(s_type_wr val)               { return val; }
static double gpu_c(d_type_wr val)              { return val; }
static cuComplex gpu_c(c_type_wr val)           { cuComplex ret; ret.x = val.r; ret.y = val.i; return ret; }
static cuDoubleComplex gpu_c(z_type_wr val)     { cuDoubleComplex ret; ret.x = val.r; ret.y = val.i; return ret; }
static float* gpu_c(s_type_wr* val)             { return val; }
static double* gpu_c(d_type_wr* val)            { return val; }
static cuComplex* gpu_c(c_type_wr* val)         { return reinterpret_cast<cuComplex*>(val); };
static cuDoubleComplex* gpu_c(z_type_wr* val)   { return reinterpret_cast<cuDoubleComplex*>(val); };

template<class T>   struct cublas_type              {};
template<>          struct cublas_type<s_type_wr>   { typedef float type; };
template<>          struct cublas_type<d_type_wr>   { typedef double type; };
template<>          struct cublas_type<c_type_wr>   { typedef cuComplex type; };
template<>          struct cublas_type<z_type_wr>   { typedef cuDoubleComplex type; };

template<class T>
cublasStatus_t set_matrix(int M, int N, int LD, const T* a, void** AA, bool set)
{
    cublasStatus_t status;

    status = cublasAlloc(M * N, sizeof(T), AA);
    if (status != CUBLAS_STATUS_SUCCESS)
        return status;

    if (set == false)
        return status;

    status = cublasSetMatrix(M, N, sizeof(T), a, LD, *AA, M);
    return status;
};

template<class V>
cublasStatus_t get_matrix(int M, int N, int ldc, const void* CC, V* c)
{
    cublasStatus_t status = cublasGetMatrix(M, N, sizeof(V), CC, M, c, ldc);
    return status;
};

//--------------------------------------------------------------------------
//                      GEMM
//--------------------------------------------------------------------------
static cublasStatus_t cublasXt_gemm(cublasXtHandle_t  handle, cublasOperation_t transa, cublasOperation_t transb,
                                    size_t m, size_t n, size_t k, const  float *alpha, const  float *A, size_t lda,
                                    const  float *B, size_t ldb, const  float *beta, float *C, size_t ldc)
{
    return cublasXtSgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
};
static cublasStatus_t cublasXt_gemm(cublasXtHandle_t  handle, cublasOperation_t transa, cublasOperation_t transb,
                                    size_t m, size_t n, size_t k, const  double *alpha, const  double *A, size_t lda,
                                    const  double *B, size_t ldb, const  double *beta, double *C, size_t ldc)
{
    return cublasXtDgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
};
static cublasStatus_t cublasXt_gemm(cublasXtHandle_t  handle, cublasOperation_t transa, cublasOperation_t transb,
                                    size_t m, size_t n, size_t k, const  cuComplex *alpha, const  cuComplex *A, size_t lda,
                                    const  cuComplex *B, size_t ldb, const  cuComplex *beta, cuComplex *C, size_t ldc)
{
    return cublasXtCgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
};
static cublasStatus_t cublasXt_gemm(cublasXtHandle_t  handle, cublasOperation_t transa, cublasOperation_t transb,
                                    size_t m, size_t n, size_t k, const  cuDoubleComplex *alpha, const  cuDoubleComplex *A, size_t lda,
                                    const  cuDoubleComplex *B, size_t ldb, const  cuDoubleComplex *beta, cuDoubleComplex *C, size_t ldc)
{
    return cublasXtZgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
};

template<class V>
static i_type_wr gemm_gpu_generic(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                              V *alpha, V *a, i_type_wr *lda, V *b, i_type_wr *ldb,
                              V *beta, V *c__, i_type_wr *ldc)
{
    typedef cublas_type<V>::type value_type;

    int A_M = 0;
    int A_N = 0;
    int B_M = 0;
    int B_N = 0;
    int C_M = *m;
    int C_N = *n;

    if (op(transa) == CUBLAS_OP_N)
    {
        A_M = *m;
        A_N = *k;
    }
    else
    {
        A_M = *k;
        A_N = *m;
    }
    if (op(transb) == CUBLAS_OP_N)
    {
        B_M = *k;
        B_N = *n;
    }
    else
    {
        B_M = *n;
        B_N = *k;        
    }

    value_type* AA = nullptr;
    value_type* BB = nullptr;
    value_type* CC = nullptr;

    cublasStatus_t status;

    status = set_matrix<V>(A_M, A_N, *lda, a, (void**)&AA, true);    
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = set_matrix<V>(B_M, B_N, *ldb, b, (void**)&BB, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = set_matrix<V>(C_M, C_N, *ldc, c__, (void**)&CC, is_zero(*beta) == false);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = cublasXt_gemm(get_cublas_handle(), op(transa), op(transb), *m, *n, *k, gpu_c(alpha), AA, A_M,
                           BB, B_M, gpu_c(beta), CC, C_M);

    if (status != CUBLAS_STATUS_SUCCESS)
    {
        release_cublas_handle();
        std::cerr << "cublass kernel error" << "\n";
        goto lab_error;
    }

    status = get_matrix<V>(C_M, C_N, *ldc, CC, c__);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "device read error" << "\n";
        goto lab_error;
    }

    cublasFree(AA);
    cublasFree(BB);
    cublasFree(CC);

    return 0;

  lab_error:
    cublasFree(AA);
    cublasFree(BB);
    cublasFree(CC);

    return 1;
};

i_type_wr dgemm_gpu_wrap_impl(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                              d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb,
                              d_type_wr *beta, d_type_wr *c__, i_type_wr *ldc)
{
    return gemm_gpu_generic<d_type_wr>(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr cgemm_gpu_wrap_impl(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                              c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb,
                              c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc)
{
    return gemm_gpu_generic<c_type_wr>(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr sgemm_gpu_wrap_impl(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                              s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb,
                              s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc)
{
    return gemm_gpu_generic<s_type_wr>(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr zgemm_gpu_wrap_impl(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                              z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb,
                              z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
{
    return gemm_gpu_generic<z_type_wr>(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
//--------------------------------------------------------------------------
//                      SYRK
//--------------------------------------------------------------------------
cublasStatus_t cublasXt_syrk(cublasXtHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans,
                           size_t n, size_t k, const float *alpha, const float *A, size_t lda, const float *beta,
                           float *C, size_t ldc)
{
    return cublasXtSsyrk(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
};
cublasStatus_t cublasXt_syrk(cublasXtHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans,
                           size_t n, size_t k, const double *alpha, const double *A, size_t lda, const double *beta,
                           double *C, size_t ldc)
{
    return cublasXtDsyrk(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
};
cublasStatus_t cublasXt_syrk(cublasXtHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans,
                           size_t n, size_t k, const cuComplex *alpha, const cuComplex *A, size_t lda, const cuComplex *beta,
                           cuComplex *C, size_t ldc)
{
    return cublasXtCsyrk(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
};
cublasStatus_t cublasXt_syrk(cublasXtHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans,
                           size_t n, size_t k, const cuDoubleComplex *alpha, const cuDoubleComplex *A, size_t lda, 
                           const cuDoubleComplex *beta, cuDoubleComplex *C, size_t ldc)
{
    return cublasXtZsyrk(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
};

template<class V>
i_type_wr syrk_gpu_generic(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, V *alpha,
                           V *a, i_type_wr *lda, V *beta, V *c__, i_type_wr *ldc)
{
    typedef cublas_type<V>::type value_type;

    int A_M = 0;
    int A_N = 0;
    int C_M = *n;
    int C_N = *n;

    if (op(trans) == CUBLAS_OP_N)
    {
        A_M = *n;
        A_N = *k;
    }
    else
    {
        A_M = *k;
        A_N = *n;
    }

    value_type* AA = nullptr;
    value_type* CC = nullptr;

    cublasStatus_t status;

    status = set_matrix<V>(A_M, A_N, *lda, a, (void**)&AA, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = set_matrix<V>(C_M, C_N, *ldc, c__, (void**)&CC, is_zero(*beta) == false);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = cublasXt_syrk(get_cublas_handle(), op_lo(uplo), op(trans), *n, *k, gpu_c(alpha), AA, A_M,
                           gpu_c(beta), CC, C_M);

    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "cublass kernel error" << "\n";
        goto lab_error;
    }

    status = get_matrix<V>(C_M, C_N, *ldc, CC, c__);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "device read error" << "\n";
        goto lab_error;
    }

    cublasFree(AA);
    cublasFree(CC);

    return 0;

  lab_error:
    cublasFree(AA);
    cublasFree(CC);

    return 1;
};

i_type_wr csyrk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha,
                              c_type_wr *a, i_type_wr *lda, c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc)
{
    return syrk_gpu_generic<c_type_wr>(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
};
i_type_wr dsyrk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, d_type_wr *a,
                              i_type_wr *lda, d_type_wr *beta, d_type_wr *c__, i_type_wr *ldc)
{
    return syrk_gpu_generic<d_type_wr>(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
};
i_type_wr ssyrk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha,
                              s_type_wr *a, i_type_wr *lda, s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc)
{
    return syrk_gpu_generic<s_type_wr>(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
};
i_type_wr zsyrk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha,
                              z_type_wr *a, i_type_wr *lda, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
{
    return syrk_gpu_generic<z_type_wr>(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
};

//--------------------------------------------------------------------------
//                      HERK
//--------------------------------------------------------------------------
cublasStatus_t cublasXt_herk(cublasXtHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans,
                             size_t n, size_t k, const float *alpha, const cuComplex *A, size_t lda,
                             const float *beta, cuComplex *C, size_t ldc)
{
    return cublasXtCherk(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
};
cublasStatus_t cublasXt_herk(cublasXtHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans,
                             size_t n, size_t k, const double *alpha, const cuDoubleComplex *A, size_t lda,
                             const double *beta, cuDoubleComplex *C, size_t ldc)
{
    return cublasXtZherk(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
};

template<class V, class VF>
i_type_wr herk_gpu_generic(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, VF *alpha,
                           V *a, i_type_wr *lda, VF *beta, V *c__, i_type_wr *ldc)
{
    typedef cublas_type<V>::type value_type;

    int A_M = 0;
    int A_N = 0;
    int C_M = *n;
    int C_N = *n;

    if (op(trans) == CUBLAS_OP_N)
    {
        A_M = *n;
        A_N = *k;
    }
    else
    {
        A_M = *k;
        A_N = *n;
    }

    value_type* AA = nullptr;
    value_type* CC = nullptr;

    cublasStatus_t status;

    status = set_matrix<V>(A_M, A_N, *lda, a, (void**)&AA, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = set_matrix<V>(C_M, C_N, *ldc, c__, (void**)&CC, is_zero(*beta) == false);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = cublasXt_herk(get_cublas_handle(), op_lo(uplo), op(trans), *n, *k, gpu_c(alpha), AA, A_M,
                           gpu_c(beta), CC, C_M);

    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "cublass kernel error" << "\n";
        goto lab_error;
    }

    status = get_matrix<V>(C_M, C_N, *ldc, CC, c__);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "device read error" << "\n";
        goto lab_error;
    }

    cublasFree(AA);
    cublasFree(CC);

    return 0;

  lab_error:
    cublasFree(AA);
    cublasFree(CC);

    return 1;
};

i_type_wr cherk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha,
                              c_type_wr *a, i_type_wr *lda, s_type_wr *beta, c_type_wr *c__, i_type_wr *ldc)
{
    return herk_gpu_generic<c_type_wr>(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
};
i_type_wr zherk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha,
                              z_type_wr *a, i_type_wr *lda, d_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
{
    return herk_gpu_generic<z_type_wr>(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
};

//--------------------------------------------------------------------------
//                      SYR2K
//--------------------------------------------------------------------------
cublasStatus_t cublasXt_syr2k(cublasXtHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans,
                              size_t n, size_t k, const float *alpha, const float *A, size_t lda,
                              const float *B, size_t ldb, const float *beta, float *C, size_t ldc)
{
    return cublasXtSsyr2k(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
cublasStatus_t cublasXt_syr2k(cublasXtHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans,
                              size_t n, size_t k, const double *alpha, const double *A, size_t lda,
                              const double *B, size_t ldb, const double *beta, double *C, size_t ldc)
{
    return cublasXtDsyr2k(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
cublasStatus_t cublasXt_syr2k(cublasXtHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans,
                              size_t n, size_t k, const cuComplex *alpha, const cuComplex *A, size_t lda,
                              const cuComplex *B, size_t ldb, const cuComplex *beta, cuComplex *C, size_t ldc)
{
    return cublasXtCsyr2k(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
cublasStatus_t cublasXt_syr2k(cublasXtHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans,
                              size_t n, size_t k, const cuDoubleComplex *alpha, const cuDoubleComplex *A, size_t lda,
                              const cuDoubleComplex *B, size_t ldb, const cuDoubleComplex *beta, cuDoubleComplex *C, size_t ldc)
{
    return cublasXtZsyr2k(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

template<class V>
i_type_wr syr2k_gpu_generic(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, V *alpha,
                            V *a, i_type_wr *lda, V *b, i_type_wr *ldb, V *beta, V *c__, i_type_wr *ldc)
{
    typedef cublas_type<V>::type value_type;

    int A_M = 0;
    int A_N = 0;
    int B_M = 0;
    int B_N = 0;
    int C_M = *n;
    int C_N = *n;

    if (op(trans) == CUBLAS_OP_N)
    {
        A_M = *n;
        A_N = *k;
        B_M = *n;
        B_N = *k;
    }
    else
    {
        A_M = *k;
        A_N = *n;
        B_M = *k;
        B_N = *n;
    }

    value_type* AA = nullptr;
    value_type* BB = nullptr;
    value_type* CC = nullptr;

    cublasStatus_t status;

    status = set_matrix<V>(A_M, A_N, *lda, a, (void**)&AA, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }
    status = set_matrix<V>(B_M, B_N, *ldb, b, (void**)&BB, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = set_matrix<V>(C_M, C_N, *ldc, c__, (void**)&CC, is_zero(*beta) == false);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = cublasXt_syr2k(get_cublas_handle(), op_lo(uplo), op(trans), *n, *k, gpu_c(alpha), AA, A_M, BB, B_M,
                           gpu_c(beta), CC, C_M);

    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "cublass kernel error" << "\n";
        goto lab_error;
    }

    status = get_matrix<V>(C_M, C_N, *ldc, CC, c__);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "device read error" << "\n";
        goto lab_error;
    }

    cublasFree(AA);
    cublasFree(BB);
    cublasFree(CC);

    return 0;

  lab_error:
    cublasFree(AA);
    cublasFree(BB);
    cublasFree(CC);

    return 1;
};

i_type_wr csyr2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha,
                               c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta,
                               c_type_wr *c__, i_type_wr *ldc)
{
    return syr2k_gpu_generic<c_type_wr>(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr dsyr2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha,
                               d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta,
                               d_type_wr *c__, i_type_wr *ldc)
{
    return syr2k_gpu_generic<d_type_wr>(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr ssyr2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha,
                               s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *beta,
                               s_type_wr *c__, i_type_wr *ldc)
{
    return syr2k_gpu_generic<s_type_wr>(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr zsyr2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha,
                               z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta,
                               z_type_wr *c__, i_type_wr *ldc)
{
    return syr2k_gpu_generic<z_type_wr>(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};

//--------------------------------------------------------------------------
//                      TRSM
//--------------------------------------------------------------------------
cublasStatus_t cublasXt_trsm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans,
                             cublasDiagType_t diag, size_t m, size_t n, const float *alpha, const float *A, size_t lda,
                             float *B, size_t ldb)
{
    return cublasXtStrsm(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
};
cublasStatus_t cublasXt_trsm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans,
                             cublasDiagType_t diag, size_t m, size_t n, const double *alpha, const double *A, size_t lda,
                             double *B, size_t ldb)
{
    return cublasXtDtrsm(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
};
cublasStatus_t cublasXt_trsm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans,
                             cublasDiagType_t diag, size_t m, size_t n, const cuComplex *alpha, const cuComplex *A, size_t lda,
                             cuComplex *B, size_t ldb)
{
    return cublasXtCtrsm(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
};
cublasStatus_t cublasXt_trsm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans,
                             cublasDiagType_t diag, size_t m, size_t n, const cuDoubleComplex *alpha, const cuDoubleComplex *A, size_t lda,
                             cuDoubleComplex *B, size_t ldb)
{
    return cublasXtZtrsm(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
};

template<class V>
i_type_wr trsm_gpu_generic(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              V *alpha, V *a, i_type_wr *lda, V *b, i_type_wr *ldb)
{
    typedef cublas_type<V>::type value_type;

    int A_M = 0;
    int A_N = 0;
    int B_M = *m;
    int B_N = *n;

    if (op_side(side) == CUBLAS_SIDE_LEFT)
    {
        A_M = *m;
        A_N = *m;
    }
    else
    {
        A_M = *n;
        A_N = *n;
    }

    value_type* AA = nullptr;
    value_type* BB = nullptr;

    cublasStatus_t status;

    status = set_matrix<V>(A_M, A_N, *lda, a, (void**)&AA, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }
    status = set_matrix<V>(B_M, B_N, *ldb, b, (void**)&BB, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = cublasXt_trsm(get_cublas_handle(), op_side(side), op_lo(uplo), op(transa), op_diag(diag), 
                           *m, *n, gpu_c(alpha), AA, A_M, BB, B_M);

    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "cublass kernel error" << "\n";
        goto lab_error;
    }

    status = get_matrix<V>(B_M, B_N, *ldb, BB, b);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "device read error" << "\n";
        goto lab_error;
    }

    cublasFree(AA);
    cublasFree(BB);

    return 0;

  lab_error:
    cublasFree(AA);
    cublasFree(BB);

    return 1;
}
i_type_wr ctrsm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb)
{
    return trsm_gpu_generic<c_type_wr>(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};
i_type_wr dtrsm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb)
{
    return trsm_gpu_generic<d_type_wr>(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};
i_type_wr strsm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb)
{
    return trsm_gpu_generic<s_type_wr>(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};
i_type_wr ztrsm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb)
{
    return trsm_gpu_generic<z_type_wr>(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};

//--------------------------------------------------------------------------
//                      SYMM
//--------------------------------------------------------------------------
cublasStatus_t cublasXt_symm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo,
                             size_t m, size_t n, const float *alpha, const float *A, size_t lda,
                             const float *B, size_t ldb, const float *beta, float *C, size_t ldc)
{
    return cublasXtSsymm(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
};
cublasStatus_t cublasXt_symm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo,
                             size_t m, size_t n, const double *alpha, const double *A, size_t lda,
                             const double *B, size_t ldb, const double *beta, double *C, size_t ldc)
{
    return cublasXtDsymm(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
};
cublasStatus_t cublasXt_symm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo,
                             size_t m, size_t n, const cuComplex *alpha, const cuComplex *A, size_t lda,
                             const cuComplex *B, size_t ldb, const cuComplex *beta, cuComplex *C, size_t ldc)
{
    return cublasXtCsymm(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
};
cublasStatus_t cublasXt_symm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo,
                             size_t m, size_t n, const cuDoubleComplex *alpha, const cuDoubleComplex *A, size_t lda,
                             const cuDoubleComplex *B, size_t ldb, const cuDoubleComplex *beta, cuDoubleComplex *C, size_t ldc)
{
    return cublasXtZsymm(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
};
template<class V>
i_type_wr symm_gpu_generic(char *side, char *uplo, i_type_wr *m, i_type_wr *n, V *alpha,
                              V *a, i_type_wr *lda, V *b, i_type_wr *ldb, V *beta,
                              V *c__, i_type_wr *ldc)
{
    typedef cublas_type<V>::type value_type;

    int A_M = 0;
    int A_N = 0;
    int B_M = *m;
    int B_N = *n;
    int C_M = *m;
    int C_N = *n;

    if (op_side(side) == CUBLAS_SIDE_LEFT)
    {
        A_M = *m;
        A_N = *m;
    }
    else
    {
        A_M = *n;
        A_N = *n;        
    }

    value_type* AA = nullptr;
    value_type* BB = nullptr;
    value_type* CC = nullptr;

    cublasStatus_t status;

    status = set_matrix<V>(A_M, A_N, *lda, a, (void**)&AA, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = set_matrix<V>(B_M, B_N, *ldb, b, (void**)&BB, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = set_matrix<V>(C_M, C_N, *ldc, c__, (void**)&CC, is_zero(*beta) == false);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = cublasXt_symm(get_cublas_handle(), op_side(side), op_lo(uplo), *m, *n, gpu_c(alpha), AA,
                           A_M, BB, B_M, gpu_c(beta), CC, C_M);

    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "cublass kernel error" << "\n";
        goto lab_error;
    }

    status = get_matrix<V>(C_M, C_N, *ldc, CC, c__);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "device read error" << "\n";
        goto lab_error;
    }

    cublasFree(AA);
    cublasFree(BB);
    cublasFree(CC);

    return 0;

  lab_error:
    cublasFree(AA);
    cublasFree(BB);
    cublasFree(CC);

    return 1;
}
i_type_wr csymm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, c_type_wr *alpha,
                              c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta,
                              c_type_wr *c__, i_type_wr *ldc)
{
    return symm_gpu_generic<c_type_wr>(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr dsymm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, d_type_wr *alpha,
                              d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta,
                              d_type_wr *c__, i_type_wr *ldc)
{
    return symm_gpu_generic<d_type_wr>(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr ssymm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, s_type_wr *alpha,
                              s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *beta,
                              s_type_wr *c__, i_type_wr *ldc)
{
    return symm_gpu_generic<s_type_wr>(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr zsymm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, z_type_wr *alpha,
                              z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta,
                              z_type_wr *c__, i_type_wr *ldc)
{
    return symm_gpu_generic<z_type_wr>(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
};

//--------------------------------------------------------------------------
//                      HEMM
//--------------------------------------------------------------------------
cublasStatus_t cublasXt_hemm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo,
                             size_t m, size_t n, const cuComplex *alpha, const cuComplex *A, size_t lda,
                             const cuComplex *B, size_t ldb, const cuComplex *beta, cuComplex *C, size_t ldc)
{
    return cublasXtChemm(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
};
cublasStatus_t cublasXt_hemm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo,
                             size_t m, size_t n, const cuDoubleComplex *alpha, const cuDoubleComplex *A, size_t lda,
                             const cuDoubleComplex *B, size_t ldb, const cuDoubleComplex *beta, cuDoubleComplex *C, size_t ldc)
{
    return cublasXtZhemm(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
};
template<class V>
i_type_wr hemm_gpu_generic(char *side, char *uplo, i_type_wr *m, i_type_wr *n, V *alpha,
                           V *a, i_type_wr *lda, V *b, i_type_wr *ldb, V *beta,
                           V *c__, i_type_wr *ldc)
{
    typedef cublas_type<V>::type value_type;

    int A_M = 0;
    int A_N = 0;
    int B_M = *m;
    int B_N = *n;
    int C_M = *m;
    int C_N = *n;

    if (op_side(side) == CUBLAS_SIDE_LEFT)
    {
        A_M = *m;
        A_N = *m;
    }
    else
    {
        A_M = *n;
        A_N = *n;
    }

    value_type* AA = nullptr;
    value_type* BB = nullptr;
    value_type* CC = nullptr;

    cublasStatus_t status;

    status = set_matrix<V>(A_M, A_N, *lda, a, (void**)&AA, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = set_matrix<V>(B_M, B_N, *ldb, b, (void**)&BB, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = set_matrix<V>(C_M, C_N, *ldc, c__, (void**)&CC, is_zero(*beta) == false);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = cublasXt_hemm(get_cublas_handle(), op_side(side), op_lo(uplo), *m, *n, gpu_c(alpha), AA,
                           A_M, BB, B_M, gpu_c(beta), CC, C_M);

    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "cublass kernel error" << "\n";
        goto lab_error;
    }

    status = get_matrix<V>(C_M, C_N, *ldc, CC, c__);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "device read error" << "\n";
        goto lab_error;
    }

    cublasFree(AA);
    cublasFree(BB);
    cublasFree(CC);

    return 0;

  lab_error:
    cublasFree(AA);
    cublasFree(BB);
    cublasFree(CC);

    return 1;
}
i_type_wr chemm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, c_type_wr *alpha,
                              c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta,
                              c_type_wr *c__, i_type_wr *ldc)
{
    return hemm_gpu_generic<c_type_wr>(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr zhemm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, z_type_wr *alpha,
                              z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta,
                              z_type_wr *c__, i_type_wr *ldc)
{
    return hemm_gpu_generic<z_type_wr>(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
};
//--------------------------------------------------------------------------
//                      HER2K
//--------------------------------------------------------------------------
cublasStatus_t cublasXt_her2k(cublasXtHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans,
                              size_t n, size_t k, const cuComplex *alpha, const cuComplex *A, size_t lda,
                              const cuComplex *B, size_t ldb, const float *beta, cuComplex *C, size_t ldc)
{
    return cublasXtCher2k(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
};
cublasStatus_t cublasXt_her2k(cublasXtHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans,
                              size_t n, size_t k, const cuDoubleComplex *alpha, const cuDoubleComplex *A, size_t lda,
                              const cuDoubleComplex *B, size_t ldb, const double *beta, cuDoubleComplex *C, size_t ldc)
{
    return cublasXtZher2k(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
};

template<class V, class VF>
i_type_wr her2k_gpu_generic(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, V *alpha,
                            V *a, i_type_wr *lda, V *b, i_type_wr *ldb, VF *beta,
                            V *c__, i_type_wr *ldc)
{
    typedef cublas_type<V>::type value_type;

    int A_M = 0;
    int A_N = 0;
    int B_M = 0;
    int B_N = 0;
    int C_M = *n;
    int C_N = *n;

    if (op(trans) == CUBLAS_OP_N)
    {
        A_M = *n;
        A_N = *k;
        B_M = *n;
        B_N = *k;
    }
    else
    {
        A_M = *k;
        A_N = *n;
        B_M = *k;
        B_N = *n;
    }

    value_type* AA = nullptr;
    value_type* BB = nullptr;
    value_type* CC = nullptr;

    cublasStatus_t status;

    status = set_matrix<V>(A_M, A_N, *lda, a, (void**)&AA, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }
    status = set_matrix<V>(B_M, B_N, *ldb, b, (void**)&BB, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = set_matrix<V>(C_M, C_N, *ldc, c__, (void**)&CC, is_zero(*beta) == false);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = cublasXt_her2k(get_cublas_handle(), op_lo(uplo), op(trans), *n, *k, gpu_c(alpha), AA, A_M,
                            BB, B_M, gpu_c(beta), CC, C_M);

    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "cublass kernel error" << "\n";
        goto lab_error;
    }

    status = get_matrix<V>(C_M, C_N, *ldc, CC, c__);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "device read error" << "\n";
        goto lab_error;
    }

    cublasFree(AA);
    cublasFree(BB);
    cublasFree(CC);

    return 0;

  lab_error:
    cublasFree(AA);
    cublasFree(BB);
    cublasFree(CC);

    return 1;
};
i_type_wr cher2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha,
                               c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, s_type_wr *beta,
                               c_type_wr *c__, i_type_wr *ldc)
{
    return her2k_gpu_generic<c_type_wr>(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr zher2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha,
                               z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, d_type_wr *beta,
                               z_type_wr *c__, i_type_wr *ldc)
{
    return her2k_gpu_generic<z_type_wr>(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};

//--------------------------------------------------------------------------
//                      TRMM
//--------------------------------------------------------------------------
cublasStatus_t cublasXt_trmm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans,
                             cublasDiagType_t diag, size_t m, size_t n, const float *alpha, const float *A, size_t lda,
                             const float *B, size_t ldb, float *C, size_t ldc)
{
    return cublasXtStrmm(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc);
};
cublasStatus_t cublasXt_trmm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans,
                             cublasDiagType_t diag, size_t m, size_t n, const double *alpha, const double *A, size_t lda,
                             const double *B, size_t ldb, double *C, size_t ldc)
{
    return cublasXtDtrmm(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc);
};
cublasStatus_t cublasXt_trmm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans,
                             cublasDiagType_t diag, size_t m, size_t n, const cuComplex *alpha, const cuComplex *A, size_t lda,
                             const cuComplex *B, size_t ldb, cuComplex *C, size_t ldc)
{
    return cublasXtCtrmm(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc);
};
cublasStatus_t cublasXt_trmm(cublasXtHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans,
                             cublasDiagType_t diag, size_t m, size_t n, const cuDoubleComplex *alpha, const cuDoubleComplex *A, size_t lda,
                             const cuDoubleComplex *B, size_t ldb, cuDoubleComplex *C, size_t ldc)
{
    return cublasXtZtrmm(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc);
};

template<class V>
i_type_wr trmm_gpu_generic(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              V *alpha, V *a, i_type_wr *lda, V *b, i_type_wr *ldb)
{
    typedef cublas_type<V>::type value_type;

    int A_M = 0;
    int A_N = 0;
    int B_M = *m;
    int B_N = *n;

    if (op_side(side) == CUBLAS_SIDE_LEFT)
    {
        A_M = *m;
        A_N = *m;
    }
    else
    {
        A_M = *n;
        A_N = *n;
    }

    value_type* AA = nullptr;
    value_type* BB = nullptr;

    cublasStatus_t status;

    status = set_matrix<V>(A_M, A_N, *lda, a, (void**)&AA, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }
    status = set_matrix<V>(B_M, B_N, *ldb, b, (void**)&BB, true);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    }

    status = cublasXt_trmm(get_cublas_handle(), op_side(side), op_lo(uplo), op(transa), op_diag(diag), 
                           *m, *n, gpu_c(alpha), AA, A_M, BB, B_M, BB, B_M);

    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "cublass kernel error" << "\n";
        goto lab_error;
    }

    status = get_matrix<V>(B_M, B_N, *ldb, BB, b);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "device read error" << "\n";
        goto lab_error;
    }

    cublasFree(AA);
    cublasFree(BB);

    return 0;

  lab_error:
    cublasFree(AA);
    cublasFree(BB);

    return 1;
};

i_type_wr ctrmm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb)
{
    return trmm_gpu_generic<c_type_wr>(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};
i_type_wr dtrmm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb)
{
    return trmm_gpu_generic<d_type_wr>(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};
i_type_wr strmm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb)
{
    return trmm_gpu_generic<s_type_wr>(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};
i_type_wr ztrmm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb)
{
    return trmm_gpu_generic<z_type_wr>(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
};

// generic stuff - the constructor which initializes all pointers to blas functions
#define CALL_SYNTAX(x) x##_gpu_wrap_impl
#include "blas_plugin_common.h"
