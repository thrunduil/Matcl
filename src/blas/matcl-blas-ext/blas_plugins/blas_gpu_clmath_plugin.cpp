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

#include "clAmdBlas.h"
#include <Windows.h>
#undef min
#undef max

//_declspec(thread) cublasXtHandle_t  cublas_handle = nullptr;
//_declspec(thread) bool              cublas_handle_initialized = false;

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
    m_threads = max_threads;
};
int global_data::get_num_threads()
{
    return m_threads;
};
void global_data::set_num_threads(int n)
{
    int max_threads = std::thread::hardware_concurrency();
    m_threads = std::min(std::max(n, 1), max_threads);
};
bool global_data::are_user_threads_allowed()
{
    return true;
};

void global_data::exit_thread()
{
    /*
    if (cublas_handle != nullptr)
        cublasXtDestroy(cublas_handle);
    */
}

/*
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

    cublas_handle_initialized = true;

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);

    return cublas_handle;
};
*/
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

clAmdBlasTranspose op(const char* type)
{
    if (type[0] == 'n' || type[0] == 'N')
        return clAmdBlasNoTrans;
    if (type[0] == 't' || type[0] == 'T')
        return clAmdBlasTrans;
    if (type[0] == 'c' || type[0] == 'C')
        return clAmdBlasConjTrans;

    return clAmdBlasNoTrans;
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
i_type_wr chemm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, c_type_wr *alpha,
                              c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta,
                              c_type_wr *c__, i_type_wr *ldc)
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
i_type_wr cher2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha,
                               c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, s_type_wr *beta,
                               c_type_wr *c__, i_type_wr *ldc)
{
    return not_implemented();
};
i_type_wr cherk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha,
                              c_type_wr *a, i_type_wr *lda, s_type_wr *beta, c_type_wr *c__, i_type_wr *ldc)
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
i_type_wr csymm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, c_type_wr *alpha,
                              c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta,
                              c_type_wr *c__, i_type_wr *ldc)
{
    return not_implemented();
};
i_type_wr csyr2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha,
                               c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta,
                               c_type_wr *c__, i_type_wr *ldc)
{
    return not_implemented();
};
i_type_wr csyrk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha,
                              c_type_wr *a, i_type_wr *lda, c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc)
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
i_type_wr ctrmm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb)
{
    return not_implemented();
};
i_type_wr ctrmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *a, i_type_wr *lda,
                              c_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr ctrsm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb)
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
i_type_wr dsymm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, d_type_wr *alpha,
                              d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta,
                              d_type_wr *c__, i_type_wr *ldc)
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
i_type_wr dsyr2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha,
                               d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta,
                               d_type_wr *c__, i_type_wr *ldc)
{
    return not_implemented();
};
i_type_wr dsyrk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, d_type_wr *a,
                              i_type_wr *lda, d_type_wr *beta, d_type_wr *c__, i_type_wr *ldc)
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
i_type_wr dtrmm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb)
{
    return not_implemented();
};
i_type_wr dtrmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *a, i_type_wr *lda,
                              d_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr dtrsm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb)
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
i_type_wr ssymm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, s_type_wr *alpha,
                              s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *beta,
                              s_type_wr *c__, i_type_wr *ldc)
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
i_type_wr ssyr2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha,
                               s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *beta,
                               s_type_wr *c__, i_type_wr *ldc)
{
    return not_implemented();
};
i_type_wr ssyrk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha,
                              s_type_wr *a, i_type_wr *lda, s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc)
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
i_type_wr strmm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb)
{
    return not_implemented();
};
i_type_wr strmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *a, i_type_wr *lda,
                              s_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr strsm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb)
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
i_type_wr zhemm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, z_type_wr *alpha,
                              z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta,
                              z_type_wr *c__, i_type_wr *ldc)
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
i_type_wr zher2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha,
                               z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, d_type_wr *beta,
                               z_type_wr *c__, i_type_wr *ldc)
{
    return not_implemented();
};
i_type_wr zherk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha,
                              z_type_wr *a, i_type_wr *lda, d_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
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
i_type_wr zsymm_gpu_wrap_impl(char *side, char *uplo, i_type_wr *m, i_type_wr *n, z_type_wr *alpha,
                              z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta,
                              z_type_wr *c__, i_type_wr *ldc)
{
    return not_implemented();
};
i_type_wr zsyr2k_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha,
                               z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta,
                               z_type_wr *c__, i_type_wr *ldc)
{
    return not_implemented();
};
i_type_wr zsyrk_gpu_wrap_impl(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha,
                              z_type_wr *a, i_type_wr *lda, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
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
i_type_wr ztrmm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb)
{
    return not_implemented();
};
i_type_wr ztrmv_gpu_wrap_impl(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *a, i_type_wr *lda,
                              z_type_wr *x, i_type_wr *incx)
{
    return not_implemented();
};
i_type_wr ztrsm_gpu_wrap_impl(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                              z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb)
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

static cl_float gpu_c(s_type_wr val)        { return val; }
static cl_double gpu_c(d_type_wr val)       { return val; }
static FloatComplex gpu_c(c_type_wr val)    { FloatComplex ret; ret.s[0] = val.r; ret.s[1] = val.i; return ret; }
static DoubleComplex gpu_c(z_type_wr val)   { DoubleComplex ret; ret.s[0] = val.r; ret.s[1] = val.i; return ret; }
static cl_float* gpu_c(s_type_wr* val)      { return val; }
static cl_double* gpu_c(d_type_wr* val)     { return val; }
static FloatComplex* gpu_c(c_type_wr* val)  { return reinterpret_cast<FloatComplex*>(val); };
static DoubleComplex* gpu_c(z_type_wr* val) { return reinterpret_cast<DoubleComplex*>(val); };

template<class T>   struct opencl_type              {};
template<>          struct opencl_type<s_type_wr>   { typedef cl_float type; };
template<>          struct opencl_type<d_type_wr>   { typedef cl_double type; };
template<>          struct opencl_type<c_type_wr>   { typedef FloatComplex type; };
template<>          struct opencl_type<z_type_wr>   { typedef DoubleComplex type; };

template<class V>
cl_int set_matrix(cl_context ctx, cl_command_queue queue, bool read_only, 
                  size_t M, size_t N, const V* a, size_t lda, cl_mem& AA, bool copy)
{
    cl_int err = 0;
    cl_mem_flags flags = (read_only) ? CL_MEM_READ_ONLY : CL_MEM_READ_WRITE;

    AA = clCreateBuffer(ctx, flags, M * N * sizeof(V), NULL, &err);

    if (err != CL_SUCCESS || AA == nullptr)
        return err;

    if (copy == false)
        return err;

    if (lda == M)
    {
        err = clEnqueueWriteBuffer(queue, AA, CL_TRUE, 0, M * N * sizeof(V), a, 0, nullptr, nullptr);
    }
    else
    {    
        for (size_t i = 0, offset = 0; i < N; ++i, offset += M * sizeof(V))
        {
            err = clEnqueueWriteBuffer(queue, AA, CL_TRUE, offset, M * sizeof(V), a, 0, nullptr, nullptr);
            a += lda;

            if (err != CL_SUCCESS)
                return err;
        }
    };
    return err;
};

template<class V>
cl_int get_matrix(cl_command_queue queue, size_t M, size_t N, V* c, size_t ld, cl_mem CC)
{
    if (ld == M)
    {
        cl_int err = clEnqueueReadBuffer(queue, CC, CL_TRUE, 0, M * N * sizeof(V), c, 0, NULL, NULL);
        return err;
    }
    else
    {
        for (size_t i = 0, offset = 0; i < N; ++i, offset += M*sizeof(V))
        {
            cl_int err = clEnqueueReadBuffer(queue, CC, CL_TRUE, offset, M * sizeof(V), c, 0, NULL, NULL);
            c += ld;

            if (err != CL_SUCCESS)
                return err;
        };

        return CL_SUCCESS;
    }
}

static void get_context(cl_context& ret_ctx, cl_command_queue& ret_queue)
{
    static cl_context st_ctx = nullptr;
    static cl_command_queue st_queue = nullptr;
    static bool initialized = false;

    if (initialized == true)
    {
        ret_ctx = st_ctx;
        ret_queue = st_queue;
        return;
    }

    static const int MAX_PLATFORM = 10;

    cl_int err;
    cl_platform_id platform[MAX_PLATFORM] = { 0 };
    cl_device_id device[MAX_PLATFORM] = { 0 };
    cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };
    cl_context ctx = 0;
    cl_command_queue queue = 0;
    int ret = 0;
    cl_uint num_platforms;
    cl_uint num_devices;

    /* Setup OpenCL environment. */
    err = clGetPlatformIDs(MAX_PLATFORM, platform, &num_platforms);
    if (err != CL_SUCCESS) {
        printf("clGetPlatformIDs() failed with %d\n", err);
        return;
    }

    /*
    for (size_t i = 0; i < num_platforms; ++i)
    {
    char info[100] = { 0 };
    size_t param_value_size_ret = 0;

    err = clGetPlatformInfo(platform[i], CL_PLATFORM_PROFILE, 100, info, &param_value_size_ret);
    err = clGetPlatformInfo(platform[i], CL_PLATFORM_VERSION, 100, info, &param_value_size_ret);
    err = clGetPlatformInfo(platform[i], CL_PLATFORM_NAME, 100, info, &param_value_size_ret);
    err = clGetPlatformInfo(platform[i], CL_PLATFORM_VENDOR, 100, info, &param_value_size_ret);
    err = clGetPlatformInfo(platform[i], CL_PLATFORM_EXTENSIONS, 100, info, &param_value_size_ret);

    err = clGetDeviceIDs(platform[i], CL_DEVICE_TYPE_GPU, MAX_PLATFORM, device, &num_devices);
    if (err != CL_SUCCESS) {
    printf("clGetDeviceIDs() failed with %d\n", err);
    }

    for (size_t j = 0; j < num_devices; ++j)
    {
    err = clGetDeviceInfo(device[j], CL_DEVICE_NAME, 100, info, &param_value_size_ret);
    err = clGetDeviceInfo(device[j], CL_DEVICE_VENDOR, 100, info, &param_value_size_ret);
    err = clGetDeviceInfo(device[j], CL_DRIVER_VERSION, 100, info, &param_value_size_ret);
    err = clGetDeviceInfo(device[j], CL_DEVICE_PROFILE, 100, info, &param_value_size_ret);
    err = clGetDeviceInfo(device[j], CL_DEVICE_VERSION, 100, info, &param_value_size_ret);
    err = clGetDeviceInfo(device[j], CL_DEVICE_EXTENSIONS, 100, info, &param_value_size_ret);

    if (err != CL_SUCCESS) {
    printf("clGetDeviceIDs() failed with %d\n", err);
    }
    };
    };
    */

    err = clGetDeviceIDs(platform[1], CL_DEVICE_TYPE_GPU, MAX_PLATFORM, device, &num_devices);
    if (err != CL_SUCCESS) {
        printf("clGetDeviceIDs() failed with %d\n", err);
        return;
    }

    props[1] = (cl_context_properties)platform[1];
    ctx = clCreateContext(props, 1, device, NULL, NULL, &err);
    if (err != CL_SUCCESS) {
        printf("clCreateContext() failed with %d\n", err);
        return;
    }

    queue = clCreateCommandQueue(ctx, device[0], 0, &err);
    if (err != CL_SUCCESS) {
        printf("clCreateCommandQueue() failed with %d\n", err);
        clReleaseContext(ctx);
        return;
    }

    /* Setup clAmdBlas. */
    err = clAmdBlasSetup();
    if (err != CL_SUCCESS) {
        printf("clAmdBlasSetup() failed with %d\n", err);
        clReleaseCommandQueue(queue);
        clReleaseContext(ctx);
        return;
    }

    st_ctx = ctx;
    st_queue = queue;
    initialized = true;

    ret_ctx = ctx;
    ret_queue = queue;
    return;
};
//--------------------------------------------------------------------------
//                      IMPLEMENTED
//--------------------------------------------------------------------------
template<class V>
i_type_wr dgemm_gpu_generic(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                            V *alpha, V *a, i_type_wr *lda, V *b, i_type_wr *ldb,
                            V *beta, V *c__, i_type_wr *ldc)
{

    typedef opencl_type<V>::type value_type;

    int A_M = 0;
    int A_N = 0;
    int B_M = 0;
    int B_N = 0;
    int C_M = *m;
    int C_N = *n;

    if (op(transa) == clAmdBlasNoTrans)
    {
        A_M = *m;
        A_N = *k;
    }
    else
    {
        A_N = *m;
        A_M = *k;
    }
    if (op(transb) == clAmdBlasNoTrans)
    {
        B_M = *k;
        B_N = *n;
    }
    else
    {
        B_N = *k;
        B_M = *n;
    }

    cl_mem AA   = nullptr;
    cl_mem BB   = nullptr;
    cl_mem CC   = nullptr;

    cl_command_queue queue = nullptr;
    cl_context ctx = nullptr;
    cl_event event = nullptr;
    get_context(ctx, queue);

    cl_int err;

    err = set_matrix<V>(ctx, queue, true, A_M, A_N, a, *lda, AA, true);
    if (err != CL_SUCCESS || AA == nullptr)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    };

    err = set_matrix<V>(ctx, queue, true, B_M, B_N, b, *ldb, BB, true);
    if (err != CL_SUCCESS || BB == nullptr)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    };

    err = set_matrix<V>(ctx, queue, false, C_M, C_N, c__, *ldc, CC, is_zero(*beta) == false);
    if (err != CL_SUCCESS || CC == nullptr)
    {
        std::cerr << "allocation error" << "\n";
        goto lab_error;
    };

    // Call clAmdBlas function
    err = clAmdBlasDgemm(clAmdBlasColumnMajor, op(transa), op(transb), *m, *n, *k, gpu_c(*alpha), AA, A_M,
                         BB, B_M, *gpu_c(beta), CC, C_M, 1, &queue, 0, NULL, &event);

    if (err != CL_SUCCESS)
    {
        std::cerr << "clblas kernel error" << "\n";
        goto lab_error;
    };

    err = clWaitForEvents(1, &event);
    if (err != CL_SUCCESS)
    {
        std::cerr << "clblas kernel error" << "\n";
        goto lab_error;
    };

    err = get_matrix<V>(queue, C_M, C_N, c__, *ldc, CC);
    if (err != CL_SUCCESS)
    {
        goto lab_error;
    };

    clReleaseMemObject(AA);
    clReleaseMemObject(BB);
    clReleaseMemObject(CC);

    return 0;

  lab_error:

    clReleaseMemObject(AA);
    clReleaseMemObject(BB);
    clReleaseMemObject(CC);

    return 1;
};

i_type_wr dgemm_gpu_wrap_impl(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                              d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb,
                              d_type_wr *beta, d_type_wr *c__, i_type_wr *ldc)
{
    return dgemm_gpu_generic<d_type_wr>(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
};
i_type_wr cgemm_gpu_wrap_impl(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                              c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb,
                              c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc)
{
    //return dgemm_gpu_generic<c_type_wr>(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    return 0;
};
i_type_wr sgemm_gpu_wrap_impl(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                              s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb,
                              s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc)
{
    //return dgemm_gpu_generic<s_type_wr>(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    return 0;
};
i_type_wr zgemm_gpu_wrap_impl(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                              z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb,
                              z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
{
    //return dgemm_gpu_generic<z_type_wr>(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    return 0;
};


// generic stuff - the constructor which initializes all pointers to blas functions
#define CALL_SYNTAX(x) x##_gpu_wrap_impl
#include "blas_plugin_common.h"
