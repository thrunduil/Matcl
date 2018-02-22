/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Paweł Kowal 2017 - 2018
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

#include <stdexcept>
#include <exception>
#include <string>

#ifndef __unix__
    #include <Windows.h>
#endif

#include "matcl-blas-lapack/blas_loader/blas_loader.h"
#include "blas_config.h"
#include <iostream>
#include <boost/smart_ptr/detail/spinlock.hpp>
#include <sstream>
#include <thread>

#ifndef __unix__ //supposedly WINDOWS version

    //TODO: get correct plugin names 
    #define BLAS_LIBRARY_MKL                "blas_mkl_plugin.dll"
    #define BLAS_LIBRARY_CLAPACK            "blas_clapack_plugin.dll"
    #define BLAS_LIBRARY_CPU_GPU_DEFAULT    "blas_plugin_cpu_gpu_default.dll"

    namespace matcl { namespace lapack
    {
        void init_env()
        {
            //it seems to be impossible to set MKL threads using omp_set_num_threads, we need
            //to set environmental variable 

            int n_threads = std::thread::hardware_concurrency();

            //TODO:
            n_threads = 1;

            std::ostringstream os;
            os << n_threads;

            std::string n_str   = os.str();

            const char lpName[] = "OMP_NUM_THREADS";
            const char* lpValue = n_str.c_str();

            SetEnvironmentVariableA(lpName, lpValue);
        };

        template<class plugins_list>
        static bool get_plugin(const plugins_list& plugins, ::blas_plugin*& plugin)
        {
            plugin          = nullptr;
            HINSTANCE mod   = nullptr;

            for (auto pos = plugins.begin(); pos != plugins.end(); ++pos)
            {
                const auto& loc_plugin = *pos;
                mod = LoadLibraryA(loc_plugin.c_str());

                if (mod != nullptr)
                    break;
            };

            if (!mod)
                return true;

            using plugin_type = blas_plugin::get_blas_plugin_type;

            plugin_type* f = (plugin_type*)GetProcAddress(mod, "get_blas_plugin");
            ::blas_plugin* ret = nullptr;

            if (f)
                ret = f();

            if (!ret)
            {
                std::cerr << "unable to init blas wrapper from blas library" << "\n";
                return false;
            }

            plugin = ret;
            return true;
        }

        template<class plugins_list>
        static bool get_plugin_mix(const plugins_list& plugins, ::blas_plugin*& plugin, 
                                    const ::blas_plugin* plugin_cpu, const ::blas_plugin* plugin_gpu)
        {
            plugin = nullptr;
            HINSTANCE mod = nullptr;

            for (auto pos = plugins.begin(); pos != plugins.end(); ++pos)
            {
                const auto& loc_plugin = *pos;
                mod = LoadLibraryA(loc_plugin.c_str());

                if (mod != nullptr)
                    break;
            };

            if (!mod)
                return true;

            using plugin_type = blas_plugin::get_blas_plugin_gpu_cpu_type;

            plugin_type* f = (plugin_type*)GetProcAddress(mod, "get_blas_plugin_gpu_cpu");
            ::blas_plugin* ret = nullptr;

            if (f)
                ret = f(plugin_cpu, plugin_gpu);

            if (!ret)
            {
                std::cerr << "unable to init blas wrapper from blas library" << "\n";
                return false;
            }

            plugin = ret;
            return true;
        }

        static ::blas_plugin* get_plugin()
        {
            init_env();

            blas_config config;

            config.add_plugin_cpu(BLAS_LIBRARY_MKL);
            config.add_plugin_cpu(BLAS_LIBRARY_CLAPACK);
            config.add_plugin_cpu_gpu(BLAS_LIBRARY_CPU_GPU_DEFAULT);

            ::blas_plugin* plugin_cpu = nullptr;
            ::blas_plugin* plugin_gpu = nullptr;

            bool ok_cpu = get_plugin(config.get_plugins_cpu(), plugin_cpu);
            bool ok_gpu = get_plugin(config.get_plugins_gpu(), plugin_gpu);

            if (ok_cpu == true && plugin_cpu == nullptr)
            {
                std::cerr << "Unable to load blas library plugin" << "\n";
                return nullptr;
            }

            if (ok_gpu == true && plugin_gpu != nullptr)
            {
                ::blas_plugin* plugin_cpu_gpu = nullptr;

                bool ok = get_plugin_mix(config.get_plugins_cpu_gpu(), plugin_cpu_gpu, 
                                         plugin_cpu, plugin_gpu);

                if (ok == true && plugin_cpu_gpu != nullptr)
                    return plugin_cpu_gpu;
                else
                    return plugin_cpu;
            }
            else
            {
                return plugin_cpu;
            };
        }

        static ::blas_plugin* m_plugin = get_plugin();

        static void init_plugin() {};
    }}

#else //supposedly LINUX version


#include <sys/types.h>
#include <dlfcn.h>

#define BLAS_LIBRARY_NAME "libblas_clapack_plugin.so"

    namespace matcl { namespace lapack
    {
        static ::blas_plugin* get_plugin()
        {
            typedef blas_plugin::get_blas_plugin plugin_type;
    
            void *mod = dlopen(BLAS_LIBRARY_NAME, RTLD_LAZY);

            if (!mod)
            {
                std::cerr << "Unable to load blas library plugin" << "\n";
                return nullptr;
            }

            plugin_type* f         = (plugin_type*)dlsym(mod, "get_blas_plugin");
            ::blas_plugin* ret     = nullptr;

	        if (f)
    		    ret = f();

	        if (!ret)
	        {
                std::cerr << "Unable to init blas wrapper from blas library" << "\n";
                return nullptr;
	        }

            return ret;
        };

        static ::blas_plugin* m_plugin = get_plugin();

        static void init_plugin(){};
    }}

#endif

namespace raw_blas_lapack
{

extern "C" 
{    
    i_type_wr get_num_threads_blas_kernel()
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->get_num_threads_fptr();
    };

    i_type_wr get_default_threads_blas_kernel()
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->get_default_threads_fptr();
    };
    
    void set_num_threads_blas_kernel(i_type_wr n)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->set_num_threads_fptr(&n);
    };
    
    bool are_user_threads_allowed()
    {
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->are_user_threads_allowed_fptr();
    };

    i_type_wr caxpy_(i_type_wr *n, c_type_wr *ca, c_type_wr *cx, i_type_wr *incx,
                     c_type_wr *cy, i_type_wr *incy)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->caxpy_fptr(n, ca, cx, incx, cy, incy);
    }

    i_type_wr ccopy_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, 
                     i_type_wr *incy)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ccopy_fptr(n, cx, incx, cy, incy);
    }

    i_type_wr cdotc_(c_type_wr *ret_val, i_type_wr *n, c_type_wr *cx, i_type_wr *incx, 
                     c_type_wr *cy, i_type_wr *incy)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->cdotc_fptr(ret_val, n, cx, incx, cy, incy);
    }

    i_type_wr cdotu_(c_type_wr *ret_val, i_type_wr *n, c_type_wr *cx, i_type_wr *incx, 
                     c_type_wr *cy, i_type_wr *incy)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->cdotu_fptr(ret_val, n, cx, incx, cy, incy);
    }

    i_type_wr cgbmv_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku,
                     c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx, 
                     c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->cgbmv_fptr(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr cgemm_(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, 
                     c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta, 
                     c_type_wr *c__, i_type_wr *ldc)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->cgemm_fptr(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr cgemv_(char *trans, i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, 
                     c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->cgemv_fptr(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr cgerc_(i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx, c_type_wr *y,
                     i_type_wr *incy, c_type_wr *a, i_type_wr *lda)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->cgerc_fptr(m, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr cgeru_(i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                     c_type_wr *y, i_type_wr *incy, c_type_wr *a, i_type_wr *lda)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->cgeru_fptr(m, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr chbmv_(char *uplo, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda,
                     c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->chbmv_fptr(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr chemm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *a, 
                     i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta, c_type_wr *c__, 
                     i_type_wr *ldc)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->chemm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr chemv_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *x, 
                     i_type_wr *incx, c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->chemv_fptr(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr cher_(char *uplo, i_type_wr *n, s_type_wr *alpha, c_type_wr *x, i_type_wr *incx, c_type_wr *a,
                    i_type_wr *lda)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->cher_fptr(uplo, n, alpha, x, incx, a, lda);
    }

    i_type_wr cher2_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx, c_type_wr *y, 
                     i_type_wr *incy, c_type_wr *a, i_type_wr *lda)
    {        
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->cher2_fptr(uplo, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr cher2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, c_type_wr *a, 
                      i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, s_type_wr *beta, c_type_wr *c__, 
                      i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->cher2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr cherk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, c_type_wr *a, 
                     i_type_wr *lda, s_type_wr *beta, c_type_wr *c__, i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->cherk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
    }

    i_type_wr chpmv_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *ap, c_type_wr *x, i_type_wr *incx,
                     c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->chpmv_fptr(uplo, n, alpha, ap, x, incx, beta, y, incy);
    }

    i_type_wr chpr_(char *uplo, i_type_wr *n, s_type_wr *alpha, c_type_wr *x, i_type_wr *incx, c_type_wr *ap)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->chpr_fptr(uplo, n, alpha, x, incx, ap);
    }
    
    i_type_wr chpr2_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx, c_type_wr *y, 
                     i_type_wr *incy, c_type_wr *ap)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->chpr2_fptr(uplo, n, alpha, x, incx, y, incy, ap);
    }

    i_type_wr crotg_(c_type_wr *ca, c_type_wr *cb, s_type_wr *c__, c_type_wr *s)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->crotg_fptr(ca, cb, c__, s);
    }
    
    i_type_wr cscal_(i_type_wr *n, c_type_wr *ca, c_type_wr *cx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->cscal_fptr(n, ca, cx, incx);
    }
    
    i_type_wr csrot_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, i_type_wr *incy, 
                     s_type_wr *c__, s_type_wr *s)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->csrot_fptr(n, cx, incx, cy, incy, c__, s);
    }

    i_type_wr csscal_(i_type_wr *n, s_type_wr *sa, c_type_wr *cx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->csscal_fptr(n, sa, cx, incx);
    }
    
    i_type_wr cswap_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->cswap_fptr(n, cx, incx, cy, incy);
    }
    
    i_type_wr csymm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *a, 
                     i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta, c_type_wr *c__, 
                     i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->csymm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr csyr2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, c_type_wr *a, 
                      i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta, c_type_wr *c__,
                      i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->csyr2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr csyrk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, c_type_wr *a, 
                     i_type_wr *lda, c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->csyrk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
    }

    i_type_wr ctbmv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, c_type_wr *a, 
                     i_type_wr *lda, c_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ctbmv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr ctbsv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, c_type_wr *a, 
                     i_type_wr *lda, c_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ctbsv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr ctpmv_(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *ap, c_type_wr *x,
                     i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ctpmv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr ctpsv_(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *ap, c_type_wr *x,
                     i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ctpsv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr ctrmm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n, 
                     c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ctrmm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr ctrmv_(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *a, i_type_wr *lda,
                     c_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ctrmv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    i_type_wr ctrsm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n, 
                     c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ctrsm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr ctrsv_(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
                     c_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ctrsv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    d_type_wr dasum_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dasum_fptr(n, dx, incx);
    }
    
    i_type_wr daxpy_(i_type_wr *n, d_type_wr *da, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->daxpy_fptr(n, da, dx, incx, dy, incy);
    }
    
    d_type_wr dcabs1_(z_type_wr *z__)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dcabs1_fptr(z__);
    }
    
    i_type_wr dcopy_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dcopy_fptr(n, dx, incx, dy, incy);
    }
    
    d_type_wr ddot_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ddot_fptr(n, dx, incx, dy, incy);
    }
    
    i_type_wr dgbmv_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
                     d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *x, i_type_wr *incx,
                     d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dgbmv_fptr(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr dgemm_(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k, d_type_wr *alpha,
                     d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta, 
                     d_type_wr *c__, i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dgemm_fptr(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr dgemv_(char *trans, i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, 
                     d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dgemv_fptr(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr dger_(i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx, d_type_wr *y, 
                    i_type_wr *incy, d_type_wr *a, i_type_wr *lda)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dger_fptr(m, n, alpha, x, incx, y, incy, a, lda);
    }

    d_type_wr dnrm2_(i_type_wr *n, d_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dnrm2_fptr(n, x, incx);
    }
    
    i_type_wr drot_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy,
                    d_type_wr *c__, d_type_wr *s)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->drot_fptr(n, dx, incx, dy, incy, c__, s);
    }

    i_type_wr drotg_(d_type_wr *da, d_type_wr *db, d_type_wr *c__, d_type_wr *s)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->drotg_fptr(da, db, c__, s);
    }
    
    i_type_wr drotm_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy, 
                     d_type_wr *dparam)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->drotm_fptr(n, dx, incx, dy, incy, dparam);
    }

    i_type_wr drotmg_(d_type_wr *dd1, d_type_wr *dd2, d_type_wr *dx1, d_type_wr *dy1, d_type_wr *dparam)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->drotmg_fptr(dd1, dd2, dx1, dy1, dparam);
    }
    
    i_type_wr dsbmv_(char *uplo, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda,
                     d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dsbmv_fptr(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr dscal_(i_type_wr *n, d_type_wr *da, d_type_wr *dx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dscal_fptr(n, da, dx, incx);
    }
    
    d_type_wr dsdot_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
    {
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dsdot_fptr(n, sx, incx, sy, incy);
    }
    
    i_type_wr dspmv_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *ap, d_type_wr *x, 
                     i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dspmv_fptr(uplo, n, alpha, ap, x, incx, beta, y, incy);
    }

    i_type_wr dspr_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx, d_type_wr *ap)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dspr_fptr(uplo, n, alpha, x, incx, ap);
    }

    i_type_wr dspr2_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx, d_type_wr *y, 
                     i_type_wr *incy, d_type_wr *ap)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dspr2_fptr(uplo, n, alpha, x, incx, y, incy, ap);
    }

    i_type_wr dswap_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dswap_fptr(n, dx, incx, dy, incy);
    }
    
    i_type_wr dsymm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *a, 
                     i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta, d_type_wr *c__, 
                     i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dsymm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr dsymv_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *x, 
                     i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dsymv_fptr(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr dsyr_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx, d_type_wr *a,
                    i_type_wr *lda)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dsyr_fptr(uplo, n, alpha, x, incx, a, lda);
    }

    i_type_wr dsyr2_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx, d_type_wr *y, 
                     i_type_wr *incy, d_type_wr *a, i_type_wr *lda)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dsyr2_fptr(uplo, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr dsyr2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, d_type_wr *a,
                      i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta, d_type_wr *c__, 
                      i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dsyr2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr dsyrk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, d_type_wr *a,
                     i_type_wr *lda, d_type_wr *beta, d_type_wr *c__, i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dsyrk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
    }

    i_type_wr dtbmv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, d_type_wr *a, 
                     i_type_wr *lda, d_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dtbmv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr dtbsv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, d_type_wr *a,
                     i_type_wr *lda, d_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dtbsv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr dtpmv_(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *ap, d_type_wr *x,
                     i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dtpmv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr dtpsv_(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *ap, d_type_wr *x, 
                     i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dtpsv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr dtrmm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n, 
                     d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dtrmm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr dtrmv_(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *a, i_type_wr *lda, 
                     d_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dtrmv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    i_type_wr dtrsm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n, 
                     d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dtrsm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr dtrsv_(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *a, i_type_wr *lda, 
                     d_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dtrsv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    d_type_wr dzasum_(i_type_wr *n, z_type_wr *zx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dzasum_fptr(n, zx, incx);
    }
    
    d_type_wr dznrm2_(i_type_wr *n, z_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->dznrm2_fptr(n, x, incx);
    }
    
    i_type_wr icamax_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->icamax_fptr(n, cx, incx);
    }
    
    i_type_wr idamax_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->idamax_fptr(n, dx, incx);
    }
    
    i_type_wr isamax_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->isamax_fptr(n, sx, incx);
    }
    
    i_type_wr izamax_(i_type_wr *n, z_type_wr *zx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->izamax_fptr(n, zx, incx);
    }

    d_type_wr sasum_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->sasum_fptr(n, sx, incx);
    }

    i_type_wr saxpy_(i_type_wr *n, s_type_wr *sa, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, 
                     i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->saxpy_fptr(n, sa, sx, incx, sy, incy);
    }

    d_type_wr scabs1_(c_type_wr *z__)
    {
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->scabs1_fptr(z__);
    }
    
    d_type_wr scasum_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->scasum_fptr(n, cx, incx);
    }
    
    d_type_wr scnrm2_(i_type_wr *n, c_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->scnrm2_fptr(n, x, incx);
    }
    
    i_type_wr scopy_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->scopy_fptr(n, sx, incx, sy, incy);
    }
    
    s_type_wr sdot_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->sdot_fptr(n, sx, incx, sy, incy);
    }
    
    d_type_wr sdsdot_(i_type_wr *n, s_type_wr *sb, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
    {
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->sdsdot_fptr(n, sb, sx, incx, sy, incy);
    }
    
    i_type_wr sgbmv_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, s_type_wr *alpha, 
                     s_type_wr *a, i_type_wr *lda, s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y,
                     i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->sgbmv_fptr(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr sgemm_(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, 
                     s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *beta, s_type_wr *c__,
                     i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->sgemm_fptr(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr sgemv_(char *trans, i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, 
                     s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->sgemv_fptr(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr sger_(i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx, s_type_wr *y, 
                    i_type_wr *incy, s_type_wr *a, i_type_wr *lda)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->sger_fptr(m, n, alpha, x, incx, y, incy, a, lda);
    }

    d_type_wr snrm2_(i_type_wr *n, s_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->snrm2_fptr(n, x, incx);
    }
    
    i_type_wr srot_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy, s_type_wr *c__,
                    s_type_wr *s)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->srot_fptr(n, sx, incx, sy, incy, c__, s);
    }

    i_type_wr srotg_(s_type_wr *sa, s_type_wr *sb, s_type_wr *c__, s_type_wr *s)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->srotg_fptr(sa, sb, c__, s);
    }
    
    i_type_wr srotm_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy, 
                     s_type_wr *sparam)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->srotm_fptr(n, sx, incx, sy, incy, sparam);
    }

    i_type_wr srotmg_(s_type_wr *sd1, s_type_wr *sd2, s_type_wr *sx1, s_type_wr *sy1, s_type_wr *sparam)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->srotmg_fptr(sd1, sd2, sx1, sy1, sparam);
    }
    
    i_type_wr ssbmv_(char *uplo, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, 
                     s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ssbmv_fptr(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr sscal_(i_type_wr *n, s_type_wr *sa, s_type_wr *sx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->sscal_fptr(n, sa, sx, incx);
    }
    
    i_type_wr sspmv_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *ap, s_type_wr *x, i_type_wr *incx, 
                     s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->sspmv_fptr(uplo, n, alpha, ap, x, incx, beta, y, incy);
    }

    i_type_wr sspr_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx, s_type_wr *ap)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->sspr_fptr(uplo, n, alpha, x, incx, ap);
    }
    
    i_type_wr sspr2_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx, s_type_wr *y, 
                     i_type_wr *incy, s_type_wr *ap)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->sspr2_fptr(uplo, n, alpha, x, incx, y, incy, ap);
    }

    i_type_wr sswap_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->sswap_fptr(n, sx, incx, sy, incy);
    }
    
    i_type_wr ssymm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *a, 
                     i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ssymm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr ssymv_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *x, 
                     i_type_wr *incx, s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ssymv_fptr(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr ssyr_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx, s_type_wr *a,
                    i_type_wr *lda)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ssyr_fptr(uplo, n, alpha, x, incx, a, lda);
    }

    i_type_wr ssyr2_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx, s_type_wr *y,
                     i_type_wr *incy, s_type_wr *a, i_type_wr *lda)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ssyr2_fptr(uplo, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr ssyr2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, s_type_wr *a,
                      i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *beta, s_type_wr *c__, 
                      i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ssyr2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr ssyrk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, s_type_wr *a, 
                     i_type_wr *lda, s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ssyrk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
    }
    
    i_type_wr stbmv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, s_type_wr *a, 
                     i_type_wr *lda, s_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->stbmv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr stbsv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, s_type_wr *a,
                     i_type_wr *lda, s_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->stbsv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr stpmv_(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *ap, s_type_wr *x, 
                     i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->stpmv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr stpsv_(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *ap, s_type_wr *x, 
                     i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->stpsv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr strmm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                     s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->strmm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr strmv_(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *a, i_type_wr *lda,
                     s_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->strmv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    i_type_wr strsm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                     s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->strsm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr strsv_(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
                     s_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->strsv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    i_type_wr zaxpy_(i_type_wr *n, z_type_wr *za, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zaxpy_fptr(n, za, zx, incx, zy, incy);
    }

    i_type_wr zcopy_(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zcopy_fptr(n, zx, incx, zy, incy);
    }
    
    i_type_wr zdotc_(z_type_wr *ret_val, i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, 
                     i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zdotc_fptr(ret_val, n, zx, incx, zy, incy);
    }

    i_type_wr zdotu_(z_type_wr *ret_val, i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, 
                     i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zdotu_fptr(ret_val, n, zx, incx, zy, incy);
    }

    i_type_wr zdrot_(i_type_wr *n, z_type_wr *cx, i_type_wr *incx, z_type_wr *cy, i_type_wr *incy, 
                     d_type_wr *c__, d_type_wr *s)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zdrot_fptr(n, cx, incx, cy, incy, c__, s);
    }

    i_type_wr zdscal_(i_type_wr *n, d_type_wr *da, z_type_wr *zx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zdscal_fptr(n, da, zx, incx);
    }
    
    i_type_wr zgbmv_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, z_type_wr *alpha,
                     z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y,
                     i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zgbmv_fptr(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr zgemm_(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, 
                     z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta, 
                     z_type_wr *c__, i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zgemm_fptr(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr zgemv_(char *trans, i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda,
                     z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zgemv_fptr(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr zgerc_(i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx, z_type_wr *y,
                     i_type_wr *incy, z_type_wr *a, i_type_wr *lda)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zgerc_fptr(m, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr zgeru_(i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx, z_type_wr *y, 
                     i_type_wr *incy, z_type_wr *a, i_type_wr *lda)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zgeru_fptr(m, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr zhbmv_(char *uplo, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, 
                     z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zhbmv_fptr(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr zhemm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *a, 
                     i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zhemm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr zhemv_(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *x,
                     i_type_wr *incx, z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zhemv_fptr(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr zher_(char *uplo, i_type_wr *n, d_type_wr *alpha, z_type_wr *x, i_type_wr *incx, z_type_wr *a,
                    i_type_wr *lda)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zher_fptr(uplo, n, alpha, x, incx, a, lda);
    }

    i_type_wr zher2_(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx, z_type_wr *y, 
                     i_type_wr *incy, z_type_wr *a, i_type_wr *lda)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zher2_fptr(uplo, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr zher2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, z_type_wr *a, 
                      i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, d_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zher2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr zherk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, z_type_wr *a,
                     i_type_wr *lda, d_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zherk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
    }

    i_type_wr zhpmv_(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *ap, z_type_wr *x, i_type_wr *incx,
                     z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zhpmv_fptr(uplo, n, alpha, ap, x, incx, beta, y, incy);
    }

    i_type_wr zhpr_(char *uplo, i_type_wr *n, d_type_wr *alpha, z_type_wr *x, i_type_wr *incx, z_type_wr *ap)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zhpr_fptr(uplo, n, alpha, x, incx, ap);
    }
    
    i_type_wr zhpr2_(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx, z_type_wr *y,
                     i_type_wr *incy, z_type_wr *ap)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zhpr2_fptr(uplo, n, alpha, x, incx, y, incy, ap);
    }

    i_type_wr zrotg_(z_type_wr *ca, z_type_wr *cb, d_type_wr *c__, z_type_wr *s)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zrotg_fptr(ca, cb, c__, s);
    }
    
    i_type_wr zscal_(i_type_wr *n, z_type_wr *za, z_type_wr *zx, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zscal_fptr(n, za, zx, incx);
    }
    
    i_type_wr zswap_(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, i_type_wr *incy)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zswap_fptr(n, zx, incx, zy, incy);
    }
    
    i_type_wr zsymm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *a,
                     i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zsymm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr zsyr2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, z_type_wr *a, 
                      i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zsyr2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr zsyrk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, z_type_wr *a, 
                     i_type_wr *lda, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->zsyrk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
    }

    i_type_wr ztbmv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, z_type_wr *a, 
                     i_type_wr *lda, z_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ztbmv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr ztbsv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, z_type_wr *a, 
                     i_type_wr *lda, z_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ztbsv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr ztpmv_(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *ap, z_type_wr *x, 
                     i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ztpmv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr ztpsv_(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *ap, z_type_wr *x, 
                     i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ztpsv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr ztrmm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n, 
                     z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ztrmm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr ztrmv_(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *a, i_type_wr *lda, 
                     z_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ztrmv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    i_type_wr ztrsm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n, 
                     z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ztrsm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr ztrsv_(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *a, i_type_wr *lda, 
                     z_type_wr *x, i_type_wr *incx)
    {   
        matcl::lapack::init_plugin();
        return matcl::lapack::m_plugin->ztrsv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }
}

}
