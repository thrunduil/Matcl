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

#include "blas_config.h"
#include "matcl-blas-lapack/blas_loader/blas_loader.h"
#include "plugin_manager.h"

#include <string>

namespace ml = matcl::lapack;

std::string raw_blas_lapack::loaded_blas_plugin_name()
{
    ml::plugin_manager::init_plugin();
    const char* name = ml::plugin_manager::plugin()->get_name_fptr();

    if (name == nullptr)
        return std::string();
    else
        return std::string(name);
}

bool raw_blas_lapack::load_blas_plugin(const std::string& path)
{
    return ml::plugin_manager::load_plugin(path);
}

void raw_blas_lapack::initialize_plugin()
{
    return ml::plugin_manager::initialize_plugin();
}

namespace raw_blas_lapack
{

extern "C" 
{    
    i_type_wr get_num_threads_blas_kernel()
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->get_num_threads_fptr();
    };

    i_type_wr get_default_threads_blas_kernel()
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->get_default_threads_fptr();
    };
    
    void set_num_threads_blas_kernel(i_type_wr n)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->set_num_threads_fptr(&n);
    };
    
    bool are_user_threads_allowed()
    {
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->are_user_threads_allowed_fptr();
    };

    i_type_wr caxpy_(i_type_wr *n, c_type_wr *ca, c_type_wr *cx, i_type_wr *incx,
                     c_type_wr *cy, i_type_wr *incy)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->caxpy_fptr(n, ca, cx, incx, cy, incy);
    }

    i_type_wr ccopy_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, 
                     i_type_wr *incy)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ccopy_fptr(n, cx, incx, cy, incy);
    }

    i_type_wr cdotc_(c_type_wr *ret_val, i_type_wr *n, c_type_wr *cx, i_type_wr *incx, 
                     c_type_wr *cy, i_type_wr *incy)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->cdotc_fptr(ret_val, n, cx, incx, cy, incy);
    }

    i_type_wr cdotu_(c_type_wr *ret_val, i_type_wr *n, c_type_wr *cx, i_type_wr *incx, 
                     c_type_wr *cy, i_type_wr *incy)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->cdotu_fptr(ret_val, n, cx, incx, cy, incy);
    }

    i_type_wr cgbmv_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku,
                     c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx, 
                     c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->cgbmv_fptr(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr cgemm_(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, 
                     c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta, 
                     c_type_wr *c__, i_type_wr *ldc)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->cgemm_fptr(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr cgemv_(char *trans, i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, 
                     c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->cgemv_fptr(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr cgerc_(i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx, c_type_wr *y,
                     i_type_wr *incy, c_type_wr *a, i_type_wr *lda)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->cgerc_fptr(m, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr cgeru_(i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx,
                     c_type_wr *y, i_type_wr *incy, c_type_wr *a, i_type_wr *lda)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->cgeru_fptr(m, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr chbmv_(char *uplo, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda,
                     c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->chbmv_fptr(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr chemm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *a, 
                     i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta, c_type_wr *c__, 
                     i_type_wr *ldc)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->chemm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr chemv_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *x, 
                     i_type_wr *incx, c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->chemv_fptr(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr cher_(char *uplo, i_type_wr *n, s_type_wr *alpha, c_type_wr *x, i_type_wr *incx, c_type_wr *a,
                    i_type_wr *lda)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->cher_fptr(uplo, n, alpha, x, incx, a, lda);
    }

    i_type_wr cher2_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx, c_type_wr *y, 
                     i_type_wr *incy, c_type_wr *a, i_type_wr *lda)
    {        
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->cher2_fptr(uplo, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr cher2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, c_type_wr *a, 
                      i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, s_type_wr *beta, c_type_wr *c__, 
                      i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->cher2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr cherk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, c_type_wr *a, 
                     i_type_wr *lda, s_type_wr *beta, c_type_wr *c__, i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->cherk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
    }

    i_type_wr chpmv_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *ap, c_type_wr *x, i_type_wr *incx,
                     c_type_wr *beta, c_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->chpmv_fptr(uplo, n, alpha, ap, x, incx, beta, y, incy);
    }

    i_type_wr chpr_(char *uplo, i_type_wr *n, s_type_wr *alpha, c_type_wr *x, i_type_wr *incx, c_type_wr *ap)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->chpr_fptr(uplo, n, alpha, x, incx, ap);
    }
    
    i_type_wr chpr2_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *incx, c_type_wr *y, 
                     i_type_wr *incy, c_type_wr *ap)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->chpr2_fptr(uplo, n, alpha, x, incx, y, incy, ap);
    }

    i_type_wr crotg_(c_type_wr *ca, c_type_wr *cb, s_type_wr *c__, c_type_wr *s)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->crotg_fptr(ca, cb, c__, s);
    }
    
    i_type_wr cscal_(i_type_wr *n, c_type_wr *ca, c_type_wr *cx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->cscal_fptr(n, ca, cx, incx);
    }
    
    i_type_wr csrot_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, i_type_wr *incy, 
                     s_type_wr *c__, s_type_wr *s)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->csrot_fptr(n, cx, incx, cy, incy, c__, s);
    }

    i_type_wr csscal_(i_type_wr *n, s_type_wr *sa, c_type_wr *cx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->csscal_fptr(n, sa, cx, incx);
    }
    
    i_type_wr cswap_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->cswap_fptr(n, cx, incx, cy, incy);
    }
    
    i_type_wr csymm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *a, 
                     i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta, c_type_wr *c__, 
                     i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->csymm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr csyr2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, c_type_wr *a, 
                      i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *beta, c_type_wr *c__,
                      i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->csyr2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr csyrk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, c_type_wr *a, 
                     i_type_wr *lda, c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->csyrk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
    }

    i_type_wr ctbmv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, c_type_wr *a, 
                     i_type_wr *lda, c_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ctbmv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr ctbsv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, c_type_wr *a, 
                     i_type_wr *lda, c_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ctbsv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr ctpmv_(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *ap, c_type_wr *x,
                     i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ctpmv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr ctpsv_(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *ap, c_type_wr *x,
                     i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ctpsv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr ctrmm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n, 
                     c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ctrmm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr ctrmv_(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *a, i_type_wr *lda,
                     c_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ctrmv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    i_type_wr ctrsm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n, 
                     c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ctrsm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr ctrsv_(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
                     c_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ctrsv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    d_type_wr dasum_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dasum_fptr(n, dx, incx);
    }
    
    i_type_wr daxpy_(i_type_wr *n, d_type_wr *da, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->daxpy_fptr(n, da, dx, incx, dy, incy);
    }
    
    d_type_wr dcabs1_(z_type_wr *z__)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dcabs1_fptr(z__);
    }
    
    i_type_wr dcopy_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dcopy_fptr(n, dx, incx, dy, incy);
    }
    
    d_type_wr ddot_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ddot_fptr(n, dx, incx, dy, incy);
    }
    
    i_type_wr dgbmv_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
                     d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *x, i_type_wr *incx,
                     d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dgbmv_fptr(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr dgemm_(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k, d_type_wr *alpha,
                     d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta, 
                     d_type_wr *c__, i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dgemm_fptr(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr dgemv_(char *trans, i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, 
                     d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dgemv_fptr(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr dger_(i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx, d_type_wr *y, 
                    i_type_wr *incy, d_type_wr *a, i_type_wr *lda)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dger_fptr(m, n, alpha, x, incx, y, incy, a, lda);
    }

    d_type_wr dnrm2_(i_type_wr *n, d_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dnrm2_fptr(n, x, incx);
    }
    
    i_type_wr drot_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy,
                    d_type_wr *c__, d_type_wr *s)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->drot_fptr(n, dx, incx, dy, incy, c__, s);
    }

    i_type_wr drotg_(d_type_wr *da, d_type_wr *db, d_type_wr *c__, d_type_wr *s)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->drotg_fptr(da, db, c__, s);
    }
    
    i_type_wr drotm_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy, 
                     d_type_wr *dparam)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->drotm_fptr(n, dx, incx, dy, incy, dparam);
    }

    i_type_wr drotmg_(d_type_wr *dd1, d_type_wr *dd2, d_type_wr *dx1, d_type_wr *dy1, d_type_wr *dparam)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->drotmg_fptr(dd1, dd2, dx1, dy1, dparam);
    }
    
    i_type_wr dsbmv_(char *uplo, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda,
                     d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dsbmv_fptr(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr dscal_(i_type_wr *n, d_type_wr *da, d_type_wr *dx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dscal_fptr(n, da, dx, incx);
    }
    
    d_type_wr dsdot_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
    {
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dsdot_fptr(n, sx, incx, sy, incy);
    }
    
    i_type_wr dspmv_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *ap, d_type_wr *x, 
                     i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dspmv_fptr(uplo, n, alpha, ap, x, incx, beta, y, incy);
    }

    i_type_wr dspr_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx, d_type_wr *ap)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dspr_fptr(uplo, n, alpha, x, incx, ap);
    }

    i_type_wr dspr2_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx, d_type_wr *y, 
                     i_type_wr *incy, d_type_wr *ap)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dspr2_fptr(uplo, n, alpha, x, incx, y, incy, ap);
    }

    i_type_wr dswap_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dswap_fptr(n, dx, incx, dy, incy);
    }
    
    i_type_wr dsymm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *a, 
                     i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta, d_type_wr *c__, 
                     i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dsymm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr dsymv_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *x, 
                     i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dsymv_fptr(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr dsyr_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx, d_type_wr *a,
                    i_type_wr *lda)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dsyr_fptr(uplo, n, alpha, x, incx, a, lda);
    }

    i_type_wr dsyr2_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx, d_type_wr *y, 
                     i_type_wr *incy, d_type_wr *a, i_type_wr *lda)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dsyr2_fptr(uplo, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr dsyr2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, d_type_wr *a,
                      i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta, d_type_wr *c__, 
                      i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dsyr2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr dsyrk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, d_type_wr *a,
                     i_type_wr *lda, d_type_wr *beta, d_type_wr *c__, i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dsyrk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
    }

    i_type_wr dtbmv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, d_type_wr *a, 
                     i_type_wr *lda, d_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dtbmv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr dtbsv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, d_type_wr *a,
                     i_type_wr *lda, d_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dtbsv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr dtpmv_(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *ap, d_type_wr *x,
                     i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dtpmv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr dtpsv_(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *ap, d_type_wr *x, 
                     i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dtpsv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr dtrmm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n, 
                     d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dtrmm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr dtrmv_(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *a, i_type_wr *lda, 
                     d_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dtrmv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    i_type_wr dtrsm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n, 
                     d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dtrsm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr dtrsv_(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *a, i_type_wr *lda, 
                     d_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dtrsv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    d_type_wr dzasum_(i_type_wr *n, z_type_wr *zx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dzasum_fptr(n, zx, incx);
    }
    
    d_type_wr dznrm2_(i_type_wr *n, z_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->dznrm2_fptr(n, x, incx);
    }
    
    i_type_wr icamax_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->icamax_fptr(n, cx, incx);
    }
    
    i_type_wr idamax_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->idamax_fptr(n, dx, incx);
    }
    
    i_type_wr isamax_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->isamax_fptr(n, sx, incx);
    }
    
    i_type_wr izamax_(i_type_wr *n, z_type_wr *zx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->izamax_fptr(n, zx, incx);
    }

    s_type_ret_wr sasum_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->sasum_fptr(n, sx, incx);
    }

    i_type_wr saxpy_(i_type_wr *n, s_type_wr *sa, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, 
                     i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->saxpy_fptr(n, sa, sx, incx, sy, incy);
    }

    s_type_ret_wr scabs1_(c_type_wr *z__)
    {
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->scabs1_fptr(z__);
    }
    
    s_type_ret_wr scasum_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->scasum_fptr(n, cx, incx);
    }
    
    s_type_ret_wr scnrm2_(i_type_wr *n, c_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->scnrm2_fptr(n, x, incx);
    }
    
    i_type_wr scopy_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->scopy_fptr(n, sx, incx, sy, incy);
    }
    
    s_type_ret_wr sdot_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->sdot_fptr(n, sx, incx, sy, incy);
    }
    
    s_type_ret_wr sdsdot_(i_type_wr *n, s_type_wr *sb, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
    {
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->sdsdot_fptr(n, sb, sx, incx, sy, incy);
    }
    
    i_type_wr sgbmv_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, s_type_wr *alpha, 
                     s_type_wr *a, i_type_wr *lda, s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y,
                     i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->sgbmv_fptr(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr sgemm_(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, 
                     s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *beta, s_type_wr *c__,
                     i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->sgemm_fptr(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr sgemv_(char *trans, i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, 
                     s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->sgemv_fptr(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr sger_(i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx, s_type_wr *y, 
                    i_type_wr *incy, s_type_wr *a, i_type_wr *lda)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->sger_fptr(m, n, alpha, x, incx, y, incy, a, lda);
    }

    s_type_ret_wr snrm2_(i_type_wr *n, s_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->snrm2_fptr(n, x, incx);
    }
    
    i_type_wr srot_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy, s_type_wr *c__,
                    s_type_wr *s)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->srot_fptr(n, sx, incx, sy, incy, c__, s);
    }

    i_type_wr srotg_(s_type_wr *sa, s_type_wr *sb, s_type_wr *c__, s_type_wr *s)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->srotg_fptr(sa, sb, c__, s);
    }
    
    i_type_wr srotm_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy, 
                     s_type_wr *sparam)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->srotm_fptr(n, sx, incx, sy, incy, sparam);
    }

    i_type_wr srotmg_(s_type_wr *sd1, s_type_wr *sd2, s_type_wr *sx1, s_type_wr *sy1, s_type_wr *sparam)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->srotmg_fptr(sd1, sd2, sx1, sy1, sparam);
    }
    
    i_type_wr ssbmv_(char *uplo, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, 
                     s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ssbmv_fptr(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr sscal_(i_type_wr *n, s_type_wr *sa, s_type_wr *sx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->sscal_fptr(n, sa, sx, incx);
    }
    
    i_type_wr sspmv_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *ap, s_type_wr *x, i_type_wr *incx, 
                     s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->sspmv_fptr(uplo, n, alpha, ap, x, incx, beta, y, incy);
    }

    i_type_wr sspr_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx, s_type_wr *ap)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->sspr_fptr(uplo, n, alpha, x, incx, ap);
    }
    
    i_type_wr sspr2_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx, s_type_wr *y, 
                     i_type_wr *incy, s_type_wr *ap)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->sspr2_fptr(uplo, n, alpha, x, incx, y, incy, ap);
    }

    i_type_wr sswap_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->sswap_fptr(n, sx, incx, sy, incy);
    }
    
    i_type_wr ssymm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *a, 
                     i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ssymm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr ssymv_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *x, 
                     i_type_wr *incx, s_type_wr *beta, s_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ssymv_fptr(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr ssyr_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx, s_type_wr *a,
                    i_type_wr *lda)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ssyr_fptr(uplo, n, alpha, x, incx, a, lda);
    }

    i_type_wr ssyr2_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx, s_type_wr *y,
                     i_type_wr *incy, s_type_wr *a, i_type_wr *lda)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ssyr2_fptr(uplo, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr ssyr2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, s_type_wr *a,
                      i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *beta, s_type_wr *c__, 
                      i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ssyr2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr ssyrk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, s_type_wr *a, 
                     i_type_wr *lda, s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ssyrk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
    }
    
    i_type_wr stbmv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, s_type_wr *a, 
                     i_type_wr *lda, s_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->stbmv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr stbsv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, s_type_wr *a,
                     i_type_wr *lda, s_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->stbsv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr stpmv_(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *ap, s_type_wr *x, 
                     i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->stpmv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr stpsv_(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *ap, s_type_wr *x, 
                     i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->stpsv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr strmm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                     s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->strmm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr strmv_(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *a, i_type_wr *lda,
                     s_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->strmv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    i_type_wr strsm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                     s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->strsm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr strsv_(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
                     s_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->strsv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    i_type_wr zaxpy_(i_type_wr *n, z_type_wr *za, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zaxpy_fptr(n, za, zx, incx, zy, incy);
    }

    i_type_wr zcopy_(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zcopy_fptr(n, zx, incx, zy, incy);
    }
    
    i_type_wr zdotc_(z_type_wr *ret_val, i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, 
                     i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zdotc_fptr(ret_val, n, zx, incx, zy, incy);
    }

    i_type_wr zdotu_(z_type_wr *ret_val, i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, 
                     i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zdotu_fptr(ret_val, n, zx, incx, zy, incy);
    }

    i_type_wr zdrot_(i_type_wr *n, z_type_wr *cx, i_type_wr *incx, z_type_wr *cy, i_type_wr *incy, 
                     d_type_wr *c__, d_type_wr *s)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zdrot_fptr(n, cx, incx, cy, incy, c__, s);
    }

    i_type_wr zdscal_(i_type_wr *n, d_type_wr *da, z_type_wr *zx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zdscal_fptr(n, da, zx, incx);
    }
    
    i_type_wr zgbmv_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, z_type_wr *alpha,
                     z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y,
                     i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zgbmv_fptr(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr zgemm_(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, 
                     z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta, 
                     z_type_wr *c__, i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zgemm_fptr(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr zgemv_(char *trans, i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda,
                     z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zgemv_fptr(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr zgerc_(i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx, z_type_wr *y,
                     i_type_wr *incy, z_type_wr *a, i_type_wr *lda)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zgerc_fptr(m, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr zgeru_(i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx, z_type_wr *y, 
                     i_type_wr *incy, z_type_wr *a, i_type_wr *lda)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zgeru_fptr(m, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr zhbmv_(char *uplo, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, 
                     z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zhbmv_fptr(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr zhemm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *a, 
                     i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zhemm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr zhemv_(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *x,
                     i_type_wr *incx, z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zhemv_fptr(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    i_type_wr zher_(char *uplo, i_type_wr *n, d_type_wr *alpha, z_type_wr *x, i_type_wr *incx, z_type_wr *a,
                    i_type_wr *lda)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zher_fptr(uplo, n, alpha, x, incx, a, lda);
    }

    i_type_wr zher2_(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx, z_type_wr *y, 
                     i_type_wr *incy, z_type_wr *a, i_type_wr *lda)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zher2_fptr(uplo, n, alpha, x, incx, y, incy, a, lda);
    }

    i_type_wr zher2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, z_type_wr *a, 
                      i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, d_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zher2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr zherk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, z_type_wr *a,
                     i_type_wr *lda, d_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zherk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
    }

    i_type_wr zhpmv_(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *ap, z_type_wr *x, i_type_wr *incx,
                     z_type_wr *beta, z_type_wr *y, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zhpmv_fptr(uplo, n, alpha, ap, x, incx, beta, y, incy);
    }

    i_type_wr zhpr_(char *uplo, i_type_wr *n, d_type_wr *alpha, z_type_wr *x, i_type_wr *incx, z_type_wr *ap)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zhpr_fptr(uplo, n, alpha, x, incx, ap);
    }
    
    i_type_wr zhpr2_(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx, z_type_wr *y,
                     i_type_wr *incy, z_type_wr *ap)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zhpr2_fptr(uplo, n, alpha, x, incx, y, incy, ap);
    }

    i_type_wr zrotg_(z_type_wr *ca, z_type_wr *cb, d_type_wr *c__, z_type_wr *s)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zrotg_fptr(ca, cb, c__, s);
    }
    
    i_type_wr zscal_(i_type_wr *n, z_type_wr *za, z_type_wr *zx, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zscal_fptr(n, za, zx, incx);
    }
    
    i_type_wr zswap_(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, i_type_wr *incy)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zswap_fptr(n, zx, incx, zy, incy);
    }
    
    i_type_wr zsymm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *a,
                     i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zsymm_fptr(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr zsyr2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, z_type_wr *a, 
                      i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zsyr2k_fptr(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
    }

    i_type_wr zsyrk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, z_type_wr *a, 
                     i_type_wr *lda, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->zsyrk_fptr(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc);
    }

    i_type_wr ztbmv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, z_type_wr *a, 
                     i_type_wr *lda, z_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ztbmv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr ztbsv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, z_type_wr *a, 
                     i_type_wr *lda, z_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ztbsv_fptr(uplo, trans, diag, n, k, a, lda, x, incx);
    }

    i_type_wr ztpmv_(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *ap, z_type_wr *x, 
                     i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ztpmv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr ztpsv_(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *ap, z_type_wr *x, 
                     i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ztpsv_fptr(uplo, trans, diag, n, ap, x, incx);
    }

    i_type_wr ztrmm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n, 
                     z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ztrmm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr ztrmv_(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *a, i_type_wr *lda, 
                     z_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ztrmv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }

    i_type_wr ztrsm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n, 
                     z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ztrsm_fptr(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }

    i_type_wr ztrsv_(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *a, i_type_wr *lda, 
                     z_type_wr *x, i_type_wr *incx)
    {   
        ml::plugin_manager::init_plugin();
        return ml::plugin_manager::plugin()->ztrsv_fptr(uplo, trans, diag, n, a, lda, x, incx);
    }
}

}
