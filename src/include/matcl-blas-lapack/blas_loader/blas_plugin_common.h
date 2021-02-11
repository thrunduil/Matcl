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

#pragma once

#include "matcl-blas-lapack/blas_loader/function_names_macros.h"

template<typename Fun>
struct function_info{};

template<class Ret, class ... Args>
struct function_info<Ret (*)(Args...)>
{
    using return_type = Ret;
};

template<class Ret, class ... Args>
struct function_info<Ret (Args...)>
{
    using return_type = Ret;
};

template<class T>
struct type_sizeof
{
    static const size_t value = sizeof(T);
};

template<>
struct type_sizeof<void>
{
    static const size_t value = sizeof(int);
};

template<typename Fun, typename F_in>
Fun blas_loader_cast(const F_in* f)
{
    using ret_1 = typename function_info<Fun>::return_type;
    using ret_2 = typename function_info<F_in>::return_type;

    static const size_t n1   = type_sizeof<ret_1>::value;
    static const size_t n2   = type_sizeof<ret_2>::value;
    static_assert(n1 == n2, "incompatible return types");
    
    return reinterpret_cast<Fun>(f);
}

#define BLAS_LOADER_ASSIGN(fun)     \
    fun##_fptr(blas_loader_cast<fun##_func_type>(&FUNCTION_NAME_##fun))

::blas_plugin::blas_plugin() 
    : 
      get_name_fptr(reinterpret_cast<get_name_type>(get_name)),
      initialize_fptr(reinterpret_cast<initialize_type>(initialize)),
      force_initialization_fptr(reinterpret_cast<initialize_type>(force_initialization)),
      get_num_threads_fptr(reinterpret_cast<get_num_threads_type>(get_num_threads)),
      get_default_threads_fptr(reinterpret_cast<get_default_threads_type>(get_default_threads)),
      set_num_threads_fptr(reinterpret_cast<set_num_threads_type>(set_num_threads)),
      are_user_threads_allowed_fptr(reinterpret_cast<user_threads_allowed_type>(are_user_threads_allowed)),

      BLAS_LOADER_ASSIGN(caxpy),
      BLAS_LOADER_ASSIGN(ccopy),
      BLAS_LOADER_ASSIGN(cdotc),
      BLAS_LOADER_ASSIGN(cdotu),
      BLAS_LOADER_ASSIGN(cgbmv),
      BLAS_LOADER_ASSIGN(cgemm),
      BLAS_LOADER_ASSIGN(cgemv),
      BLAS_LOADER_ASSIGN(cgerc),
      BLAS_LOADER_ASSIGN(cgeru),
      BLAS_LOADER_ASSIGN(chbmv),
      BLAS_LOADER_ASSIGN(chemm),
      BLAS_LOADER_ASSIGN(chemv),
      BLAS_LOADER_ASSIGN(cher),
      BLAS_LOADER_ASSIGN(cher2),
      BLAS_LOADER_ASSIGN(cher2k),
      BLAS_LOADER_ASSIGN(cherk),
      BLAS_LOADER_ASSIGN(chpmv),
      BLAS_LOADER_ASSIGN(chpr),
      BLAS_LOADER_ASSIGN(chpr2),
      BLAS_LOADER_ASSIGN(crotg),
      BLAS_LOADER_ASSIGN(cscal),
      BLAS_LOADER_ASSIGN(csrot),
      BLAS_LOADER_ASSIGN(csscal),
      BLAS_LOADER_ASSIGN(cswap),
      BLAS_LOADER_ASSIGN(csymm),
      BLAS_LOADER_ASSIGN(csyr2k),
      BLAS_LOADER_ASSIGN(csyrk),
      BLAS_LOADER_ASSIGN(ctbmv),
      BLAS_LOADER_ASSIGN(ctbsv),
      BLAS_LOADER_ASSIGN(ctpmv),
      BLAS_LOADER_ASSIGN(ctpsv),
      BLAS_LOADER_ASSIGN(ctrmm),
      BLAS_LOADER_ASSIGN(ctrmv),
      BLAS_LOADER_ASSIGN(ctrsm),
      BLAS_LOADER_ASSIGN(ctrsv),
      BLAS_LOADER_ASSIGN(dasum),
      BLAS_LOADER_ASSIGN(daxpy),
      BLAS_LOADER_ASSIGN(dcabs1),
      BLAS_LOADER_ASSIGN(dcopy),
      BLAS_LOADER_ASSIGN(ddot),
      BLAS_LOADER_ASSIGN(dgbmv),
      BLAS_LOADER_ASSIGN(dgemm),
      BLAS_LOADER_ASSIGN(dgemv),
      BLAS_LOADER_ASSIGN(dger),
      BLAS_LOADER_ASSIGN(dnrm2),
      BLAS_LOADER_ASSIGN(drot),
      BLAS_LOADER_ASSIGN(drotg),
      BLAS_LOADER_ASSIGN(drotm),
      BLAS_LOADER_ASSIGN(drotmg),
      BLAS_LOADER_ASSIGN(dsbmv),
      BLAS_LOADER_ASSIGN(dscal),
      BLAS_LOADER_ASSIGN(dsdot),
      BLAS_LOADER_ASSIGN(dspmv),
      BLAS_LOADER_ASSIGN(dspr),
      BLAS_LOADER_ASSIGN(dspr2),
      BLAS_LOADER_ASSIGN(dswap),
      BLAS_LOADER_ASSIGN(dsymm),
      BLAS_LOADER_ASSIGN(dsymv),
      BLAS_LOADER_ASSIGN(dsyr),
      BLAS_LOADER_ASSIGN(dsyr2),
      BLAS_LOADER_ASSIGN(dsyr2k),
      BLAS_LOADER_ASSIGN(dsyrk),
      BLAS_LOADER_ASSIGN(dtbmv),
      BLAS_LOADER_ASSIGN(dtbsv),
      BLAS_LOADER_ASSIGN(dtpmv),
      BLAS_LOADER_ASSIGN(dtpsv),
      BLAS_LOADER_ASSIGN(dtrmm),
      BLAS_LOADER_ASSIGN(dtrmv),
      BLAS_LOADER_ASSIGN(dtrsm),
      BLAS_LOADER_ASSIGN(dtrsv),
      BLAS_LOADER_ASSIGN(dzasum),
      BLAS_LOADER_ASSIGN(dznrm2),
      BLAS_LOADER_ASSIGN(icamax),
      BLAS_LOADER_ASSIGN(idamax),
      BLAS_LOADER_ASSIGN(isamax),
      BLAS_LOADER_ASSIGN(izamax),
      BLAS_LOADER_ASSIGN(sasum),
      BLAS_LOADER_ASSIGN(saxpy),
      BLAS_LOADER_ASSIGN(scabs1),
      BLAS_LOADER_ASSIGN(scasum),
      BLAS_LOADER_ASSIGN(scnrm2),
      BLAS_LOADER_ASSIGN(scopy),
      BLAS_LOADER_ASSIGN(sdot),
      BLAS_LOADER_ASSIGN(sdsdot),
      BLAS_LOADER_ASSIGN(sgbmv),
      BLAS_LOADER_ASSIGN(sgemm),
      BLAS_LOADER_ASSIGN(sgemv),
      BLAS_LOADER_ASSIGN(sger),
      BLAS_LOADER_ASSIGN(snrm2),
      BLAS_LOADER_ASSIGN(srot),
      BLAS_LOADER_ASSIGN(srotg),
      BLAS_LOADER_ASSIGN(srotm),
      BLAS_LOADER_ASSIGN(srotmg),
      BLAS_LOADER_ASSIGN(ssbmv),
      BLAS_LOADER_ASSIGN(sscal),
      BLAS_LOADER_ASSIGN(sspmv),
      BLAS_LOADER_ASSIGN(sspr),
      BLAS_LOADER_ASSIGN(sspr2),
      BLAS_LOADER_ASSIGN(sswap),
      BLAS_LOADER_ASSIGN(ssymm),
      BLAS_LOADER_ASSIGN(ssymv),
      BLAS_LOADER_ASSIGN(ssyr),
      BLAS_LOADER_ASSIGN(ssyr2),
      BLAS_LOADER_ASSIGN(ssyr2k),
      BLAS_LOADER_ASSIGN(ssyrk),
      BLAS_LOADER_ASSIGN(stbmv),
      BLAS_LOADER_ASSIGN(stbsv),
      BLAS_LOADER_ASSIGN(stpmv),
      BLAS_LOADER_ASSIGN(stpsv),
      BLAS_LOADER_ASSIGN(strmm),
      BLAS_LOADER_ASSIGN(strmv),
      BLAS_LOADER_ASSIGN(strsm),
      BLAS_LOADER_ASSIGN(strsv),
      BLAS_LOADER_ASSIGN(zaxpy),
      BLAS_LOADER_ASSIGN(zcopy),
      BLAS_LOADER_ASSIGN(zdotc),
      BLAS_LOADER_ASSIGN(zdotu),
      BLAS_LOADER_ASSIGN(zdrot),
      BLAS_LOADER_ASSIGN(zdscal),
      BLAS_LOADER_ASSIGN(zgbmv),
      BLAS_LOADER_ASSIGN(zgemm),
      BLAS_LOADER_ASSIGN(zgemv),
      BLAS_LOADER_ASSIGN(zgerc),
      BLAS_LOADER_ASSIGN(zgeru),
      BLAS_LOADER_ASSIGN(zhbmv),
      BLAS_LOADER_ASSIGN(zhemm),
      BLAS_LOADER_ASSIGN(zhemv),
      BLAS_LOADER_ASSIGN(zher),
      BLAS_LOADER_ASSIGN(zher2),
      BLAS_LOADER_ASSIGN(zher2k),
      BLAS_LOADER_ASSIGN(zherk),
      BLAS_LOADER_ASSIGN(zhpmv),
      BLAS_LOADER_ASSIGN(zhpr),
      BLAS_LOADER_ASSIGN(zhpr2),
      BLAS_LOADER_ASSIGN(zrotg),
      BLAS_LOADER_ASSIGN(zscal),
      BLAS_LOADER_ASSIGN(zswap),
      BLAS_LOADER_ASSIGN(zsymm),
      BLAS_LOADER_ASSIGN(zsyr2k),
      BLAS_LOADER_ASSIGN(zsyrk),
      BLAS_LOADER_ASSIGN(ztbmv),
      BLAS_LOADER_ASSIGN(ztbsv),
      BLAS_LOADER_ASSIGN(ztpmv),
      BLAS_LOADER_ASSIGN(ztpsv),
      BLAS_LOADER_ASSIGN(ztrmm),
      BLAS_LOADER_ASSIGN(ztrmv),
      BLAS_LOADER_ASSIGN(ztrsm),
      BLAS_LOADER_ASSIGN(ztrsv)
{}
