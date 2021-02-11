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

// define names of exporting functions
#ifndef FUNCTION_NAME_caxpy
    #define FUNCTION_NAME_caxpy CALL_SYNTAX(caxpy)
#endif

#ifndef FUNCTION_NAME_ccopy
    #define FUNCTION_NAME_ccopy CALL_SYNTAX(ccopy)
#endif

#ifndef FUNCTION_NAME_cdotc
    #define FUNCTION_NAME_cdotc CALL_SYNTAX(cdotc)
#endif

#ifndef FUNCTION_NAME_cdotu
    #define FUNCTION_NAME_cdotu CALL_SYNTAX(cdotu)
#endif

#ifndef FUNCTION_NAME_cgbmv
    #define FUNCTION_NAME_cgbmv CALL_SYNTAX(cgbmv)
#endif

#ifndef FUNCTION_NAME_cgemm
    #define FUNCTION_NAME_cgemm CALL_SYNTAX(cgemm)
#endif

#ifndef FUNCTION_NAME_cgemv
    #define FUNCTION_NAME_cgemv CALL_SYNTAX(cgemv)
#endif

#ifndef FUNCTION_NAME_cgerc
    #define FUNCTION_NAME_cgerc CALL_SYNTAX(cgerc)
#endif

#ifndef FUNCTION_NAME_cgeru
    #define FUNCTION_NAME_cgeru CALL_SYNTAX(cgeru)
#endif

#ifndef FUNCTION_NAME_chbmv
    #define FUNCTION_NAME_chbmv CALL_SYNTAX(chbmv)
#endif

#ifndef FUNCTION_NAME_chemm
    #define FUNCTION_NAME_chemm CALL_SYNTAX(chemm)
#endif

#ifndef FUNCTION_NAME_chemv
    #define FUNCTION_NAME_chemv CALL_SYNTAX(chemv)
#endif

#ifndef FUNCTION_NAME_cher
    #define FUNCTION_NAME_cher CALL_SYNTAX(cher)
#endif

#ifndef FUNCTION_NAME_cher2
    #define FUNCTION_NAME_cher2 CALL_SYNTAX(cher2)
#endif

#ifndef FUNCTION_NAME_cher2k
    #define FUNCTION_NAME_cher2k CALL_SYNTAX(cher2k)
#endif

#ifndef FUNCTION_NAME_cherk
    #define FUNCTION_NAME_cherk CALL_SYNTAX(cherk)
#endif

#ifndef FUNCTION_NAME_chpmv
    #define FUNCTION_NAME_chpmv CALL_SYNTAX(chpmv)
#endif

#ifndef FUNCTION_NAME_chpr
    #define FUNCTION_NAME_chpr CALL_SYNTAX(chpr)
#endif

#ifndef FUNCTION_NAME_chpr2
    #define FUNCTION_NAME_chpr2 CALL_SYNTAX(chpr2)
#endif

#ifndef FUNCTION_NAME_crotg
    #define FUNCTION_NAME_crotg CALL_SYNTAX(crotg)
#endif

#ifndef FUNCTION_NAME_cscal
    #define FUNCTION_NAME_cscal CALL_SYNTAX(cscal)
#endif

#ifndef FUNCTION_NAME_csrot
    #define FUNCTION_NAME_csrot CALL_SYNTAX(csrot)
#endif

#ifndef FUNCTION_NAME_csscal
    #define FUNCTION_NAME_csscal CALL_SYNTAX(csscal)
#endif

#ifndef FUNCTION_NAME_cswap
    #define FUNCTION_NAME_cswap CALL_SYNTAX(cswap)
#endif

#ifndef FUNCTION_NAME_csymm
    #define FUNCTION_NAME_csymm CALL_SYNTAX(csymm)
#endif

#ifndef FUNCTION_NAME_csyr2k
    #define FUNCTION_NAME_csyr2k CALL_SYNTAX(csyr2k)
#endif

#ifndef FUNCTION_NAME_csyrk
    #define FUNCTION_NAME_csyrk CALL_SYNTAX(csyrk)
#endif

#ifndef FUNCTION_NAME_ctbmv
    #define FUNCTION_NAME_ctbmv CALL_SYNTAX(ctbmv)
#endif

#ifndef FUNCTION_NAME_ctbsv
    #define FUNCTION_NAME_ctbsv CALL_SYNTAX(ctbsv)
#endif

#ifndef FUNCTION_NAME_ctpmv
    #define FUNCTION_NAME_ctpmv CALL_SYNTAX(ctpmv)
#endif

#ifndef FUNCTION_NAME_ctpsv
    #define FUNCTION_NAME_ctpsv CALL_SYNTAX(ctpsv)
#endif

#ifndef FUNCTION_NAME_ctrmm
    #define FUNCTION_NAME_ctrmm CALL_SYNTAX(ctrmm)
#endif

#ifndef FUNCTION_NAME_ctrmv
    #define FUNCTION_NAME_ctrmv CALL_SYNTAX(ctrmv)
#endif

#ifndef FUNCTION_NAME_ctrsm
    #define FUNCTION_NAME_ctrsm CALL_SYNTAX(ctrsm)
#endif

#ifndef FUNCTION_NAME_ctrsv
    #define FUNCTION_NAME_ctrsv CALL_SYNTAX(ctrsv)
#endif

#ifndef FUNCTION_NAME_dasum
    #define FUNCTION_NAME_dasum CALL_SYNTAX(dasum)
#endif

#ifndef FUNCTION_NAME_daxpy
    #define FUNCTION_NAME_daxpy CALL_SYNTAX(daxpy)
#endif

#ifndef FUNCTION_NAME_dcabs1
    #define FUNCTION_NAME_dcabs1 CALL_SYNTAX(dcabs1)
#endif

#ifndef FUNCTION_NAME_dcopy
    #define FUNCTION_NAME_dcopy CALL_SYNTAX(dcopy)
#endif

#ifndef FUNCTION_NAME_ddot
    #define FUNCTION_NAME_ddot CALL_SYNTAX(ddot)
#endif

#ifndef FUNCTION_NAME_dgbmv
    #define FUNCTION_NAME_dgbmv CALL_SYNTAX(dgbmv)
#endif

#ifndef FUNCTION_NAME_dgemm
    #define FUNCTION_NAME_dgemm CALL_SYNTAX(dgemm)
#endif

#ifndef FUNCTION_NAME_dgemv
    #define FUNCTION_NAME_dgemv CALL_SYNTAX(dgemv)
#endif

#ifndef FUNCTION_NAME_dger
    #define FUNCTION_NAME_dger CALL_SYNTAX(dger)
#endif

#ifndef FUNCTION_NAME_dnrm2
    #define FUNCTION_NAME_dnrm2 CALL_SYNTAX(dnrm2)
#endif

#ifndef FUNCTION_NAME_drot
    #define FUNCTION_NAME_drot CALL_SYNTAX(drot)
#endif

#ifndef FUNCTION_NAME_drotg
    #define FUNCTION_NAME_drotg CALL_SYNTAX(drotg)
#endif

#ifndef FUNCTION_NAME_drotm
    #define FUNCTION_NAME_drotm CALL_SYNTAX(drotm)
#endif

#ifndef FUNCTION_NAME_drotmg
    #define FUNCTION_NAME_drotmg CALL_SYNTAX(drotmg)
#endif

#ifndef FUNCTION_NAME_dsbmv
    #define FUNCTION_NAME_dsbmv CALL_SYNTAX(dsbmv)
#endif

#ifndef FUNCTION_NAME_dscal
    #define FUNCTION_NAME_dscal CALL_SYNTAX(dscal)
#endif

#ifndef FUNCTION_NAME_dsdot
    #define FUNCTION_NAME_dsdot CALL_SYNTAX(dsdot)
#endif

#ifndef FUNCTION_NAME_dspmv
    #define FUNCTION_NAME_dspmv CALL_SYNTAX(dspmv)
#endif

#ifndef FUNCTION_NAME_dspr
    #define FUNCTION_NAME_dspr CALL_SYNTAX(dspr)
#endif

#ifndef FUNCTION_NAME_dspr2
    #define FUNCTION_NAME_dspr2 CALL_SYNTAX(dspr2)
#endif

#ifndef FUNCTION_NAME_dswap
    #define FUNCTION_NAME_dswap CALL_SYNTAX(dswap)
#endif

#ifndef FUNCTION_NAME_dsymm
    #define FUNCTION_NAME_dsymm CALL_SYNTAX(dsymm)
#endif

#ifndef FUNCTION_NAME_dsymv
    #define FUNCTION_NAME_dsymv CALL_SYNTAX(dsymv)
#endif

#ifndef FUNCTION_NAME_dsyr
    #define FUNCTION_NAME_dsyr CALL_SYNTAX(dsyr)
#endif

#ifndef FUNCTION_NAME_dsyr2
    #define FUNCTION_NAME_dsyr2 CALL_SYNTAX(dsyr2)
#endif

#ifndef FUNCTION_NAME_dsyr2k
    #define FUNCTION_NAME_dsyr2k CALL_SYNTAX(dsyr2k)
#endif

#ifndef FUNCTION_NAME_dsyrk
    #define FUNCTION_NAME_dsyrk CALL_SYNTAX(dsyrk)
#endif

#ifndef FUNCTION_NAME_dtbmv
    #define FUNCTION_NAME_dtbmv CALL_SYNTAX(dtbmv)
#endif

#ifndef FUNCTION_NAME_dtbsv
    #define FUNCTION_NAME_dtbsv CALL_SYNTAX(dtbsv)
#endif

#ifndef FUNCTION_NAME_dtpmv
    #define FUNCTION_NAME_dtpmv CALL_SYNTAX(dtpmv)
#endif

#ifndef FUNCTION_NAME_dtpsv
    #define FUNCTION_NAME_dtpsv CALL_SYNTAX(dtpsv)
#endif

#ifndef FUNCTION_NAME_dtrmm
    #define FUNCTION_NAME_dtrmm CALL_SYNTAX(dtrmm)
#endif

#ifndef FUNCTION_NAME_dtrmv
    #define FUNCTION_NAME_dtrmv CALL_SYNTAX(dtrmv)
#endif

#ifndef FUNCTION_NAME_dtrsm
    #define FUNCTION_NAME_dtrsm CALL_SYNTAX(dtrsm)
#endif

#ifndef FUNCTION_NAME_dtrsv
    #define FUNCTION_NAME_dtrsv CALL_SYNTAX(dtrsv)
#endif

#ifndef FUNCTION_NAME_dzasum
    #define FUNCTION_NAME_dzasum CALL_SYNTAX(dzasum)
#endif

#ifndef FUNCTION_NAME_dznrm2
    #define FUNCTION_NAME_dznrm2 CALL_SYNTAX(dznrm2)
#endif

#ifndef FUNCTION_NAME_icamax
    #define FUNCTION_NAME_icamax CALL_SYNTAX(icamax)
#endif

#ifndef FUNCTION_NAME_idamax
    #define FUNCTION_NAME_idamax CALL_SYNTAX(idamax)
#endif

#ifndef FUNCTION_NAME_isamax
    #define FUNCTION_NAME_isamax CALL_SYNTAX(isamax)
#endif

#ifndef FUNCTION_NAME_izamax
    #define FUNCTION_NAME_izamax CALL_SYNTAX(izamax)
#endif

#ifndef FUNCTION_NAME_sasum
    #define FUNCTION_NAME_sasum CALL_SYNTAX(sasum)
#endif

#ifndef FUNCTION_NAME_saxpy
    #define FUNCTION_NAME_saxpy CALL_SYNTAX(saxpy)
#endif

#ifndef FUNCTION_NAME_scabs1
    #define FUNCTION_NAME_scabs1 CALL_SYNTAX(scabs1)
#endif

#ifndef FUNCTION_NAME_scasum
    #define FUNCTION_NAME_scasum CALL_SYNTAX(scasum)
#endif

#ifndef FUNCTION_NAME_scnrm2
    #define FUNCTION_NAME_scnrm2 CALL_SYNTAX(scnrm2)
#endif

#ifndef FUNCTION_NAME_scopy
    #define FUNCTION_NAME_scopy CALL_SYNTAX(scopy)
#endif

#ifndef FUNCTION_NAME_sdot
    #define FUNCTION_NAME_sdot CALL_SYNTAX(sdot)
#endif

#ifndef FUNCTION_NAME_sdsdot
    #define FUNCTION_NAME_sdsdot CALL_SYNTAX(sdsdot)
#endif

#ifndef FUNCTION_NAME_sgbmv
    #define FUNCTION_NAME_sgbmv CALL_SYNTAX(sgbmv)
#endif

#ifndef FUNCTION_NAME_sgemm
    #define FUNCTION_NAME_sgemm CALL_SYNTAX(sgemm)
#endif

#ifndef FUNCTION_NAME_sgemv
    #define FUNCTION_NAME_sgemv CALL_SYNTAX(sgemv)
#endif

#ifndef FUNCTION_NAME_sger
    #define FUNCTION_NAME_sger CALL_SYNTAX(sger)
#endif

#ifndef FUNCTION_NAME_snrm2
    #define FUNCTION_NAME_snrm2 CALL_SYNTAX(snrm2)
#endif

#ifndef FUNCTION_NAME_srot
    #define FUNCTION_NAME_srot CALL_SYNTAX(srot)
#endif

#ifndef FUNCTION_NAME_srotg
    #define FUNCTION_NAME_srotg CALL_SYNTAX(srotg)
#endif

#ifndef FUNCTION_NAME_srotm
    #define FUNCTION_NAME_srotm CALL_SYNTAX(srotm)
#endif

#ifndef FUNCTION_NAME_srotmg
    #define FUNCTION_NAME_srotmg CALL_SYNTAX(srotmg)
#endif

#ifndef FUNCTION_NAME_ssbmv
    #define FUNCTION_NAME_ssbmv CALL_SYNTAX(ssbmv)
#endif

#ifndef FUNCTION_NAME_sscal
    #define FUNCTION_NAME_sscal CALL_SYNTAX(sscal)
#endif

#ifndef FUNCTION_NAME_sspmv
    #define FUNCTION_NAME_sspmv CALL_SYNTAX(sspmv)
#endif

#ifndef FUNCTION_NAME_sspr
    #define FUNCTION_NAME_sspr CALL_SYNTAX(sspr)
#endif

#ifndef FUNCTION_NAME_sspr2
    #define FUNCTION_NAME_sspr2 CALL_SYNTAX(sspr2)
#endif

#ifndef FUNCTION_NAME_sswap
    #define FUNCTION_NAME_sswap CALL_SYNTAX(sswap)
#endif

#ifndef FUNCTION_NAME_ssymm
    #define FUNCTION_NAME_ssymm CALL_SYNTAX(ssymm)
#endif

#ifndef FUNCTION_NAME_ssymv
    #define FUNCTION_NAME_ssymv CALL_SYNTAX(ssymv)
#endif

#ifndef FUNCTION_NAME_ssyr
    #define FUNCTION_NAME_ssyr CALL_SYNTAX(ssyr)
#endif

#ifndef FUNCTION_NAME_ssyr2
    #define FUNCTION_NAME_ssyr2 CALL_SYNTAX(ssyr2)
#endif

#ifndef FUNCTION_NAME_ssyr2k
    #define FUNCTION_NAME_ssyr2k CALL_SYNTAX(ssyr2k)
#endif

#ifndef FUNCTION_NAME_ssyrk
    #define FUNCTION_NAME_ssyrk CALL_SYNTAX(ssyrk)
#endif

#ifndef FUNCTION_NAME_stbmv
    #define FUNCTION_NAME_stbmv CALL_SYNTAX(stbmv)
#endif

#ifndef FUNCTION_NAME_stbsv
    #define FUNCTION_NAME_stbsv CALL_SYNTAX(stbsv)
#endif

#ifndef FUNCTION_NAME_stpmv
    #define FUNCTION_NAME_stpmv CALL_SYNTAX(stpmv)
#endif

#ifndef FUNCTION_NAME_stpsv
    #define FUNCTION_NAME_stpsv CALL_SYNTAX(stpsv)
#endif

#ifndef FUNCTION_NAME_strmm
    #define FUNCTION_NAME_strmm CALL_SYNTAX(strmm)
#endif

#ifndef FUNCTION_NAME_strmv
    #define FUNCTION_NAME_strmv CALL_SYNTAX(strmv)
#endif

#ifndef FUNCTION_NAME_strsm
    #define FUNCTION_NAME_strsm CALL_SYNTAX(strsm)
#endif

#ifndef FUNCTION_NAME_strsv
    #define FUNCTION_NAME_strsv CALL_SYNTAX(strsv)
#endif

#ifndef FUNCTION_NAME_zaxpy
    #define FUNCTION_NAME_zaxpy CALL_SYNTAX(zaxpy)
#endif

#ifndef FUNCTION_NAME_zcopy
    #define FUNCTION_NAME_zcopy CALL_SYNTAX(zcopy)
#endif

#ifndef FUNCTION_NAME_zdotc
    #define FUNCTION_NAME_zdotc CALL_SYNTAX(zdotc)
#endif

#ifndef FUNCTION_NAME_zdotu
    #define FUNCTION_NAME_zdotu CALL_SYNTAX(zdotu)
#endif

#ifndef FUNCTION_NAME_zdrot
    #define FUNCTION_NAME_zdrot CALL_SYNTAX(zdrot)
#endif

#ifndef FUNCTION_NAME_zdscal
    #define FUNCTION_NAME_zdscal CALL_SYNTAX(zdscal)
#endif

#ifndef FUNCTION_NAME_zgbmv
    #define FUNCTION_NAME_zgbmv CALL_SYNTAX(zgbmv)
#endif

#ifndef FUNCTION_NAME_zgemm
    #define FUNCTION_NAME_zgemm CALL_SYNTAX(zgemm)
#endif

#ifndef FUNCTION_NAME_zgemv
    #define FUNCTION_NAME_zgemv CALL_SYNTAX(zgemv)
#endif

#ifndef FUNCTION_NAME_zgerc
    #define FUNCTION_NAME_zgerc CALL_SYNTAX(zgerc)
#endif

#ifndef FUNCTION_NAME_zgeru
    #define FUNCTION_NAME_zgeru CALL_SYNTAX(zgeru)
#endif

#ifndef FUNCTION_NAME_zhbmv
    #define FUNCTION_NAME_zhbmv CALL_SYNTAX(zhbmv)
#endif

#ifndef FUNCTION_NAME_zhemm
    #define FUNCTION_NAME_zhemm CALL_SYNTAX(zhemm)
#endif

#ifndef FUNCTION_NAME_zhemv
    #define FUNCTION_NAME_zhemv CALL_SYNTAX(zhemv)
#endif

#ifndef FUNCTION_NAME_zher
    #define FUNCTION_NAME_zher CALL_SYNTAX(zher)
#endif

#ifndef FUNCTION_NAME_zher2
    #define FUNCTION_NAME_zher2 CALL_SYNTAX(zher2)
#endif

#ifndef FUNCTION_NAME_zher2k
    #define FUNCTION_NAME_zher2k CALL_SYNTAX(zher2k)
#endif

#ifndef FUNCTION_NAME_zherk
    #define FUNCTION_NAME_zherk CALL_SYNTAX(zherk)
#endif

#ifndef FUNCTION_NAME_zhpmv
    #define FUNCTION_NAME_zhpmv CALL_SYNTAX(zhpmv)
#endif

#ifndef FUNCTION_NAME_zhpr
    #define FUNCTION_NAME_zhpr CALL_SYNTAX(zhpr)
#endif

#ifndef FUNCTION_NAME_zhpr2
    #define FUNCTION_NAME_zhpr2 CALL_SYNTAX(zhpr2)
#endif

#ifndef FUNCTION_NAME_zrotg
    #define FUNCTION_NAME_zrotg CALL_SYNTAX(zrotg)
#endif

#ifndef FUNCTION_NAME_zscal
    #define FUNCTION_NAME_zscal CALL_SYNTAX(zscal)
#endif

#ifndef FUNCTION_NAME_zswap
    #define FUNCTION_NAME_zswap CALL_SYNTAX(zswap)
#endif

#ifndef FUNCTION_NAME_zsymm
    #define FUNCTION_NAME_zsymm CALL_SYNTAX(zsymm)
#endif

#ifndef FUNCTION_NAME_zsyr2k
    #define FUNCTION_NAME_zsyr2k CALL_SYNTAX(zsyr2k)
#endif

#ifndef FUNCTION_NAME_zsyrk
    #define FUNCTION_NAME_zsyrk CALL_SYNTAX(zsyrk)
#endif

#ifndef FUNCTION_NAME_ztbmv
    #define FUNCTION_NAME_ztbmv CALL_SYNTAX(ztbmv)
#endif

#ifndef FUNCTION_NAME_ztbsv
    #define FUNCTION_NAME_ztbsv CALL_SYNTAX(ztbsv)
#endif

#ifndef FUNCTION_NAME_ztpmv
    #define FUNCTION_NAME_ztpmv CALL_SYNTAX(ztpmv)
#endif

#ifndef FUNCTION_NAME_ztpsv
    #define FUNCTION_NAME_ztpsv CALL_SYNTAX(ztpsv)
#endif

#ifndef FUNCTION_NAME_ztrmm
    #define FUNCTION_NAME_ztrmm CALL_SYNTAX(ztrmm)
#endif

#ifndef FUNCTION_NAME_ztrmv
    #define FUNCTION_NAME_ztrmv CALL_SYNTAX(ztrmv)
#endif

#ifndef FUNCTION_NAME_ztrsm
    #define FUNCTION_NAME_ztrsm CALL_SYNTAX(ztrsm)
#endif

#ifndef FUNCTION_NAME_ztrsv
    #define FUNCTION_NAME_ztrsv CALL_SYNTAX(ztrsv)
#endif
