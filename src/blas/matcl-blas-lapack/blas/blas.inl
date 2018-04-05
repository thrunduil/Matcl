/*
 * This file is a part of Matrix Computation Library (MATCL)
 *
 * Copyright (c) Pawe³ Kowal 2017 - 2018
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/blas/details/config_blas_lib.h"
#include <float.h>
#include <iostream>
#include "matcl-blas-lapack/blas/config_blas.h"
#include "matcl-blas-lapack/blas_ext/rot.inl"

#ifndef INLINE_TYPE
#define INLINE_TYPE inline
#endif

namespace matcl { namespace lapack
{
    
//-----------------------------------------------------------------
//                          UTILS
//-----------------------------------------------------------------

BLAS_EXPORT INLINE_TYPE
void xerbla(const char* srname, int info )
{
    std::cerr<<"blas.dll XERBLA: " << srname << " returned info: " << info << std::endl;
}

//-----------------------------------------------------------------
//                          ISNAN
//-----------------------------------------------------------------
template<> BLAS_EXPORT INLINE_TYPE 
bool isnan<s_type>(s_type a)
{
    return std::isnan(a) ? true : false;
};

template<> BLAS_EXPORT INLINE_TYPE 
bool isnan<d_type>(d_type a)
{
    return std::isnan(a) ? true : false;
};

template<> BLAS_EXPORT INLINE_TYPE 
bool isnan<c_type>(c_type a)
{
    return isnan(real(a)) || isnan(imag(a));
};

template<> BLAS_EXPORT INLINE_TYPE 
bool isnan<z_type>(z_type a)
{
    return isnan(real(a)) || isnan(imag(a));
};

// BLAS Level1
//-----------------------------------------------------------------
//                          SCAL
//-----------------------------------------------------------------
template<> BLAS_EXPORT INLINE_TYPE 
void scal<s_type,s_type>(const i_type n, const s_type a, s_type *x, const i_type incx)
{
    BLAS_NAME(sscal)(_rc(&n),_rc(&a),_rc(x),_rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE 
void scal<d_type,d_type>(const i_type n, const d_type a, d_type *x, const i_type incx)
{
    BLAS_NAME(dscal)(_rc(&n),_rc(&a),_rc(x),_rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE 
void scal<c_type,c_type>(const i_type n, const c_type a, c_type *x, const i_type incx)
{
    BLAS_NAME(cscal)(_rc(&n),_rc(&a),_rc(x),_rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE 
void scal<z_type,z_type>(const i_type n, const z_type a, z_type *x, const i_type incx)
{
    BLAS_NAME(zscal)(_rc(&n),_rc(&a),_rc(x),_rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE 
void scal<s_type,c_type>(const i_type n, const s_type a, c_type *x, const i_type incx)
{
    BLAS_NAME(csscal)(_rc(&n),_rc(&a),_rc(x),_rc(&incx));
}

template<> BLAS_EXPORT INLINE_TYPE 
void scal<d_type,z_type>(const i_type n, const d_type a, z_type *x, const i_type incx)
{
    BLAS_NAME(zdscal)(_rc(&n),_rc(&a),_rc(x),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE 
void sscal (i_type n, s_type a, s_type *x, i_type incx)
{
    BLAS_NAME(sscal)(_rc(&n),_rc(&a),_rc(x),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE 
void cscal (i_type n, c_type a, c_type *x, i_type incx)
{
    BLAS_NAME(cscal)(_rc(&n),_rc(&a),_rc(x),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE 
void csscal(i_type n, s_type a, c_type *x, i_type incx)
{
    BLAS_NAME(csscal)(_rc(&n),_rc(&a),_rc(x),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE 
void dscal (i_type n, d_type a, d_type *x, i_type incx)
{
    BLAS_NAME(dscal)(_rc(&n),_rc(&a),_rc(x),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE 
void zdscal(i_type n, d_type a, z_type *x, i_type incx)
{
    BLAS_NAME(zdscal)(_rc(&n),_rc(&a),_rc(x),_rc(&incx));
}; 

BLAS_EXPORT INLINE_TYPE 
void zscal (i_type n, z_type a, z_type *x, i_type incx)
{
    BLAS_NAME(zscal)(_rc(&n),_rc(&a),_rc(x),_rc(&incx));
}; 

//-----------------------------------------------------------------
//                          AXPY
//-----------------------------------------------------------------
template<> BLAS_EXPORT INLINE_TYPE 
void axpy<s_type>(i_type n, s_type alpha, const s_type *x, i_type incx, s_type *y, 
                  i_type incy)
{
    BLAS_NAME(saxpy)(_rc(&n), _rc(&alpha), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void axpy<d_type>(i_type n, d_type alpha, const d_type *x, i_type incx, d_type *y, 
                  i_type incy)
{
    BLAS_NAME(daxpy)(_rc(&n), _rc(&alpha), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void axpy<c_type>(i_type n, c_type alpha, const c_type *x, i_type incx, c_type *y, 
                  i_type incy)
{
    BLAS_NAME(caxpy)(_rc(&n), _rc(&alpha), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void axpy<z_type>(i_type n, z_type alpha, const z_type *x, i_type incx, z_type *y,
                  i_type incy)
{
    BLAS_NAME(zaxpy)(_rc(&n), _rc(&alpha), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
};

BLAS_EXPORT INLINE_TYPE 
void saxpy(i_type n, s_type alpha, const s_type *x, i_type incx, s_type *y, i_type incy)
{
    BLAS_NAME(saxpy)(_rc(&n), _rc(&alpha), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
};

BLAS_EXPORT INLINE_TYPE 
void caxpy(i_type n, c_type alpha, const c_type *x, i_type incx, c_type *y, i_type incy)
{
    BLAS_NAME(caxpy)(_rc(&n), _rc(&alpha), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
}; 

BLAS_EXPORT INLINE_TYPE
void daxpy(i_type n, d_type alpha, const d_type *x, i_type incx, d_type *y, i_type incy)
{
    BLAS_NAME(daxpy)(_rc(&n), _rc(&alpha), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
void zaxpy(i_type n, z_type alpha, const z_type *x, i_type incx, z_type *y, i_type incy)
{
    BLAS_NAME(zaxpy)(_rc(&n), _rc(&alpha), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
}; 

//-----------------------------------------------------------------
//                          AMAX
//-----------------------------------------------------------------
template<> BLAS_EXPORT INLINE_TYPE 
i_type amax<s_type>(i_type n, const s_type *x, i_type incx)
{
    return BLAS_NAME(isamax)(_rc(&n), _rc(x), _rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE 
i_type amax<d_type>(i_type n, const d_type *x, i_type incx)
{
    i_type ret = BLAS_NAME(idamax)(_rc(&n), _rc(x), _rc(&incx));
    return ret;
};

template<> BLAS_EXPORT INLINE_TYPE 
i_type amax<c_type>(i_type n, const c_type *x, i_type incx)
{
    return BLAS_NAME(icamax)(_rc(&n), _rc(x), _rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE 
i_type amax<z_type>(i_type n, const z_type *x, i_type incx)
{
    return BLAS_NAME(izamax)(_rc(&n), _rc(x), _rc(&incx));
};

BLAS_EXPORT INLINE_TYPE
i_type  idamax(i_type n, const d_type *x, i_type incx)
{
    return BLAS_NAME(idamax)(_rc(&n),_rc(x),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE
i_type  isamax(i_type n, const s_type *x, i_type incx)
{
    return BLAS_NAME(isamax)(_rc(&n),_rc(x),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE
i_type  icamax(i_type n, const c_type *x, i_type incx)
{
    return BLAS_NAME(icamax)(_rc(&n),_rc(x),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE
i_type  izamax(i_type n, const z_type *x, i_type incx)
{
    return BLAS_NAME(izamax)(_rc(&n),_rc(x),_rc(&incx));
};

//-----------------------------------------------------------------
//                          SWAP
//-----------------------------------------------------------------
template<> BLAS_EXPORT INLINE_TYPE 
void swap<s_type>(i_type n, s_type *x, i_type incx, s_type *y, i_type incy)
{
    BLAS_NAME(sswap)(_rc(&n), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void swap<d_type>(i_type n, d_type *x, i_type incx, d_type *y, i_type incy)
{
    BLAS_NAME(dswap)(_rc(&n), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void swap<c_type>(i_type n, c_type *x, i_type incx, c_type *y, i_type incy)
{
    BLAS_NAME(cswap)(_rc(&n), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void swap<z_type>(i_type n, z_type *x, i_type incx, z_type *y, i_type incy)
{
    BLAS_NAME(zswap)(_rc(&n), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
void   sswap(i_type n, s_type *x, i_type incx, s_type *y, i_type incy)
{
    BLAS_NAME(sswap)(_rc(&n), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
void   cswap(i_type n, c_type *x, i_type incx, c_type *y, i_type incy)
{
    BLAS_NAME(cswap)(_rc(&n), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
}; 

BLAS_EXPORT INLINE_TYPE
void   dswap(i_type n, d_type *x, i_type incx, d_type *y, i_type incy)
{
    BLAS_NAME(dswap)(_rc(&n), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
void   zswap(i_type n, z_type *x, i_type incx, z_type *y, i_type incy)
{
    BLAS_NAME(zswap)(_rc(&n), _rc(x), _rc(&incx), _rc(y), _rc(&incy));
}; 

//-----------------------------------------------------------------
//                          COPY
//-----------------------------------------------------------------
template<> BLAS_EXPORT INLINE_TYPE 
void copy<s_type>(i_type n, const s_type *x, i_type incx, s_type *y, i_type incy)
{
    BLAS_NAME(scopy)(_rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void copy<d_type>(i_type n, const d_type *x, i_type incx, d_type *y, i_type incy)
{
    BLAS_NAME(dcopy)(_rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void copy<c_type>(i_type n, const c_type *x, i_type incx, c_type *y, i_type incy)
{
    BLAS_NAME(ccopy)(_rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void copy<z_type>(i_type n, const z_type *x, i_type incx, z_type *y, i_type incy)
{
    BLAS_NAME(zcopy)(_rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
void   scopy(i_type n, const s_type *x, i_type incx, s_type *y, i_type incy)
{
    BLAS_NAME(scopy)(_rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
void   ccopy(i_type n, const c_type *x, i_type incx, c_type *y, i_type incy)
{
    BLAS_NAME(ccopy)(_rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
void   dcopy(i_type n, const d_type *x, i_type incx, d_type *y, i_type incy)
{
    BLAS_NAME(dcopy)(_rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
void   zcopy(i_type n, const z_type *x, i_type incx, z_type *y, i_type incy)
{
    BLAS_NAME(zcopy)(_rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
};

//-----------------------------------------------------------------
//                          DOT
//-----------------------------------------------------------------
template<> BLAS_EXPORT INLINE_TYPE
s_type dot<s_type>(i_type n, const s_type *x, i_type incx, const s_type *y, i_type incy)
{
    return sdot(n,x,incx,y,incy);
};

template<> BLAS_EXPORT INLINE_TYPE
d_type dot<d_type>(i_type n, const d_type *x, i_type incx, const d_type *y, i_type incy)
{
    return ddot(n,x,incx,y,incy);
};

template<> BLAS_EXPORT INLINE_TYPE
c_type dot<c_type>(i_type n, const c_type *x, i_type incx, const c_type *y, i_type incy)
{
    return cdotu(n,x,incx,y,incy);
};

template<> BLAS_EXPORT INLINE_TYPE
z_type dot<z_type>(i_type n, const z_type *x, i_type incx, const z_type *y, i_type incy)
{
    return zdotu(n,x,incx,y,incy);
};

BLAS_EXPORT INLINE_TYPE
s_type sdot(i_type n, const s_type *x, i_type incx, const s_type *y, i_type incy)
{
    return BLAS_NAME(sdot)(_rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
d_type ddot(i_type n, const d_type *x, i_type incx, const d_type *y, i_type incy)
{
    return BLAS_NAME(ddot)(_rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
c_type cdotu(i_type n, const c_type *x, i_type incx, const c_type *y, i_type incy)
{
    c_type out;
    BLAS_NAME(cdotu)(_rc(&out), _rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
    return out;
};

BLAS_EXPORT INLINE_TYPE
z_type zdotu(i_type n, const z_type *x, i_type incx, const z_type *y, i_type incy)
{
    z_type out;
    BLAS_NAME(zdotu)(_rc(&out),_rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
    return out;
};

//-----------------------------------------------------------------
//                          ROT
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void lapack::rot<s_type>(i_type n, s_type *x, i_type incx, s_type *y, i_type incy, 
                         s_type C, s_type S)
{
    rot2<s_type>::eval(n, x, incx, y, incy, C, S);    
};

template<> BLAS_EXPORT INLINE_TYPE
void lapack::rot<d_type>(i_type n, d_type *x, i_type incx, d_type *y, i_type incy, 
                         d_type C, d_type S)
{
    rot2<d_type>::eval(n, x, incx, y, incy, C, S);
    //BLAS_NAME(drot)(&n,x,&incx,y,&incy,&C,&S);
};

template<> BLAS_EXPORT INLINE_TYPE
void lapack::rot<c_type>(i_type n, c_type *x, i_type incx, c_type *y, i_type incy, 
                         s_type C, c_type S)
{
    BLAS_NAME(crot)(&n,_rc(x),&incx,_rc(y),&incy,&C,_rc(&S));
};

template<> BLAS_EXPORT INLINE_TYPE
void lapack::rot<z_type>(i_type n, z_type *x, i_type incx, z_type *y, i_type incy, 
                         d_type C, z_type S)
{
    BLAS_NAME(zrot)(&n,_rc(x),&incx,_rc(y),&incy,&C,_rc(&S));
};

BLAS_EXPORT INLINE_TYPE
void lapack::srot(i_type n, s_type *x, i_type incx, s_type *y, i_type incy, s_type C,
                  s_type S)
{
    rot2<s_type>::eval(n, x, incx, y, incy, C, S);
    //BLAS_NAME(srot)(&n,x,&incx,y,&incy,&C,&S);
};

BLAS_EXPORT INLINE_TYPE
void lapack::drot(i_type n, d_type *x, i_type incx, d_type *y, i_type incy, d_type C,
                  d_type S)
{
    rot2<d_type>::eval(n, x, incx, y, incy, C, S);
    //BLAS_NAME(drot)(&n,x,&incx,y,&incy,&C,&S);
};

BLAS_EXPORT INLINE_TYPE
void lapack::crot(i_type n, c_type *x, i_type incx, c_type *y, i_type incy, s_type C, 
                  c_type S)
{
    BLAS_NAME(crot)(&n,_rc(x),&incx,_rc(y),&incy,_rc(&C),_rc(&S));
};

BLAS_EXPORT INLINE_TYPE
void lapack::zrot(i_type n, z_type *x, i_type incx, z_type *y, i_type incy, d_type C, 
                  z_type S)
{
    BLAS_NAME(zrot)(&n,_rc(x),&incx,_rc(y),&incy,&C,_rc(&S));
};

//-----------------------------------------------------------------
//                          DOT
//-----------------------------------------------------------------
template<> BLAS_EXPORT INLINE_TYPE
s_type dotc<s_type>(i_type n, const s_type *x, i_type incx, const s_type *y, i_type incy)
{
    return sdot(n,x,incx,y,incy);
};

template<> BLAS_EXPORT INLINE_TYPE
d_type dotc<d_type>(i_type n, const d_type *x, i_type incx, const d_type *y, i_type incy)
{
    return ddot(n,x,incx,y,incy);
};

template<> BLAS_EXPORT INLINE_TYPE
c_type dotc<c_type>(i_type n, const c_type *x, i_type incx, const c_type *y, i_type incy)
{
    return cdotc(n,x,incx,y,incy);
};

template<> BLAS_EXPORT INLINE_TYPE
z_type dotc<z_type>(i_type n, const z_type *x, i_type incx, const z_type *y, i_type incy)
{
    return zdotc(n,x,incx,y,incy);
};

BLAS_EXPORT INLINE_TYPE
c_type cdotc(i_type n, const c_type *x, i_type incx, const c_type *y, i_type incy)
{
    c_type out;
    BLAS_NAME(cdotc)(_rc(&out), _rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
    return out;
};

BLAS_EXPORT INLINE_TYPE
z_type zdotc(i_type n, const z_type *x, i_type incx, const z_type *y, i_type incy)
{
    z_type out;
    BLAS_NAME(zdotc)(_rc(&out),_rc(&n),_rc(x),_rc(&incx),_rc(y),_rc(&incy));
    return out;
};

//-----------------------------------------------------------------
//                          GER
//-----------------------------------------------------------------
template<> BLAS_EXPORT INLINE_TYPE 
void geru<s_type>(i_type m, i_type n, s_type alpha, const s_type *x, i_type incx, 
                  const s_type *y, i_type incy, s_type *a, i_type lda)
{
    BLAS_NAME(sger)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                    _rc(a),_rc(&lda));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gerc<s_type>(i_type m, i_type n, s_type alpha, const s_type *x, i_type incx, 
                  const s_type *y, i_type incy, s_type *a, i_type lda)
{
    BLAS_NAME(sger)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                    _rc(a),_rc(&lda));
};

template<> BLAS_EXPORT INLINE_TYPE 
void geru<d_type>(i_type m, i_type n, d_type alpha, const d_type *x, i_type incx, 
                  const d_type *y, i_type incy, d_type *a, i_type lda)
{
    BLAS_NAME(dger)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                    _rc(a),_rc(&lda));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gerc<d_type>(i_type m, i_type n, d_type alpha, const d_type *x, i_type incx, 
                  const d_type *y, i_type incy, d_type *a, i_type lda)
{
    BLAS_NAME(dger)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                    _rc(a),_rc(&lda));
};

template<> BLAS_EXPORT INLINE_TYPE 
void geru<c_type>(i_type m, i_type n, c_type alpha, const c_type *x, i_type incx,
                  const c_type *y, i_type incy, c_type *a, i_type lda)
{
    BLAS_NAME(cgeru)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                     _rc(a),_rc(&lda));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gerc<c_type>(i_type m, i_type n, c_type alpha, const c_type *x, i_type incx, 
                  const c_type *y, i_type incy, c_type *a, i_type lda)
{
    BLAS_NAME(cgerc)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                     _rc(a),_rc(&lda));
};

template<> BLAS_EXPORT INLINE_TYPE 
void geru<z_type>(i_type m, i_type n, z_type alpha, const z_type *x, i_type incx, 
                  const z_type *y, i_type incy, z_type *a, i_type lda)
{
    BLAS_NAME(zgeru)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                     _rc(a),_rc(&lda));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gerc<z_type>(i_type m, i_type n, z_type alpha, const z_type *x, i_type incx, 
                  const z_type *y, i_type incy, z_type *a, i_type lda)
{
    BLAS_NAME(zgerc)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                     _rc(a),_rc(&lda));
};

BLAS_EXPORT INLINE_TYPE 
void sger(i_type m, i_type n, s_type alpha, const s_type *x, i_type incx, 
          const s_type *y, i_type incy, s_type *a, i_type lda)
{
    BLAS_NAME(sger)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                    _rc(a),_rc(&lda));
};

BLAS_EXPORT INLINE_TYPE
void dger(i_type m, i_type n, d_type alpha, const d_type *x, i_type incx, 
          const d_type *y, i_type incy, d_type *a, i_type lda)
{
    BLAS_NAME(dger)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                    _rc(a),_rc(&lda));
};

BLAS_EXPORT INLINE_TYPE
void cgerc(i_type m, i_type n, c_type alpha, const c_type *x, i_type incx, 
           const c_type *y, i_type incy, c_type *a, i_type lda)
{
    BLAS_NAME(cgerc)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                     _rc(a),_rc(&lda));
};

BLAS_EXPORT INLINE_TYPE
void cgeru(i_type m, i_type n, c_type alpha, const c_type *x, i_type incx, 
           const c_type *y, i_type incy, c_type *a, i_type lda)
{
    BLAS_NAME(cgeru)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                     _rc(a),_rc(&lda));
};

BLAS_EXPORT INLINE_TYPE
void zgerc(i_type m, i_type n, z_type alpha, const z_type *x, i_type incx, 
           const z_type *y, i_type incy, z_type *a, i_type lda)
{
    BLAS_NAME(zgerc)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                     _rc(a),_rc(&lda));
};

BLAS_EXPORT INLINE_TYPE
void zgeru(i_type m, i_type n, z_type alpha, const z_type *x, i_type incx, 
           const z_type *y, i_type incy, z_type *a, i_type lda)
{
    BLAS_NAME(zgeru)(_rc(&m),_rc(&n),_rc(&alpha),_rc(x),_rc(&incx),_rc(y),_rc(&incy),
                     _rc(a),_rc(&lda));
};

//-----------------------------------------------------------------
//                          GEMV
//-----------------------------------------------------------------
template<> BLAS_EXPORT INLINE_TYPE
void gemv<s_type>(const char *trans, i_type m, i_type n, s_type alpha, const s_type *a, 
                  i_type lda, const s_type *x, i_type incx, s_type beta, s_type *y, 
                  i_type incy)
{
    BLAS_NAME(sgemv)(_rc(trans), _rc(&m), _rc(&n), _rc(&alpha), _rc(a), _rc(&lda), 
                     _rc(x), _rc(&incx), _rc(&beta), _rc(y), _rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE
void gemv<d_type>(const char *trans, i_type m, i_type n, d_type alpha, const d_type *a,
                  i_type lda, const d_type *x, i_type incx, d_type beta, d_type *y, 
                  i_type incy)
{
    BLAS_NAME(dgemv)(_rc(trans), _rc(&m), _rc(&n), _rc(&alpha), _rc(a), _rc(&lda), 
                     _rc(x), _rc(&incx), _rc(&beta), _rc(y), _rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE
void gemv<c_type>(const char *trans, i_type m, i_type n, c_type alpha, const c_type *a,
                  i_type lda, const c_type *x, i_type incx, c_type beta, c_type *y, 
                  i_type incy)
{
    BLAS_NAME(cgemv)(_rc(trans), _rc(&m), _rc(&n), _rc(&alpha), _rc(a), _rc(&lda), 
                     _rc(x), _rc(&incx), _rc(&beta), _rc(y), _rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE
void gemv<z_type>(const char *trans, i_type m, i_type n, z_type alpha, const z_type *a, 
                  i_type lda, const z_type *x, i_type incx, z_type beta, z_type *y, 
                  i_type incy)
{
    BLAS_NAME(zgemv)(_rc(trans), _rc(&m), _rc(&n), _rc(&alpha), _rc(a), _rc(&lda), 
                     _rc(x), _rc(&incx), _rc(&beta), _rc(y), _rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
void sgemv(const char *trans, i_type m, i_type n, s_type alpha, const s_type *a, 
           i_type lda, const s_type *x, i_type incx, s_type beta, s_type *y, 
           i_type incy)
{
    BLAS_NAME(sgemv)(_rc(trans), _rc(&m), _rc(&n), _rc(&alpha), _rc(a), _rc(&lda), 
                     _rc(x), _rc(&incx), _rc(&beta), _rc(y), _rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
void cgemv(const char *trans, i_type m, i_type n, c_type alpha, const c_type *a, 
           i_type lda, const c_type *x, i_type incx, c_type beta, c_type *y, 
           i_type incy)
{
    BLAS_NAME(cgemv)(_rc(trans), _rc(&m), _rc(&n), _rc(&alpha), _rc(a), _rc(&lda), 
                     _rc(x), _rc(&incx), _rc(&beta), _rc(y), _rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
void dgemv(const char *trans, i_type m, i_type n, d_type alpha, const d_type *a, 
           i_type lda, const d_type *x, i_type incx, d_type beta, d_type *y, 
           i_type incy)
{
    BLAS_NAME(dgemv)(_rc(trans), _rc(&m), _rc(&n), _rc(&alpha), _rc(a), _rc(&lda), 
                     _rc(x), _rc(&incx), _rc(&beta), _rc(y), _rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
void zgemv(const char *trans, i_type m, i_type n, z_type alpha, const z_type *a, 
           i_type lda, const z_type *x, i_type incx, z_type beta, z_type *y,
           i_type incy)
{
    BLAS_NAME(zgemv)(_rc(trans), _rc(&m), _rc(&n), _rc(&alpha), _rc(a), _rc(&lda), 
                     _rc(x), _rc(&incx), _rc(&beta), _rc(y), _rc(&incy));
};

//-----------------------------------------------------------------
//                          TRMV
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void trmv<s_type>(const char *uplo, const char *transa, const char *diag, i_type n, 
                          const s_type *a, i_type lda, s_type *b, i_type incx)
{
    BLAS_NAME(strmv)(_rc(uplo),_rc(transa),_rc(diag),_rc(&n),_rc(a),_rc(&lda),
                     _rc(b),_rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE
void trmv<d_type>(const char *uplo, const char *transa, const char *diag, i_type n, 
                          const d_type *a, i_type lda, d_type *b, i_type incx)
{
    BLAS_NAME(dtrmv)(_rc(uplo),_rc(transa),_rc(diag),_rc(&n),_rc(a),_rc(&lda),
                     _rc(b),_rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE
void trmv<c_type>(const char *uplo, const char *transa, const char *diag, i_type n, 
                          const c_type *a, i_type lda, c_type *b, i_type incx)
{
    BLAS_NAME(ctrmv)(_rc(uplo),_rc(transa),_rc(diag),_rc(&n),_rc(a),_rc(&lda),
                     _rc(b),_rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE
void trmv<z_type>(const char *uplo, const char *transa, const char *diag, i_type n, 
                          const z_type *a, i_type lda, z_type *b, i_type incx)
{
    BLAS_NAME(ztrmv)(_rc(uplo),_rc(transa),_rc(diag),_rc(&n),_rc(a),_rc(&lda),
                     _rc(b),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE 
void strmv(const char *uplo, const char *transa, const char *diag, i_type n, 
           const s_type *a, i_type lda, s_type *b, i_type incx)
{
    BLAS_NAME(strmv)(_rc(uplo),_rc(transa),_rc(diag),_rc(&n),_rc(a),_rc(&lda),
                     _rc(b),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE 
void ctrmv(const char *uplo, const char *transa, const char *diag, i_type n, 
           const c_type *a, i_type lda, c_type *b, i_type incx)
{
    BLAS_NAME(ctrmv)(_rc(uplo),_rc(transa),_rc(diag),_rc(&n),_rc(a),_rc(&lda),
                     _rc(b),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE 
void dtrmv(const char *uplo, const char *transa, const char *diag, i_type n, 
           const d_type *a, i_type lda, d_type *b, i_type incx)
{
    BLAS_NAME(dtrmv)(_rc(uplo),_rc(transa),_rc(diag),_rc(&n),_rc(a),_rc(&lda),
                     _rc(b),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE 
void ztrmv(const char *uplo, const char *transa, const char *diag, i_type n, 
           const z_type *a, i_type lda, z_type *b, i_type incx)
{
    BLAS_NAME(ztrmv)(_rc(uplo),_rc(transa),_rc(diag),_rc(&n),_rc(a),_rc(&lda),
                     _rc(b),_rc(&incx));
};

//-----------------------------------------------------------------
//                          TRSM
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void trsm<s_type>(const char *side, const char *uplo, const char *transa, 
                  const char *diag, i_type m, i_type n, s_type alpha, 
                  const s_type *a, i_type lda, s_type *b, i_type ldb)
{
    BLAS_NAME(strsm)(_rc(side), _rc(uplo), _rc(transa), _rc(diag), _rc(&m), 
                     _rc(&n), _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb));
};

template<> BLAS_EXPORT INLINE_TYPE 
void trsm<d_type>(const char *side, const char *uplo, const char *transa, 
                  const char *diag, i_type m, i_type n, d_type alpha, const d_type *a,
                  i_type lda, d_type *b, i_type ldb)
{
    BLAS_NAME(dtrsm)(_rc(side), _rc(uplo), _rc(transa), _rc(diag), _rc(&m), _rc(&n),
                     _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb));
};

template<> BLAS_EXPORT INLINE_TYPE 
void trsm<c_type>(const char *side, const char *uplo, const char *transa, 
                  const char *diag, i_type m, i_type n, c_type alpha, const c_type *a,
                  i_type lda, c_type *b, i_type ldb)
{
    BLAS_NAME(ctrsm)(_rc(side), _rc(uplo), _rc(transa), _rc(diag), _rc(&m), _rc(&n),
                     _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb));
};

template<> BLAS_EXPORT INLINE_TYPE 
void trsm<z_type>(const char *side, const char *uplo, const char *transa, 
                  const char *diag, i_type m, i_type n, z_type alpha, const z_type *a, 
                  i_type lda, z_type *b, i_type ldb)
{
    BLAS_NAME(ztrsm)(_rc(side), _rc(uplo), _rc(transa), _rc(diag), _rc(&m), _rc(&n), 
                     _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb));
};

BLAS_EXPORT INLINE_TYPE
void strsm(const char *side, const char *uplo, const char *transa, const char *diag,
            i_type m, i_type n, s_type alpha, const s_type *a, i_type lda, s_type *b, 
           i_type ldb)
{
    BLAS_NAME(strsm)(_rc(side), _rc(uplo), _rc(transa), _rc(diag), _rc(&m), _rc(&n), 
                     _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb));
};

BLAS_EXPORT INLINE_TYPE
void ctrsm(const char *side, const char *uplo, const char *transa, const char *diag,
            i_type m, i_type n, c_type alpha, const c_type *a, i_type lda, c_type *b, 
           i_type ldb)
{
    BLAS_NAME(ctrsm)(_rc(side), _rc(uplo), _rc(transa), _rc(diag), _rc(&m), _rc(&n), 
                     _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb));
};

BLAS_EXPORT INLINE_TYPE
void dtrsm(const char *side, const char *uplo, const char *transa, const char *diag,
            i_type m, i_type n, d_type alpha, const d_type *a, i_type lda, d_type *b, 
           i_type ldb)
{
    BLAS_NAME(dtrsm)(_rc(side), _rc(uplo), _rc(transa), _rc(diag), _rc(&m), _rc(&n), 
                     _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb));
};

BLAS_EXPORT INLINE_TYPE
void ztrsm(const char *side, const char *uplo, const char *transa, const char *diag,
            i_type m, i_type n, z_type alpha, const z_type *a, i_type lda, z_type *b, 
           i_type ldb)
{
    BLAS_NAME(ztrsm)(_rc(side), _rc(uplo), _rc(transa), _rc(diag), _rc(&m), _rc(&n), 
                     _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb));
};

//-----------------------------------------------------------------
//                          GEMM
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void gemm<s_type>(const char *transa, const char *transb, i_type m, i_type n, i_type k,
                  s_type alpha, const s_type *a, i_type lda, const s_type *b, i_type ldb, 
                  s_type beta, s_type *c, i_type ldc)
{
    BLAS_NAME(sgemm)(_rc(transa), _rc(transb), _rc(&m), _rc(&n), _rc(&k),
           _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb), _rc(&beta), _rc(c), 
                     _rc(&ldc));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gemm<d_type>(const char *transa, const char *transb, i_type m, i_type n, i_type k,
                  d_type alpha, const d_type *a, i_type lda, const d_type *b, i_type ldb, 
                  d_type beta, d_type *c, i_type ldc)
{
    BLAS_NAME(dgemm)(_rc(transa), _rc(transb), _rc(&m), _rc(&n), _rc(&k),
           _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb), _rc(&beta), _rc(c), _rc(&ldc));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gemm<c_type>(const char *transa, const char *transb, i_type m, i_type n, i_type k,
                  c_type alpha, const c_type *a, i_type lda, const c_type *b, i_type ldb, 
                  c_type beta, c_type *c, i_type ldc)
{
    BLAS_NAME(cgemm)(_rc(transa), _rc(transb), _rc(&m), _rc(&n), _rc(&k),
           _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb), _rc(&beta), _rc(c), _rc(&ldc));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gemm<z_type>(const char *transa, const char *transb, i_type m, i_type n, i_type k,
                  z_type alpha, const z_type *a, i_type lda, const z_type *b, i_type ldb,
                  z_type beta, z_type *c, i_type ldc)
{
    BLAS_NAME(zgemm)(_rc(transa), _rc(transb), _rc(&m), _rc(&n), _rc(&k),
           _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb), _rc(&beta), _rc(c), _rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE
void sgemm(const char *transa, const char *transb, i_type m, i_type n, i_type k,
           s_type alpha, const s_type *a, i_type lda, const s_type *b, i_type ldb,
           s_type beta, s_type *c, i_type ldc)
{
    BLAS_NAME(sgemm)(_rc(transa), _rc(transb), _rc(&m), _rc(&n), _rc(&k),
           _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb), _rc(&beta), _rc(c), _rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE
void cgemm(const char *transa, const char *transb, i_type m, i_type n, i_type k,
           c_type alpha, const c_type *a, i_type lda, const c_type *b, i_type ldb, 
           c_type beta, c_type *c, i_type ldc)
{
    BLAS_NAME(cgemm)(_rc(transa), _rc(transb), _rc(&m), _rc(&n), _rc(&k),
           _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb), _rc(&beta), _rc(c), _rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE
void dgemm(const char *transa, const char *transb, i_type m, i_type n, i_type k,
           d_type alpha, const d_type *a, i_type lda, const d_type *b, i_type ldb,
           d_type beta, d_type *c, i_type ldc)
{
    BLAS_NAME(dgemm)(_rc(transa), _rc(transb), _rc(&m), _rc(&n), _rc(&k),
           _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb), _rc(&beta), _rc(c), _rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE
void zgemm(const char *transa, const char *transb, i_type m, i_type n, i_type k,
           z_type alpha, const z_type *a, i_type lda, const z_type *b, i_type ldb, 
           z_type beta, z_type *c, i_type ldc)
{
    BLAS_NAME(zgemm)(_rc(transa), _rc(transb), _rc(&m), _rc(&n), _rc(&k),
           _rc(&alpha), _rc(a), _rc(&lda), _rc(b), _rc(&ldb), _rc(&beta), _rc(c), _rc(&ldc));
};

//-----------------------------------------------------------------
//                          HERK
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void herk<s_type>(const char *uplo, const char *trans, i_type n, i_type k, s_type alpha, 
                  const s_type *a, i_type lda, s_type beta, s_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'C' || trans[0] == 'c')
    {
        trans2 = 'T';
    };

    BLAS_NAME(ssyrk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                _rc(&beta), _rc(c), _rc(&ldc));
};

template<> BLAS_EXPORT INLINE_TYPE 
void herk<d_type>(const char *uplo, const char *trans, i_type n, i_type k, d_type alpha, 
                  const d_type *a, i_type lda, d_type beta, d_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'C' || trans[0] == 'c')
    {
        trans2 = 'T';
    };

    BLAS_NAME(dsyrk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                _rc(&beta), _rc(c), _rc(&ldc));
};

template<> BLAS_EXPORT INLINE_TYPE 
void herk<c_type>(const char *uplo, const char *trans, i_type n, i_type k, s_type alpha,
                  const c_type *a, i_type lda, s_type beta, c_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'T' || trans[0] == 't')
    {
        trans2 = 'C';
    };

    BLAS_NAME(cherk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                    _rc(&beta), _rc(c), _rc(&ldc));
};

template<> BLAS_EXPORT INLINE_TYPE 
void herk<z_type>(const char *uplo, const char *trans, i_type n, i_type k, d_type alpha, 
                  const z_type *a, i_type lda, d_type beta, z_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'T' || trans[0] == 't')
    {
        trans2 = 'C';
    };

    BLAS_NAME(zherk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                    _rc(&beta), _rc(c), _rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE
void cherk(const char *uplo, const char *trans, i_type n, i_type k, s_type alpha, 
           const c_type *a, i_type lda, s_type beta, c_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'T' || trans[0] == 't')
    {
        trans2 = 'C';
    };

    BLAS_NAME(cherk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                    _rc(&beta), _rc(c), _rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE
void zherk(const char *uplo, const char *trans, i_type n, i_type k, d_type alpha, 
           const z_type *a, i_type lda, d_type beta, z_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'T' || trans[0] == 't')
    {
        trans2 = 'C';
    };

    BLAS_NAME(zherk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                    _rc(&beta), _rc(c), _rc(&ldc));
};

//-----------------------------------------------------------------
//                          syrk
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void syrk<s_type>(const char *uplo, const char *trans, i_type n, i_type k, 
                  const s_type& alpha, const s_type *a, i_type lda, const s_type& beta,
                  s_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'C' || trans[0] == 'c')
    {
        trans2 = 'T';
    };

    BLAS_NAME(ssyrk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                _rc(&beta), _rc(c), _rc(&ldc));
};

template<> BLAS_EXPORT INLINE_TYPE 
void syrk<d_type>(const char *uplo, const char *trans, i_type n, i_type k, 
                  const d_type& alpha, const d_type *a, i_type lda, const d_type& beta, 
                  d_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'C' || trans[0] == 'c')
    {
        trans2 = 'T';
    };

    BLAS_NAME(dsyrk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                _rc(&beta), _rc(c), _rc(&ldc));
};

template<> BLAS_EXPORT INLINE_TYPE 
void syrk<c_type>(const char *uplo, const char *trans, i_type n, i_type k, 
                  const c_type& alpha, const c_type *a, i_type lda, const c_type& beta, 
                  c_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'C' || trans[0] == 'c')
    {
        trans2 = 'T';
    };

    BLAS_NAME(csyrk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                _rc(&beta), _rc(c), _rc(&ldc));
};

template<> BLAS_EXPORT INLINE_TYPE 
void syrk<z_type>(const char *uplo, const char *trans, i_type n, i_type k, 
                  const z_type& alpha, const z_type *a, i_type lda, const z_type& beta, 
                  z_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'C' || trans[0] == 'c')
    {
        trans2 = 'T';
    };

    BLAS_NAME(zsyrk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                _rc(&beta), _rc(c), _rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE
void ssyrk(const char *uplo, const char *trans, i_type n, i_type k, 
           const s_type& alpha, const s_type *a, i_type lda, const s_type& beta,
           s_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'C' || trans[0] == 'c')
    {
        trans2 = 'T';
    };

    BLAS_NAME(ssyrk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                    _rc(&beta), _rc(c), _rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE
void dsyrk(const char *uplo, const char *trans, i_type n, i_type k, 
           const d_type& alpha, const d_type *a, i_type lda, const d_type& beta, 
           d_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'C' || trans[0] == 'c')
    {
        trans2 = 'T';
    };

    BLAS_NAME(dsyrk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                    _rc(&beta), _rc(c), _rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE
void csyrk(const char *uplo, const char *trans, i_type n, i_type k, 
           const c_type& alpha, const c_type *a, i_type lda, const c_type& beta, 
           c_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'C' || trans[0] == 'c')
    {
        trans2 = 'T';
    };

    BLAS_NAME(csyrk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                    _rc(&beta), _rc(c), _rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE
void zsyrk(const char *uplo, const char *trans, i_type n, i_type k, 
           const z_type& alpha, const z_type *a, i_type lda, const z_type& beta,
           z_type* c, i_type ldc)
{
    char trans2 = trans[0];
    if (trans[0] == 'C' || trans[0] == 'c')
    {
        trans2 = 'T';
    };

    BLAS_NAME(zsyrk)(_rc(uplo), _rc(&trans2), _rc(&n), _rc(&k), _rc(&alpha), _rc(a), _rc(&lda), 
                    _rc(&beta), _rc(c), _rc(&ldc));
};

//-----------------------------------------------------------------
//                          TRMM
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void trmm<s_type>(const char *side, const char *uplo, const char *transa, 
                  const char *diag, i_type m, i_type n, s_type alpha, const s_type *a, 
                  i_type lda, s_type *b, i_type ldb)
{
    BLAS_NAME(strmm)(_rc(side),_rc(uplo),_rc(transa),_rc(diag),_rc(&m),_rc(&n),
                     _rc(&alpha),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

template<> BLAS_EXPORT INLINE_TYPE 
void trmm<d_type>(const char *side, const char *uplo, const char *transa, 
                  const char *diag, i_type m, i_type n, d_type alpha, const d_type *a,
                  i_type lda, d_type *b, i_type ldb)
{
    BLAS_NAME(dtrmm)(_rc(side),_rc(uplo),_rc(transa),_rc(diag),_rc(&m),_rc(&n),
                     _rc(&alpha),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

template<> BLAS_EXPORT INLINE_TYPE 
void trmm<c_type>(const char *side, const char *uplo, const char *transa, 
                  const char *diag, i_type m, i_type n, c_type alpha, const c_type *a,
                  i_type lda, c_type *b, i_type ldb)
{
    BLAS_NAME(ctrmm)(_rc(side),_rc(uplo),_rc(transa),_rc(diag),_rc(&m),_rc(&n),
                     _rc(&alpha),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

template<> BLAS_EXPORT INLINE_TYPE 
void trmm<z_type>(const char *side, const char *uplo, const char *transa, 
                  const char *diag, i_type m, i_type n, z_type alpha, const z_type *a,
                  i_type lda, z_type *b, i_type ldb)
{
    BLAS_NAME(ztrmm)(_rc(side),_rc(uplo),_rc(transa),_rc(diag),_rc(&m),_rc(&n),
                     _rc(&alpha),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

BLAS_EXPORT INLINE_TYPE 
void strmm(const char *side, const char *uplo, const char *transa,
           const char *diag, i_type m, i_type n, s_type alpha, const s_type *a, 
           i_type lda, s_type *b, i_type ldb)
{
    BLAS_NAME(strmm)(_rc(side),_rc(uplo),_rc(transa),_rc(diag),_rc(&m),_rc(&n),
                     _rc(&alpha),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

BLAS_EXPORT INLINE_TYPE 
void ctrmm(const char *side, const char *uplo, const char *transa, 
           const char *diag, i_type m, i_type n, c_type alpha, const c_type *a, 
           i_type lda, c_type *b, i_type ldb)
{
    BLAS_NAME(ctrmm)(_rc(side),_rc(uplo),_rc(transa),_rc(diag),_rc(&m),_rc(&n),
                     _rc(&alpha),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
}; 

BLAS_EXPORT INLINE_TYPE 
void dtrmm(const char *side, const char *uplo, const char *transa, 
           const char *diag, i_type m, i_type n, d_type alpha, const d_type *a, 
           i_type lda, d_type *b, i_type ldb)
{
    BLAS_NAME(dtrmm)(_rc(side),_rc(uplo),_rc(transa),_rc(diag),_rc(&m),_rc(&n),
                     _rc(&alpha),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

BLAS_EXPORT INLINE_TYPE 
void ztrmm(const char *side, const char *uplo, const char *transa,
           const char *diag, i_type m, i_type n, z_type alpha, const z_type *a, 
           i_type lda, z_type *b, i_type ldb)
{
    BLAS_NAME(ztrmm)(_rc(side),_rc(uplo),_rc(transa),_rc(diag),_rc(&m),_rc(&n),
                     _rc(&alpha),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
}; 

//-----------------------------------------------------------------
//                          SYMM
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void symm<s_type>(const char *side, const char *uplo, i_type m, i_type n, 
                  s_type alpha, const s_type *a, i_type lda, const s_type *b, 
                  i_type ldb, s_type beta, s_type *c, i_type ldc)
{
    BLAS_NAME(ssymm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),
                     _rc(&lda),_rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
};

template<> BLAS_EXPORT INLINE_TYPE 
void symm<d_type>(const char *side, const char *uplo, i_type m, i_type n,
                  d_type alpha, const d_type *a, i_type lda, const d_type *b, 
                  i_type ldb, d_type beta, d_type *c, i_type ldc)
{
    BLAS_NAME(dsymm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),
                     _rc(&lda),_rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
};

template<> BLAS_EXPORT INLINE_TYPE 
void symm<c_type>(const char *side, const char *uplo, i_type m, i_type n,
                  c_type alpha, const c_type *a, i_type lda, const c_type *b, 
                  i_type ldb, c_type beta, c_type *c, i_type ldc)
{
    BLAS_NAME(csymm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),
                     _rc(&lda),_rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
};

template<> BLAS_EXPORT INLINE_TYPE 
void symm<z_type>(const char *side, const char *uplo, i_type m, i_type n, 
                  z_type alpha, const z_type *a, i_type lda, const z_type *b,
                  i_type ldb, z_type beta, z_type *c, i_type ldc)
{
    BLAS_NAME(zsymm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),
                     _rc(&lda),_rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE 
void ssymm(const char *side, const char *uplo, i_type m, i_type n, s_type alpha,
           const s_type *a, i_type lda, const s_type *b, i_type ldb, s_type beta,
           s_type *c, i_type ldc)
{
    BLAS_NAME(ssymm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),
                     _rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE 
void csymm(const char *side, const char *uplo, i_type m, i_type n, c_type alpha, 
           const c_type *a, i_type lda,const c_type *b, i_type ldb, c_type beta, 
           c_type *c, i_type ldc)
{
    BLAS_NAME(csymm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),
                     _rc(&lda),_rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE 
void dsymm(const char *side, const char *uplo, i_type m, i_type n, d_type alpha,
           const d_type *a, i_type lda, const d_type *b, i_type ldb, d_type beta,
           d_type *c, i_type ldc)
{
    BLAS_NAME(dsymm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),
                     _rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
};

BLAS_EXPORT INLINE_TYPE 
void zsymm(const char *side, const char *uplo, i_type m, i_type n, z_type alpha,
           const z_type *a, i_type lda, const z_type *b, i_type ldb, z_type beta, 
           z_type *c, i_type ldc)
{
    BLAS_NAME(zsymm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),
                     _rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
};

//-----------------------------------------------------------------
//                          SYMV
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void symv<s_type>(const char *uplo, i_type n, s_type alpha, const s_type *a, i_type lda, 
                  const s_type *x, i_type incx, s_type beta, s_type *y, i_type incy)
{
    BLAS_NAME(ssymv)(_rc(uplo),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),_rc(&incx),
                     _rc(&beta),_rc(y),_rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void symv<d_type>(const char *uplo, i_type n, d_type alpha, const d_type *a, i_type lda, 
                  const d_type *x, i_type incx, d_type beta, d_type *y, i_type incy)
{
    BLAS_NAME(dsymv)(_rc(uplo),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),_rc(&incx),
                     _rc(&beta),_rc(y),_rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void symv<c_type>(const char *uplo, i_type n, c_type alpha, const c_type *a, i_type lda, 
                  const c_type *x, i_type incx, c_type beta, c_type *y, i_type incy)
{
    LAPACK_NAME(csymv)(_rc(uplo),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),_rc(&incx),
                       _rc(&beta),_rc(y),_rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void symv<z_type>(const char *uplo, i_type n, z_type alpha, const z_type *a, i_type lda, 
                  const z_type *x, i_type incx, z_type beta, z_type *y, i_type incy)
{
    LAPACK_NAME(zsymv)(_rc(uplo),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),_rc(&incx),
                       _rc(&beta),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE
void ssymv(const char *uplo, i_type n, s_type alpha, const s_type *a, i_type lda,
           const s_type *x, i_type incx, s_type beta, s_type *y, i_type incy)
{
    BLAS_NAME(ssymv)(_rc(uplo),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),_rc(&incx),
                     _rc(&beta),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE 
void dsymv(const char *uplo, i_type n, d_type alpha, const d_type *a, i_type lda,
           const d_type *x, i_type incx, d_type beta, d_type *y, i_type incy)
{
    BLAS_NAME(dsymv)(_rc(uplo),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),_rc(&incx),
                     _rc(&beta),_rc(y),_rc(&incy));
};

//-----------------------------------------------------------------
//                          HEMM
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void hemm<c_type>(const char *side, const char *uplo, i_type m, i_type n,
                  c_type alpha, const c_type *a, i_type lda, const c_type *b, 
                  i_type ldb, c_type beta, c_type *c, i_type ldc)
{
    BLAS_NAME(chemm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),
                     _rc(&lda),_rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
}; 

template<> BLAS_EXPORT INLINE_TYPE 
void hemm<z_type>(const char *side, const char *uplo, i_type m, i_type n, 
                  z_type alpha, const z_type *a, i_type lda, const z_type *b,
                  i_type ldb, z_type beta, z_type *c, i_type ldc)
{
    BLAS_NAME(zhemm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),
                     _rc(&lda),_rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
}; 

template<> BLAS_EXPORT INLINE_TYPE 
void hemm<s_type>(const char *side, const char *uplo, i_type m, i_type n,
                  s_type alpha, const s_type *a, i_type lda, const s_type *b, 
                  i_type ldb, s_type beta, s_type *c, i_type ldc)
{
    BLAS_NAME(ssymm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),
                     _rc(&lda),_rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
}; 

template<> BLAS_EXPORT INLINE_TYPE 
void hemm<d_type>(const char *side, const char *uplo, i_type m, i_type n,
                  d_type alpha, const d_type *a, i_type lda, const d_type *b,
                  i_type ldb, d_type beta, d_type *c, i_type ldc)
{
    BLAS_NAME(dsymm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),
                     _rc(&lda),_rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
}; 

BLAS_EXPORT INLINE_TYPE  
void chemm(const char *side, const char *uplo, i_type m, i_type n, c_type alpha, 
           const c_type *a, i_type lda, const c_type *b, i_type ldb, c_type beta,
           c_type *c, i_type ldc)
{
    BLAS_NAME(chemm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),
                     _rc(&lda),_rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
}; 

BLAS_EXPORT INLINE_TYPE 
void zhemm(const char *side, const char *uplo, i_type m, i_type n, z_type alpha,
           const z_type *a, i_type lda, const z_type *b, i_type ldb, z_type beta, 
           z_type *c, i_type ldc)
{
    BLAS_NAME(zhemm)(_rc(side),_rc(uplo),_rc(&m),_rc(&n),_rc(&alpha),_rc(a),
                     _rc(&lda),_rc(b),_rc(&ldb), _rc(&beta),_rc(c),_rc(&ldc));
}; 

//-----------------------------------------------------------------
//                          GMBV
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void gbmv<s_type>(const char *trans, i_type m, i_type n, i_type kl, i_type ku,
                  s_type alpha, const s_type *a, i_type lda, const s_type *x,
                  i_type incx, s_type beta, s_type *y, i_type incy)
{
    BLAS_NAME(sgbmv)(_rc(trans),_rc(&m),_rc(&n),_rc(&kl),_rc(&ku),_rc(&alpha),
                     _rc(a),_rc(&lda),_rc(x),_rc(&incx), _rc(&beta),_rc(y),_rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gbmv<d_type>(const char *trans, i_type m, i_type n, i_type kl, i_type ku, 
                  d_type alpha, const d_type *a, i_type lda, const d_type *x, 
                  i_type incx, d_type beta, d_type *y, i_type incy)
{
    BLAS_NAME(dgbmv)(_rc(trans),_rc(&m),_rc(&n),_rc(&kl),_rc(&ku),_rc(&alpha),
                     _rc(a),_rc(&lda),_rc(x),_rc(&incx), _rc(&beta),_rc(y),_rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gbmv<c_type>(const char *trans, i_type m, i_type n, i_type kl, i_type ku, 
                  c_type alpha, const c_type *a, i_type lda, const c_type *x, i_type incx, 
                  c_type beta, c_type *y, i_type incy)
{
    BLAS_NAME(cgbmv)(_rc(trans),_rc(&m),_rc(&n),_rc(&kl),_rc(&ku),_rc(&alpha),_rc(a),
                     _rc(&lda),_rc(x),_rc(&incx), _rc(&beta),_rc(y),_rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gbmv<z_type>(const char *trans, i_type m, i_type n, i_type kl, i_type ku, 
                  z_type alpha, const z_type *a, i_type lda, const z_type *x, 
                  i_type incx, z_type beta, z_type *y, i_type incy)
{
    BLAS_NAME(zgbmv)(_rc(trans),_rc(&m),_rc(&n),_rc(&kl),_rc(&ku),_rc(&alpha),
                     _rc(a),_rc(&lda),_rc(x),_rc(&incx), _rc(&beta),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE 
void sgbmv(const char *trans, i_type m, i_type n, i_type kl, i_type ku, 
           s_type alpha, const s_type *a, i_type lda, const s_type *x, i_type incx, 
           s_type beta, s_type *y, i_type incy)
{
    BLAS_NAME(sgbmv)(_rc(trans),_rc(&m),_rc(&n),_rc(&kl),_rc(&ku),_rc(&alpha),
                     _rc(a),_rc(&lda),_rc(x),_rc(&incx), _rc(&beta),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE 
void cgbmv(const char *trans, i_type m, i_type n, i_type kl, i_type ku, 
           c_type alpha, const c_type *a, i_type lda, const c_type *x, i_type incx,
           c_type beta, c_type *y, i_type incy)
{
     BLAS_NAME(cgbmv)(_rc(trans),_rc(&m),_rc(&n),_rc(&kl),_rc(&ku),_rc(&alpha),
                      _rc(a),_rc(&lda),_rc(x),_rc(&incx), _rc(&beta),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE 
void dgbmv(const char *trans, i_type m, i_type n, i_type kl, i_type ku, d_type alpha, 
           const d_type *a, i_type lda, const d_type *x, i_type incx, d_type beta, 
           d_type *y, i_type incy)
{
    BLAS_NAME(dgbmv)(_rc(trans),_rc(&m),_rc(&n),_rc(&kl),_rc(&ku),_rc(&alpha),_rc(a),
                     _rc(&lda),_rc(x),_rc(&incx), _rc(&beta),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE 
void zgbmv(const char *trans, i_type m, i_type n, i_type kl, i_type ku, z_type alpha, 
           const z_type *a, i_type lda, const z_type *x, i_type incx, z_type beta, 
           z_type *y, i_type incy)
{
    BLAS_NAME(zgbmv)(_rc(trans),_rc(&m),_rc(&n),_rc(&kl),_rc(&ku),_rc(&alpha),
                     _rc(a),_rc(&lda),_rc(x),_rc(&incx), _rc(&beta),_rc(y),_rc(&incy));
};

//-----------------------------------------------------------------
//                          HEMV
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void hemv<c_type>(const char *uplo, i_type n, c_type alpha, const c_type *a, i_type lda, 
                  const c_type *x, i_type incx, c_type beta, c_type *y, i_type incy)
{
    BLAS_NAME(chemv)(_rc(uplo),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),_rc(&incx),
                     _rc(&beta),_rc(y),_rc(&incy));
}; 

template<> BLAS_EXPORT INLINE_TYPE 
void hemv<z_type>(const char *uplo, i_type n, z_type alpha, const z_type *a, i_type lda, 
                  const z_type *x, i_type incx, z_type beta, z_type *y, i_type incy)
{
    BLAS_NAME(zhemv)(_rc(uplo),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),_rc(&incx),
                     _rc(&beta),_rc(y),_rc(&incy));
}; 

template<> BLAS_EXPORT INLINE_TYPE 
void hemv<s_type>(const char *uplo, i_type n, s_type alpha, const s_type *a, i_type lda, 
                  const s_type *x, i_type incx, s_type beta, s_type *y, i_type incy)
{
    BLAS_NAME(ssymv)(_rc(uplo),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),_rc(&incx),
                     _rc(&beta),_rc(y),_rc(&incy));
}; 

template<> BLAS_EXPORT INLINE_TYPE 
void hemv<d_type>(const char *uplo, i_type n, d_type alpha, const d_type *a, i_type lda, 
                  const d_type *x, i_type incx, d_type beta, d_type *y, i_type incy)
{
    BLAS_NAME(dsymv)(_rc(uplo),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),_rc(&incx),
                     _rc(&beta),_rc(y),_rc(&incy));
}; 


BLAS_EXPORT INLINE_TYPE
void chemv(const char *uplo, i_type n, c_type alpha, const c_type *a, i_type lda, 
           const c_type *x, i_type incx, c_type beta, c_type *y, i_type incy)
{
    BLAS_NAME(chemv)(_rc(uplo),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),_rc(&incx),
                     _rc(&beta),_rc(y),_rc(&incy));
}; 

BLAS_EXPORT INLINE_TYPE
void zhemv(const char *uplo, i_type n, z_type alpha, const z_type *a, i_type lda, 
           const z_type *x, i_type incx, z_type beta, z_type *y, i_type incy)
{
    BLAS_NAME(zhemv)(_rc(uplo),_rc(&n),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),_rc(&incx),
                     _rc(&beta),_rc(y),_rc(&incy));
}; 

//-----------------------------------------------------------------
//                          SBMV
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void sbmv<s_type>(const char *uplo, i_type n, i_type k, s_type alpha, const s_type *a,
                  i_type lda, const s_type *x, i_type incx, s_type beta, s_type *y, 
                  i_type incy)
{
    BLAS_NAME(ssbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                     _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
}; 

template<> BLAS_EXPORT INLINE_TYPE 
void sbmv<d_type>(const char *uplo, i_type n, i_type k, d_type alpha, const d_type *a, 
                  i_type lda, const d_type *x, i_type incx, d_type beta, d_type *y,
                  i_type incy)
{
    BLAS_NAME(dsbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                     _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void sbmv<c_type>(const char *uplo, i_type n, i_type k, c_type alpha, const c_type *a,
                  i_type lda, const c_type *x, i_type incx, c_type beta, c_type *y, 
                  i_type incy)
{
    LAPACK_NAME(csbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                       _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
};

template<> BLAS_EXPORT INLINE_TYPE 
void sbmv<z_type>(const char *uplo, i_type n, i_type k, z_type alpha, const z_type *a, 
                  i_type lda, const z_type *x, i_type incx, z_type beta, z_type *y,
                  i_type incy)
{
    LAPACK_NAME(zsbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                       _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE 
void ssbmv(const char *uplo, i_type n, i_type k, s_type alpha, const s_type *a, 
           i_type lda, const s_type *x, i_type incx, s_type beta, s_type *y, i_type incy)
{
    BLAS_NAME(ssbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                     _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE 
void dsbmv(const char *uplo, i_type n, i_type k, d_type alpha, const d_type *a, 
           i_type lda, const d_type *x, i_type incx, d_type beta, d_type *y, i_type incy)
{
    BLAS_NAME(dsbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                     _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE 
void csbmv(const char *uplo, i_type n, i_type k, c_type alpha, const c_type *a, 
           i_type lda, const c_type *x, i_type incx, c_type beta, c_type *y, i_type incy)
{
    LAPACK_NAME(csbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                       _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
};

BLAS_EXPORT INLINE_TYPE 
void zsbmv(const char *uplo, i_type n, i_type k, z_type alpha, const z_type *a, 
           i_type lda, const z_type *x, i_type incx, z_type beta, z_type *y, i_type incy)
{
    LAPACK_NAME(zsbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                       _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
};

//-----------------------------------------------------------------
//                          TBMV
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void tbmv<s_type>(const char *uplo, const char *trans, const char *diag, i_type n, i_type k, 
                  const s_type *a, i_type lda, s_type *x, i_type incx)
{
    BLAS_NAME(stbmv)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&k),_rc(a),
                     _rc(&lda),_rc(x),_rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE 
void tbmv<d_type>(const char *uplo, const char *trans, const char *diag, i_type n,
                  i_type k, const d_type *a, i_type lda, d_type *x, i_type incx)
{
    BLAS_NAME(dtbmv)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&k),_rc(a),_rc(&lda),
                     _rc(x),_rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE 
void tbmv<c_type>(const char *uplo, const char *trans, const char *diag, i_type n, 
                  i_type k, const c_type *a, i_type lda, c_type *x, i_type incx)
{
    BLAS_NAME(ctbmv)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&k),_rc(a),
                     _rc(&lda),_rc(x),_rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE 
void tbmv<z_type>(const char *uplo, const char *trans, const char *diag, i_type n,
                  i_type k, const z_type *a, i_type lda, z_type *x, i_type incx)
{
    BLAS_NAME(ztbmv)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&k),_rc(a),
                     _rc(&lda),_rc(x),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE 
void stbmv(const char *uplo, const char *trans, const char *diag, i_type n, i_type k,
           const s_type *a, i_type lda, s_type *x, i_type incx)
{
    BLAS_NAME(stbmv)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&k),_rc(a),
                     _rc(&lda),_rc(x),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE 
void ctbmv(const char *uplo, const char *trans, const char *diag, i_type n, i_type k,
           const c_type *a, i_type lda, c_type *x, i_type incx)
{
    BLAS_NAME(ctbmv)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&k),_rc(a),_rc(&lda),
                     _rc(x),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE 
void dtbmv(const char *uplo, const char *trans, const char *diag, i_type n, i_type k,
           const d_type *a, i_type lda, d_type *x, i_type incx)
{
    BLAS_NAME(dtbmv)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&k),_rc(a),_rc(&lda),
                     _rc(x),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE 
void ztbmv(const char *uplo, const char *trans, const char *diag, i_type n, i_type k,
           const z_type *a, i_type lda, z_type *x, i_type incx)
{
    BLAS_NAME(ztbmv)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&k),_rc(a),_rc(&lda),
                     _rc(x),_rc(&incx));
};

//-----------------------------------------------------------------
//                          HBMV
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void hbmv<c_type>(const char *uplo, i_type n, i_type k, c_type alpha, const c_type *a, 
                  i_type lda, const c_type *x, i_type incx, c_type beta, c_type *y,
                  i_type incy)
{
    BLAS_NAME(chbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                     _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
}; 

template<> BLAS_EXPORT INLINE_TYPE 
void hbmv<z_type>(const char *uplo, i_type n, i_type k, z_type alpha, const z_type *a, 
                  i_type lda, const z_type *x, i_type incx, z_type beta, z_type *y,
                  i_type incy)
{
    BLAS_NAME(zhbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                     _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
}; 

template<> BLAS_EXPORT INLINE_TYPE 
void hbmv<d_type>(const char *uplo, i_type n, i_type k, d_type alpha, const d_type *a, 
                  i_type lda, const d_type *x, i_type incx, d_type beta, d_type *y, 
                  i_type incy)
{
    BLAS_NAME(dsbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                     _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
}; 

template<> BLAS_EXPORT INLINE_TYPE 
void hbmv<s_type>(const char *uplo, i_type n, i_type k, s_type alpha, const s_type *a,
                  i_type lda, const s_type *x, i_type incx, s_type beta, s_type *y, 
                  i_type incy)
{
    BLAS_NAME(ssbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                     _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
}; 

BLAS_EXPORT INLINE_TYPE 
void chbmv(const char *uplo, i_type n, i_type k, c_type alpha, const c_type *a, 
           i_type lda, const c_type *x, i_type incx, c_type beta, c_type *y, i_type incy)
{
    BLAS_NAME(chbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                     _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
}; 

BLAS_EXPORT INLINE_TYPE 
void zhbmv(const char *uplo, i_type n, i_type k, z_type alpha, const z_type *a, 
           i_type lda, const z_type *x, i_type incx, z_type beta, z_type *y, i_type incy)
{
    BLAS_NAME(zhbmv)(_rc(uplo),_rc(&n),_rc(&k),_rc(&alpha),_rc(a),_rc(&lda),_rc(x),
                     _rc(&incx),_rc(&beta),_rc(y),_rc(&incy));
}; 

//-----------------------------------------------------------------
//                          tbsv
//-----------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void tbsv<s_type>(const char *uplo, const char *trans, const char *diag, i_type n, i_type k, 
                  const s_type *a, i_type lda, s_type *x, i_type incx)
{
    BLAS_NAME(stbsv)(_rc(uplo), _rc(trans), _rc(diag), _rc(&n), _rc(&k), _rc(a),
                     _rc(&lda), _rc(x), _rc(&incx));
}; 

template<> BLAS_EXPORT INLINE_TYPE 
void tbsv<d_type>(const char *uplo, const char *trans, const char *diag, i_type n, 
                  i_type k, const d_type *a,  i_type lda, d_type *x, i_type incx)
{
    BLAS_NAME(dtbsv)(_rc(uplo), _rc(trans), _rc(diag), _rc(&n), _rc(&k), _rc(a), 
                     _rc(&lda), _rc(x), _rc(&incx));
}; 

template<> BLAS_EXPORT INLINE_TYPE 
void tbsv<c_type>(const char *uplo, const char *trans, const char *diag, i_type n,
                  i_type k, const c_type *a, i_type lda, c_type *x, i_type incx)
{
    BLAS_NAME(ctbsv)(_rc(uplo), _rc(trans), _rc(diag), _rc(&n), _rc(&k), _rc(a), 
                     _rc(&lda), _rc(x), _rc(&incx));
}; 

template<> BLAS_EXPORT INLINE_TYPE 
void tbsv<z_type>(const char *uplo, const char *trans, const char *diag, i_type n, 
                  i_type k, const z_type *a, i_type lda, z_type *x, i_type incx)
{
    BLAS_NAME(ztbsv)(_rc(uplo), _rc(trans), _rc(diag), _rc(&n), _rc(&k), _rc(a),
                     _rc(&lda), _rc(x), _rc(&incx));
}; 

BLAS_EXPORT INLINE_TYPE
void stbsv(const char *uplo, const char *trans, const char *diag, i_type n,
           i_type k, const s_type *a, i_type lda, s_type *x, i_type incx)
{
    BLAS_NAME(stbsv)(_rc(uplo), _rc(trans), _rc(diag), _rc(&n), _rc(&k), _rc(a),
                     _rc(&lda), _rc(x), _rc(&incx));
};

BLAS_EXPORT INLINE_TYPE
void dtbsv(const char *uplo, const char *trans, const char *diag, i_type n, 
           i_type k, const d_type *a, i_type lda, d_type *x, i_type incx)
{
    BLAS_NAME(dtbsv)(_rc(uplo), _rc(trans), _rc(diag), _rc(&n), _rc(&k), _rc(a),
                     _rc(&lda), _rc(x), _rc(&incx));
};

BLAS_EXPORT INLINE_TYPE
void ctbsv(const char *uplo, const char *trans, const char *diag, i_type n, 
           i_type k, const c_type *a, i_type lda, c_type *x, i_type incx)
{
    BLAS_NAME(ctbsv)(_rc(uplo), _rc(trans), _rc(diag), _rc(&n), _rc(&k), _rc(a),
                     _rc(&lda), _rc(x), _rc(&incx));
};

BLAS_EXPORT INLINE_TYPE
void ztbsv(const char *uplo, const char *trans, const char *diag, i_type n, 
           i_type k, const z_type *a, i_type lda, z_type *x, i_type incx)
{
    BLAS_NAME(ztbsv)(_rc(uplo), _rc(trans), _rc(diag), _rc(&n), _rc(&k), 
                     _rc(a), _rc(&lda), _rc(x), _rc(&incx));
};

};};

#undef INLINE_TYPE
