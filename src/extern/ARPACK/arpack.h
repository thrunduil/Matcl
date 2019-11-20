/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2011-2016
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

#include "matcl-blas-lapack/lapack/lapack.h"

namespace matcl 
{ 
namespace arpack
{

using matcl::lapack::i_type;
using matcl::lapack::d_type;
using matcl::lapack::s_type;
using matcl::lapack::z_type;
using matcl::lapack::c_type;
using matcl::lapack::l_type;
    
extern "C" 
{

void dnaupd_(   i_type *ido, const char *bmat, const i_type *n, const char *which,
                i_type *nev, const d_type *tol, d_type *resid,
                i_type *ncv, d_type *V, i_type *ldv,
                i_type *iparam, i_type *ipntr, d_type *workd,
                d_type *workl, i_type *lworkl, i_type *info);

void dneupd_(   l_type *rvec, const char *HowMny, l_type *select,
                d_type *dr, d_type *di, d_type *Z,
                i_type *ldz, d_type *sigmar,
                d_type *sigmai, d_type *workev,
                const char *bmat, const i_type *n, const char *which,
                i_type *nev, const d_type *tol, d_type *resid,
                i_type *ncv, d_type *V, i_type *ldv,
                i_type *iparam, i_type *ipntr,
                d_type *workd, d_type *workl,
                i_type *lworkl, i_type *info);

void dsaupd_(   i_type *ido, const char *bmat, const i_type *n, const char *which,
                i_type *nev, const d_type *tol, d_type *resid,
                i_type *ncv, d_type *V, i_type *ldv,
                i_type *iparam, i_type *ipntr, d_type *workd,
                d_type *workl, i_type *lworkl, i_type *info);

void dseupd_(   l_type *rvec, const char *HowMny, l_type *select,
                d_type *d, d_type *Z,
                i_type *ldz, d_type *sigma,
                const char *bmat, const i_type *n, const char *which,
                i_type *nev, const d_type *tol, d_type *resid,
                i_type *ncv, d_type *V, i_type *ldv,
                i_type *iparam, i_type *ipntr,
                d_type *workd, d_type *workl,
                i_type *lworkl, i_type *info);

void snaupd_(   i_type *ido, const char *bmat, const i_type *n, const char *which,
                i_type *nev, const s_type *tol, s_type *resid,
                i_type *ncv, s_type *V, i_type *ldv,
                i_type *iparam, i_type *ipntr, s_type *workd,
                s_type *workl, i_type *lworkl, i_type *info);

void sneupd_(   l_type *rvec, const char *HowMny, l_type *select,
                s_type *dr, s_type *di, s_type *Z,
                i_type *ldz, s_type *sigmar,
                s_type *sigmai, s_type *workev,
                const char *bmat, const i_type *n, const char *which,
                i_type *nev, const s_type *tol, s_type *resid,
                i_type *ncv, s_type *V, i_type *ldv,
                i_type *iparam, i_type *ipntr,
                s_type *workd, s_type *workl,
                i_type *lworkl, i_type *info);

void ssaupd_(   i_type *ido, const char *bmat, const i_type *n, const char *which,
                i_type *nev, const s_type *tol, s_type *resid,
                i_type *ncv, s_type *V, i_type *ldv,
                i_type *iparam, i_type *ipntr, s_type *workd,
                s_type *workl, i_type *lworkl, i_type *info);

void sseupd_(   l_type *rvec, const char *HowMny, l_type *select,
                s_type *d, s_type *Z,
                i_type *ldz, s_type *sigma,
                const char *bmat, const i_type *n, const char *which,
                i_type *nev, const s_type *tol, s_type *resid,
                i_type *ncv, s_type *V, i_type *ldv,
                i_type *iparam, i_type *ipntr,
                s_type *workd, s_type *workl,
                i_type *lworkl, i_type *info);

void znaupd_(   i_type *ido, const char *bmat, const i_type *n, const char *which,
                i_type *nev, const d_type *tol, z_type *resid,
                i_type *ncv, z_type *V, i_type *ldv,
                i_type *iparam, i_type *ipntr, z_type *workd,
                z_type *workl, i_type *lworkl, d_type *rwork, i_type *info);

void zneupd_(   l_type *rvec, const char *HowMny, l_type *select,
                z_type *d, z_type *Z,
                i_type *ldz, z_type *sigma, z_type *workev,
                const char *bmat, const i_type *n, const char *which,
                i_type *nev, const d_type *tol, z_type *resid,
                i_type *ncv, z_type *V, i_type *ldv,
                i_type *iparam, i_type *ipntr,
                z_type *workd, z_type *workl,
                i_type *lworkl, d_type *rwork, i_type *info);

void cnaupd_(   i_type *ido, const char *bmat, const i_type *n, const char *which,
                i_type *nev, const s_type *tol, c_type *resid,
                i_type *ncv, c_type *V, i_type *ldv,
                i_type *iparam, i_type *ipntr, c_type *workd,
                c_type *workl, i_type *lworkl, s_type *rwork, i_type *info);

void cneupd_(   l_type *rvec, const char *HowMny, l_type *select,
                c_type *d, c_type *Z,
                i_type *ldz, c_type *sigma, c_type *workev,
                const char *bmat, const i_type *n, const char *which,
                i_type *nev, const s_type *tol, c_type *resid,
                i_type *ncv, c_type *V, i_type *ldv,
                i_type *iparam, i_type *ipntr,
                c_type *workd, c_type *workl,
                i_type *lworkl, s_type *rwork, i_type *info);

void dnaitr_(   i_type *ido, const char *bmat, const i_type *n, const i_type *k, const i_type *np, 
                const i_type *nb, d_type* resid, d_type* rnorm, d_type* v, const i_type *ldv, d_type* h, 
                const i_type *ldh, i_type *ipntr, d_type* workd, i_type* info);
void snaitr_(   i_type *ido, const char *bmat, const i_type *n, const i_type *k, const i_type *np, 
                const i_type *nb, s_type* resid, s_type* rnorm, s_type* v, const i_type *ldv, s_type* h, 
                const i_type *ldh, i_type *ipntr, s_type* workd, i_type* info);
void znaitr_(   i_type *ido, const char *bmat, const i_type *n, const i_type *k, const i_type *np, 
                const i_type *nb, z_type* resid, d_type* rnorm, z_type* v, const i_type *ldv, z_type* h, 
                const i_type *ldh, i_type *ipntr, z_type* workd, i_type* info);
void cnaitr_(   i_type *ido, const char *bmat, const i_type *n, const i_type *k, const i_type *np, 
                const i_type *nb, c_type* resid, s_type* rnorm, c_type* v, const i_type *ldv, c_type* h, 
                const i_type *ldh, i_type *ipntr, c_type* workd, i_type* info);

void dsaitr_(   i_type* ido, const char* bmat, const i_type* n, const i_type* k, const i_type* np, 
                const i_type* mode, d_type* resid, d_type* rnorm, d_type* v, const i_type* ldv, d_type* h, 
                const i_type* ldh, i_type* ipntr, d_type* workd, i_type* info);
void ssaitr_(   i_type* ido, const char* bmat, const i_type* n, const i_type* k, const i_type* np, 
                const i_type* mode, s_type* resid, s_type* rnorm, s_type* v, const i_type* ldv, s_type* h, 
                const i_type* ldh, i_type* ipntr, s_type* workd, i_type* info);

}

}}