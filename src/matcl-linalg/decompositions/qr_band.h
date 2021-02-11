/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
namespace matcl { namespace lapack
{

// from band_qr DGbBQRF QR factors an M by N band matrix in banded format, with blocking. 
template<class Val>
void gbbqrf(i_type NB, i_type M, i_type N, i_type ML, i_type MU, Val* A, 
                i_type LDA, Val* tau, Val* WORK, i_type& INFO );

//dgbbqr2
// from band_qr DGBBQR2 QR factors an M by N band matrix in banded format, with no blocking. 
template<class Val>
void gbbqr2(i_type M, i_type N, i_type ML, i_type MU, Val* A, 
                i_type LDA, Val* tau, Val* WORK, i_type& INFO);

// from band_qr DGEBQR2 QR factors an M by N band matrix in GE format, with no blocking. 
template<class Val>
void gebqr2(i_type M, i_type N, i_type ML, i_type MU, Val* A, i_type LDA, Val* tau, 
            Val* WORK, i_type& INFO );

// from band_qr DGEBQRF QR factors an M by N band matrix in GE format, with blocking. 
template<class Val>
void gebqrf(i_type NB, i_type M, i_type N, i_type ML, i_type MU, Val *A, i_type LDA, 
            Val* tau, Val* WORK, i_type& INFO );

};};