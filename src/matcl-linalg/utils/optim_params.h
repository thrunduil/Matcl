/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "matcl-linalg/utils/utils.h"

namespace matcl { namespace details
{

struct linalg_optim_params
{
    static const Integer max_block_size_ORMQR   = 32;

    static Integer  block_size_ORMQR()              { return 32;  };
    static Real     max_row_density_chol_update()   { return 0.2; };
    static Integer  qrband_block_size()             { return 32;  };
    static Integer  qrband_block_crossover()        { return 96;  };
    static Integer  qrband_band_crossover()         { return 10;  };

    // if number of subdiagonal or superdiagonals lu is less than
    // threashold * N, then select givens version of qr instead of
    // the standard qr
    static Real     qr_givens_qr_threashold()       { return 0.1; };

    // if number of subdiagonal or superdiagonals lu is less than
    // threashold * N, then select qr instead of lu factorization
    static Real     dense_qr_lu_threashold()        { return 0.1; };

    // if number of subdiagonal or superdiagonals lu is less than
    // threashold * N, then select qr instead of lu factorization
    static Real     band_qr_lu_threashold()         { return 0.1; };

    // maximum bandwidth allowed for band symmetric Hessenberg
    static Real     crossover_bd_hess_sym()         { return 0.2; };

    // if number of subdiagonal or superdiagonals lu is less than
    // threashold * N, then select qr instead of lu factorization
    static Real     sparse_qr_lu_threashold()       { return 0.001; };    
};

};}
