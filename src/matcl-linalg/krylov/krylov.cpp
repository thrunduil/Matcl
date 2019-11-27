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

#include "matcl-linalg/krylov/krylov.h"
#include "matcl-linalg/krylov/block_arnoldi.h"
#include "matcl-linalg/krylov/arnoldi.h"

namespace matcl
{

ret_arnoldi matcl::blk_arnoldi(const linear_operator& A, const Matrix& v, Integer k, Real tol)
{
    Integer max_k   = k;
    Integer max_kb  = v.cols();

    block_arnoldi_iteration bli(A, max_k, max_kb);
    bli.run(v, k, tol);

    Matrix V    = bli.get_V();
    Matrix H    = bli.get_H();
    Matrix r    = bli.get_resid();
    Real rn     = bli.get_norm_resid();

    return ret_arnoldi(V,H,r,rn);
};

ret_arnoldi matcl::arnoldi(const linear_operator& A, const Matrix& v, Integer k, Real tol)
{
    Integer max_k   = k;

    arnoldi_iteration bli(A, max_k);
    bli.run(v, k, tol);

    Matrix V    = bli.get_V();
    Matrix H    = bli.get_H();
    Matrix r    = bli.get_resid();
    Real rn     = bli.get_norm_resid();

    return ret_arnoldi(V,H,r,rn);
};

ret_arnoldi matcl::arnoldi(const linear_operator& A, const linear_operator& B, 
                                        const Matrix& v, Integer k, Real tol)
{
    Integer max_k   = k;

    arnoldi_b_iteration bli(A, B, false, max_k);
    bli.run(v, k, tol);

    Matrix V    = bli.get_V();
    Matrix H    = bli.get_H();
    Matrix r    = bli.get_resid();
    Real rn     = bli.get_norm_resid();

    return ret_arnoldi(V,H,r,rn);
};
ret_arnoldi matcl::lanczos(const linear_operator& A, const linear_operator& B, 
                                        const Matrix& v, Integer k, Real tol)
{
    Integer max_k   = k;

    arnoldi_b_iteration bli(A, B, true, max_k);
    bli.run(v, k, tol);

    Matrix V    = bli.get_V();
    Matrix H    = bli.get_H();
    Matrix r    = bli.get_resid();
    Real rn     = bli.get_norm_resid();

    return ret_arnoldi(V,H,r,rn);
};

}