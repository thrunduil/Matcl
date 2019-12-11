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

#include "matcl-matrep/matrix/matrix.h"
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-matrep/matrix/permvec.h"
#include "matcl-linalg/special_matrices/unitary_matrix.h"

namespace matcl { namespace details
{
    class quern_solver_impl;
}};

namespace matcl
{

// class for sparse QR decomposition using Givens rotations
class quern_solver
{
    private:
        using impl_type = details::quern_solver_impl;
        using impl_ptr  = std::shared_ptr<impl_type>;

    private:
        impl_ptr    m_impl;
        bool        m_economy;

    public:
        // create sparse QR decomposition. QR = A or QR = A' if trans == true
        // trans = true is cheaper, so use this if QR=A' is required
        quern_solver(const Matrix& A, bool with_q, bool economy, bool trans);
        quern_solver(Matrix&& A, bool with_q, bool economy, bool trans);

        // get the unitary matrix Q of size M x K, where K = M if economy = false,
        // and K = min(M,N) if economy = true, where A has size M x N; returned matrix 
        // is not created explicitly
        unitary_matrix  get_unitary_matrix() const;

        // get the sparse upper-triangular factor R. If R' is required 
        // use get_r_trans which is cheaper
        Matrix          get_r() const;

        // get the sparse lower-triangular factor R'
        Matrix          get_r_trans() const;
};

};
