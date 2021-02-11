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

#include "matcl-linalg/general/linalg_exception.h"

#pragma warning(push)
#pragma warning(disable:4251)

namespace matcl { namespace error
{
 
/// this class is throws if computation of eigenvectors in functions in hess_eig
/// family failed, due to convergence failure; members of this class contain 
/// information about current computation status and convergence status
class MATCL_LINALG_EXPORT hess_eig_failed : public matcl_exception_linalg
{
    public:
        /// current left eigenvectors if computed; matrix of size N x M, where N
        /// is the size of hessenberg matrix and M is the number of selected eigenvalues
        Matrix  m_VL;

        /// current right eigenvectors if computed; matrix of size N x M, where N
        /// is the size of hessenberg matrix and M is the number of selected eigenvalues
        Matrix  m_VR;

        /// integer matrix of size M; if i-th element is 0, then i-th left eigenvector
        /// converged otherwise not converged
        Matrix  m_fail_L;

        /// integer matrix of size M; if i-th element is 0, then i-th right eigenvector
        /// converged otherwise not converged
        Matrix  m_fail_R;

        /// true if left eigenvectors were computed
        bool    m_eval_L;

        /// true if right eigenvectors were computed
        bool    m_eval_R;

        /// number of eigenvectors which failed to converge
        Integer m_failed_eig;

    public:
        /// standard constructor
        hess_eig_failed(const Matrix& VL, const Matrix& VR, const Matrix& fail_L, const Matrix& fail_R, 
                        bool eval_L, bool eval_R, Integer failed_eig)
            : m_VL(VL), m_VR(VR), m_fail_L(fail_L), m_fail_R(fail_R), m_eval_L(eval_L), m_eval_R(eval_R)
            , m_failed_eig(failed_eig)
        {};

    private:
        virtual const char* what(exception_message_linalg& em) const;
};

};};

#pragma warning(pop)