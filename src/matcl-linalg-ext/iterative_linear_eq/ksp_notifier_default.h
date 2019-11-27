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

#if 0
TODO

#include "matcl-linalg/iterative_linear_eq/ksp_solver_impl.h"
#include "matcl-matrep/matcl_matrep.h"

namespace matcl { namespace details
{

//---------------------------------------------------------------
//                  default_notifier
//---------------------------------------------------------------

class default_notifier : public ksp_solver_notifier
{
    private:
        enum solve_state 
        { 
            not_initialized, st_begin_solve, st_begin_vector, st_iterations, st_end_vector, st_end_solve
        };

    private:
        Integer                 m_it;
        Real                    m_norm;
        Real                    m_anorm;
        std::string             m_id;
        const ksp_solver_impl*  m_owner;
        solve_state             m_solve_state;

	public:
        default_notifier(const ksp_solver_impl* owner);

        virtual ~default_notifier();

        virtual void begin_solve(const ksp_solver& solver, const linear_operator& op, Integer nvec) override;
        virtual void end_solve() override;
        virtual void begin_solve_vector(Integer i) override;
        virtual void end_solve_vector() override;
        virtual void begin_iterations(const Matrix& RHS) override;
        virtual void end_iterations(const linear_operator& op_A, const Matrix& RHS, const Matrix& sol) override;
		virtual void report_iteration(Integer itc, matcl::Real rnorm, KSP ksp) override;
		virtual void report_iteration_ls(Integer itc, matcl::Real rnorm, matcl::Real arnorm, KSP ksp) override;
        virtual void report_error(const std::string& msg) override;
        virtual void report_divergence(const std::string& str) override;
        virtual void report_convergence(const std::string& str) override;
};

}};

#endif