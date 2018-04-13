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
#include "matcl-linalg/iterative_linear_eq/ksp_notifier_default.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/utils/linalg_utils.h"

#pragma warning (push)
#pragma warning (disable:4127)
#pragma warning (disable:4101)
#pragma warning (disable:4100)

#include "petscmat.h"
#include "petsc/private/kspimpl.h"
#include "petsc/private/matimpl.h"
#include "petscksp.h" 

#pragma warning (pop)

#include <iostream>

namespace matcl
{
//---------------------------------------------------------------
//                  ksp_solver_notifier
//---------------------------------------------------------------

ksp_solver_notifier::ksp_solver_notifier()
    :m_least_squares(false)
{};

bool ksp_solver_notifier::is_least_squares() const
{
    return m_least_squares;
};

void ksp_solver_notifier::set_least_squares(bool is_ls)
{
    m_least_squares = is_ls;
};

Matrix ksp_solver_notifier::get_residuals(KSP ksp) const
{
    PetscErrorCode ierr;
    Vec            resid;

    ierr        = KSPBuildResidual(ksp,NULL,NULL,&resid);

    if (ierr)
    {
        Integer N = 1;
        VecGetSize(ksp->vec_rhs, &N);

        return details::make_nan_matrix(N, 1, value_code::v_real);
    }

    Matrix res  = petsc::petsc_to_mmlib(resid);

    VecDestroy(&resid);
    return res;
};

Matrix ksp_solver_notifier::get_solution(KSP ksp) const
{
    PetscErrorCode ierr;
    Vec            sol;

    ierr        = KSPBuildSolution(ksp,NULL,&sol);

    if (ierr)
    {
        Integer N = 1;
        VecGetSize(ksp->vec_rhs, &N);

        return details::make_nan_matrix(N, 1, value_code::v_real);
    }

    Matrix res  = petsc::petsc_to_mmlib(sol);

    //VecDestroy(&sol);
    return res;
};

Matrix ksp_solver_notifier::get_rhs(KSP ksp) const
{
    Matrix res  = petsc::petsc_to_mmlib(ksp->vec_rhs);
    return res;
};

/// calculate norm of true residuals
Real ksp_solver_notifier::true_resid_norm(KSP ksp, basic_vector_norm norm_type) const
{
    PetscErrorCode ierr;
    Vec            resid;
    PetscReal      truenorm = 0.0;

    NormType nt;

    switch(norm_type)
    {
        case basic_vector_norm::norm_1:     nt = NORM_1; break;
        case basic_vector_norm::norm_2:     nt = NORM_2; break;
        case basic_vector_norm::norm_inf:   nt = NORM_INFINITY; break;

        default:
            nt = NORM_2;
            break;
    }

    ierr    = KSPBuildResidual(ksp,NULL,NULL,&resid);
    
    if (ierr)
        return constants::nan();

    ierr    = VecNorm(resid,nt,&truenorm);

    if (ierr)
        return constants::nan();

    VecDestroy(&resid);

    return truenorm;
};

/// calculate norm of RHS vector
Real ksp_solver_notifier::rhs_norm(KSP ksp, basic_vector_norm norm_type) const
{
    Real bnorm = 0.0;
    NormType nt;

    switch(norm_type)
    {
        case basic_vector_norm::norm_1:     nt = NORM_1; break;
        case basic_vector_norm::norm_2:     nt = NORM_2; break;
        case basic_vector_norm::norm_inf:   nt = NORM_INFINITY; break;

        default:
            nt = NORM_2;
            break;
    }

    PetscErrorCode ierr = ::VecNorm(ksp->vec_rhs, nt, &bnorm);
    
    if (ierr)
        return constants::nan();

    return bnorm;
}

};

namespace matcl { namespace details
{

//---------------------------------------------------------------
//                  default_notifier
//---------------------------------------------------------------

//TODO: remove this
default_notifier::default_notifier(const ksp_solver_impl* owner)
    : m_it(0), m_norm(0.0), m_anorm(0.0), m_owner(owner)
    , m_solve_state(not_initialized)
{};

default_notifier::~default_notifier()
{}

void default_notifier::begin_solve(const ksp_solver& solver, const linear_operator& op, Integer nvec)
{
    (void)nvec;
    (void)op;

    m_solve_state       = st_begin_solve;
    std::string ps      = solver.solver_name();
    std::string pc      = solver.precond_name();
    m_id                = solver.get_id();

    std::cout << "solver " << m_id << ": " << ps << " " << pc << "\n"; 
}
void default_notifier::end_solve()
{
    m_solve_state       = st_end_solve;
    //std::cout << "end solve " << "\n";
}
void default_notifier::begin_solve_vector(Integer i)
{
    m_solve_state       = st_begin_vector;
    (void)i;
    //std::cout << "begin solve vector " << i << "\n";
}
void default_notifier::end_solve_vector()
{
    m_solve_state       = st_end_vector;
    //std::cout << "end solve vector" << "\n";
}
void default_notifier::begin_iterations(const Matrix& RHS)
{
    m_solve_state       = st_iterations;

    (void)RHS;
    //std::cout << "begin iterations " << "\n";
}
void default_notifier::end_iterations(const linear_operator& op_A, const Matrix& RHS, const Matrix& sol)
{    
    Matrix resid = op_A.mmul_right(sol, trans_type::no_trans) - RHS;

    Real true_norm  = norm(resid, 2.0);

    if (this->is_least_squares())
    {
        std::cout << m_id << " " << "iterations: " << m_it << ", rnorm: " << m_norm << ", true rnorm: " 
                << true_norm << ", normal rnorm: " << m_anorm << "\n";
    }
    else
    {
        std::cout << m_id << " " << "iterations: " << m_it << ", rnorm: " << m_norm << ", true rnorm: " 
                << true_norm << "\n";
    }
}
                     
void default_notifier::report_iteration(Integer itc, matcl::Real rnorm, KSP ksp)
{
    if (m_solve_state < st_iterations && itc == 0 && false)
    {
        std::string ps      = m_owner->get_solver();
        std::string pc      = m_owner->get_precond();
        m_id                = m_owner->get_id();

        std::cout << "inner iteration" << "; ";
        std::cout << "solver " << m_id << ": " << ps << " " << pc << "\n";         
    };

    m_it    = itc;
    m_norm  = rnorm;
    (void)ksp;
    //std::cout << "report iteration: " << itc << " " << rnorm << "\n";
}
void default_notifier::report_iteration_ls(Integer itc, matcl::Real rnorm, matcl::Real arnorm, KSP ksp)
{
    m_it    = itc;
    m_norm  = rnorm;
    m_anorm = arnorm;
    (void)ksp;
    //std::cout << "report iteration: " << itc << " " << rnorm << "\n";
}

void default_notifier::report_error(const std::string& msg)
{            
    std::cout << "error: " << msg << "\n";
}

void default_notifier::report_divergence(const std::string& str)
{
    std::string ps      = m_owner->get_solver();
    std::cout << "solver " << m_id << "; divergence: " << str << "\n";
}
void default_notifier::report_convergence(const std::string& str)
{
    std::string ps      = m_owner->get_solver();
    std::cout << "solver " << m_id << "; convergence: " << str << "\n";
}

}};

#endif