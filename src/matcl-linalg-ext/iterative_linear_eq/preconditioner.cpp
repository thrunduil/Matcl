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

#include "matcl-linalg/linear_eq/preconditioner.h"
#include "matcl-linalg/iterative_linear_eq/ksp_solver_impl.h"
#include "matcl-linalg/petsc_utils/petsc_option.h"
#include "matcl-linalg/petsc_utils/petsc_algs.h"
#include "matcl-linalg/linear_eq/linsolve_object.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-linalg/special_matrices/struct_flag_linalg.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-linalg/options/options_petsc.h"

namespace matcl
{

static options make_options(const options& opts, petsc_precond type)
{
    options ret = opts;
    ret.set(matcl::opt::petsc::preconditioner(type));
    ret.set(matcl::opt::petsc::solver(petsc_solver::preonly));
    return ret;
};
static options make_options(const options& opts)
{
    options ret = opts;

    ret.set(matcl::opt::petsc::preconditioner(petsc_precond::none));
    ret.set(matcl::opt::petsc::solver(petsc_solver::preonly));
    return ret;
};
static options make_options()
{
    options ret;
    ret.set(matcl::opt::petsc::preconditioner(petsc_precond::none));
    ret.set(matcl::opt::petsc::solver(petsc_solver::preonly));
    return ret;
};

//----------------------------------------------------------------------
//                      preconditioner
//----------------------------------------------------------------------
precond::precond()
{};

precond::precond(const Matrix& P, petsc_precond type, const options& opts)
    :ksp_solver(P,P,make_options(opts,type))
{};

precond::precond(const linsolve_obj& P)
    :ksp_solver(P.base_matrix(), make_options(options(), petsc_precond::none))
{
    ksp_solver::set_precond(P);
};

precond::~precond()
{};

std::string precond::precond_name(petsc_precond precond)
{
    return petsc::converter_precond().to_string((Integer)precond);
};

petsc_precond precond::precond_code(const std::string& precond_name)
{
    return petsc::converter_precond().to_code(precond_name);
};

//--------------------------------------------------------------------
//                   block preconditioner
//--------------------------------------------------------------------

block_precond::block_precond()
{};

block_precond::block_precond(const Matrix& P, Integer n_blocks)
    :ksp_solver(P,P,make_options())
{
    m_impl->divide_on_blocks(n_blocks);
};

block_precond::block_precond(const Matrix& P, const Matrix& block_sizes)
    :ksp_solver(P,P,make_options())
{
    m_impl->divide_on_blocks(block_sizes);
};

block_precond::block_precond(const Matrix& P, Integer n, Integer overlap,
    asm_composition comp, asm_restrict restr)
    :ksp_solver(P,P,make_options())
{
    m_impl->divide_on_blocks_asm(n, overlap, comp, restr);
};

block_precond::block_precond(const Matrix& P, const std::vector<Matrix>& partition, Integer overlap, 
    asm_composition comp, asm_restrict restr)
    :ksp_solver(P,P,make_options())
{
    m_impl->divide_on_blocks_asm(partition, overlap, comp, restr);
};

block_precond::~block_precond()
{};

Matrix block_precond::get_block_precond_matrix(Integer i) const
{
    return m_impl->get_block_preconditioner_matrix(i);
};

Matrix block_precond::get_partitioning(Integer i) const
{
    return m_impl->get_partitioning(i);
}

void block_precond::set_block_solver(Integer i, const linsolve_obj& precond)
{
    m_impl->set_block_solver(i, precond);
};
void block_precond::set_block_solver(Integer i, const ksp_solver& precond)
{
    m_impl->set_block_solver(i, precond);
};
void block_precond::build_block_solver(Integer i, const options& opts)
{
    m_impl->build_block_solver(i, opts);
};

Integer block_precond::number_of_blocks() const
{
    return m_impl->number_of_blocks();
};

//--------------------------------------------------------------------
//                   composite preconditioner
//--------------------------------------------------------------------

composite_precond::composite_precond()
{};

/// build a preconditioner by composing together several preconditioners
composite_precond::composite_precond(const Matrix& P, petsc_composite_type type,
    const std::vector<ksp_solver>& precond_list)
    :ksp_solver(P, make_options(options()))
{
    m_impl->set_composite_preconditioner(type, precond_list);
};
composite_precond::composite_precond(const Matrix& P, petsc_composite_type type,
    const std::vector<linsolve_obj>& precond_list)
    :ksp_solver(P, make_options(options()))
{
    m_impl->set_composite_preconditioner(type, precond_list);
};

composite_precond::composite_precond(const Matrix& P, const ksp_solver& M1, const ksp_solver& M2)
    :ksp_solver(P, make_options(options()))
{
    m_impl->set_mult_preconditioner(M1, M2);
};

composite_precond::composite_precond(const Matrix& P, const linsolve_obj& M1, const linsolve_obj& M2)
    :ksp_solver(P, make_options(options()))
{
    m_impl->set_mult_preconditioner(M1, M2);
};

composite_precond::~composite_precond()
{};

//--------------------------------------------------------------------
//                   Galerkin preconditioner
//--------------------------------------------------------------------
galerkin_precond::galerkin_precond()
{};

galerkin_precond::galerkin_precond(const Matrix& B, const Matrix& PR, Matrix& S)
    :ksp_solver(B, make_options(options()))
{
    S = m_impl->set_galerkin_RP(PR);
};

galerkin_precond::galerkin_precond(const Matrix& B, const Matrix& P, const Matrix& R, Matrix& S)
    :ksp_solver(B, make_options(options()))
{
    S = m_impl->set_galerkin_RP(R, P);
};

void galerkin_precond::set_galerkin_solver(const ksp_solver& S_solver)
{
    m_impl->set_galerkin_solver(S_solver);
};

galerkin_precond::~galerkin_precond()
{};

//--------------------------------------------------------------------
//                   Field Split Preconditioner
//--------------------------------------------------------------------

field_split_precond::field_split_precond()
{};

field_split_precond::field_split_precond(const Matrix& P, const std::vector<Matrix>& partitions,
    field_split_type fst)
    :ksp_solver(P, make_options())
{
    m_impl->field_split(partitions, fst);
};

field_split_precond::~field_split_precond()
{};

Matrix field_split_precond::get_block_precond_matrix(Integer i) const
{
    return m_impl->fieldsplit_get_block_preconditioner_matrix(i);
};
Matrix field_split_precond::get_partitioning(Integer i) const
{
    return m_impl->fieldsplit_get_partitioning(i);
};

void field_split_precond::set_block_solver(Integer i, const linsolve_obj& precond)
{
    m_impl->fieldsplit_set_block_solver(i, precond);
};

void field_split_precond::set_block_solver(Integer i, const ksp_solver& precond)
{
    m_impl->fieldsplit_set_block_solver(i, precond);
};

void field_split_precond::build_block_solver(Integer i, const options& opts)
{
    m_impl->fieldsplit_build_block_solver(i, opts);
};

/// return number of blocks of the preconditioner (0 if is_partitioned == false)
Integer field_split_precond::number_of_blocks() const
{
    return m_impl->fieldsplit_number_of_blocks();
};

//--------------------------------------------------------------------
//                   Multigrid
//--------------------------------------------------------------------

multigrid::multigrid()
{};

multigrid::multigrid(const Matrix& B, Integer lev, mg_type type, mg_cycle_type cycle, Integer mc)
    :ksp_solver(B, make_options(options()))
{
    m_impl->set_multigrid(lev, type, cycle, mc);
};

multigrid::~multigrid()
{};

Matrix multigrid::set_interp_restr(const Matrix& PR)
{
    return m_impl->mg_set_interp_restr(PR);
};

Matrix multigrid::set_interp_restr(const Matrix& P, const Matrix& R)
{
    return m_impl->mg_set_interp_restr(P, R);
};
Matrix multigrid::set_interp_restr_lev(const Matrix& S, const Matrix& PR)
{
    return m_impl->mg_set_interp_restr_lev(S, PR);
}
Matrix multigrid::set_interp_restr_lev(const Matrix& S, const Matrix& P, const Matrix& R)
{
    return m_impl->mg_set_interp_restr_lev(S, P, R);
}

void multigrid::set_post_smoother(const ksp_solver& smoother, Integer n_steps)
{
    return m_impl->mg_set_smoother_post(smoother, n_steps);
};

void multigrid::set_pre_smoother(const ksp_solver& smoother, Integer n_steps)
{
    return m_impl->mg_set_smoother_pre(smoother, n_steps);
};

void multigrid::coarser_level()
{
    return m_impl->coarser_level();
};
void multigrid::set_solver(const ksp_solver& S_solver)
{
    return m_impl->mg_set_solver(S_solver);
};

ksp_solver multigrid::richardson_smoother(const Matrix& S, Real omega, bool jacobi, bool auto_scale)
{
    namespace op = opt::petsc;

    options opts;
    opts.set(op::solver(petsc_solver::richardson));

    if (jacobi)
        opts.set(op::preconditioner(petsc_precond::jacobi));
    else
        opts.set(op::preconditioner(petsc_precond::none));

    opts.set(op::richardson_damping(omega));

    if (auto_scale == true)
        opts.set(op::richardson_auto(true));
    else
        opts.set(op::richardson_auto(false));

    ksp_solver ret(S, opts);
    return ret;
};


#if 0

// fieldsplit schur is removed; configuration this preconditioner
// from source seems to be too hard (if possible at all)
class MATCL_LINALG_EXPORT field_split_schur_precond : public ksp_solver
{
    public:
        /// uninitialized preconditioner
        field_split_schur_precond();

        /// divide preconditioner matrix P on two blocks; later preconditioners must
        /// be set for each block
        field_split_schur_precond(const Matrix& P, const Matrix& part_1, const Matrix& part_2,
            schur_fact_type sft);

        /// standard destructor
        ~field_split_schur_precond();

        /// return first diagonal block of the matrix used to create preconditioner;
        Matrix                  get_block_precond_matrix() const;

        Matrix                  get_schur_matrix() const;

        /// return indices of elements that are assigned to i-th block (i = 1,2)
        Matrix                  get_partitioning(Integer i) const;

        /// set preconditioner forfirst block of partitioned matrix
        void                    set_block_precond(const linear_operator& precond);
        void                    set_block_precond(const linsolve_obj& precond);
        void                    set_block_precond(const ksp_solver& precond);

        /// build preconditioner for first block of partitioned matrix from options
        /// preconditioner is given by general solver, one should probably select preonly solver
        void                    build_block_precond(const options& opts);

        void                    set_schur_precond_type(schur_precond_type t);
        void                    set_schur_precond_lsc(const options& opts = options());
};

field_split_schur_precond::field_split_schur_precond()
{};

field_split_schur_precond::field_split_schur_precond(const Matrix& P, const Matrix& part_1, const Matrix& part_2,
    schur_fact_type sft)
    :ksp_solver(P, make_options_schur())
{
    m_impl->field_split_schur(part_1, part_2, sft);
};

field_split_schur_precond::~field_split_schur_precond()
{};

void field_split_schur_precond::set_schur_precond_lsc(const options& opts)
{
    return m_impl->set_schur_precond_lsc(opts);
};

void field_split_schur_precond::set_schur_precond_type(schur_precond_type t)
{
    return m_impl->set_schur_precond_type(t);
};

Matrix field_split_schur_precond::get_partitioning(Integer i) const
{
    return m_impl->fieldsplit_get_partitioning(i);
};
Matrix field_split_schur_precond::get_block_precond_matrix() const
{
    return m_impl->fieldsplit_get_block_preconditioner_matrix(1);
};
Matrix field_split_schur_precond::get_schur_matrix() const
{
    return m_impl->fieldsplit_get_schur_matrix();
};

void field_split_schur_precond::set_block_precond(const linear_operator& precond)
{
    m_impl->fieldsplit_set_block_preconditioner(1, precond);
};

void field_split_schur_precond::set_block_precond(const linsolve_obj& precond)
{
    m_impl->fieldsplit_set_block_preconditioner(1, precond);
};

void field_split_schur_precond::set_block_precond(const ksp_solver& precond)
{
    m_impl->fieldsplit_set_block_preconditioner(1, precond.to_linsolve_object());
};

void field_split_schur_precond::build_block_precond(const options& opts)
{
    m_impl->fieldsplit_build_block_preconditioner(1, opts);
};


#endif
}

#endif