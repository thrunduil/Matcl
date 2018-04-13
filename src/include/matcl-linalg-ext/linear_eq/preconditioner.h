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

#include "matcl-linalg/general/config_linalg.h"
#include "matcl-linalg/details/linalg_fwd.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-core/options/matcl_options.h"
#include "matcl-linalg/linear_eq/petsc_enums.h"
#include "matcl-linalg/special_matrices/matrix_functors.h"
#include "matcl-linalg/norms_error/norm.h"
#include "matcl-matrep/matrix/permvec.h"

#include "matcl-linalg/linear_eq/ksp_solver.h"

#pragma warning( push )
#pragma warning( disable: 4251 ) // class 'shared_ptr<>' needs to have dll-interface to be used by clients

namespace matcl
{

/// build preconditioner to use with ksp_solver
class MATCL_LINALG_EXPORT precond : public ksp_solver
{
    public:
        /// uninitialized preconditioner
        precond();

        /// build standard preconditioner for a matrix P of given type
        precond(const Matrix& P, petsc_precond type, const options& opts = options());

        /// build preconditioner from a linsolve object
        precond(const linsolve_obj& P);

        /// standard destructor
        ~precond();

        /// string name of petsc preconditioner from code
        static std::string      precond_name(petsc_precond precond);

        /// code of petsc preconditioner from string name
        static petsc_precond    precond_code(const std::string& precond_name);
};

/// note that this is unsafe to use iterative method as preconditioner unless solution
/// given by such solver is linear function of RHS (this is not true for Krylov space
/// methods), if this condition is not satisfied then flexible Krylov solvers should be
/// used

//--------------------------------------------------------------------
//                   block preconditioner
//--------------------------------------------------------------------

/// build block preconditioner that uses block Jacobi method or additive Schwarz method
class MATCL_LINALG_EXPORT block_precond : public ksp_solver
{
    public:
        /// uninitialized preconditioner
        block_precond();

        /// divide preconditioner matrix P on blocks and use block Jacobi method;
        /// one may specify total number of equal blocks or size of each individual block;
        /// later preconditioner must be set for each block 
        block_precond(const Matrix& P, Integer n_blocks);
        block_precond(const Matrix& P, const Matrix& block_sizes);

        /// divide preconditioner matrix P on n blocks possibly overlaped and use
        /// additive Schwarz method; later preconditioner must be set for each block 
        block_precond(const Matrix& P, Integer n, Integer overlap,
            asm_composition comp = asm_composition::additive, asm_restrict restr = asm_restrict::restricted);

        /// divide preconditioner matrix P on n blocks possibly overlaped and use
        /// additive Schwarz method; number of blocks is given by size of the vector
        /// 'partition', each element of this vector gives assignment of matrix elements
        /// to given block and must be integer matrix with indices (1-based), all elements
        /// must be assigned, but may be assigned to many blocks; later preconditioner must
        /// be set for each block 
        block_precond(const Matrix& P, const std::vector<Matrix>& partition, Integer overlap, 
            asm_composition comp = asm_composition::additive, asm_restrict restr = asm_restrict::restricted);

        /// standard destructor
        ~block_precond();

        /// return i-th diagonal block of the matrix used to create preconditioner;
        Matrix                  get_block_precond_matrix(Integer i) const;

        /// return indices of elements that are assigned to i-th block
        Matrix                  get_partitioning(Integer i) const;

        /// set solver for i-th block of partitioned matrix
        void                    set_block_solver(Integer i, const linsolve_obj& solver);
        void                    set_block_solver(Integer i, const ksp_solver& solver);

        /// build solver for i-th block of partitioned matrix from options
        void                    build_block_solver(Integer i, const options& opts);

        /// return number of blocks of the preconditioner (0 if is_partitioned == false)
        Integer                 number_of_blocks() const;
};

//--------------------------------------------------------------------
//                   Field Split Preconditioner
//--------------------------------------------------------------------
/// divide preconditioner on blocks and assign preconditioners to
/// diagonal blocks
class MATCL_LINALG_EXPORT field_split_precond : public ksp_solver
{
    public:
        /// uninitialized preconditioner
        field_split_precond();

        /// divide preconditioner matrix P on blocks; blocks may overlap; later
        /// preconditioners must be set for each block
        field_split_precond(const Matrix& P, const std::vector<Matrix>& partitions,
            field_split_type fst);

        /// standard destructor
        ~field_split_precond();

        /// return i-th diagonal block of the matrix used to create preconditioner;
        Matrix                  get_block_precond_matrix(Integer i) const;

        /// return indices of elements that are assigned to i-th block
        Matrix                  get_partitioning(Integer i) const;

        /// set solver for i-th block of partitioned matrix
        void                    set_block_solver(Integer i, const linsolve_obj& solver);
        void                    set_block_solver(Integer i, const ksp_solver& solver);

        /// build solver for i-th block of partitioned matrix from options
        void                    build_block_solver(Integer i, const options& opts);

        /// return number of blocks of the preconditioner (0 if is_partitioned == false)
        Integer                 number_of_blocks() const;
};

//--------------------------------------------------------------------
//                   composite preconditioner
//--------------------------------------------------------------------

/// add or multitly preconditioners
class MATCL_LINALG_EXPORT composite_precond : public ksp_solver
{
    public:
        /// uninitialized preconditioner
        composite_precond();

        /// build a preconditioner by composing together several preconditioners
        composite_precond(const Matrix& P, petsc_composite_type type,
            const std::vector<ksp_solver>& precond_list);
        composite_precond(const Matrix& P, petsc_composite_type type,
            const std::vector<linsolve_obj>& precond_list);

        /// multiply two procenditioners, i.e. apply first M1 and then apply M2 on
        /// obtained result; see the note about setting Krylov solvers as preconditioners
        composite_precond(const Matrix& P, const ksp_solver& M1, const ksp_solver& M2);
        composite_precond(const Matrix& P, const linsolve_obj& M1, const linsolve_obj& M2);

        /// standard destructor
        ~composite_precond();
};

//--------------------------------------------------------------------
//                   Galerkin preconditioner
//--------------------------------------------------------------------

/// the preconditioner matrix is given by P*S*R, where P is N x k interpolation
/// matrix S is k x k operator matrix and R is k x N restriction matrix, k < N;
/// note that this preconditioner should not be used directly but only as a 
/// building block of other preconditioner

class MATCL_LINALG_EXPORT galerkin_precond : public ksp_solver
{
    public:
        /// uninitialized preconditioner
        galerkin_precond();

        /// set interpolation and restriction matrix, where P = R'; the PR argument
        /// can be the P matrix or the R matrix; function returns the matrix S constructed
        /// as S = PR' * B * PR (assuming PR has size N x k), where B is the preconditioner
        /// matrix assigned to this object; this matrix (or any other of size k x k)
        /// can be used to build ksp_solver, which must be set using set_galerkin_solver
        galerkin_precond(const Matrix& B, const Matrix& PR, Matrix& S);

        /// set interpolation and restriction matrix; return the matrix S constructed as
        /// S = P' * B * R', where B is the preconditioner matrix assigned to this object; 
        /// this matrix (or any other of size k x k) can be used to build ksp_solver, which
        /// must be set using set_galerkin_solver
        galerkin_precond(const Matrix& B, const Matrix& P, const Matrix& R, Matrix& S);

        /// set solver for the matrix S
        void                    set_galerkin_solver(const ksp_solver& S_solver);

        /// standard destructor
        ~galerkin_precond();
};

//--------------------------------------------------------------------
//                   Multigrid
//--------------------------------------------------------------------

/// multigrid preconditioner
class MATCL_LINALG_EXPORT multigrid : public ksp_solver
{
    public:
        /// uninitialized preconditioner
        multigrid();

        /// standard destructor
        ~multigrid();

        /// create multigrid preconditioner for a matrix B with lev levels (at least
        /// two); later interpolation and restriction operators as well as solver for
        /// coarser grid must be set for lev - 1 levels;
        /// one can set type of multigrid and type of cycle; (see mg_type and
        /// mg_cycle_type for details); for multiplicative mc sets the number of
        /// cycles to use for each preconditioner step of multigrid
        multigrid(const Matrix& B, Integer lev, mg_type type = mg_type::multiplicative, 
            mg_cycle_type cycle = mg_cycle_type::v, Integer mc = 1);

        //-----------------------------------------------------------------------------
        //                          BUILD GRID
        //-----------------------------------------------------------------------------
        /// one must call set_interp_restr and set_post_smoother, set_pre_smoother
        /// on every level except the coarsest and set_solver on the coarsest grid;
        /// calling sequence:
        ///     mg = multigrid(B);
        ///
        ///     p1 = make_smoother(B);   for some function make_smoother
        ///     S1 = mg.set_interp_restr(PR1);        
        ///     mg.set_post_smoother(p1); mg.set_pre_smoother(p1);
        ///     mg.coarser_level();
        ///
        ///     p2 = make_smoother(S1);   for some function make_smoother
        ///     S2 = mg.set_interp_restr_lev(S1, PR2);        
        ///     mg.set_post_smoother(p2); mg.set_pre_smoother(p2);
        ///     mg.coarser_level();
        ///
        ///     [...]
        ///
        ///     pk = make_solver(Sk); for some function make_solver
        ///     mg.set_solver(pk)

        /// set interpolation and restriction matrix for the finest level, where P = R'; 
        /// the PR argument can be the P matrix or the R matrix; function returns the matrix
        /// on coarser grid S constructed as S = PR' * B * PR (assuming PR has size N x k), 
        /// where B is the preconditioner matrix; this matrix S (or any other of size k x k) 
        /// can be used to build ksp_solver
        Matrix                  set_interp_restr(const Matrix& PR);

        /// set interpolation and restriction matrix for the finest level; return the matrix
        /// on coarser grid S constructed as S = P' * B * R', where B is the preconditioner 
        /// matrix; this matrix S (or any other of size k x k) can be used to build ksp_solver
        Matrix                  set_interp_restr(const Matrix& P, const Matrix& R);

        /// set interpolation and restriction matrix for current level; the matrix S is the
        /// matrix returned by set_interp_restr called previously and is used to construct
        /// a matrix S on coarser level;
        Matrix                  set_interp_restr_lev(const Matrix& S, const Matrix& PR);
        Matrix                  set_interp_restr_lev(const Matrix& S, const Matrix& P, const Matrix& R);

        /// set solver that is applied after interpolation from coarser grid (post-smoother)
        /// and number of performed steps; on the finest level the smoother should be 
        /// constructed from the matrix B supplied to the constructor; on coarser levels
        /// should be constructed for matrix passed to set_interp_restr on this level;
        /// smoother may be constructed for other matrix of the same size;
        /// set_post_smoother must be caller right after setting interpolation/restriction
        /// matrices; this function cannot be called on the coarsest grid
        /// generally smoother should have good smoothing property but can be given by any
        /// solver; if smoother is a Krylov solver, then one should consider using flexible 
        /// solvers (like fgremes, fcg) at top level
        void                    set_post_smoother(const ksp_solver& smoother, Integer n_steps = 2);

        /// set solver that is applied before restriction to coarser grid (pre-smoother)
        /// see set_post_smoother for details; usually this is the same smoother as supplied
        /// to set_post_smoother; this function cannot be called on the coarsest grid
        void                    set_pre_smoother(const ksp_solver& smoother, Integer n_steps = 2);

        /// go to next (coarser) level
        void                    coarser_level();

        /// set solver for the coarsest operator; matrix of this operator (based on Galerkin
        /// process) is returned by set_interp_restr, but can be differenct;
        /// usually S_solver is a direct solver
        void                    set_solver(const ksp_solver& S_solver);

        //-----------------------------------------------------------------------------
        //                          HELPERS
        //-----------------------------------------------------------------------------

        /// build Richardson smoother for given matrix using Jacobi preconditioner if
        /// jacobi = true or unpreconditioned version; use auto_scale = true to select
        /// damping factor omega optimally at each iteration (at small cost)
        static ksp_solver       richardson_smoother(const Matrix& S, Real omega = 2.0/3.0, 
                                    bool jacobi = false, bool auto_scale = false);
};

}

#pragma warning( pop )

