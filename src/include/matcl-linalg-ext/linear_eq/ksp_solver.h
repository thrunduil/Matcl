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

#pragma warning( push )
#pragma warning( disable: 4251 ) // class 'shared_ptr<>' needs to have dll-interface to be used by clients

namespace matcl
{

/// abstact class for inspecting convergense of petsc solvers for
class MATCL_LINALG_EXPORT ksp_solver_notifier
{
    private:
        bool            m_least_squares;

	public:
        ksp_solver_notifier();

        virtual ~ksp_solver_notifier() {}; 

        /// return true if least squares problem is solved
        bool            is_least_squares() const;

        /// order of functions defines order of events        

        /// begin solve for n_vec right-hand-side vectors
        virtual void    begin_solve(const ksp_solver& solver, const linear_operator& op, Integer n_vec) = 0;

        /// begin solve for i-th RHS vector
        virtual void    begin_solve_vector(Integer i) = 0;

        /// begin iterations for given RHS vector
        virtual void    begin_iterations(const Matrix& RHS) = 0;        
            
        /// report iteration itc and (estimated) 2-norm of (preconditioned) residual
		virtual void	report_iteration(Integer itc, matcl::Real resid_norm, KSP ksp) = 0;

        /// report iteration itc for least square problems with (estimated) 2-norm of
        /// (preconditioned) residual resid_norm and of residual of normal equations 
        /// resid_anorm
		virtual void	report_iteration_ls(Integer itc, matcl::Real resid_norm, matcl::Real resid_anorm, 
                                            KSP ksp) = 0;

        /// solver did not converged
        virtual void    report_divergence(const std::string& str) = 0;

        /// solver converged
        virtual void    report_convergence(const std::string& str) = 0;

        /// error occured
        virtual void    report_error(const std::string& msg) = 0;

        /// iterations for op(A) * sol = RHS finished; operator op_A represents op(A)
        virtual void    end_iterations(const linear_operator& op_A, const Matrix& RHS, const Matrix& sol) = 0;

        /// solve for given RHS vector finished
        virtual void    end_solve_vector() = 0;

        /// begin solve for all RHS vectors
        virtual void    end_solve() = 0;

    protected:
        /// these functions can only be called when KSP object is available;
        /// KSP object cannot be stored

        /// calculare true residuals
        Matrix          get_residuals(KSP ksp) const;

        /// get right-hand-side vector
        Matrix          get_rhs(KSP ksp) const;

        /// build solution at given iteration
        Matrix          get_solution(KSP ksp) const;

        /// calculate norm of true residuals
        Real            true_resid_norm(KSP ksp, basic_vector_norm norm_type) const;

        /// calculate norm of RHS vector
        Real            rhs_norm(KSP ksp, basic_vector_norm norm_type) const;

    private:
        void            set_least_squares(bool is_ls);

        friend details::ksp_solver_impl;
};

/// wrapper for the Petsc's KSP algorithms for solving nonsingular linear equations:
///         op(A) * X = B,                  (1)
/// where op(A) = A, trans(A), or ctrans(A).
/// Petsc is not thread safe, only one thread can use ksp_solver otherwise exception
/// will be throw; ksp_solver is not copy-on-wite object as most of mmlib objects - 
/// any change in one object will be visible in all references to that object
class MATCL_LINALG_EXPORT ksp_solver
{
    private:
        using impl_type         = details::ksp_solver_impl;
        using impl_ptr          = std::shared_ptr<impl_type>;

        friend details::ksp_solver_impl;

    public:
        /// type of shared poiner to notifier
        using notifier          = std::shared_ptr<ksp_solver_notifier>;

    protected:
        impl_ptr                m_impl;

    public: 
        //--------------------------------------------------------------------
        //                      building the solver
        //--------------------------------------------------------------------

        /// prepare iterative solver for a matrix A = 1.0
        ksp_solver();

        /// prepare iterative solver for a matrix or linear operator A; matrix A is
        /// not internally balanced; any balancing should be performed externally; 
        /// selects solvers and preconditioners according to options, initialize 
        /// internal data structure and possibly make factorizations if required; 
        /// see opt::petsc namespace for details, the most important option are
        /// solver and preconditioner; see petsc_solver and petsc_precond enumerations
        /// for list of available solvers and preconditioners
        ksp_solver(const Matrix& A, const options& opts = options());
        explicit ksp_solver(const linear_operator& A, const options& opts = options());

        /// convert linsolve object to ksp_solver
        ksp_solver(const linsolve_obj& A);

        /// prepare iterative solver for a matrix or linear operator A;  preconditioner
        /// is constructed based on the matrix B; the preconditioner can also be set
        /// later using set_preconditioner function
        ksp_solver(const Matrix& A, const Matrix& B, const options& opts = options());
        ksp_solver(const linear_operator& A, const Matrix& B, const options& opts = options());

        /// standard destructor
        ~ksp_solver();

        /// prepare iterative solver for new matrix or linear operator; see 
        /// constructors for details; Petsc KSP objects returned by get_KSP function
        /// are invalidated
        ksp_solver&              operator()(const Matrix& A, const options& opts);
        ksp_solver&              operator()(const linear_operator& A, const options& opts);
        ksp_solver&              operator()(const Matrix& A, const Matrix& B, const options& opts);
        ksp_solver&              operator()(const linear_operator& A, const Matrix& B, const options& opts);

        /// prepare iterative solver for new matrix or linear operator; use
        /// options associated with this object; see constructors for details;
        /// Petsc KSP objects returned by get_KSP function are invalidated
        ksp_solver&              operator()(const Matrix& A);
        ksp_solver&              operator()(const linear_operator& A);
        ksp_solver&              operator()(const Matrix& A, const Matrix& B);
        ksp_solver&              operator()(const linear_operator& A, const Matrix& B);

        /// change operator A but do not change preconditioner; Petsc KSP objects
        /// returned by get_KSP function are still valid; any null spaces set
        /// are removed; size of new_A must be the same as size of the operator A
        void                    reset_operator(const Matrix& new_A);
        void                    reset_operator(const linear_operator& new_A);

        /// build new preconditioner from the matrix B keeping operator A
        /// unchanged; if options are not given, then options from this object
        /// will be used; this function is logically equivalent to operator()(A,B)
        /// but does not create new Petsc matrix for the operator A and doest not
        /// remove null spaces; Petsc KSP objects returned by get_KSP function are
        /// invalidated
        void                    rebuild_precond(const Matrix& B);
        void                    rebuild_precond(const Matrix& B, const options& opts);

        /// change options; this may require complete reinitialization or possibly
        /// new making factorizations, threfore can be very costly;
        /// changing options governing convergence checks is fast
        void                    set_options(const options& opts);

        /// set convergence tolerances; maxit: maximum number of iterations
        /// rtol, atol, dtol: see option atol, rtol and dtol; use -2 to set
        /// default value
        void                    set_tolerances(Integer maxit = -2, Real rtol = -2.0, 
                                    Real atol = -2.0, Real dtol = -2.0);

        /// change notifier; notifier can be null
        void                    set_notifier(const notifier& notifier);

        /// inform the solver about null space; Nr must form orthonormal basis
        /// of the right null space of A (i.e. A * Nr = 0, Nr' * Nr = I); if 
        /// with_const is true, then additionally A * e = 0, where e_i = 1 forall i;
        /// setting null space may improve convergence, solution is searched only
        /// is space orthogonal to the null space; notice that preconditioners may
        /// still break down for singular problems
        void                    set_nullspace_right(const Matrix& Nr, bool with_const = false);

        /// inform the solver left null space; Nl' must form orthonormal basis
        /// of the right null space of A (i.e. Nl' * A = 0, Nl' * Nl = I); if 
        /// with_const is true, then additionally e' * A = 0, where e_i = 1 forall i
        void                    set_nullspace_left(const Matrix& Nr, bool with_const = false);

        /// remove associated right null space
        void                    remove_nullspace_right();

        /// remove associated left null space
        void                    remove_nullspace_left();

        /// set preconditioner given as a linear operator (or a matrix), i.e. the 
        /// preconditioning matrix should approximate the inverse of the operator A
        void                    set_precond(const linear_operator& precond);

        /// set preconditioner given by linsolve object
        void                    set_precond(const linsolve_obj& precond);

        /// set preconditioner given by iterative solver; note that this is unsafe to
        /// use iterative method as preconditioner unless solution given by such solver
        /// is linear function of RHS (this is not true for Krylov space methods), if
        /// this condition is not satisfied then flexible Krylov solvers should be used
        void                    set_precond(const ksp_solver& precond);

        /// string identifier passed to notifier
        void                    set_id(const std::string& id);

        /// get associated identifier
        const std::string&      get_id() const;

        //--------------------------------------------------------------------
        //          const functions defined for all matrix types
        //--------------------------------------------------------------------

        /// value code of elements stored in matrix A
        value_code              get_value_code() const;

        /// number of rows of the matrix A
        Integer                 rows() const;

        /// number of columns of the matrix A
        Integer                 cols() const;

        /// check if all elements are finite
        bool                    all_finite() const;

        //--------------------------------------------------------------------
        //                  solve functions
        //--------------------------------------------------------------------

        /// solve the problem (1), the rhs matrix b can have many columns; return
        /// current estimation of solution and flag, whether solver converged for all
        /// vectors; more detailed convergence report can be constructed from information
        /// sent to notifier
        tuple<Matrix,bool>      solve(const matcl::Matrix& b, trans_type t = trans_type::no_trans) const; 

        /// version of solve function; initial guess is supplied
        tuple<Matrix,bool>      solve(const matcl::Matrix& b,const Matrix& start, 
                                    trans_type t = trans_type::no_trans) const;

        //------------------------------------------------------------------
        //                      UTILITIES
        //------------------------------------------------------------------

        /// get solver type
        std::string             solver_name() const;

        /// ret preconditioner
        std::string             precond_name() const;

        /// string name of petsc solver from code
        static std::string      solver_name(petsc_solver sol);

        /// code of petsc solver from string name
        static petsc_solver     solver_code(const std::string& sol_name);

        /// return options that take default values set by petsc
        options                 options_missing() const;

        /// return options that was used by solve and have defined values
        options                 options_set() const;

        /// return stored petsc's KSP object
        KSP                     get_KSP() const;

        /// get linear operator associated with this object
        linear_operator         get_linear_operator() const;

        /// get matrix associated with this object; exceptions is thrown if
        /// ksp_solver is initialized with a linear operator
        Matrix                  get_matrix() const;

        /// return true if this solver was constructed for a matrix
        bool                    is_from_matrix() const;

        /// get preconditioner matrix associated with this object; exceptions is
        /// thrown if precondioner is not given by a matrix
        Matrix                  get_precond_matrix() const;

        /// return true if preconditioner is given by a matrix
        bool                    is_precond_from_matrix() const;

        /// construct linsolve object from this solver
        linsolve_obj            to_linsolve_object2() const;

        void                    set_option(const std::string& opt, const std::string& val);
};

/// solve least squares problem 
///     min |op(A) * X - B|_2 
/// where A is square or rectangular (overdetermined or underdetermined) and may 
/// have any rank; such problem can also be solved by ksp_solver
class MATCL_LINALG_EXPORT ksp_solver_ls : public ksp_solver
{
    public: 
        /// prepare iterative solver for a matrix A = 1.0
        ksp_solver_ls();

        /// prepare iterative solver for a matrix or linear operator A; see
        /// constructors in ksp_solver for details; (least squares solvers
        /// do not support preconditioning)
        ksp_solver_ls(const Matrix& A, const options& opts = options());
        ksp_solver_ls(Matrix&& A, const options& opts = options());
        ksp_solver_ls(const linear_operator& A, const options& opts = options());

        /// standard destructor
        ~ksp_solver_ls();

        /// get estimated standard error of solution if requested (approximation to 
        /// unbiased estimator given by std(Y)_i = s * sqrt( (A'*A)^-1_ii ), and
        /// s = r'*r / (N - K), r = A * Y - B ); see opt::petsc::lsqr_std_est for 
        /// details; resulting vector has size K x 1, where K is number of columns 
        /// of A (for standard problem) or number of rows of A (for transposed problems)
        Matrix                  get_solution_std() const;
};

}

#pragma warning( pop )

