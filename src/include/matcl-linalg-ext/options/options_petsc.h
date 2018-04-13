/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2015-2016
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
#include "matcl-core/options/matcl_options.h"
#include "matcl-linalg/linear_eq/petsc_enums.h"

namespace matcl { namespace opt { namespace petsc
{

//---------------------------------------------------------------
//                  CONVERGENCE TEST
//---------------------------------------------------------------
/// maximum number of iterations; 
/// corresponding string option: "-ksp_max_it"
class max_it : public option_base<Integer, max_it>
{
    private:
        using base_type         = option_base<Integer, max_it>;
        using opt_type          = optional<Integer>;

    public:
        max_it()                : base_type() {};
        max_it(opt_type x)      : base_type(x) {};

        static void config()
        {
            m_description       = "maximum number of iterations allowed";
            m_default_value     = -1;
            m_validator         = validator_positive_or_val<Integer>(-1);
        };
};

/// convergence test based on absolute value of residual norm;
/// by default for left preconditioning it is the 2-norm of the preconditioned
/// residual, and the 2-norm of the residual for right preconditioning
/// corresponding string option: "-ksp_atol"
class atol : public option_base<Real, atol>
{
    private:
        using base_type         = option_base<Real, atol>;

    public:
        atol()                  : base_type() {};
        atol(opt_type x)        : base_type(x) {};

        static void config()
        {
            m_description       = "convergence test based on absolute value of residual norm";
            m_default_value     = Real(-1.0);
            m_validator         = validator_nonnegative_or_val<Real>(-1.0);
        };
};

/// convergence test based on relative value of residual norm; iterations 
/// are stopped if norm(r) <= rtol*norm(b), where r are residuals and b is 
/// rhs vector; by default for left preconditioning it is the 2-norm of the
/// preconditioned residual, and the 2-norm of the residual for right 
/// preconditioning
/// corresponding string option: "-ksp_rtol"
class rtol : public option_base<Real, rtol>
{
    private:
        using base_type         = option_base<Real, rtol>;

    public:
        rtol()                  : base_type() {};
        rtol(opt_type x)        : base_type(x) {};

        static void config()
        {
            m_description       = "convergence test based on relative value of residual norm";
            m_default_value     = Real(-1.0);
            m_validator         = validator_nonnegative_or_val<Real>(-1.0);
        };
};

/// if residuals norm satisfies norm(r) >= divtol*norm(b) then divergense is
/// declared, where r are residuals and b is rhs vector; by default for left 
/// preconditioning it is the 2-norm of the preconditioned residual, and the 
/// 2-norm of the residual for right preconditioning
/// corresponding string option: "-ksp_divtol"
class divtol : public option_base<Real, divtol>
{
    private:
        using base_type         = option_base<Real, divtol>;

    public:
        divtol()                : base_type() {};
        divtol(opt_type x)      : base_type(x) {};

        static void config()
        {
            m_description       = "divergence test based on increase of residual norm";
            m_default_value     = Real(-1.0);
            m_validator         = validator_nonnegative_or_val<Real>(-1.0);
        };
};

/// norm that is passed in the Krylov convergence test routines
/// corresponding string option: "-ksp_norm_type"
class norm_type : public option_base<Integer, norm_type>
{
    private:
        using base_type         = option_base<Integer, norm_type>;
        using opt_type          = optional<Integer>;

    public:
        norm_type()                 : base_type() {};
        norm_type(opt_type x)       : base_type(x) {};
        norm_type(petsc_conv_norm x): base_type((Integer)x) {};

        static void config()
        {
            m_description       = "norm that is passed in the Krylov convergence test routines";
            m_default_value     = (Integer)petsc_conv_norm::default_val;
            m_validator         = validator_enum((Integer)petsc_conv_norm::last);
        };
};

//---------------------------------------------------------------
//                  SOLVERT CONFIGURATION
//---------------------------------------------------------------
/// solver type
/// corresponding string option: "-ksp_type"
class solver : public option_base<Integer, solver>
{
    private:
        using base_type         = option_base<Integer, solver>;
        using opt_type          = optional<Integer>;

    public:
        solver()                : base_type() {};
        solver(opt_type x)      : base_type(x) {};
        solver(petsc_solver x)  : base_type((Integer)x) {};

        static void config()
        {
            m_description       = "petsc ksp solver";
            m_default_value     = (Integer)petsc_solver::default_val;
            m_validator         = validator_enum((Integer)petsc_solver::last);
        };
};

//---------------------------------------------------------------
//                      SMOOTHERS
//---------------------------------------------------------------
/// set approximations to the smallest and largest eigenvalues of the
/// preconditioned operator; if these are accurate you will get much 
/// faster convergence.
/// corresponding string option: "-ksp_chebyshev_eigenvalues"
class cheb_eig : public option_base<std::vector<Real>, cheb_eig>
{
    private:
        using base_type         = option_base<std::vector<Real>, cheb_eig>;
        using opt_type          = optional<std::vector<Real>>;

    public:
        cheb_eig()              : base_type() {};
        cheb_eig(opt_type x)    : base_type(x) {};
        cheb_eig(std::initializer_list<Real> init)
                                : base_type(std::vector<Real>(init)){};

        static void config()
        {
            m_description       = "set approximations to the smallest and largest eigenvalues"
                                  " of the preconditioned operator";
            m_default_value     = std::vector<Real>();
            m_validator         = [](opt_type x)->opt_type
                                {
                                    if (x && x.value().size()!= 0 
                                        && (x.value().size() != 2 || x.value()[0] > x.value()[1]))
                                    {
                                        std::ostringstream os;
                                        os << "invalid option: invalid minimum and maximum eigenvalues";
                                        throw std::runtime_error(os.str());
                                    }
                                    return x;
                                };

        };
};

/// automatically estimate the eigenvalues to use for Chebyshev; 
/// the Chebyshev bounds are set using 
///     minbound = v(1)*minest + v(2)*maxest
///     maxbound = v(4)*minest + v(4)*maxest
/// where minest, maxest are minimum and maximum eigenvalues estimated
/// using Lanczos or Arnoldi iterations
/// for example transform (0, 0.1; 0, 1.1) will targets the "upper" part of the spectrum, 
/// as desirable for use with multigrid as a smoother.
/// corresponding string option: "-ksp_chebyshev_esteig"
class cheb_eig_est : public option_base<std::vector<Real>, cheb_eig_est>
{
    private:
        using base_type         = option_base<std::vector<Real>, cheb_eig_est>;
        using opt_type          = optional<std::vector<Real>>;

    public:
        cheb_eig_est()              : base_type() {};
        cheb_eig_est(opt_type x)    : base_type(x) {};
        cheb_eig_est(std::initializer_list<Real> init)
                                    : base_type(std::vector<Real>(init)){};

        static void config()
        {
            m_description       = "automatically estimate the eigenvalues to use for Chebyshev";
            m_default_value     = std::vector<Real>();
            m_validator         = [](opt_type x)->opt_type
                                {
                                    if (x && x.value().size()!= 0 && (x.value().size() != 4))
                                    {
                                        std::ostringstream os;
                                        os << "invalid option: invalid eigenvalue filtering parameters";
                                        throw std::runtime_error(os.str());
                                    }
                                    return x;
                                };

        };
};

/// number of iterations in eigenvalue estimation for chebyshev
/// method
/// corresponding string option: "-ksp_chebyshev_esteig_steps"
class cheb_eig_steps : public option_base<Integer, cheb_eig_steps>
{
    private:
        using base_type         = option_base<Integer, cheb_eig_steps>;

    public:
        cheb_eig_steps()            : base_type() {};
        cheb_eig_steps(opt_type x)  : base_type(x) {};

        static void config()
        {
            m_description       = "number of iterations in eigenvalue estimation";
            m_default_value     = -1;
            m_validator         = validator_positive_or_val<Integer>(-1);
        };
};

/// set the damping factor
/// corresponding string option: "-ksp_richardson_scale"
class richardson_damping : public option_base<Real, richardson_damping>
{
    private:
        using base_type         = option_base<Real, richardson_damping>;

    public:
        richardson_damping()            : base_type() {};
        richardson_damping(opt_type x)  : base_type(x) {};

        static void config()
        {
            m_description       = "set the damping factor";
            m_default_value     = -1.0;
            m_validator         = validator_nonnegative_or_val<Real>(-1.0);
        };
};

/// sets Richardson to automatically determine optimal damping at each iteration to minimize
/// the 2-norm of the preconditioned residual
/// "-ksp_richardson_self_scale"
class richardson_auto : public option_base<Integer, richardson_auto>
{
    private:
        using base_type         = option_base<Integer, richardson_auto>;

    public:
        richardson_auto()           : base_type() {};
        richardson_auto(opt_type x) : base_type(x) {};
        richardson_auto(bool x)     : base_type(x ? 1 : 0) {};
        richardson_auto(Integer x)  : base_type(x == -1 ? -1 : (x == 0 ? 0 : 1)) {};

        static void config()
        {
            m_description       = "sets Richardson to automatically determine optimal damping factor";
            m_default_value     = -1;
            m_validator         = validator_range<Integer>(-1,1);
        };
};

//---------------------------------------------------------------
//                  GMRES FAMILLY
//---------------------------------------------------------------

/// the number of Krylov directions to orthogonalize against;
/// corresponding string option: "-ksp_gmres_restart", "-ksp_gcr_restart"
class gmres_restart : public option_base<Integer, gmres_restart>
{
    private:
        using base_type         = option_base<Integer, gmres_restart>;
        using opt_type          = optional<Integer>;

    public:
        gmres_restart()            : base_type() {};
        gmres_restart(opt_type x)  : base_type(x) {};

        static void config()
        {
            m_description       = "the number of Krylov directions to orthogonalize against";
            m_default_value     = -1;
            m_validator         = validator_positive_or_val<Integer>(-1);
        };
};

/// set tolerance for determining happy breakdown;
/// happy breakdown is the rare case in GMRES where an 'exact' solution is obtained
/// after a certain number of iterations; if you attempt more iterations after this
/// point unstable things can happen hence very occasionally you may need to set this
/// value to detect this condition 
/// corresponding string option: "-ksp_gmres_haptol", "-ksp_lcd_haptol"
class gmres_haptol : public option_base<Real, gmres_haptol>
{
    private:
        using base_type         = option_base<Real, gmres_haptol>;
        using opt_type          = optional<Real>;

    public:
        gmres_haptol()              : base_type() {};
        gmres_haptol(opt_type x)    : base_type(x) {};

        static void config()
        {
            m_description       = "set tolerance for determining happy breakdown";
            m_default_value     = -1.0;
            m_validator         = validator_nonnegative_or_val<Real>(-1.0);
        };
};

/// use classical or modified Gram-Schmidt to orthogonalize against the Krylov space
/// if true, then unmodifed Gram-Schmidt procedure is used, otherwise modifed Gram-Schmidt
/// is used, which is more stable but slower
/// corresponding string options: "-ksp_gmres_classicalgramschmidt", 
/// "-ksp_gmres_modifiedgramschmidt"
class gmres_gs : public option_base<Integer, gmres_gs>
{
    private:
        using base_type         = option_base<Integer, gmres_gs>;
        using opt_type          = optional<Integer>;

    public:
        gmres_gs()              : base_type() {};
        gmres_gs(opt_type x)    : base_type(x) {};
        gmres_gs(bool x)        : base_type(x ? 1 : 0) {};
        gmres_gs(Integer x)     : base_type(x == -1 ? -1 : (x == 0 ? 0 : 1)) {};

        static void config()
        {
            m_description       = "use classical or modified Gram-Schmidt to orthogonalize against the Krylov space";
            m_default_value     = -1;
            m_validator         = validator_range<Integer>(-1,1);
        };
};

/// determine if iterative refinement is used to increase the stability 
/// of the classical Gram-Schmidt orthogonalization.
/// corresponding string option: "-ksp_gmres_cgs_refinement_type"
class gmres_refine : public option_base<Integer, gmres_refine>
{
    private:
        using base_type         = option_base<Integer, gmres_refine>;
        using opt_type          = optional<Integer>;

    public:
        gmres_refine()                  : base_type() {};
        gmres_refine(opt_type x)        : base_type(x) {};
        gmres_refine(orth_refine_type x): base_type((Integer)x) {};

        static void config()
        {
            m_description       = "iterative refinement used to increase the stability"
                                  " of the classical Gram-Schmidt orthogonalization";
            m_default_value     = (Integer)orth_refine_type::default_val;
            m_validator         = validator_enum((Integer)orth_refine_type::last);
        };
};

/// dgmres option
/// maximum number of eigenvalues that can be extracted during the iterative process
/// corresponding string option: "-ksp_dgmres_max_eigen"
class dgmres_max_eigen : public option_base<Integer, dgmres_max_eigen>
{
    private:
        using base_type         = option_base<Integer, dgmres_max_eigen>;

    public:
        dgmres_max_eigen()                  : base_type() {};
        dgmres_max_eigen(opt_type x)        : base_type(x) {};

        static void config()
        {
            m_description       = "maximum number of eigenvalues that can be extracted during the iterative process";
            m_default_value     = -1;
            m_validator         = validator_nonnegative_or_val<Integer>(-1);
        };
};

/// dgmres option
/// number of smallest eigenvalues to extract at each restart
/// corresponding string option: "-ksp_dgmres_eigen"
class dgmres_eigen : public option_base<Integer, dgmres_eigen>
{
    private:
        using base_type         = option_base<Integer, dgmres_eigen>;

    public:
        dgmres_eigen()              : base_type() {};
        dgmres_eigen(opt_type x)    : base_type(x) {};

        static void config()
        {
            m_description       = "number of smallest eigenvalues to extract at each restart";
            m_default_value     = -1;
            m_validator         = validator_nonnegative_or_val<Integer>(-1);
        };
};

/// lgmres option
/// number of error approximations to augment the Krylov space with
/// corresponding string option: "-ksp_lgmres_augment"
class lgmres_augment : public option_base<Integer, lgmres_augment>
{
    private:
        using base_type         = option_base<Integer, lgmres_augment>;

    public:
        lgmres_augment()            : base_type() {};
        lgmres_augment(opt_type x)  : base_type(x) {};

        static void config()
        {
            m_description       = "number of error approximations to augment the Krylov space with";
            m_default_value     = -1;
            m_validator         = validator_nonnegative_or_val<Integer>(-1);
        };
};

//---------------------------------------------------------------
//                   CG AND BICG FAMILY
//---------------------------------------------------------------

/// the number of Krylov directions to orthogonalize against;
/// corresponding string option: "-ksp_lcd_restart", "-ksp_fcg_mmax"
using cg_restart = gmres_restart;

/// number of Krylov search directions for bcgsl solver;
/// for large num_dir it is common for the polynomial update problem to become singular
/// (due to happy breakdown for smallish test problems, but also for larger problems). 
/// consequently, by default, the system is solved by pseudoinverse, which allows the 
/// iteration to complete successfully
/// corresponding string option: "-ksp_bcgsl_ell"
class bcgsl_num_dir : public option_base<Integer, bcgsl_num_dir>
{
    private:
        using base_type         = option_base<Integer, bcgsl_num_dir>;

    public:
        bcgsl_num_dir()            : base_type() {};
        bcgsl_num_dir(opt_type x)  : base_type(x) {};

        static void config()
        {
            m_description       = "number of Krylov search directions";
            m_default_value     = -1;
            m_validator         = validator_positive_or_val<Integer>(-1);
        };
};

/// strategy of orthogonalization a new direction against previous directions
/// corresponding string option: "-ksp_fcg_truncation_type"
class fcg_orthog : public option_base<Integer, fcg_orthog>
{
    private:
        using base_type         = option_base<Integer, fcg_orthog>;
        using opt_type          = optional<Integer>;

    public:
        fcg_orthog()                : base_type() {};
        fcg_orthog(opt_type x)      : base_type(x) {};
        fcg_orthog(fcg_orth_type x) : base_type((Integer)x) {};

        static void config()
        {
            m_description       = "orthogonalization strategy";
            m_default_value     = (Integer)fcg_orth_type::default_val;
            m_validator         = validator_enum((Integer)fcg_orth_type::last);
        };
};

//---------------------------------------------------------------
//                      LEAST SQUARES
//---------------------------------------------------------------
/// if true then estimate standard error of solution and estimated
/// standard error can be accessed from the solver
/// corresponding string options: "-ksp_lsqr_set_standard_error", 
class lsqr_std_est : public option_base<Integer, lsqr_std_est>
{
    private:
        using base_type         = option_base<Integer, lsqr_std_est>;
        using opt_type          = optional<Integer>;

    public:
        lsqr_std_est()              : base_type() {};
        lsqr_std_est(opt_type x)    : base_type(x) {};
        lsqr_std_est(bool x)        : base_type(x ? 1 : 0) {};
        lsqr_std_est(Integer x)     : base_type(x == -1 ? -1 : (x == 0 ? 0 : 1)) {};

        static void config()
        {
            m_description       = "if true then estimate standard error of solution";
            m_default_value     = -1;
            m_validator         = validator_range<Integer>(-1,1);
        };
};

//---------------------------------------------------------------
//                      TRUST REGION SOLVERS
//---------------------------------------------------------------
/// Trust Region radius
/// corresponding string option: "-ksp_gltr_radius", "-ksp_nash_radius",
/// "-ksp_stcg_radius"
class trust_region_radius : public option_base<Real, trust_region_radius>
{
    private:
        using base_type         = option_base<Real, trust_region_radius>;

    public:
        trust_region_radius()           : base_type() {};
        trust_region_radius(opt_type x) : base_type(x) {};

        static void config()
        {
            m_description       = "Trust Region radius";
            m_default_value     = -1.0;
            m_validator         = validator_nonnegative_or_val<Real>(-1.0);
        };
};

//---------------------------------------------------------------
//                      PRECONDITIONERS
//---------------------------------------------------------------

/// preconditioner
/// corresponding string option: "-pc_type"
class preconditioner : public option_base<Integer, preconditioner>
{
    private:
        using base_type         = option_base<Integer, preconditioner>;
        using opt_type          = optional<Integer>;

    public:
        preconditioner()                : base_type() {};
        preconditioner(opt_type x)      : base_type(x) {};
        preconditioner(petsc_precond x) : base_type((Integer)x) {};

        static void config()
        {
            m_description       = "petsc ksp preconditioner";
            m_default_value     = (Integer)petsc_precond::default_val;
            m_validator         = validator_enum((Integer)petsc_precond::last);
        };
};

/// set the preconditioning side (left, right or symmetric); note that side
/// must be valid for given solver
/// corresponding string option: "-ksp_pc_side"
class precond_side : public option_base<Integer, precond_side>
{
    private:
        using base_type         = option_base<Integer, precond_side>;

    public:
        precond_side()                      : base_type() {};
        precond_side(opt_type x)            : base_type(x) {};
        precond_side(petsc_precond_side x)  : base_type((Integer)x) {};

        static void config()
        {
            m_description       = "set the preconditioning side";
            m_default_value     = (Integer)petsc_precond_side::default_val;
            m_validator         = validator_enum((Integer)petsc_precond_side::last);
        };
};

/// relaxation factor for successive over-relaxation algorithms
/// corresponding string option: "-pc_sor_omega", "-pc_eisenstat_omega"
class sor_omega : public option_base<Real, sor_omega>
{
    private:
        using base_type         = option_base<Real, sor_omega>;

    public:
        sor_omega()             : base_type() {};
        sor_omega(opt_type x)   : base_type(x) {};

        static void config()
        {
            m_description       = "relaxation factor";
            m_default_value     = -1.0;
            m_validator         = validator_range_or_val<Real>(-1.0, 0.0, 2.0, true, true);
        };
};

/// type of relaxation (backward, forward, symmetric)
/// corresponding string option: "-pc_sor_symmetric", "-pc_sor_backward",
/// "-pc_sor_forward"
class sor_type : public option_base<Integer, sor_type>
{
    private:
        using base_type         = option_base<Integer, sor_type>;

    public:
        sor_type()                  : base_type() {};
        sor_type(opt_type x)        : base_type(x) {};
        sor_type(petsc_sor_type x)  : base_type((Integer)x) {};

        static void config()
        {
            m_description       = "type of relaxation (backward, forward, symmetric)";
            m_default_value     = (Integer)petsc_sor_type::default_val;
            m_validator         = validator_enum((Integer)petsc_sor_type::last);
        };
};

/// number of inner iterations to be used by SOR-type preconditioners
/// corresponding string options: "-pc_sor_its", "-pc_sor_lits"
class sor_iter : public option_base<Integer, sor_iter>
{
    private:
        using base_type         = option_base<Integer, sor_iter>;

    public:
        sor_iter()                 : base_type() {};
        sor_iter(opt_type x)       : base_type(x) {};

        static void config()
        {
            m_description       = "number of inner iterations to be used by SOR-type preconditioners";
            m_default_value     = -1;
            m_validator         = validator_positive_or_val<Integer>(-1);
        };
};

/// type of Jacobi preconditioner (based on diagonal entries, maximums in
/// each row or sum of rows)
/// corresponding string option: "-pc_jacobi_abs", "-pc_jacobi_rowmax",
/// "-pc_jacobi_rowsum"
class jacobi_type : public option_base<Integer, jacobi_type>
{
    private:
        using base_type         = option_base<Integer, jacobi_type>;

    public:
        jacobi_type()                   : base_type() {};
        jacobi_type(opt_type x)         : base_type(x) {};
        jacobi_type(petsc_jacobi_type x): base_type((Integer)x) {};

        static void config()
        {
            m_description       = "type of Jacobi preconditioner";
            m_default_value     = (Integer)petsc_jacobi_type::default_val;
            m_validator         = validator_enum((Integer)petsc_jacobi_type::last);
        };
};

/// use the absolute value of entry defined according to jacobi_type
/// corresponding string options: "-pc_jacobi_abs", 
class jacobi_abs : public option_base<Integer, jacobi_abs>
{
    private:
        using base_type         = option_base<Integer, jacobi_abs>;
        using opt_type          = optional<Integer>;

    public:
        jacobi_abs()            : base_type() {};
        jacobi_abs(opt_type x)  : base_type(x) {};
        jacobi_abs(bool x)      : base_type(x ? 1 : 0) {};
        jacobi_abs(Integer x)   : base_type(x == -1 ? -1 : (x == 0 ? 0 : 1)) {};

        static void config()
        {
            m_description       = "use the absolute value of entry (diagonal, maxrow, or sumrows)";
            m_default_value     = -1;
            m_validator         = validator_range<Integer>(-1,1);
        };
};

/// relaxation factor for kaczmarz_omega algorithm
/// corresponding string option: "-pc_kaczmarz_lambda"
class kaczmarz_omega : public option_base<Real, kaczmarz_omega>
{
    private:
        using base_type         = option_base<Real, kaczmarz_omega>;

    public:
        kaczmarz_omega()             : base_type() {};
        kaczmarz_omega(opt_type x)   : base_type(x) {};

        static void config()
        {
            m_description       = "relaxation factor";
            m_default_value     = -1.0;
            m_validator         = validator_positive_or_val<Real>(-1.0);
        };
};

/// apply row projections symmetrically in Kaczmarz algorithm
/// corresponding string options: "-pc_kaczmarz_symmetric", 
class kaczmarz_sym : public option_base<Integer, kaczmarz_sym>
{
    private:
        using base_type         = option_base<Integer, kaczmarz_sym>;

    public:
        kaczmarz_sym()            : base_type() {};
        kaczmarz_sym(opt_type x)  : base_type(x) {};
        kaczmarz_sym(bool x)      : base_type(x ? 1 : 0) {};
        kaczmarz_sym(Integer x)   : base_type(x == -1 ? -1 : (x == 0 ? 0 : 1)) {};

        static void config()
        {
            m_description       = "apply row projections symmetrically";
            m_default_value     = -1;
            m_validator         = validator_range<Integer>(-1,1);
        };
};

//---------------------------------------------------------------
//                      FACTORIZATIONS
//---------------------------------------------------------------

/// parameter k in ILU(k) (incomplete LU) factorization - nonzero
/// structure of factorized matrix A is determined based on nonzero
/// structure of A^(k+1)
/// corresponding string options: "-pc_factor_levels", 
class ilu_k : public option_base<Integer, ilu_k>
{
    private:
        using base_type         = option_base<Integer, ilu_k>;

    public:
        ilu_k()                 : base_type() {};
        ilu_k(opt_type x)       : base_type(x) {};

        static void config()
        {
            m_description       = "level of fills for ilu - nonzero structure based on A^(k+1)";
            m_default_value     = -1;
            m_validator         = validator_nonnegative_or_val<Integer>(-1);
        };
};

/// allow fill into empty diagonal entry
/// corresponding string options: "-pc_factor_diagonal_fill", 
class allow_diagonal_fill : public option_base<Integer, allow_diagonal_fill>
{
    private:
        using base_type         = option_base<Integer, allow_diagonal_fill>;

    public:
        allow_diagonal_fill()           : base_type() {};
        allow_diagonal_fill(opt_type x) : base_type(x) {};
        allow_diagonal_fill(bool x)     : base_type(x ? 1 : 0) {};
        allow_diagonal_fill(Integer x)  : base_type(x == -1 ? -1 : (x == 0 ? 0 : 1)) {};

        static void config()
        {
            m_description       = "allow fill into empty diagonal entry";
            m_default_value     = -1;
            m_validator         = validator_range<Integer>(-1,1);
        };
};

/// sparse matrix ordering algorithm
/// corresponding string option: "-pc_factor_mat_ordering_type"
class ordering : public option_base<Integer, ordering>
{
    private:
        using base_type         = option_base<Integer, ordering>;

    public:
        ordering()                  : base_type() {};
        ordering(opt_type x)        : base_type(x) {};
        ordering(petsc_ordering x)  : base_type((Integer)x) {};

        static void config()
        {
            m_description       = "sparse matrix ordering algorithm";
            m_default_value     = (Integer)petsc_ordering::default_val;
            m_validator         = validator_enum((Integer)petsc_ordering::last);
        };
};

/// reorder to remove zeros from diagonal; if this option is set
/// then column pivoing is applied to remove zeroes on main diagonal;
/// value of this option is use to test for (effective) zero values
/// corresponding string option: "-pc_factor_nonzeros_along_diagonal"
class reororder_zero_diagonal : public option_base<Real, reororder_zero_diagonal>
{
    private:
        using base_type         = option_base<Real, reororder_zero_diagonal>;

    public:
        reororder_zero_diagonal()           : base_type() {};
        reororder_zero_diagonal(opt_type x) : base_type(x) {};

        static void config()
        {
            m_description       = "reorder to remove zeros from diagonal";
            m_default_value     = -1.0;
            m_validator         = validator_nonnegative_or_val<Real>(-1.0);
        };
};

/// pivot is considered zero if less than tol given by this option
/// corresponding string option: "-pc_factor_zeropivot"
class zero_pivot : public option_base<Real, zero_pivot>
{
    private:
        using base_type         = option_base<Real, zero_pivot>;

    public:
        zero_pivot()           : base_type() {};
        zero_pivot(opt_type x) : base_type(x) {};

        static void config()
        {
            m_description       = "pivot is considered zero if less than tol";
            m_default_value     = -1.0;
            m_validator         = validator_nonnegative_or_val<Real>(-1.0);
        };
};

/// adds a quantity to the diagonal of the matrix during numerical
/// factorization, thus the matrix has nonzero pivots (zero pivot is
/// detected using option zero_pivot); exact interpretation of this
/// option depends on shift type
/// corresponding string option: "-pc_factor_shift_amount", "-pc_sor_diagonal_shift"
class shift_amount : public option_base<Real, shift_amount>
{
    private:
        using base_type         = option_base<Real, shift_amount>;

    public:
        shift_amount()           : base_type() {};
        shift_amount(opt_type x) : base_type(x) {};

        static void config()
        {
            m_description       = "adds a quantity to the diagonal of the matrix during numerical factorization";
            m_default_value     = -1.0;
            m_validator         = validator_nonnegative_or_val<Real>(-1.0);
        };
};

/// type of shift to add to diagonal; see petsc_shift_type for details
/// corresponding string option: "-pc_factor_shift_type"
class shift_type : public option_base<Integer, shift_type>
{
    private:
        using base_type         = option_base<Integer, shift_type>;

    public:
        shift_type()                      : base_type() {};
        shift_type(opt_type x)            : base_type(x) {};
        shift_type(petsc_shift_type x)    : base_type((Integer)x) {};

        static void config()
        {
            m_description       = "type of shift to add to diagonal";
            m_default_value     = (Integer)petsc_shift_type::default_val;
            m_validator         = validator_enum((Integer)petsc_shift_type::last);
        };
};

/// automatically finds rows with zero or negative diagonal and uses Schur complement
/// with no preconditioner as the solver; used in field_split_precond
/// corresponding string option: "-pc_fieldsplit_detect_saddle_point"
class detect_saddle : public option_base<Integer, detect_saddle>
{
    private:
        using base_type         = option_base<Integer, detect_saddle>;

    public:
        detect_saddle()             : base_type() {};
        detect_saddle(opt_type x)   : base_type(x) {};

        static void config()
        {
            m_description       = "automatically finds rows with zero or negative diagonal"
                                  " in field_split_precond";
            m_default_value     = -1;
            m_validator         = validator_nonnegative_or_val<Integer>(-1);
        };
};

}}};
