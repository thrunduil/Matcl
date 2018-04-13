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

namespace matcl
{

/// PETSC iterative solver
/// list of available solvers and more detailed descriptions
/// can be found at http://www.mcs.anl.gov/petsc/index.html
enum class petsc_solver
{
    default_val,/// default PETSC solver
    preonly,    /// this implements a stub method that applies ONLY the preconditioner

    /// general problems

    gmres,      /// Generalized Minimal Residual method, (Saad and Schultz, 1986), with restart
                /// Left and right preconditioning are supported, but not symmetric preconditioning
    fgmres,     /// flexible Generalized Minimal Residual method developed by Saad with restart;
                /// only right preconditioning is supported.
    lgmres,     /// the GMRES method with augmented approximation space with approximations to the
                /// error from previous restart cycles; supports both left and right preconditioning, 
                /// but not symmetric
    dgmres,     /// the deflated GMRES as defined
	            /// J. Erhel, K. Burrage and B. Pohl, Restarted GMRES preconditioned by deflation,
                /// J. Computational and Applied Mathematics, 69(1996).
                /// D. NUENTSA WAKAM and F. PACULL, Memory Efficient Hybrid Algebraic Solvers for Linear
                /// Systems Arising from Compressible Flows, Computers and Fluids, In Press, 
                /// http://dx.doi.org/10.1016/j.compfluid.2012.03.023
                /// the adaptive strategy allows to switch to the deflated GMRES when the stagnation occurs;
                /// left and right preconditioning are supported, but not symmetric preconditioning;
                /// complex arithmetic is not yet supported
    cgne,       /// preconditioned conjugate gradient method to the A' * A, without explicitly 
                /// forming A' * A; If you intend to solve least squares problems use LSQR;
                /// it merely uses PCG algorithm with the matrix defined by A' * A and preconditioner
                /// defined by B' * B where B is the preconditioner for A; A must be square;
                /// only left preconditioning is supported
    tfqmr,      /// a transpose free QMR (quasi minimal residual), supports left and right preconditioning,
                /// but not symmetric; the "residual norm" computed in this algorithm is actually just an 
                /// upper bound on the actual residual norm; that is for left preconditioning it is a bound
                /// on the preconditioned residual and for right preconditioning it is a bound on the true
                /// residual
    tcqmr,      /// a variant of QMR (quasi minimal residual) developed by Tony Chan
                /// Tony F. Chan, Lisette de Pillis, and Henk van der Vorst, Transpose free formulations of 
                /// Lanczos type methods for nonsymmetric linear systems, Numerical Algorithms, Volume 17, 1998. 
                /// supports either left or right preconditioning, but not symmetric; the "residual norm" 
                /// computed in this algorithm is actually just an upper bound on the actual residual norm;
                /// that is for left preconditioning it is a bound on the preconditioned residual and for right 
                /// preconditioning it is a bound on the true residual
    bicg,       /// implements the Biconjugate gradient method (similar to running the conjugate gradient on 
                /// the normal equations); supports only left preconditioning
    bcgs,       /// implements the BiCGStab (Stabilized version of BiConjugate Gradient Squared) method;
                /// supports left and right preconditioning but not symmetric
    bcgsl,      /// implements a slight variant of the Enhanced BiCGStab(L) algorithm:
                /// G.L.G. Sleijpen, H.A. van der Vorst, D.R. Fokkema, "BiCGStab(L) and other hybrid BiCG methods", 
                ///         Numerical Algorithms, 7, 1994.
                /// D.R. Fokkema, "Enhanced implementation of BiCGStab(L) for solving linear systems of equations", 
                ///         preprint from www.citeseer.com.in
                /// supports left and right preconditioning but not symmetric
    ibcgs,      /// the IBiCGStab (Improved Stabilized version of BiConjugate Gradient Squared) method in an 
                /// alternative form to have only a single global reduction operation instead of the usual 3 (or 4)
                /// supports left and right preconditioning; not supported for complex numbers.
    fbcgs,      /// the flexible BiCGStab method; only allow right preconditioning
    fbcgsr,     /// mathematically equivalent variant of FBiCGSTab; only allow right preconditioning
    cgs,        /// implements the CGS (Conjugate Gradient Squared) method; does not apply transpose of the matrix;
                /// supports left and right preconditioning, but not symmetric.
    gcr,        /// the preconditioned Generalized Conjugate Residual method; the GCR method permits the use 
                /// of a preconditioner which may vary from one iteration to the next; support only for right
                /// preconditioning.

    /// positive definite problems
    
    cg,         /// the preconditioned conjugate gradient (PCG) method; The PCG method requires 
                /// both the matrix and preconditioner to be symmetric positive (or negative) (semi)
                /// definite; only left preconditioning is supported    
    cr,         /// the (preconditioned) conjugate residuals method; the operator and the preconditioner 
                /// must be symmetric for this method. The preconditioner must be positive definite and
                /// the operator positive semi definite; support only for left preconditioning.
    fcg,        /// flexible Conjugate Gradient method; support left preconditioning only
    lcd,        /// the LCD (left conjugate direction) method; support only for left preconditioning

    /// symmetric problems
    
    minres,     /// the MINRES (Minimum Residual) method; The operator and the preconditioner must be 
                /// symmetric and the preconditioner must be positive definite for this method. 
                /// Supports only left preconditioning.
    symmlq,     /// the SYMMLQ method; the operator and the preconditioner must be symmetric; the 
                /// preconditioner must be positive definite; supports only left preconditioning.

    /// least squares problems

    lsqr,       /// the LSQR method - a conjugate gradient type method for solving linear equations
                /// and least-squares problems: AX = b, or min |AX-B|_2 or 
                /// min sqrt(|AX-b|_2) + sqrt(l|X|_2), where A is square or rectangular (overdetermined
                /// or underdetermined) and may have any rank, l is a damping parameter; the 
                /// preconditioner needs to work for the normal equations A'*A; supports only left 
                /// preconditioning

    /// constrained positive definite problems

    nash,       /// conjugate gradient method subject to a constraint on the solution norm: min |A*X-B|_2
                /// s.t. |X|_2 <= delta; the preconditioner supplied should be symmetric and positive 
                /// definite; and A should be symmetric semi-positive definite
    stcg,       /// see nash
    gltr,       /// see nash
    qcg,        /// conjugate gradient method subject to a constraint on the solution norm: min |A*X-B|_2
                /// s.t. |D*X|_2 <= delta, where D is scaling matrix, the preconditioner supplied should 
                /// be symmetric and positive definite; and A should be symmetric semi-positive definite
                /// only symmetric preconditioning is allowe (jacobi, icc)

    /// smoothers

    richardson, /// the preconditioned Richardson iterative method with iteration for general problems
                /// x^{n+1} = x^{n} + scale*B(b - A x^{n}), where B is the application of preconditioner;
                /// this method usually will not converge unless scale is very small;
                /// For some preconditioners, currently SOR, the convergence test is skipped to improve
                /// speed, thus it always iterates the maximum number of iterations you've selected. 
                /// When any monitor is turned on, the norm is computed at each iteration and so the 
                /// convergence test is run; supports only left preconditioning
    chebyshev,  /// the preconditioned Chebyshev iterative method; both the matrix and preconditioner
                /// must be to be symmetric positive (semi) definite; this method required good estimates
                /// of the smallest and largest eigenvalues of the preconditioned operator;
                /// only support for left preconditioning

    last
};

/// norm that is passed in the Krylov convergence test routines
enum class petsc_conv_norm
{
    default_val,    /// default PETSC norm, i.e. preconditioned
    none,           /// do not compute a norm during the Krylov process, iterations will be finished
                    /// after maximum number of iterations
    preconditioned, /// compute the norm of the preconditioned residual B*(b - A*x), if left preconditioning
    true_resid,     /// compute the norm of the true residual (b - A*x), may be costly to compute
    normal,         /// compute the norm of residual sqrt((b - A*x)*B*(b - A*x))

    last
};

/// PETSC preconditioner
enum class petsc_precond
{
    default_val,

    none,       /// no preconditioning

    /// factorizations
    ilu,        /// incomplete factorization preconditioners
    icc,        /// incomplete Cholesky factorization preconditioners    
    lu,         /// uses a direct solver, based on LU factorization, as a preconditioner    
    svd,        /// use pseudo inverse defined by SVD of operator
    cholesky,   /// uses a direct solver, based on Cholesky factorization, as a preconditioner

    /// Gauss–Seidel type methods
    sor,        /// (S)SOR (successive over relaxation, Gauss-Seidel) preconditioning
    eisenstat,  /// An implementation of SSOR (symmetric successive over relaxation, symmetric
                /// Gauss-Seidel) preconditioning that incorporates Eisenstat's trick to reduce
                /// the amount of computation needed.    
    jacobi,     /// Jacobi (i.e. diagonal scaling preconditioning); susing symmetric side preconditioning
                /// can scale each side of the matrix by the square root of the diagonal entries;
                /// zero entries along the diagonal are replaced with the value 1.0
    pbjacobi,   /// Point block Jacobi preconditioner
    kaczmarz,   /// Kaczmarz iteration; if the linear system is consistent, then Kaczmarz iteration
                /// converges to the minimum-norm solution, provided that the iterations start with
                /// the zero vector

    last
};

/// iterative refinement procedure in Gramm-Schmidt orthogonalization
enum class orth_refine_type
{
    default_val,        /// default PETSC procedure
    never,              /// never refine Gramm-Schmidt orthogonalization
    if_needed,          /// do the classical Gram-Schmidt process and one step of iterative refinement if
                        /// an estimate of the orthogonality of the resulting vectors indicates poor orthogonality
    always,             /// do two steps of the classical Gram-Schmidt process

    last
};

/// type of orthogonalization of a new search direction against all
/// previous directions in fcg solver
enum class fcg_orth_type
{
    default_val,        /// default PETSC procedure
    full,               /// orthogonalize against all stored directions
    partial,            /// orthogonalize against k last directions, where k changes from 1 to m
                        /// cyclicaly (m is number of stored search directions)

    last,
};

/// ordering algorithm
enum class petsc_ordering
{
    default_val,        /// default PETSC algorithm
    natural,            /// no ordering
    nd,                 /// nested dissection (symmetric ordering)
    w1d,                /// one-way dissection (symmetric ordering)
    rcm,                /// reverse Cuthill-McKee (symmetric ordering)
    qmd,                /// quotient minimum degree (symmetric ordering)
    rowlength,          /// order by number of nonzeros (symmetric ordering)
    nonzero_diagonal,
    last
};

/// type of shift to add to diagonal during numerical phase of factorization
enum class petsc_shift_type
{
    default_val,        /// default PETSC algorithm
    none,               /// never add shift to zero diagonal
    nonzero,            /// force |diag| > zeropivot*rs, where zeropivot is
                        /// the tolerance used to detect small pivot (option zero_pivot)
                        /// and rs is active row sum of abs(offdiagonals)
    positive_definite,  /// force factored matrix to be positive definite and diagonally
                        /// dominant
    inblocks,           /// if |diag| <= zeropivot, then change diag to diag + shift,
                        /// (this may cause that diag is even smaller)

    last
};

/// type of relaxation for successive over-relaxation algorithms
enum class petsc_sor_type
{
    default_val,        /// default PETSC algorithm
    forward,            /// forward version
    backward,           /// backward version
    symmetric,          /// symmetric version

    last
};

/// side of preconditioner
enum class petsc_precond_side
{
    default_val,        /// default PETSC preconditioner side
    left,               /// left preconditioning
    right,              /// right preconditioning
    symmetric,          /// symmetric preconditioning

    last
};

/// type of jacobi preconditioner
enum class petsc_jacobi_type
{
    default_val,        /// default PETSC jacobi version
    diag,               /// use values of entries on main diagonal
    max_abs_row,        /// use maximum absolute values in each row
    sumrows,            /// use sum of entries in each row

    last
};

/// type of local solver composition for additive Schwarz method preconditioner
enum class asm_composition
{
    additive,           /// local additive combination
    multiplicative      /// local multiplicative combination; more accurate but slower
};

/// type of restriction/extension for additive Schwarz method preconditioner
enum class asm_restrict
{
    basic,              /// symmetric version where residuals from the ghost points are used
                        /// and computed values in ghost regions are added together;
                        /// classical standard additive Schwarz.
    restricted,         /// residuals from ghost points are used but computed values in ghost
                        /// region are discarded;
    interpolate,        /// residuals from ghost points are not used, computed values in ghost
                        /// region are added back in.
    none                /// Residuals from ghost points are not used, computed ghost values are
                        /// discarded. Not very good.
};

/// type of composition of preconditioners
enum class petsc_composite_type
{
    mult,               /// preconditioners are applied sequentially to the residual
                        /// freshly computed after the previous preconditioner application
    sum,                /// results from application of all preconditioners are added together
    sym_mult,           /// preconditioners are applied sequentially to the residual freshly
                        /// computed from first preconditioner to last and then back (use only
                        /// for symmetric matrices and preconditions)
};

/// type of preconditioner build for a split A = [A00, A012; A10 A11;] (or more blocks)
enum class field_split_type
{
    jacobi,             /// block Jacobi: preconditioner is constructed from diagonal blocks: 
                        /// P = [inv(A00), 0; 0, inv(A11)]
    gauss_seidel,       /// forward block Gauss-Seidel constructed from diagonal and subdiagonals
                        /// [I 0; 0 inv(A11)] x [I 0; -A10 I] x [inv(A00) 0; 0 I]
    sym_gauss_seidel,   /// symmetric block Gauss-Seidel: preconditioner is applied as:
                        /// [inv(A00) 0; 0 I] x [I -A10; 0 I] x [A00 0; 0 inv(A11)] x
                        /// [I 0, -A10 I] x [inv(A00) 0; 0 I]
};

/// multigrid cycle type
enum class mg_cycle_type
{
    v,      /// V-cycle
    w,      /// W-cycle
};

/// multigrid type
enum class mg_type
{
    multiplicative, /// traditional V or W cycle
    additive,       /// the additive multigrid preconditioner where all levels are smoothed
                    /// before updating the residual. This only uses the down smoother, in the
                    /// preconditioner the upper smoother is ignored
    full,           /// same as multiplicative except one also performs grid sequencing, that is 
                    /// starts on the coarsest grid, performs a cycle, interpolates to the next,
                    /// performs a cycle etc.
    kaskade         /// like full multigrid except one never goes back to a coarser level from a finer
};

}
