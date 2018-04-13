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

#include "matcl-linalg/general/config_linalg.h"
#include "matcl-core/options/matcl_options.h"
#include "matcl-core/lib_functions/constants.h"

namespace matcl { namespace opt { namespace linsolve
{

/// pivoting strategy
enum class pivot_type
{
    partial = 0,    /// maximum element in given row
    rook,           /// maximum element in row and column
    complete,       /// maximum element in a matrix

    last
};

/// ordering method for Cholesky factorization
enum class chol_ordering_type
{
    default_val,    /// default method used by given solver, usually this is the
                    /// best ordering
    natural,        /// no reordering

    amd,            /// use minimum degree (AMD)
    cholmod_nested, /// use cholmod's nested dissection

    last
};

/// ordering method for lu (or ilu) factorization
/// other ordering methods must be used explicitly by permuting the matrix and 
/// selecting natural ordering; see matcl_graph.h for available reordering
/// methods, usually colamd gives the best ordering
enum class lu_ordering_type
{
    natural,        /// no reordering
    colamd,         /// use order_colamnd

    last
};

/// incomplete lu factorization type
enum class ilu_method
{
    silu,           /// standard pivoted ilu
    milu,           /// modified ilu, diagonal entries of U factor are modified as
                    /// U(i,i) := U(i,i) + sign(U(i,i)) * sum(dropped entries)
    milu_abs,       /// modified ilu, diagonal entries of U factor are modified as
                    /// U(i,i) := U(i,i) + sign(U(i,i)) * sum(|dropped entries|)

    last
};

/// LU solver used for LU factorization with partial pivoting
enum class lu_solver_type
{
    lusol,          /// use lusol, notice that lusol uses Markowitz pivoting strategy
                    /// and ignores user provider ordering; lusol usually generates 
                    /// sparser factors than superlu but is slower
    superlu,        /// use superlu (default)
    last
};

//---------------------------------------------------------
//              options for lu decomposition
//---------------------------------------------------------

/// pivoting strategy
class pivot : public option_base<Integer, pivot>
{
    private:
        using base_type         = option_base<Integer, pivot>;
        using opt_type          = optional<Integer>;

    public:
        pivot()                 : base_type() {};
        pivot(opt_type x)       : base_type(x) {};
        pivot(pivot_type x)     : base_type((Integer)x) {};

        static void config()
        {
            m_description       = "pivoting strategy";
            m_default_value     = (Integer)pivot_type::partial;
            m_validator         = validator_enum((Integer)pivot_type::last);
        };
};

/// pivoting strategy
class lu_solver : public option_base<Integer, lu_solver>
{
    private:
        using base_type         = option_base<Integer, lu_solver>;
        using opt_type          = optional<Integer>;

    public:
        lu_solver()                 : base_type() {};
        lu_solver(opt_type x)       : base_type(x) {};
        lu_solver(lu_solver_type x) : base_type((Integer)x) {};

        static void config()
        {
            m_description       = "lu sparse solver";
            m_default_value     = (Integer)lu_solver_type::superlu;
            m_validator         = validator_enum((Integer)lu_solver_type::last);
        };
};

/// tolerance used to test singularity of diagonal elements of the 
/// U factor; element v is considered as zero if |v| < tol
/// if tol is negative, then default tolerance is used,
/// tol is set to 10.0 * eps * |A|_F / sqrt(min(M,N)), where
/// |A|_F is the Frobenius norm of a MxN matrix
class tol : public option_base<Real, tol>
{
    private:
        using base_type         = option_base<Real, tol>;
        using opt_type          = optional<Real>;

    public:
        tol()                   : base_type() {};
        tol(opt_type x)         : base_type(x) {};

        static void config()
        {
            m_description       = "tolerance used to test singularity";
            m_default_value     = Real(-1.0);
        };
};

/// column pivot tolerance; given column is selected as a pivot if element 
/// is greater than 1/tol_c * v, where v is maximul element in given row, 
class tol_c : public option_base<Real, tol_c>
{
    private:
        using base_type         = option_base<Real, tol_c>;
        using opt_type          = optional<Real>;

    public:
        tol_c()                 : base_type() {};
        tol_c(opt_type x)       : base_type(x) {};

        static void config()
        {
            m_description       = "column pivot tolerance";
            m_default_value     = Real(0.1);
            m_validator         = validator_range<Real>(0.0, 1.0);
        };
};

/// row pivot tolerance; given row is selected as a pivot if element 
/// is greater than 1/tol_v * v, where v is maximul element in given column, 
class tol_r : public option_base<Real, tol_r>
{
    private:
        using base_type         = option_base<Real, tol_r>;
        using opt_type          = optional<Real>;

    public:
        tol_r()                 : base_type() {};
        tol_r(opt_type x)       : base_type(x) {};

        static void config()
        {
            m_description       = "row pivot tolerance";
            m_default_value     = Real(0.1);
            m_validator         = validator_range<Real>(0.0, 1.0);
        };
};

/// if this option is true, then the L factor can be
/// any lower triangular matrix; then LU factorization
/// of a lower triangular matrix A gives L = A, U = I
class allow_nonunit_L : public option_base<bool, allow_nonunit_L>
{
    private:
        using base_type        = option_base<bool, allow_nonunit_L>;
        using opt_type         = optional<bool>;

    public:
        allow_nonunit_L()              : base_type() {};
        allow_nonunit_L(opt_type x)    : base_type(x) {};

        static void config()
        {
            m_description       = "if true, then the L factor can be any lower tringular matrix";
            m_default_value     = true;
        };
};

//---------------------------------------------------------
//          options controling factorization used
//          when creating linsolve object
//---------------------------------------------------------
/// if this option is true, then rank revealing factorizations will be
/// used, i.e. rook pivoted lu or permuted Cholesky factorization
class use_rr : public option_base<bool, use_rr>
{
    private:
        using base_type        = option_base<bool, use_rr>;
        using opt_type         = optional<bool>;

    public:
        use_rr()            : base_type() {};
        use_rr(opt_type x)  : base_type(x) {};

        static void config()
        {
            m_description       = "if true, then rank revealing factorizations will be used";
            m_default_value     = false;
        };
};

/// if this option is true, then LDL factorization will be used
/// for dense symmetric/hermitian indefinite matrices
class can_use_ldl : public option_base<bool, can_use_ldl>
{
    private:
        using base_type        = option_base<bool, can_use_ldl>;
        using opt_type         = optional<bool>;

    public:
        can_use_ldl()          : base_type() {};
        can_use_ldl(opt_type x): base_type(x) {};

        static void config()
        {
            m_description       = "if true, then use LDL factorization for symmetric indefinite problems";
            m_default_value     = true;
        };
};

/// if this option is true, then matrix is balanced first before
/// factorization is performed; this option is not used by rank
/// revealing factorizations
class do_balancing : public option_base<bool, do_balancing>
{
    private:
        using base_type        = option_base<bool, do_balancing>;
        using opt_type         = optional<bool>;

    public:
        do_balancing()          : base_type() {};
        do_balancing(opt_type x): base_type(x) {};

        static void config()
        {
            m_description       = "if true, then balancing is performed when factorization is required";
            m_default_value     = true;
        };
};

/// if this option is true, then matrix is balanced first before
/// factorization is performed; balancing functions use small value
/// tolerances based on tol_sing option; this option is used only 
/// by rank revealing factorizations (i.e. lu with rook or complete
/// pivoting, pivoted cholesky, pivoted qr)
class do_balancing_rr : public option_base<bool, do_balancing_rr>
{
    private:
        using base_type        = option_base<bool, do_balancing_rr>;
        using opt_type         = optional<bool>;

    public:
        do_balancing_rr()          : base_type() {};
        do_balancing_rr(opt_type x): base_type(x) {};

        static void config()
        {
            m_description       = "if true, then balancing is performed when rank revealing"
                                  " factorization is required";
            m_default_value     = false;
        };
};

/// control small diagonals corrections added to triangular factors T
/// created when constructing linsolve object; if |d_ii| < tol, then
/// i-th diagonal element is modified to sign(d_ii) * tol; if 
/// min_diag_pivot = 0 then no correction is applied; if min_diag_pivot < 0
/// then tol = -min_diag_pivot; if min_diag_pivot > 0, then 
/// tol = eps^min_diag_pivot * max(|T|), where eps is epsilon for given 
/// value type
class tol_sing : public option_base<Real, tol_sing>
{
    private:
        using base_type         = option_base<Real, tol_sing>;
        using opt_type          = optional<Real>;

    public:
        tol_sing()              : base_type() {};
        tol_sing(opt_type x)    : base_type(x) {};

        static void config()
        {
            m_description       = "tolerance used to correct small diagonal elements in triangular factors";
            m_default_value     = 0.0;
        };
};

//---------------------------------------------------------------------
//      options controling solution phase in linsolve objects
//---------------------------------------------------------------------

/// if true then iterative refinement is performed; iterative refinement
/// is performed inly for direct solvers when no diagonal corrections
/// was added
class use_ir : public option_base<bool, use_ir>
{
    private:
        using base_type         = option_base<bool, use_ir>;
        using opt_type          = optional<bool>;

    public:
        use_ir()                : base_type() {};
        use_ir(opt_type x)      : base_type(x) {};

        static void config()
        {
            m_description       = "use iterative refinement?";
            m_default_value     = true;
        };
};

/// check componentwise forward error bound of computed solution, i.e.
/// |Y - Y*|_inf / |Y|_inf <= ferr, where Y is computed solution and Y*
/// is true solution; this options set required number of correct digits
/// of large entries of the solution; if 0 then quality of the solution 
/// is not tested; value +-k requires k correct digits; if estimated number
/// of correct digits is less than k, then warning is issued if sign is positive
/// and error is thrown is sign is negative
class test_sol : public option_base<Integer, test_sol>
{
    private:
        using base_type         = option_base<Integer, test_sol>;
        using opt_type          = optional<Integer>;

    public:
        test_sol()             : base_type() {};
        test_sol(opt_type x)   : base_type(x) {};

        static void config()
        {
            m_description       = "check correctness of solution?";
            m_default_value     = 2;
        };
};

//---------------------------------------------------------
//          iterative refinement
//---------------------------------------------------------
/// iterative refinement is performed by explicitly calling 
/// refinement functions from linsolve_obj

/// maximum number of iterative refinement steps, values less 
/// than 1 implies no iterative refinement; in case of 
/// highly ill-conditioned problem value 100 is recommended
class refinement_iter : public option_base<Integer, refinement_iter>
{
    private:
        using base_type         = option_base<Integer, refinement_iter>;
        using opt_type          = optional<Integer>;

    public:
        refinement_iter()           : base_type() {};
        refinement_iter(opt_type x) : base_type(x) {};

        static void config()
        {
            m_description       = "maximum number of iterative refinement steps";
            m_default_value     = 10;
        };
};

/// minimum error decrease allowed, iterative refinement is stopped if
/// |dx_i|_inf / |dx_{i-1}|_inf >= refinement_rho, where |dx_i| is the 
/// solution correction in step i; for highly ill-conditioned problems
/// value 0.9 is recommended
class refinement_rho : public option_base<Real, refinement_rho>
{
    private:
        using base_type         = option_base<Real, refinement_rho>;
        using opt_type          = optional<Real>;

    public:
        refinement_rho()            : base_type() {};
        refinement_rho(opt_type x)  : base_type(x) {};

        static void config()
        {
            m_description       = "minimum error decrease allowed";
            m_default_value     = Real(0.5);
            m_validator         = validator_range<Real>(0.0, 1.0, true);
        };
};

//---------------------------------------------------------
//          balancing
//---------------------------------------------------------
/// balancing is only used when creating linsolve objects; explicit call
/// to factorization routines never do balancing

/// maximum number of balancing iterations, values less 
/// than 1 implies no balancing; 
class balancing_iter : public option_base<Integer, balancing_iter>
{
    private:
        using base_type         = option_base<Integer, balancing_iter>;
        using opt_type          = optional<Integer>;

    public:
        balancing_iter()            : base_type() {};
        balancing_iter(opt_type x)  : base_type(x) {};

        static void config()
        {
            m_description       = "maximum number of balancing iterations";
            m_default_value     = 20;
        };
};

/// tolenance used to check convergense of diagonal scale factors;
/// negative value implies default tolerance used by algorithm, that use
/// this option
class balancing_tol : public option_base<Real, balancing_tol>
{
    private:
        using base_type         = option_base<Real, balancing_tol>;
        using opt_type          = optional<Real>;

    public:
        balancing_tol()             : base_type() {};
        balancing_tol(opt_type x)   : base_type(x) {};

        static void config()
        {
            m_description       = "tolenance used to check convergense of diagonal scale factors";
            m_default_value     = -1.0;
        };
};

//---------------------------------------------------------
//          reordering
//---------------------------------------------------------
/// ordering method used for Cholesky factorization, see
/// chol_ordering_type for details
class chol_ordering : public option_base<Integer, chol_ordering>
{
    private:
        using base_type         = option_base<Integer, chol_ordering>;
        using opt_type          = optional<Integer>;

    public:
        chol_ordering()                 : base_type() {};
        chol_ordering(opt_type x)       : base_type(x) {};
        chol_ordering(chol_ordering_type x)  : base_type((Integer)x) {};

        static void config()
        {
            m_description       = "reordering method for Cholesky factorization";
            m_default_value     = (Integer)chol_ordering_type::default_val;
            m_validator         = validator_enum((Integer)chol_ordering_type::last);
        };
};

/// ordering method used for LU factorization, see
/// lu_ordering_type for details
class lu_ordering : public option_base<Integer, lu_ordering>
{
    private:
        using base_type         = option_base<Integer, lu_ordering>;
        using opt_type          = optional<Integer>;

    public:
        lu_ordering()                 : base_type() {};
        lu_ordering(opt_type x)       : base_type(x) {};
        lu_ordering(lu_ordering_type x)  : base_type((Integer)x) {};

        static void config()
        {
            m_description       = "reordering method for LU factorization";
            m_validator         = validator_enum((Integer)lu_ordering_type::last);
        };
};

//---------------------------------------------------------
//          incomplete factorizations
//---------------------------------------------------------
/// elements lower that drop_tol are treated as zero; this
/// value is used directly for L factor, for U factor tolerance
/// |A(:,i)|_oo * drop_tol is used
class drop_tol : public option_base<Real, drop_tol>
{
    private:
        using base_type         = option_base<Real, drop_tol>;
        using opt_type          = optional<Real>;

    public:
        drop_tol()              : base_type() {};
        drop_tol(opt_type x)    : base_type(x) {};

        static void config()
        {
            m_description       = "threshold for dropping";
            m_default_value     = 1e-4;
            m_validator         = validator_range<Real>(0.0, 1.0);
        };
};

/// if some U(i,i) = 0, so that U is exactly singular, then a small number is 
/// added to the diagonal entry, that is U(i,i) = |A(:,i)|_oo * diag_fill_tol ^ pow
/// where pow is a factor depending on density of given column
class diag_fill_tol : public option_base<Real, diag_fill_tol>
{
    private:
        using base_type         = option_base<Real, diag_fill_tol>;
        using opt_type          = optional<Real>;

    public:
        diag_fill_tol()              : base_type() {};
        diag_fill_tol(opt_type x)    : base_type(x) {};

        static void config()
        {
            m_description       = "threshold for zero pivot perturbation";
            m_default_value     = 1e-2;
            m_validator         = validator_range<Real>(0.0, 1.0);
        };
};

/// expected upper bound on fill ratio, i.e. nnz(L+U-I)/nnz(A); it is not 
/// guaranteed that actual fill ration will be not greater that fill_factor
class fill_factor : public option_base<Real, fill_factor>
{
    private:
        using base_type         = option_base<Real, fill_factor>;
        using opt_type          = optional<Real>;

    public:
        fill_factor()              : base_type() {};
        fill_factor(opt_type x)    : base_type(x) {};

        static void config()
        {
            m_description       = "upper bound on fill ratio, i.e. nnz(L+U-I)/nnz(A)";
            m_default_value     = 10.0;
            m_validator         = validator_range<Real>(1.0, constants::inf<Real>());
        };
};

class ilu_type : public option_base<Integer, ilu_type>
{
    private:
        using base_type         = option_base<Integer, ilu_type>;
        using opt_type          = optional<Integer>;

    public:
        ilu_type()              : base_type() {};
        ilu_type(opt_type x)    : base_type(x) {};
        ilu_type(ilu_method x)  : base_type((Integer)x) {};

        static void config()
        {
            m_description       = "version of incomplete lu to use";
            m_default_value     = (Integer)ilu_method::silu;
            m_validator         = validator_enum((Integer)ilu_method::last);
        };
};

}}};
