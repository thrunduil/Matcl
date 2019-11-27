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

#include "matcl-linalg/petsc_utils/petsc_option.h"
#include "matcl-linalg/petsc_utils/petsc_objects.h"
#include "matcl-linalg/petsc_utils/petsc_ext.h"

#include "matcl-linalg/general/mmlib_petsc_exception.h"
#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/linear_eq/preconditioner.h"

#pragma warning (push)
#pragma warning (disable:4127)
#pragma warning (disable:4101)
#pragma warning (disable:4100)

#include "petscmat.h"
#include "petsc/private/kspimpl.h"
#include "petsc/private/matimpl.h"
#include "petscksp.h" 

#pragma warning (pop)

#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/options/options_petsc.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-linalg/special_matrices/struct_flag_linalg.h"

//---------------------------------------------------------------
//                  typedefs from PETSC
//---------------------------------------------------------------

typedef struct {
  PetscInt  nwork_n,nwork_m;
  Vec       *vwork_m;   /* work vectors of length m, where the system is size m x n */
  Vec       *vwork_n;   /* work vectors of length n */
  Vec       se;         /* Optional standard error vector */
  PetscBool se_flg;     /* flag for -ksp_lsqr_set_standard_error */
  PetscReal arnorm;     /* Norm of the vector A.r */
  PetscReal anorm;      /* Frobenius norm of the matrix A */
  PetscReal rhs_norm;   /* Norm of the right hand side */
} KSP_LSQR;

namespace matcl { namespace details
{

static PetscErrorCode KSP_monitor_definition(KSP ksp, PetscInt its, PetscReal rnorm, void* notifier_data)
{
    using notifier_ptr      = std::shared_ptr<ksp_solver_notifier>;
    notifier_ptr* not_void  = (notifier_ptr*)notifier_data;

    notifier_ptr notifier   = *not_void;

    if (notifier == nullptr)
        return 0;

    bool is_ls              = notifier->is_least_squares();    

    if (is_ls)
    {
        ::KSP_LSQR* lsqr    = (::KSP_LSQR*)ksp->data;

        notifier->report_iteration_ls(its, rnorm, lsqr->arnorm, ksp);
    }
    else
    {
        notifier->report_iteration(its, rnorm, ksp);
    };

    return 0;
}

mg_data::mg_data()
    :m_levels(0), m_current_level(0), m_solver_set(false), m_smoother_set(false)
{};

void mg_data::set_levels(Integer nlev)
{
    m_levels        = nlev;
    m_current_level = nlev-1;

    m_smooth_post.resize(m_levels);
    m_smooth_pre.resize(m_levels);
};

void mg_data::coarser_level()
{
    --m_current_level;

    if (m_smoother_set == false)
        throw error::multigrid_smoother_not_set();

    if (m_current_level < 0)
        throw error::invalid_mg_level(m_current_level);

    m_smoother_set = false;
}
void mg_data::set_solver()
{
    m_solver_set = true;
};
void mg_data::set_smoother_post(const ksp_solver& s)
{
    m_smoother_set = true;
    m_smooth_post[m_current_level]  = ksp_bool(s,true);
};
void mg_data::set_smoother_pre(const ksp_solver& s)
{
    m_smoother_set = true;
    m_smooth_pre[m_current_level]  = ksp_bool(s, true);
};
mg_data::ksp_vec&& mg_data::get_post()
{
    return std::move(m_smooth_post);
}
mg_data::ksp_vec&& mg_data::get_pre()
{
    return std::move(m_smooth_pre);
}

//---------------------------------------------------------------
//                  ksp_solver_impl
//---------------------------------------------------------------
ksp_solver_impl::ksp_solver_impl(const Matrix& A_in, const options& opts)
    : m_ksp_lock(new matcl::petsc::petsc_lock_helper())
    , m_mat_A(preprocess_mat(A_in))    
    , m_petsc_mat_A(details::petsc_matrix::create(m_mat_A))
{
    m_mat_P         = m_mat_A;
    m_petsc_mat_P   = m_petsc_mat_A;
    
    initialize();
    reset_options(opts);
}

ksp_solver_impl::ksp_solver_impl(const Matrix& A_in, const Matrix& B_in, const options& opts)
    : m_ksp_lock(new matcl::petsc::petsc_lock_helper())
    , m_mat_A(preprocess_mat(A_in)), m_mat_P(preprocess_mat(B_in))
    , m_petsc_mat_A(details::petsc_matrix::create(m_mat_A))
    , m_petsc_mat_P(details::petsc_matrix::create(m_mat_P))
{
    check_preconditioner(B_in.rows(), B_in.cols(), A_in.rows(), A_in.cols());
    initialize();
    reset_options(opts);
}

ksp_solver_impl::ksp_solver_impl(const linear_operator& A_in, const options& opts)
    : m_ksp_lock(new matcl::petsc::petsc_lock_helper())
    , m_petsc_mat_A(details::petsc_matrix::create(A_in))
{
    m_petsc_mat_P = m_petsc_mat_A;

    initialize();
    reset_options(opts);
}

ksp_solver_impl::ksp_solver_impl(const linear_operator& A_in, const Matrix& B_in, const options& opts)
    : m_ksp_lock(new matcl::petsc::petsc_lock_helper()), m_mat_P(preprocess_mat(B_in))
    , m_petsc_mat_A(details::petsc_matrix::create(A_in))
    , m_petsc_mat_P(details::petsc_matrix::create(m_mat_P))
{
    check_preconditioner(B_in.rows(), B_in.cols(), A_in.rows(), A_in.cols());
    initialize();
    reset_options(opts);
}

void ksp_solver_impl::set_shell_ksp(KSP* ksp_ref, const linear_operator& solv, KSP& old)
{        
    ksp_solver ksp_solv(solv);
    old         = *ksp_ref;
    KSP new_ksp = ksp_solv.get_KSP();

    m_ksp_vec.push_back(ksp_solv);

    PetscObjectReference((PetscObject)new_ksp);

    *ksp_ref = new_ksp;
}
void ksp_solver_impl::set_shell_ksp(KSP* ksp_ref, const linsolve_obj& solv, KSP& old)
{        
    ksp_solver ksp_solv(solv);
    old         = *ksp_ref;
    KSP new_ksp = ksp_solv.get_KSP();

    m_ksp_vec.push_back(ksp_solv);

    PetscObjectReference((PetscObject)new_ksp);

    *ksp_ref = new_ksp;
}

void ksp_solver_impl::set_shell_ksp(KSP* ksp_ref, const ksp_solver& solv, KSP& old)
{
    old         = *ksp_ref;
    KSP new_ksp = solv.get_KSP();

    m_ksp_vec.push_back(solv);

    PetscObjectReference((PetscObject)new_ksp);

    *ksp_ref = new_ksp;
}

void ksp_solver_impl::check_preconditioner(Integer rows, Integer cols, Integer mat_rows, Integer mat_cols) const
{
    (void)mat_rows;

    if (rows != cols)
        throw error::square_matrix_required(rows, cols);

    if (cols != mat_cols)
        throw error::invalid_size2(rows, cols, this->cols(), this->cols());
};

void ksp_solver_impl::set_id(const std::string& id)
{
    m_id = id;
}
const std::string& ksp_solver_impl::get_id() const
{
    return m_id;
};

void ksp_solver_impl::reset_operator(const Matrix& A_in)
{
    if (A_in.rows() != this->rows() || A_in.cols() != this->cols())
        throw error::invalid_size2(A_in.rows(), A_in.cols(), this->rows(), this->cols());

    // just for precausion increase refcounts;
    Matrix A                = preprocess_mat(A_in);
    Matrix A_old            = m_mat_A;
    m_mat_A                 = A;

    petsc_matrix_ptr pm_A   = details::petsc_matrix::create(A);

    reinitialize(pm_A);
    reset_options();
};

void ksp_solver_impl::reset_operator(const linear_operator& A_in)
{
    if (A_in.rows() != this->rows() || A_in.cols() != this->cols())
        throw error::invalid_size2(A_in.rows(), A_in.cols(), this->rows(), this->cols());

    petsc_matrix_ptr pm_A   = details::petsc_matrix::create(A_in);
    
    reinitialize(pm_A);
    reset_options();
};

void ksp_solver_impl::rebuild_preconditioner(const Matrix& B_in, const options& opts)
{
    Matrix B                = preprocess_mat(B_in);
    petsc_matrix_ptr pm_P   = details::petsc_matrix::create(B);

    check_preconditioner(B_in.rows(), B_in.cols(), m_mat_A.rows(), m_mat_A.cols());    

    reuse_preconditioner(false);
    
    KSPSetOperators(*m_ksp, m_petsc_mat_A->get_Aksp(), pm_P->get_Aksp());
    reuse_preconditioner(true);

    m_petsc_mat_P       = pm_P;
    m_mat_P             = B;

    m_pc_shell_vec.clear();
    m_block_data.clear();
    m_setup_called      = false;    
    m_matrix_vec.clear();
    m_ksp_vec.clear();
    m_mg_data           = mg_data();

    reset_options(opts);
};

void ksp_solver_impl::restore_null_spaces(const ksp_solver_impl& base_ksp)
{
    m_has_right_null    = base_ksp.m_has_right_null;
    m_has_left_null     = base_ksp.m_has_left_null;
};

Matrix ksp_solver_impl::preprocess_mat(const Matrix& A)
{
    Matrix B = A;
    B.diag(0).add_sparse();
    return B;
}

void ksp_solver_impl::initialize()
{
    m_ksp       = ksp_solver_solver(new mp::smart_KSP(mp::smart_KSP::create(PETSC_COMM_SELF)));
    m_notifier  = notifier(new default_notifier(this));
    
    KSPSetOperators(*m_ksp, m_petsc_mat_A->get_Aksp(), m_petsc_mat_P->get_Aksp());

    m_has_right_null    = false;
    m_has_left_null     = false;
    m_ls_problem        = false;
    m_setup_called      = false;
}

void ksp_solver_impl::reinitialize(const petsc_matrix_ptr& pm_A)
{
    reuse_preconditioner(true);

    KSPSetOperators(*m_ksp, pm_A->get_Aksp(), m_petsc_mat_P->get_Aksp());
    m_petsc_mat_A   = pm_A;

    m_has_right_null    = false;
    m_has_left_null     = false;
    m_setup_called      = false;

    remove_nullspace_right();
    remove_nullspace_left();
};

void ksp_solver_impl::set_tolerances(Integer maxit, Real rtol, Real atol, Real dtol)
{
    if (maxit != -2)
        m_set_options.set(opt::petsc::max_it(maxit));
    if (rtol != -2.0)
        m_set_options.set(opt::petsc::rtol(rtol));
    if (atol != -2.0)
        m_set_options.set(opt::petsc::atol(atol));
    if (dtol != -2.0)
        m_set_options.set(opt::petsc::atol(dtol));

    petsc::option_stack st = m_options_handler->install();
    ::KSPSetTolerances(get_KSP(),rtol,atol,dtol, maxit);
};

void ksp_solver_impl::reset_options(const options& opts0)
{
    options opts            = correct_options(opts0);
    bool nonsquare_allowed  = check_nonsquare_allowed(opts);

    if (nonsquare_allowed == false && (m_petsc_mat_A->rows() != m_petsc_mat_A->cols()))
        throw error::square_matrix_required(m_petsc_mat_A->rows(), m_petsc_mat_A->cols());

    if (m_notifier)
        m_notifier->set_least_squares(nonsquare_allowed);

    opts.make_unique();
    m_set_options           = opts;
    m_missing_opts          = options();  
    m_options_handler       = opt_handler::make(&m_set_options, &m_missing_opts);      

    petsc::option_stack st  = m_options_handler->install();

    KSPMonitorCancel(*m_ksp);
    KSPSetFromOptions(*m_ksp);    
    KSPMonitorSet(*m_ksp, &details::KSP_monitor_definition, &m_notifier, PETSC_NULL);
    reuse_preconditioner(true);

    m_setup_called          = false;
};

void ksp_solver_impl::reset_options()
{
    reuse_preconditioner(true);
    m_setup_called  = false;
};

void ksp_solver_impl::finalize_build()
{
    if (m_mg_data.m_levels > 0 && m_mg_data.m_solver_set == false)
        throw error::multigrid_solver_not_set();

    setup();
};
void ksp_solver_impl::setup() const
{
    if (m_setup_called == true)
        return;

    petsc::option_stack st  = m_options_handler->install();

    KSPSetUp(*m_ksp);
    m_setup_called = true;
};

bool ksp_solver_impl::is_partitioned(partition_type& pt) const
{
    PC pc       = get_PC();
    PCType type;
    ::PCGetType(pc, &type);

    std::string stype   = type ? std::string(type) : "";

    if (stype == std::string(PCBJACOBI))
    {
        pt = BJAC;
        return true;
    }
    else if (stype == std::string(PCASM))
    {
        pt = ASM;
        return true;
    }
    else
    {
        return false;
    };
};

void ksp_solver_impl::set_notifier(const notifier& notif)
{
    m_notifier = notif;
    KSPMonitorCancel(*m_ksp);
    KSPMonitorSet(*m_ksp, &details::KSP_monitor_definition, &m_notifier, PETSC_NULL);
}

void ksp_solver_impl::set_ls()
{
    m_ls_problem        = true;
}

Matrix ksp_solver_impl::get_solution_std() const
{
    if (m_ls_problem == false)
        throw error::petsc_solution_std_not_available();

    ::KSP_LSQR* lsqr    = (::KSP_LSQR*)get_KSP()->data;

    if (lsqr->se == nullptr)
        throw error::petsc_solution_std_not_available();

    return petsc::petsc_to_mmlib(lsqr->se);
};

ksp_solver_impl::~ksp_solver_impl()
{}

void ksp_solver_impl::set_nullspace_right(const Matrix& Nr, bool with_const)
{
    m_has_right_null    = false;    
    
    if (Nr.rows() != m_petsc_mat_A->cols())
        throw error::invalid_size2(Nr.rows(), Nr.cols(), m_petsc_mat_A->cols(), Nr.cols());
    if (Nr.cols() > Nr.rows())
        throw error::invalid_size2(Nr.rows(), Nr.cols(), Nr.rows(), Nr.rows());

    if (Nr.all_finite() == false)
        throw error::invalid_unitary_matrix();

    if (Nr.cols() == 0 && with_const == false)
        return;

    //TODO: non real types
    Matrix Nrc  = convert(Nr, mat_code::real_dense);

    petsc::mmlib_to_petsc_null(Nrc, m_null_right, with_const);
    ::MatSetNullSpace(m_petsc_mat_A->get_Aksp(), *m_null_right);

    m_has_right_null    = true;
};

void ksp_solver_impl::set_nullspace_left(const Matrix& Nl, bool with_const)
{
    m_has_left_null     = false;    
    
    if (Nl.rows() != m_petsc_mat_A->rows())
        throw error::invalid_size2(Nl.rows(), Nl.cols(), m_petsc_mat_A->rows(), Nl.cols());
    if (Nl.cols() > Nl.rows())
        throw error::invalid_size2(Nl.rows(), Nl.cols(), Nl.rows(), Nl.rows());

    if (Nl.all_finite() == false)
        throw error::invalid_unitary_matrix();

    if (Nl.cols() == 0 && with_const == false)
        return;

    //TODO: non real types
    Matrix Nlc  = convert(Nl, mat_code::real_dense);

    petsc::mmlib_to_petsc_null(Nlc, m_null_left, with_const);
    ::MatSetTransposeNullSpace(m_petsc_mat_A->get_Aksp(), *m_null_left);

    m_has_left_null = true;
};

void ksp_solver_impl::set_preconditioner(const linsolve_obj& ls)
{
    check_preconditioner(ls.rows(), ls.cols(), this->rows(), this->cols());

    petsc::option_stack st  = m_options_handler->install();

    PC old_pc;
    KSPGetPC(get_KSP(), &old_pc);

    pcshell_ptr pc;
    petsc::linsolve_to_pcshell(ls, old_pc, pc);
    ::KSPSetPC(get_KSP(), *pc);

    m_pc_shell_vec.push_back(pc);
    m_setup_called  = false;
    reuse_preconditioner(true);
};

void ksp_solver_impl::set_preconditioner(const linear_operator& ls)
{
    check_preconditioner(ls.rows(), ls.cols(), this->rows(), this->cols());

    petsc::option_stack st  = m_options_handler->install();

    PC old_pc;
    KSPGetPC(get_KSP(), &old_pc);

    pcshell_ptr pc;
    petsc::linop_to_pcshell(ls, old_pc, pc);
    ::KSPSetPC(get_KSP(), *pc);

    m_pc_shell_vec.push_back(pc);
    m_setup_called  = false;
    reuse_preconditioner(true);
};
void ksp_solver_impl::set_preconditioner(const ksp_solver& ls)
{
    check_preconditioner(ls.rows(), ls.cols(), this->rows(), this->cols());

    ls.m_impl->finalize_build();

    petsc::option_stack st  = m_options_handler->install();

    PC old_pc;
    KSPGetPC(get_KSP(), &old_pc);

    PCSetType(old_pc, PCKSP);
    KSP* ksp_ref;
    PCKSPGetKSPRef(old_pc,&ksp_ref);

    KSP old;
    set_shell_ksp(ksp_ref, ls, old);

    // PCSHELL has proper refcount (?)
    ::KSPDestroy(&old);

    m_setup_called  = false;
    reuse_preconditioner(true);
};

void ksp_solver_impl::reuse_preconditioner(bool val) const
{
    ::KSPSetReusePreconditioner(get_KSP(), val == true ? PETSC_TRUE : PETSC_FALSE);
};

void ksp_solver_impl::divide_on_blocks(Integer n_blocks)
{
    if (this->rows() != this->cols())
        throw error::square_matrix_required(this->rows(), this->cols());

    if (m_pc_shell_vec.size() > 0)
        throw error::unable_partition_external_preconditioner();

    petsc::option_stack st  = m_options_handler->install();

    ::KSPSetType(get_KSP(), KSPPREONLY);
    ::PCSetType(get_PC(), PCBJACOBI);
    ::PCBJacobiSetTotalBlocks(get_PC(), n_blocks, nullptr);    

    options old_opts    = this->options_set();
    old_opts.set(opt::petsc::solver(petsc_solver::preonly));
    old_opts.set(opt::petsc::preconditioner(petsc_precond::none));

    petsc::scoped_options sopts(*m_options_handler, &old_opts);
    ::KSPSetUp(get_KSP());

    m_block_data.resize(n_blocks);
}

void ksp_solver_impl::divide_on_blocks_asm(Integer n_blocks, Integer overlap, asm_composition comp, 
                                asm_restrict restr)
{
    if (this->rows() != this->cols())
        throw error::square_matrix_required(this->rows(), this->cols());

    if (m_pc_shell_vec.size() > 0)
        throw error::unable_partition_external_preconditioner();

    petsc::option_stack st  = m_options_handler->install();

    ::KSPSetType(get_KSP(), KSPPREONLY);
    ::PCSetType(get_PC(), PCASM);
    ::PCASMSetLocalSubdomains(get_PC(), n_blocks, nullptr, nullptr);    
    ::PCASMSetOverlap(get_PC(),overlap);

    switch (comp)
    {
        case asm_composition::additive:
            ::PCASMSetLocalType(get_PC(), PC_COMPOSITE_ADDITIVE);
            break;
        case asm_composition::multiplicative:
            ::PCASMSetLocalType(get_PC(), PC_COMPOSITE_MULTIPLICATIVE);
            break;
    };

    switch (restr)
    {
        case asm_restrict::basic:
            ::PCASMSetType(get_PC(), PC_ASM_BASIC);
            break;
        case asm_restrict::restricted:
            ::PCASMSetType(get_PC(), PC_ASM_RESTRICT);
            break;
        case asm_restrict::interpolate:
            ::PCASMSetType(get_PC(), PC_ASM_INTERPOLATE);
            break;
        case asm_restrict::none:
            ::PCASMSetType(get_PC(), PC_ASM_NONE);
            break;
    };

    options old_opts    = this->options_set();
    old_opts.set(opt::petsc::solver(petsc_solver::preonly));
    old_opts.set(opt::petsc::preconditioner(petsc_precond::none));

    petsc::scoped_options sopts(*m_options_handler, &old_opts);
    ::KSPSetUp(get_KSP());

    m_block_data.resize(n_blocks);
}

void ksp_solver_impl::divide_on_blocks_asm(const std::vector<Matrix>& partition, Integer overlap, 
                        asm_composition comp, asm_restrict restr)
{
    if (this->rows() != this->cols())
        throw error::square_matrix_required(this->rows(), this->cols());

    if (m_pc_shell_vec.size() > 0)
        throw error::unable_partition_external_preconditioner();

    petsc::option_stack st  = m_options_handler->install();

    Integer n_blocks    = partition.size();
    ::KSPSetType(get_KSP(), KSPPREONLY);
    ::PCSetType(get_PC(), PCASM);

    using IS_ptr        = std::shared_ptr<petsc::petsc_IS>;

    std::vector<IS_ptr> vec_is_ptr(n_blocks);
    std::vector<IS>     vec_is(n_blocks);    

    for (Integer i = 0; i < n_blocks; ++i)
    {
        Matrix part_l   = partition[i];

        if (part_l.is_vector() == false)
            throw error::vector_required(part_l.rows(), part_l.cols());

        IS_ptr loc_is;
        const raw::integer_dense& mat_d = part_l.impl<raw::integer_dense>();
        petsc::mmlib_to_IS(mat_d, this->cols(), loc_is);

        vec_is_ptr[i]   = loc_is;
        vec_is[i]       = *loc_is;
    };

    IS* ptr_is = vec_is.data();
    ::PCASMSetLocalSubdomains(get_PC(), n_blocks, ptr_is, ptr_is);
    ::PCASMSetOverlap(get_PC(),overlap);

    switch (comp)
    {
        case asm_composition::additive:
            ::PCASMSetLocalType(get_PC(), PC_COMPOSITE_ADDITIVE);
            break;
        case asm_composition::multiplicative:
            ::PCASMSetLocalType(get_PC(), PC_COMPOSITE_MULTIPLICATIVE);
            break;
    };

    switch (restr)
    {
        case asm_restrict::basic:
            ::PCASMSetType(get_PC(), PC_ASM_BASIC);
            break;
        case asm_restrict::restricted:
            ::PCASMSetType(get_PC(), PC_ASM_RESTRICT);
            break;
        case asm_restrict::interpolate:
            ::PCASMSetType(get_PC(), PC_ASM_INTERPOLATE);
            break;
        case asm_restrict::none:
            ::PCASMSetType(get_PC(), PC_ASM_NONE);
            break;
    };

    options old_opts    = this->options_set();
    old_opts.set(opt::petsc::solver(petsc_solver::preonly));
    old_opts.set(opt::petsc::preconditioner(petsc_precond::none));

    petsc::scoped_options sopts(*m_options_handler, &old_opts);
    ::KSPSetUp(get_KSP());

    m_block_data.resize(n_blocks);
};

void ksp_solver_impl::field_split(const std::vector<Matrix>& partition, field_split_type fst)
{
    if (this->rows() != this->cols())
        throw error::square_matrix_required(this->rows(), this->cols());

    if (m_pc_shell_vec.size() > 0)
        throw error::unable_partition_external_preconditioner();

    petsc::option_stack st  = m_options_handler->install();

    ::KSPSetType(get_KSP(), KSPPREONLY);

    Integer n_blocks    = partition.size();
    ::PCSetType(get_PC(), PCFIELDSPLIT);

    using IS_ptr        = std::shared_ptr<petsc::petsc_IS>;

    std::vector<IS_ptr> vec_is_ptr(n_blocks);

    for (Integer i = 0; i < n_blocks; ++i)
    {
        Matrix part_l   = partition[i];

        if (part_l.is_vector() == false)
            throw error::vector_required(part_l.rows(), part_l.cols());

        IS_ptr loc_is;
        const raw::integer_dense& mat_d = part_l.impl<raw::integer_dense>();
        petsc::mmlib_to_IS(mat_d, this->cols(), loc_is);

        vec_is_ptr[i]   = loc_is;
    };

    for (Integer i = 0; i < n_blocks; ++i)
        ::PCFieldSplitSetIS(get_PC(), nullptr, *vec_is_ptr[i]);

    switch (fst)
    {
        case field_split_type::gauss_seidel:
            ::PCFieldSplitSetType(get_PC(), PC_COMPOSITE_MULTIPLICATIVE);
            break;
        case field_split_type::jacobi:
            ::PCFieldSplitSetType(get_PC(), PC_COMPOSITE_ADDITIVE);
            break;
        case field_split_type::sym_gauss_seidel:
            ::PCFieldSplitSetType(get_PC(), PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE);
            break;
    };

    options old_opts    = this->options_set();
    old_opts.set(opt::petsc::solver(petsc_solver::preonly));
    old_opts.set(opt::petsc::preconditioner(petsc_precond::none));

    petsc::scoped_options sopts(*m_options_handler, &old_opts);
    ::KSPSetUp(get_KSP());

    m_block_data.resize(n_blocks);
};

void ksp_solver_impl::field_split_schur(const Matrix& part1, const Matrix& part2, schur_fact_type sft)
{
    if (this->rows() != this->cols())
        throw error::square_matrix_required(this->rows(), this->cols());

    if (m_pc_shell_vec.size() > 0)
        throw error::unable_partition_external_preconditioner();

    if (part1.length() + part2.length() != this->cols())
        throw error::invalid_partition(part1.length() + part2.length(), this->cols());

    petsc::option_stack st  = m_options_handler->install();

    ::KSPSetType(get_KSP(), KSPPREONLY);
    ::PCSetType(get_PC(), PCFIELDSPLIT);

    using IS_ptr        = std::shared_ptr<petsc::petsc_IS>;

    std::vector<IS_ptr> vec_is_ptr(2);

    {
        Matrix part_l   = part1;

        if (part_l.is_vector() == false)
            throw error::vector_required(part_l.rows(), part_l.cols());

        IS_ptr loc_is;
        const raw::integer_dense& mat_d = part_l.impl<raw::integer_dense>();
        petsc::mmlib_to_IS(mat_d, this->cols(), loc_is);

        vec_is_ptr[0]   = loc_is;
    };
    {
        Matrix part_l   = part2;

        if (part_l.is_vector() == false)
            throw error::vector_required(part_l.rows(), part_l.cols());

        IS_ptr loc_is;
        const raw::integer_dense& mat_d = part_l.impl<raw::integer_dense>();
        petsc::mmlib_to_IS(mat_d, this->cols(), loc_is);

        vec_is_ptr[1]   = loc_is;
    };

    for (Integer i = 0; i < 2; ++i)
        ::PCFieldSplitSetIS(get_PC(), nullptr, *vec_is_ptr[i]);

    ::PCFieldSplitSetType(get_PC(), PC_COMPOSITE_SCHUR);

    switch (sft)
    {
        case schur_fact_type::diag:
            ::PCFieldSplitSetSchurFactType(get_PC(), PC_FIELDSPLIT_SCHUR_FACT_DIAG);
            break;
        case schur_fact_type::lower:
            ::PCFieldSplitSetSchurFactType(get_PC(), PC_FIELDSPLIT_SCHUR_FACT_LOWER);
            break;
        case schur_fact_type::upper:
            ::PCFieldSplitSetSchurFactType(get_PC(), PC_FIELDSPLIT_SCHUR_FACT_UPPER);
            break;
        case schur_fact_type::full:
            ::PCFieldSplitSetSchurFactType(get_PC(), PC_FIELDSPLIT_SCHUR_FACT_FULL);
            break;
    };
    //PCFieldSplitSchurPrecondition

    options old_opts    = this->options_set();
    old_opts.set(opt::petsc::solver(petsc_solver::preonly));
    old_opts.set(opt::petsc::preconditioner(petsc_precond::none));

    petsc::scoped_options sopts(*m_options_handler, &old_opts);
    ::KSPSetUp(get_KSP());

    m_block_data.resize(2);
};

void ksp_solver_impl::set_schur_precond_type(schur_precond_type t)
{
    petsc::option_stack st  = m_options_handler->install();

    switch (t)
    {
        case schur_precond_type::a11:
            PCFieldSplitSetSchurPre(get_PC(), PC_FIELDSPLIT_SCHUR_PRE_A11, nullptr);
            break;
        case schur_precond_type::full:
            PCFieldSplitSetSchurPre(get_PC(), PC_FIELDSPLIT_SCHUR_PRE_FULL, nullptr);
            break;
    }

    //PC_FIELDSPLIT_SCHUR_PRE_SELF,
    //PC_FIELDSPLIT_SCHUR_PRE_SELFP,
    //PC_FIELDSPLIT_SCHUR_PRE_USER,
};

void ksp_solver_impl::set_schur_precond_lsc(const options& opts)
{
    PetscInt N;
    KSP* sub_ksp;

    petsc::option_stack st  = m_options_handler->install();
    
    KSPSetUp(get_KSP());

    PCFieldSplitGetSubKSP(get_PC(), &N, &sub_ksp);

    if (N != 2)
        throw error::invalid_number_of_blocks(N, 2);

    KSP ksp_schur   = sub_ksp[1];
    PC pc_schur;
    
    KSPGetPC(ksp_schur, &pc_schur);
    PCSetType(pc_schur, PCLSC);

    petsc::scoped_options sopts(*m_options_handler, &opts);
    KSPSetUp(ksp_schur);
};

void ksp_solver_impl::divide_on_blocks(const Matrix& block_sizes)
{
    if (this->rows() != this->cols())
        throw error::square_matrix_required(this->rows(), this->cols());

    if (block_sizes.is_vector() == false)
        throw error::vector_required(block_sizes.rows(), block_sizes.cols());

    if (m_pc_shell_vec.size() > 0)
        throw error::unable_partition_external_preconditioner();

    raw::integer_dense mat_d = block_sizes.impl<raw::integer_dense>().make_explicit();
    check_partition(this->rows(), mat_d);

    Integer n_blocks        = mat_d.size();
    const Integer* blocks   = mat_d.ptr();

    petsc::option_stack st  = m_options_handler->install();

    ::KSPSetType(get_KSP(), KSPPREONLY);
    ::PCSetType(get_PC(), PCBJACOBI);
    ::PCBJacobiSetTotalBlocks(get_PC(), n_blocks, blocks);

    options old_opts    = this->options_set();
    old_opts.set(opt::petsc::solver(petsc_solver::preonly));
    old_opts.set(opt::petsc::preconditioner(petsc_precond::none));

    petsc::scoped_options sopts(*m_options_handler, &old_opts);
    ::KSPSetUp(get_KSP());

    m_block_data.resize(n_blocks);
};

KSP ksp_solver_impl::get_block_KSP(Integer block) const
{
    partition_type pt;

    if (is_partitioned(pt) == false)
        throw error::preconditioner_is_not_partitioned();

    ::PetscInt  n_blocks;
    ::PetscInt  first_local;

    petsc::option_stack st  = m_options_handler->install();

    KSP* sub_ksp;

    if (pt == BJAC)
        ::PCBJacobiGetSubKSP(get_PC(), &n_blocks, &first_local, &sub_ksp);
    else
        ::PCASMGetSubKSP(get_PC(), &n_blocks, &first_local, &sub_ksp);
    
    if (block < 1 || block > n_blocks)
        throw error::invalid_block_index(block, n_blocks);

    return sub_ksp[block-1];
};
KSP* ksp_solver_impl::get_block_KSP_ref(Integer block) const
{
    partition_type pt;

    if (is_partitioned(pt) == false)
        throw error::preconditioner_is_not_partitioned();

    ::PetscInt  n_blocks;
    ::PetscInt  first_local;

    petsc::option_stack st  = m_options_handler->install();

    KSP** sub_ksp_ref;

    if (pt == BJAC)
        ::PCBJacobiGetSubKSPRef(get_PC(), &n_blocks, &first_local, &sub_ksp_ref);
    else
        ::PCASMGetSubKSPRef(get_PC(), &n_blocks, &first_local, &sub_ksp_ref);
    
    if (block < 1 || block > n_blocks)
        throw error::invalid_block_index(block, n_blocks);

    return sub_ksp_ref[block-1];
};

KSP ksp_solver_impl::get_block_KSP_fieldsplit(Integer block) const
{
    ::PetscInt  n_blocks;

    petsc::option_stack st  = m_options_handler->install();

    KSP* sub_ksp;

    ::PCFieldSplitGetSubKSP(get_PC(), &n_blocks, &sub_ksp);
    
    if (block < 1 || block > n_blocks)
        throw error::invalid_block_index(block, n_blocks);

    return sub_ksp[block-1];
};
KSP* ksp_solver_impl::get_block_KSP_fieldsplit_ref(Integer block) const
{
    ::PetscInt  n_blocks;

    petsc::option_stack st  = m_options_handler->install();

    KSP** sub_ksp_ref;

    ::PCFieldSplitGetSubKSPRef(get_PC(), &n_blocks, &sub_ksp_ref);
    
    if (block < 1 || block > n_blocks)
        throw error::invalid_block_index(block, n_blocks);

    return sub_ksp_ref[block-1];
};

Matrix ksp_solver_impl::get_partitioning(Integer block) const
{
    partition_type pt;

    if (is_partitioned(pt) == false)
        throw error::preconditioner_is_not_partitioned();

    if (pt != ASM)
        throw error::preconditioner_is_not_partitioned();

    petsc::option_stack st  = m_options_handler->install();

    PetscInt n_blocks;
    IS* is;
    IS* is_loc; 
    ::PCASMGetLocalSubdomains(get_PC(), &n_blocks, &is, &is_loc);

    if (block < 1 || block > n_blocks)
        throw error::invalid_block_index(block, n_blocks);

    Matrix p;    
    petsc::IS_to_int_mat(is[block-1], p);

    return p;
};

Matrix ksp_solver_impl::get_block_preconditioner_matrix(Integer i) const
{
    KSP sub_ksp = get_block_KSP(i);

    Mat Amat, Pmat;
    ::KSPGetOperators(sub_ksp, &Amat, &Pmat);

    return petsc::petsc_mat_to_mmlib(Amat);
};

void ksp_solver_impl::build_block_solver(Integer i, const options& opts0)
{
    KSP sub_ksp = get_block_KSP(i);

    options opts    = correct_options(opts0);
    
    petsc::option_stack st  = m_options_handler->install();

    petsc::scoped_options sopts(*m_options_handler, &opts);
    ::KSPSetFromOptions(sub_ksp);
};

void ksp_solver_impl::set_block_solver(Integer i, const linsolve_obj& precond)
{
    petsc::option_stack st  = m_options_handler->install();

    KSP* sub_ksp_ref = get_block_KSP_ref(i);

    Mat Amat, Pmat;
    ::KSPGetOperators(*sub_ksp_ref, &Amat, &Pmat);
    
    ::PetscInt m, n;
    ::MatGetSize(Amat, &m, &n);

    check_preconditioner(precond.rows(), precond.cols(), m, n);

    KSP old;
    set_shell_ksp(sub_ksp_ref, precond, old);    

    // bjacobi and asm have proper refcount (?)
    ::KSPDestroy(&old);
};
void ksp_solver_impl::set_block_solver(Integer i, const ksp_solver& precond)
{
    petsc::option_stack st  = m_options_handler->install();

    KSP* sub_ksp_ref = get_block_KSP_ref(i);

    Mat Amat, Pmat;
    ::KSPGetOperators(*sub_ksp_ref, &Amat, &Pmat);
    
    ::PetscInt m, n;
    ::MatGetSize(Amat, &m, &n);

    check_preconditioner(precond.rows(), precond.cols(), m, n);

    KSP old;
    set_shell_ksp(sub_ksp_ref, precond, old);

    // bjacobi and asm have proper refcount (?)
    ::KSPDestroy(&old);
};

void ksp_solver_impl::set_composite_preconditioner(petsc_composite_type type,
                            const std::vector<ksp_solver>& precond_list)
{
    if (this->rows() != this->cols())
        throw error::square_matrix_required(this->rows(), this->cols());

    petsc::option_stack st  = m_options_handler->install();

    ::KSPSetType(get_KSP(), KSPPREONLY);
    ::PCSetType(get_PC(), PCCOMPOSITE);

    switch (type)
    {
        case petsc_composite_type::mult:
            ::PCCompositeSetType(get_PC(),PC_COMPOSITE_MULTIPLICATIVE);
            break;
        case petsc_composite_type::sum:
            ::PCCompositeSetType(get_PC(),PC_COMPOSITE_ADDITIVE);
            break;
        case petsc_composite_type::sym_mult:
            ::PCCompositeSetType(get_PC(),PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE);
            break;
    }

    for (size_t i = 0; i < precond_list.size(); ++i)
        ::PCCompositeAddPC(get_PC(), PCSHELL);

    for (size_t i = 0; i < precond_list.size(); ++i)
    {
        const ksp_solver& loc_pc = precond_list[i];
        check_preconditioner(loc_pc.rows(), loc_pc.cols(), this->rows(), this->cols());

        PC sub_pc;
        ::PCCompositeGetPC(get_PC(), i, &sub_pc);

        PCSetType(sub_pc, PCKSP);
        KSP* ksp_ref;
        PCKSPGetKSPRef(sub_pc,&ksp_ref);

        KSP old;
        set_shell_ksp(ksp_ref, loc_pc, old);

        // PCSHELL has proper refcount (?)
        ::KSPDestroy(&old);
    };
};

void ksp_solver_impl::set_composite_preconditioner(petsc_composite_type type,
                            const std::vector<linsolve_obj>& precond_list)
{
    if (this->rows() != this->cols())
        throw error::square_matrix_required(this->rows(), this->cols());

    petsc::option_stack st  = m_options_handler->install();

    ::KSPSetType(get_KSP(), KSPPREONLY);
    ::PCSetType(get_PC(), PCCOMPOSITE);

    switch (type)
    {
        case petsc_composite_type::mult:
            ::PCCompositeSetType(get_PC(),PC_COMPOSITE_MULTIPLICATIVE);
            break;
        case petsc_composite_type::sum:
            ::PCCompositeSetType(get_PC(),PC_COMPOSITE_ADDITIVE);
            break;
        case petsc_composite_type::sym_mult:
            ::PCCompositeSetType(get_PC(),PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE);
            break;
    }

    for (size_t i = 0; i < precond_list.size(); ++i)
        ::PCCompositeAddPC(get_PC(), PCSHELL);

    for (size_t i = 0; i < precond_list.size(); ++i)
    {
        const linsolve_obj& loc_pc = precond_list[i];
        check_preconditioner(loc_pc.rows(), loc_pc.cols(), this->rows(), this->cols());

        PC sub_pc;
        ::PCCompositeGetPC(get_PC(), i, &sub_pc);

        pcshell_ptr pc;
        petsc::linsolve_to_pcshell(loc_pc, sub_pc, pc);        

        m_pc_shell_vec.push_back(pc);
    };
};

void ksp_solver_impl::set_mult_preconditioner(const ksp_solver& R, const ksp_solver& S)
{
    //TODO:
    //this implements PC_COMPOSITE_SPECIAL; but PC_COMPOSITE_SPECIAL does not allow for
    //transposed solves
    linsolve_obj lo = linsolve_seq(R.to_linsolve_object2(), S.to_linsolve_object2());
    return set_preconditioner(lo);
};

void ksp_solver_impl::set_mult_preconditioner(const linsolve_obj& R, const linsolve_obj& S)
{
    //this implements PC_COMPOSITE_SPECIAL; but PC_COMPOSITE_SPECIAL does not allow for
    //transposed solves
    linsolve_obj lo = linsolve_seq(R, S);
    return set_preconditioner(lo);
};

Matrix ksp_solver_impl::set_galerkin_RP(const Matrix& RP)
{
    if (this->rows() != this->cols())
        throw error::square_matrix_required(this->rows(), this->cols());

    if (m_pc_shell_vec.size() > 0)
        throw error::unable_partition_external_preconditioner();

    petsc::option_stack st  = m_options_handler->install();

    ::KSPSetType(get_KSP(), KSPPREONLY);
    ::PCSetType(get_PC(), PCGALERKIN);

    petsc_matrix_ptr pm_RP  = details::petsc_matrix::create(RP);

    Integer N       = this->rows();
    Integer RP_M    = RP.rows();
    Integer RP_N    = RP.cols();

    if (RP_M != N && RP_N != N)
        throw error::invalid_size(RP_M, RP_N);

    if (std::max(RP_M, RP_N) > N)
        throw error::invalid_size(RP_M, RP_N);

    PC pc           = get_PC();
    bool restr;

    if (RP_M < RP_N)
    {
        ::PCGalerkinSetRestriction(pc, pm_RP->get_Aksp());
        restr           = true;
    }
    else
    {
        ::PCGalerkinSetInterpolation(pc, pm_RP->get_Aksp());
        restr           = false;
    };

    Matrix mat_P        = get_preconditioner_matrix();
    Matrix mat_S;

    if (restr == false)
    {
        mat_S           = mmul(RP, mat_P, trans_type::conj_trans);
        mat_S           = mmul(std::move(mat_S), RP, trans_type::no_trans, trans_type::no_trans);
    }
    else
    {
        mat_S           = mmul(RP, mat_P, trans_type::no_trans);
        mat_S           = mmul(std::move(mat_S), RP, trans_type::no_trans, trans_type::conj_trans);
    };

    m_matrix_vec.push_back(pm_RP);
    return mat_S;
};

void ksp_solver_impl::set_galerkin_solver(const ksp_solver& S_solver)
{
    KSP* S_ref;
    ::PCGalerkinGetKSPRef(get_PC(), &S_ref);

    KSP old;
    set_shell_ksp(S_ref, S_solver, old);

    // galerkin has proper refcount (?)
    ::KSPDestroy(&old);
};

Matrix ksp_solver_impl::set_galerkin_RP(const Matrix& R, const Matrix& P)
{
    if (this->rows() != this->cols())
        throw error::square_matrix_required(this->rows(), this->cols());

    if (m_pc_shell_vec.size() > 0)
        throw error::unable_partition_external_preconditioner();

    petsc::option_stack st  = m_options_handler->install();

    ::KSPSetType(get_KSP(), KSPPREONLY);
    ::PCSetType(get_PC(), PCGALERKIN);

    Integer N       = this->rows();
    Integer R_M     = R.rows();
    Integer R_N     = R.cols();
    Integer P_M     = P.rows();
    Integer P_N     = P.cols();

    if (R_N != N)
        throw error::invalid_size(R_M, R_N);
    if (P_M != N)
        throw error::invalid_size(P_M, P_N);

    if (R_M > N)
        throw error::invalid_size(R_M, R_N);
    if (P_N > N)
        throw error::invalid_size(P_M, P_N);
    if (P_N != R_M)
        throw error::invalid_size(P_N, R_M);

    PC pc                   = get_PC();

    petsc_matrix_ptr pm_R  = details::petsc_matrix::create(R);
    petsc_matrix_ptr pm_P  = details::petsc_matrix::create(P);

    ::PCGalerkinSetRestriction(pc, pm_P->get_Aksp());
    ::PCGalerkinSetInterpolation(pc, pm_R->get_Aksp());

    Matrix mat_P            = get_preconditioner_matrix();
    Matrix mat_S            = mmul(P, mat_P, trans_type::conj_trans);
    mat_S                   = mmul(std::move(mat_S), R, trans_type::no_trans, trans_type::conj_trans);

    m_matrix_vec.push_back(pm_R);
    m_matrix_vec.push_back(pm_P);

    return mat_S;
};

void ksp_solver_impl::set_multigrid(Integer n_lev, mg_type type, mg_cycle_type cycle, Integer mc)
{
    if (this->rows() != this->cols())
        throw error::square_matrix_required(this->rows(), this->cols());

    if (m_pc_shell_vec.size() > 0)
        throw error::unable_partition_external_preconditioner();

    if (n_lev < 2)
        throw error::invalid_multigrid_levels(n_lev);

    m_mg_data.set_levels(n_lev);
    petsc::option_stack st  = m_options_handler->install();

    ::KSPSetType(get_KSP(), KSPPREONLY);

    ::PCSetType(get_PC(), PCMG);
    ::PCMGSetLevels(get_PC(), n_lev, nullptr);

    PCMGType form;
    
    switch (type)
    {
        case mg_type::multiplicative:
            form = PC_MG_MULTIPLICATIVE;
            break;
        case mg_type::additive:
            form = PC_MG_ADDITIVE;
            break;
        case mg_type::full:
            form = PC_MG_FULL;
            break;
        case mg_type::kaskade:
            form = PC_MG_KASKADE;
            break;
        default:
            form = PC_MG_MULTIPLICATIVE;
            break;
    }
        
    ::PCMGSetType(get_PC(), form);
    
    PCMGCycleType cycle_type;

    switch(cycle)
    {
        case mg_cycle_type::v:
            cycle_type = PC_MG_CYCLE_V;
            break;
        case mg_cycle_type::w:
            cycle_type = PC_MG_CYCLE_W;
            break;
        default:
            cycle_type = PC_MG_CYCLE_V;
            break;
    }    

    PCMGSetCycleType(get_PC(), cycle_type);

    mc = std::max(mc, 1);
    PCMGMultiplicativeSetCycles(get_PC(), mc);

    mg_create_vectors(n_lev, n_lev - 1, this->rows());
}


Matrix ksp_solver_impl::mg_set_interp_restr_lev(const Matrix& SP, const Matrix& RP)
{
    if (SP.is_square() == false)
        throw error::square_matrix_required(SP.rows(), SP.cols());

    petsc::option_stack st  = m_options_handler->install();

    petsc_matrix_ptr pm_RP  = details::petsc_matrix::create(RP);

    Integer N       = SP.rows();
    Integer RP_M    = RP.rows();
    Integer RP_N    = RP.cols();

    if (RP_M != N && RP_N != N)
        throw error::invalid_size(RP_M, RP_N);

    if (std::max(RP_M, RP_N) > N)
        throw error::invalid_size(RP_M, RP_N);

    PC pc           = get_PC();
    bool restr;

    Integer level   = m_mg_data.m_current_level;

    if (RP_M < RP_N)
    {
        ::PCMGSetRestriction(pc, level, pm_RP->get_Aksp());
        restr           = true;
    }
    else
    {
        ::PCMGSetInterpolation(pc, level, pm_RP->get_Aksp());
        restr           = false;
    };

    ::PCMGSetGalerkin(pc, PETSC_TRUE);

    Matrix mat_P        = SP;
    Matrix mat_S;

    if (restr == false)
    {
        mat_S           = mmul(RP, mat_P, trans_type::conj_trans);
        mat_S           = mmul(std::move(mat_S), RP, trans_type::no_trans, trans_type::no_trans);
    }
    else
    {
        mat_S           = mmul(RP, mat_P, trans_type::no_trans);
        mat_S           = mmul(std::move(mat_S), RP, trans_type::no_trans, trans_type::conj_trans);
    };

    m_matrix_vec.push_back(pm_RP);
    Integer max_levels  = m_mg_data.m_levels;

    mg_create_vectors(max_levels, level - 1, mat_S.rows());

    return mat_S;
};

Matrix ksp_solver_impl::mg_set_interp_restr_lev(const Matrix& SP, const Matrix& R, const Matrix& P)
{
    if (SP.is_square() == false)
        throw error::square_matrix_required(SP.rows(), SP.cols());

    petsc::option_stack st  = m_options_handler->install();

    Integer N       = SP.rows();
    Integer R_M     = R.rows();
    Integer R_N     = R.cols();
    Integer P_M     = P.rows();
    Integer P_N     = P.cols();

    if (R_N != N)
        throw error::invalid_size(R_M, R_N);
    if (P_M != N)
        throw error::invalid_size(P_M, P_N);

    if (R_M > N)
        throw error::invalid_size(R_M, R_N);
    if (P_N > N)
        throw error::invalid_size(P_M, P_N);
    if (P_N != R_M)
        throw error::invalid_size(P_N, R_M);

    PC pc                   = get_PC();

    petsc_matrix_ptr pm_R  = details::petsc_matrix::create(R);
    petsc_matrix_ptr pm_P  = details::petsc_matrix::create(P);

    Integer level           = m_mg_data.m_current_level;

    ::PCMGSetRestriction(pc, level, pm_P->get_Aksp());
    ::PCMGSetInterpolation(pc, level, pm_R->get_Aksp());
    ::PCMGSetGalerkin(pc, PETSC_TRUE);

    Matrix mat_P            = SP;
    Matrix mat_S            = mmul(P, mat_P, trans_type::conj_trans);
    mat_S                   = mmul(std::move(mat_S), R, trans_type::no_trans, trans_type::conj_trans);

    m_matrix_vec.push_back(pm_R);
    m_matrix_vec.push_back(pm_P);
    Integer max_levels      = m_mg_data.m_levels;

    mg_create_vectors(max_levels, level - 1, mat_S.rows());

    return mat_S;
};

Matrix ksp_solver_impl::mg_set_interp_restr(const Matrix& RP)
{
    return mg_set_interp_restr_lev(get_preconditioner_matrix(), RP);
}
Matrix ksp_solver_impl::mg_set_interp_restr(const Matrix& P, const Matrix& R)
{
    return mg_set_interp_restr_lev(get_preconditioner_matrix(), P, R);
}

void ksp_solver_impl::coarser_level()
{
    m_mg_data.coarser_level();
};

void ksp_solver_impl::mg_create_vectors(Integer num_levels, Integer level, Integer size)
{
    Vec v_R, v_X, v_B;
        
    if (level < num_levels - 1)
    {
        //rhs and sol cannot be set on the finest level
        VecCreateSeq(PETSC_COMM_SELF, size, &v_B);
        PCMGSetRhs(get_PC(), level, v_B);

        VecCreateSeq(PETSC_COMM_SELF, size, &v_X);
        PCMGSetX(get_PC(), level, v_X);
    };

    if (level > 0)
    {
        //residuals cannot be set for the coarsest grid
        VecCreateSeq(PETSC_COMM_SELF, size, &v_R);
        PCMGSetR(get_PC(), level, v_R);
    };
};

void ksp_solver_impl::mg_set_solver(const ksp_solver& S_solver)
{       
    if (m_mg_data.m_current_level != 0)
        throw error::invalid_mg_set_solver(m_mg_data.m_current_level);

    // setup must be called before installing solvers; setup clears
    // all solvers
    options old_opts    = this->options_set();
    old_opts.set(opt::petsc::solver(petsc_solver::preonly));
    old_opts.set(opt::petsc::preconditioner(petsc_precond::none));

    petsc::scoped_options sopts(*m_options_handler, &old_opts);
    ::KSPSetUp(get_KSP());

    // initialize smoothers
    Integer n_lev       = m_mg_data.m_levels;

    using ksp_bool      = std::pair<ksp_solver,bool>;
    using ksp_solv_vec  = std::vector<ksp_bool>;

    ksp_solv_vec vec_post    = m_mg_data.get_post();
    ksp_solv_vec vec_pre     = m_mg_data.get_pre();

    for (Integer i = 1; i < n_lev; ++i)
    {
        bool set_post       = vec_post[i].second;
        bool set_pre        = vec_pre[i].second;

        if (set_post)
        {            
            ksp_solver s_post   = vec_post[i].first;

            KSP* S_ref;
            ::PCMGGetSolverUpRef(get_PC(), i, &S_ref);

            KSP old;
            set_shell_ksp(S_ref, s_post, old);

            ::PCMGDestroySolverUpKSP(get_PC(), i, &old);
        };

        if (set_pre)
        {
            ksp_solver s_pre   = vec_pre[i].first;

            KSP* S_ref;
            ::PCMGGetSolverDownRef(get_PC(), i, &S_ref);

            KSP old;
            set_shell_ksp(S_ref, s_pre, old);

            ::PCMGDestroySolverDownKSP(get_PC(), i, &old);
        };
    };

    KSP* S_ref;
    ::PCMGGetSolverDownRef(get_PC(), 0, &S_ref);

    KSP old;
    set_shell_ksp(S_ref, S_solver, old);

    ::PCMGDestroySolverDownKSP(get_PC(), 0, &old);

    m_mg_data.set_solver();
};

void ksp_solver_impl::mg_set_smoother_post(const ksp_solver& smoother, Integer n_steps)
{    
    Integer level   = m_mg_data.m_current_level;

    if (level <= 0)
        throw error::invalid_mg_set_smoother(level);

    n_steps = std::max(n_steps, 0);
    KSPSetTolerances(smoother.get_KSP(), PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, n_steps);

    m_mg_data.set_smoother_post(smoother);
};

void ksp_solver_impl::mg_set_smoother_pre(const ksp_solver& smoother, Integer n_steps)
{    
    Integer level   = m_mg_data.m_current_level;

    if (level <= 0)
        throw error::invalid_mg_set_smoother(level);

    n_steps = std::max(n_steps, 0);
    KSPSetTolerances(smoother.get_KSP(), PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, n_steps);
    m_mg_data.set_smoother_pre(smoother);
};

Integer ksp_solver_impl::number_of_blocks() const
{
    ::PetscInt  n_block;
    ::PetscInt  first_local;

    KSP* sub_ksp;

    ::PCBJacobiGetSubKSP(get_PC(), &n_block, &first_local, &sub_ksp);

    return n_block;
};

void ksp_solver_impl::check_partition(Integer size, const mr::integer_dense& part) const
{
    Integer s           = part.size();
    const Integer* ptr  = part.ptr();
    Integer sum         = 0;

    for (Integer i = 0; i < s; ++i)
    {
        Integer val     = ptr[i];
        sum             += val;

        if (val <= 0)
            throw error::invalid_partition_block(val);
    };

    if (sum != size)
        throw error::invalid_partition_size(size, sum);
};

void ksp_solver_impl::remove_nullspace_right()
{
    if(m_has_right_null == true)
    {
        ::MatSetNullSpace(m_petsc_mat_A->get_Aksp(),nullptr);
        m_null_right    = petsc_mat_null();
    }

    m_has_right_null    = false;
};

void ksp_solver_impl::remove_nullspace_left()
{
    if(m_has_left_null == true)
    {
        ::MatSetTransposeNullSpace(m_petsc_mat_A->get_Aksp(),nullptr);
        m_null_left     = petsc_mat_null();
    }

    m_has_left_null     = false;
};

KSP ksp_solver_impl::get_KSP() const
{
    return *m_ksp;
};

PC ksp_solver_impl::get_PC() const
{
    PC pc;
    ::KSPGetPC(get_KSP(), &pc);
    return pc;
};

std::string ksp_solver_impl::get_solver() const
{
    KSPType type;
    ::KSPGetType(get_KSP(), &type);

    return type? std::string(type) : "";
};

std::string ksp_solver_impl::get_precond() const
{
    PCType type;
    ::PCGetType(get_PC(), &type);

    return type? std::string(type) : "";
}

options ksp_solver_impl::options_missing() const
{
    return m_missing_opts;
};

options ksp_solver_impl::options_set() const
{
    return m_set_options;
};

void ksp_solver_impl::set_option(const std::string& opt, const std::string& val)
{
    petsc::option_stack st  = m_options_handler->install();
    ::PetscOptionsSetValue(opt.c_str(), val.c_str());
};

void ksp_solver_impl::add_missing(const options& miss)
{
    m_missing_opts.set(miss);
};

value_code ksp_solver_impl::get_value_code() const
{
    return m_petsc_mat_A->get_value_code();
};

Integer ksp_solver_impl::rows() const
{
    return m_petsc_mat_A->rows();
};

Integer ksp_solver_impl::cols() const
{
    return m_petsc_mat_A->cols();
};

bool ksp_solver_impl::all_finite() const
{
    return m_petsc_mat_A->all_finite();
};

bool ksp_solver_impl::is_from_matrix() const
{
    return m_petsc_mat_A->is_shell() == false;
};

bool ksp_solver_impl::is_preconditioner_from_matrix() const
{
    if (m_pc_shell_vec.size() > 0)
        return false;

    PCType type;
    ::PCGetType(get_PC(), &type);

    if (type == nullptr)
        return false;

    std::string stype   = std::string(type);

    if (stype == PCSHELL || stype == PCCOMPOSITE)
        return false;

    return true;
};

Matrix ksp_solver_impl::get_preconditioner_matrix() const
{
    Mat Amat, Pmat;
    ::KSPGetOperators(get_KSP(), &Amat, &Pmat);

    return petsc::petsc_mat_to_mmlib(Pmat);
};

linear_operator ksp_solver_impl::get_linear_operator() const
{
    return m_petsc_mat_A->get_linop();
};

Matrix ksp_solver_impl::get_matrix() const
{
    return m_petsc_mat_A->get_matrix();
};

void ksp_solver_impl::fieldsplit_build_block_solver(Integer i, const options& opts0)
{
    KSP sub_ksp = get_block_KSP_fieldsplit(i);

    options opts            = correct_options(opts0);    
    petsc::option_stack st  = m_options_handler->install();

    petsc::scoped_options sopts(*m_options_handler, &opts);
    ::KSPSetFromOptions(sub_ksp);
};

void ksp_solver_impl::fieldsplit_set_block_solver(Integer i, const linsolve_obj& precond)
{
    petsc::option_stack st  = m_options_handler->install();

    KSP* sub_ksp_ref = get_block_KSP_fieldsplit_ref(i);

    Mat Amat, Pmat;
    ::KSPGetOperators(*sub_ksp_ref, &Amat, &Pmat);
    
    ::PetscInt m, n;
    ::MatGetSize(Amat, &m, &n);

    check_preconditioner(precond.rows(), precond.cols(), m, n);

    KSP old;
    set_shell_ksp(sub_ksp_ref, precond, old);

    // fieldsplit has proper refcount (?)
    ::KSPDestroy(&old);
}
void ksp_solver_impl::fieldsplit_set_block_solver(Integer i, const ksp_solver& precond)
{
    petsc::option_stack st  = m_options_handler->install();

    KSP* sub_ksp_ref = get_block_KSP_fieldsplit_ref(i);

    Mat Amat, Pmat;
    ::KSPGetOperators(*sub_ksp_ref, &Amat, &Pmat);
    
    ::PetscInt m, n;
    ::MatGetSize(Amat, &m, &n);

    check_preconditioner(precond.rows(), precond.cols(), m, n);

    KSP old;
    set_shell_ksp(sub_ksp_ref, precond, old);

    // fieldsplit has proper refcount (?)
    ::KSPDestroy(&old);
}

Matrix ksp_solver_impl::fieldsplit_get_block_preconditioner_matrix(Integer i) const
{    
    KSP sub_ksp = get_block_KSP_fieldsplit(i);

    Mat Amat, Pmat;
    ::KSPGetOperators(sub_ksp, &Amat, &Pmat);

    return petsc::petsc_mat_to_mmlib(Amat);
}

Matrix ksp_solver_impl::fieldsplit_get_schur_matrix() const
{   
    Mat S0;
    ::PCFieldSplitSchurGetS(get_PC(), &S0);

    Mat S;
    MatSchurComplementGetPmat(S0, MAT_INITIAL_MATRIX, &S);
    return petsc::petsc_mat_to_mmlib(S);
}


Matrix ksp_solver_impl::fieldsplit_get_partitioning(Integer block) const
{    
    Integer n_blocks = fieldsplit_number_of_blocks();

    if (block < 1 || block > n_blocks)
        throw error::invalid_block_index(block, n_blocks);

    std::ostringstream is_str;
    is_str << block - 1;
    std::string is_name = is_str.str();

    petsc::option_stack st  = m_options_handler->install();

    IS is;
    PCFieldSplitGetIS(get_PC(), is_name.c_str(), &is);

    if (is == nullptr)
        throw error::invalid_block_index(block, n_blocks);

    Matrix p;    
    petsc::IS_to_int_mat(is, p);

    return p;
};

Integer ksp_solver_impl::fieldsplit_number_of_blocks() const
{
    ::PetscInt  n_block;
    KSP* sub_ksp;

    ::PCFieldSplitGetSubKSP(get_PC(), &n_block, &sub_ksp);

    return n_block;
};

tuple<Matrix, bool> ksp_solver_impl::solve(const ksp_solver& owner, const matcl::Matrix& b, const Matrix& start0, 
                                          trans_type tA, bool has_start)
{
    petsc::option_stack st  = m_options_handler->install();

    Integer M, N;

    if (tA == trans_type::no_trans)
    {
        M   = m_petsc_mat_A->rows();
        N   = m_petsc_mat_A->cols();
    }
    else
    {
        M   = m_petsc_mat_A->cols();
        N   = m_petsc_mat_A->rows();        
    };

    if (b.rows() != M)
        throw error::error_lsolve(M, N, b.rows(), b.cols());

    if (has_start == true)
    {
        if(N != start0.rows() || start0.cols() != b.cols())
            throw error::error_lsolve(M, N, start0.rows(), start0.cols());
    };

    //TODO: complex?, integer?
    bool is_compl   = matrix_traits::is_float_complex(b.get_value_code());

    if (has_start == true)
        is_compl    = is_compl || matrix_traits::is_float_complex(start0.get_value_code());

    if (is_compl)
        throw error::error_complex_value_type_not_allowed();

    if (b.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("ksp_solver_impl");

    if (has_start == true && start0.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("ksp_solver_impl");

    using ret_type          = tuple<Matrix,bool>;

    bool isv                = m_petsc_mat_A->all_finite();
    isv                     = isv && b.all_finite();

    if (has_start == true && start0.all_finite() == false)
        has_start           = false;

    value_code vc           = matrix_traits::unify_value_types(m_petsc_mat_A->get_value_code(), 
                                                            b.get_value_code());
    vc                      = matrix_traits::unify_value_types(vc, value_code::v_float);
    vc                      = matrix_traits::real_value_type(vc);

    if (isv == false)
    {
        Matrix solution     = details::make_nan_matrix(M, b.cols(), vc);
        return ret_type(solution,true);
    };   

    //lock petsc
    petsc::petsc_lock_helper lock;

    using Mat_D             = raw::real_dense;

    //TODO: non real matrix
    Matrix b_full_real      = convert(b, mat_code::real_dense);
    raw::real_dense b_rep   = b_full_real.get_impl_unique<raw::real_dense>();
    
    if (b_rep.cols() == 0) 
    {
        Matrix solution     = spzeros(M, b_rep.cols(), 0, vc);
        return ret_type(solution,true);
    }

    setup();

    if (has_start == true)
        ::KSPSetInitialGuessNonzero(*m_ksp, PETSC_TRUE);
    else
        ::KSPSetInitialGuessNonzero(*m_ksp, PETSC_FALSE);

    bool has_failed = false;

    if (m_notifier)
        m_notifier->begin_solve(owner, m_petsc_mat_A->get_linop(), b_rep.cols());

    Real* b_ptr     = b_rep.ptr();
    Integer b_ld    = b_rep.ld();

    Matrix x0       = has_start ? start0 : make_dense_noinit(N, b.cols(), vc);
    Mat_D x_rep     = x0.get_impl_unique<Mat_D>();

    Real* x_ptr     = x_rep.ptr();
    Integer x_ld    = x_rep.ld();

    for (int i = 0; i < b_rep.cols(); i++)
    {
        if (m_notifier)
        {
            Matrix RHS  = Matrix(b_rep,false)(colon(), i + 1);
            m_notifier->begin_solve_vector(i+1);
            m_notifier->begin_iterations(RHS);
        };

        mp::smart_vec bksp(mp::smart_vec::vec_create_seq_with_array(PETSC_COMM_SELF,1, M, b_ptr + i * b_ld));
        mp::smart_vec xksp(mp::smart_vec::vec_create_seq_with_array(PETSC_COMM_SELF, 1, N, x_ptr + i * x_ld));

        try
        {
            if (tA == trans_type::no_trans)
                KSPSolve(*m_ksp, bksp, xksp);
            else
                KSPSolveTranspose(*m_ksp, bksp, xksp);

            if (m_notifier) 
            {
                linear_operator op_A    = linop_trans(this->get_linear_operator(), tA);
                Matrix RHS              = Matrix(b_rep,false)(colon(), i + 1);
                Matrix sol              = Matrix(x_rep,false)(colon(), i + 1);

                m_notifier->end_iterations(op_A, RHS, sol);
            };
        }
        catch(matcl::error::mmlib_petsc_error & ex)
        {
            if (m_notifier)
            {
                m_notifier->report_error(ex.what());
                m_notifier->end_solve_vector();
                m_notifier->end_solve();
            }

            has_failed = true;

            // do not use rethrow here
            // exception propagation through PETSC is suspicious
            throw error::internal_petsc_error(ex.what());
        }

        KSPConvergedReason reason;
        KSPGetConvergedReason(*m_ksp, &reason);               

        if (reason < 0)
        {
            if(m_notifier)
                m_notifier->report_divergence(KSPConvergedReasons[reason]);

            has_failed = true;
        }
        else if (reason > 0)
        {
            if(m_notifier)
                m_notifier->report_convergence(KSPConvergedReasons[reason]);
        }

        if(m_notifier)
            m_notifier->end_solve_vector();
    }
    
    if (m_notifier) 
        m_notifier->end_solve();

    return ret_type(Matrix(x_rep,false), has_failed == false);
}

petsc_solver ksp_solver_impl::get_solver_type(const options& opts)
{
    options miss;
    petsc::petsc_options_handler oh(&opts, &miss);

    std::string name    = "-ksp_type";
    std::string val;
    bool has = oh.get_option<std::string>("", name, name, val);

    if (has == false)
        return petsc_solver::default_val;

    return petsc::converter_solver().to_code(val);
};

bool ksp_solver_impl::check_nonsquare_allowed(const options& opts)
{
    petsc_solver solver = get_solver_type(opts);

    //least square solvers
    if (solver == petsc_solver::lsqr)
        return true;

    return false;
}

options ksp_solver_impl::correct_options(const options& opts)  const
{
    //remove options that crush petsc
    options opt = opts;

    if (m_ls_problem == true)
    {
        //lsqr is the only least squares solver
        opt.set(opt::petsc::solver(petsc_solver::lsqr));
    };

    opt.remove(string_option<bool>("-ksp_monitor_lg_range"));
    opt.remove(string_option<bool>("-ksp_monitor_lg_residualnorm"));
    opt.remove(string_option<bool>("-ksp_monitor_lg_true_residualnorm"));
    opt.remove(string_option<std::string>("-ksp_monitor_python"));
    opt.remove(matcl::string_option<bool>("-ksp_monitor_solution"));
    opt.remove(matcl::string_option<bool>("-ksp_plot_eigenvalues"));    
    opt.remove(matcl::string_option<std::string>("-ksp_gmres_krylov_monitor"));
    opt.remove(matcl::string_option<std::string>("-pc_svd_monitor"));

    //options ignored or handled in other way
    opt.remove(matcl::string_option<std::string>("-ksp_convergence_test"));
    opt.remove(matcl::string_option<bool>("-ksp_error_if_not_converged"));    
    opt.remove(matcl::string_option<bool>("-ksp_diagonal_scale_fix"));        
    //diagonal scaling is removed; PETSC scaling is useless
    opt.remove(matcl::string_option<bool>("-ksp_diagonal_scale"));
    //use ARPACK instead PETSC; not working for all solvers; for some gives
    //wrong estimates (LSQR for example)
    opt.remove(matcl::string_option<bool>("-ksp_compute_eigenvalues"));
    opt.remove(matcl::string_option<bool>("-ksp_compute_singularvalues"));
    opt.remove(matcl::string_option<std::string>("-ksp_monitor_singular_value"));
    opt.remove(matcl::string_option<bool>("-ksp_constant_null_space"));
    opt.remove(matcl::string_option<std::string>("-ksp_monitor"));    
    opt.remove(matcl::string_option<std::string>("-ksp_monitor_short"));        
    opt.remove(matcl::string_option<bool>("-ksp_test_null_space"));
    opt.remove(matcl::string_option<std::string>("-ksp_monitor_true_residual"));
    opt.remove(matcl::string_option<std::string>("-ksp_monitor_max"));
    opt.remove(matcl::string_option<std::string>("-ksp_monitor_range"));    
    opt.remove(matcl::string_option<std::string>("-ksp_lsqr_monitor"));
    //opt.remove(matcl::string_option<std::string>("-pc_factor_mat_solver_package"));

    opt.remove(matcl::string_option<bool>("-pc_factor_in_place"));
    	
    opt.remove(matcl::string_option<bool>("-ksp_fischer_guess"));
    opt.remove(matcl::string_option<bool>("-ksp_knoll"));
    opt.remove(string_option<bool>("-ksp_initial_guess_nonzero"));
    opt.remove(string_option<bool>("-ksp_reuse_preconditioner"));
    opt.remove(string_option<bool>("-ksp_diagonal_scale"));
    opt.remove(string_option<bool>("-ksp_diagonal_scale_fix"));
    opt.remove(string_option<bool>("-pc_mg_dump_matlab"));
    opt.remove(string_option<bool>("-pc_mg_dump_binary"));

    //ptscotch is the only available partitioner
    opt.set(matcl::string_option<std::string>("-mat_partitioning_type","ptscotch"));

    return opt;
};

}};

#endif