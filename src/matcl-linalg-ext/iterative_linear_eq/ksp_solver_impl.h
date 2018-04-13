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

#if 0
TODO

#include "matcl-linalg/general/config_linalg.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-linalg/linear_eq/ksp_solver.h"

#include <memory>

#include "matcl-linalg/petsc_utils/petsc_objects.h"
#include "matcl-linalg/petsc_utils/petsc_option.h"
#include "matcl-linalg/petsc_utils/petsc_library.h"
#include "matcl-linalg/petsc_utils/petsc_matrix.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

namespace matcl { namespace details
{

namespace mp    = matcl::petsc;
namespace mr    = matcl::raw;

// this enums should not be public; PETSC schur is not ready

// type of preconditioner build for a split A = [A00, A012; A10 A11;] using Schur preconditioner
// the Schur complement is S = A11 - A10 x inv(A00) x A01
enum class schur_fact_type
{
    diag,               /// preconditioner is given as [inv(A00) 0; 0 -inv(S)]
    lower,              /// preconditioner is given by [A00 0; A10 S]^-1
    upper,              /// preconditioner is given by [A00 A01; 0 S]^-1
    full                /// preconditioner is given by [I - inv(A00)xA01; 0 I] x [inv(A00) 0; 0 inv(S)]
                        /// x [I 0; -A10xinv(A00) I]
};

enum class schur_precond_type
{
    a11,
    full
};

struct mg_data
{
    using ksp_bool  = std::pair<ksp_solver,bool>;
    using ksp_vec   = std::vector<ksp_bool>;

    Integer     m_levels;
    Integer     m_current_level;
    bool        m_solver_set;
    bool        m_smoother_set;
    ksp_vec     m_smooth_post;
    ksp_vec     m_smooth_pre;

    void        set_levels(Integer nlev);
    void        coarser_level();
    void        set_solver();
    void        set_smoother_post(const ksp_solver& s);
    void        set_smoother_pre(const ksp_solver& s);

    ksp_vec&&   get_post();
    ksp_vec&&   get_pre();

    mg_data();
};

//---------------------------------------------------------------
//                  ksp_solver_impl
//---------------------------------------------------------------
class ksp_solver_impl
{
    public:
        enum partition_type { BJAC, ASM };

    private:
        using petsc_lock        = std::shared_ptr<matcl::petsc::petsc_lock_helper>;
        using petsc_matrix_ptr  = std::shared_ptr<details::petsc_matrix>;
        using ksp_solver_solver = std::shared_ptr<matcl::petsc::smart_KSP>;
        using petsc_mat_null    = std::shared_ptr<petsc::petsc_null_mat>;
        using notifier          = std::shared_ptr<ksp_solver_notifier>;
        using reorder_ret       = tuple<permvec,permvec>;
        using pcshell_ptr       = std::shared_ptr<petsc::petsc_pcshell>;
        using pcshell_vec       = std::vector<pcshell_ptr>;
        using ptr_type          = std::shared_ptr<ksp_solver_impl>;
        using block_data        = std::vector<pcshell_ptr>;
        using options_ptr       = std::shared_ptr<options>;
        using opt_handler       = mp::petsc_options_handler;
        using opt_handler_ptr   = opt_handler::sptr_type;
        using matrix_vec        = std::vector<petsc_matrix_ptr>;
        using ksp_vec           = std::vector<ksp_solver>;

    public: 
        ksp_solver_impl(const Matrix& A_in, const options& opts);
        ksp_solver_impl(const linear_operator& A_in, const options& opts);

        ksp_solver_impl(const Matrix& A_in, const Matrix& P_in, const options& opts);
        ksp_solver_impl(const linear_operator& A_in, const Matrix& P_in, const options& opts);

        ~ksp_solver_impl();

        void                set_id(const std::string& id);
        const std::string&  get_id() const;

        void                reset_operator(const Matrix& A);
        void                reset_operator(const linear_operator& A);
        void                rebuild_preconditioner(const Matrix& B, const options& opts);

        void                set_notifier(const notifier& notifier);
        options             options_missing() const;    
        options             options_set() const;    
        tuple<Matrix,bool>  solve(const ksp_solver& owner, const matcl::Matrix& b,const Matrix& start, 
                                trans_type tA,  bool has_start);
        void                set_option(const std::string& opt, const std::string& val);
        void                add_missing(const options& miss);

        void                set_nullspace_right(const Matrix& Nr, bool with_const);
        void                set_nullspace_left(const Matrix& Nr, bool with_const);
        void                remove_nullspace_right();
        void                remove_nullspace_left();
        void                set_preconditioner(const linsolve_obj& ls);
        void                set_preconditioner(const linear_operator& ls);
        void                set_preconditioner(const ksp_solver& ls);
        void                reuse_preconditioner(bool val) const;

        void                divide_on_blocks(Integer n_blocks);
        void                divide_on_blocks(const Matrix& block_sizes);
        void                divide_on_blocks_asm(Integer n_blocks, Integer overlap, asm_composition comp, 
                                asm_restrict restr);
        void                divide_on_blocks_asm(const std::vector<Matrix>& partition, Integer overlap, 
                                asm_composition comp, asm_restrict restr);

        void                build_block_solver(Integer i, const options& opts);
        Matrix              get_block_preconditioner_matrix(Integer i) const;
        void                set_block_solver(Integer i, const linsolve_obj& precond);
        void                set_block_solver(Integer i, const ksp_solver& precond);
        Matrix              get_partitioning(Integer i) const;
        Integer             number_of_blocks() const;

        void                set_composite_preconditioner(petsc_composite_type type,
                                    const std::vector<ksp_solver>& precond_list);
        void                set_composite_preconditioner(petsc_composite_type type,
                                    const std::vector<linsolve_obj>& precond_list);
        void                set_mult_preconditioner(const ksp_solver& R, const ksp_solver& S);
        void                set_mult_preconditioner(const linsolve_obj& R, const linsolve_obj& S);
        Matrix              set_galerkin_RP(const Matrix& P);        
        Matrix              set_galerkin_RP(const Matrix& R, const Matrix& P);
        void                set_galerkin_solver(const ksp_solver& S_solver);

        void                field_split(const std::vector<Matrix>& partition, field_split_type fst);
        void                field_split_schur(const Matrix& p1, const Matrix& p2, schur_fact_type sft);
        void                set_schur_precond_lsc(const options& opts);
        void                set_schur_precond_type(schur_precond_type t);
        void                set_schur_precond_type(const linsolve_obj& precond);

        void                fieldsplit_build_block_solver(Integer i, const options& opts);
        Matrix              fieldsplit_get_block_preconditioner_matrix(Integer i) const;
        Matrix              fieldsplit_get_schur_matrix() const;
        void                fieldsplit_set_block_solver(Integer i, const ksp_solver& precond);
        void                fieldsplit_set_block_solver(Integer i, const linsolve_obj& precond);
        Matrix              fieldsplit_get_partitioning(Integer i) const;
        Integer             fieldsplit_number_of_blocks() const;

        void                set_multigrid(Integer lev, mg_type type, mg_cycle_type cycle, Integer mc);
        Matrix              mg_set_interp_restr(const Matrix& RP);
        Matrix              mg_set_interp_restr(const Matrix& P, const Matrix& R);
        Matrix              mg_set_interp_restr_lev(const Matrix& S, const Matrix& RP);
        Matrix              mg_set_interp_restr_lev(const Matrix& S, const Matrix& P, const Matrix& R);
        void                mg_create_vectors(Integer num_levels, Integer level, Integer size);
        void                mg_set_solver(const ksp_solver& S_solver);        
        void                mg_set_smoother_pre(const ksp_solver& smoother, Integer n_steps);
        void                mg_set_smoother_post(const ksp_solver& smoother, Integer n_steps);
        void                coarser_level();

        void                set_ls();
        KSP                 get_KSP() const;
        PC                  get_PC() const;
        Matrix              get_solution_std() const;
        void                reset_options(const options& opts);
        void                reset_options();
        void                set_tolerances(Integer maxit, Real rtol, Real atol, Real dtol);
        linear_operator     get_linear_operator() const;
        Matrix              get_matrix() const;
        value_code          get_value_code() const;
        Matrix              get_preconditioner_matrix() const;
        bool                is_preconditioner_from_matrix() const;

        Integer             rows() const;
        Integer             cols() const;
        bool                all_finite() const;
        bool                is_from_matrix() const;

        bool                is_partitioned(partition_type& pt) const;        
        std::string         get_solver() const;
        std::string         get_precond() const;
        void                finalize_build();        

    private:
        void                initialize();
        void                reinitialize(const petsc_matrix_ptr& pm_A);
        void                restore_null_spaces(const ksp_solver_impl& base_ksp);
        static Matrix       preprocess_mat(const Matrix& A);
        static bool         check_nonsquare_allowed(const options& opts);
        static petsc_solver get_solver_type(const options& opts);
        options             correct_options(const options& opts)  const;
        void                setup() const;
        void                check_preconditioner(Integer rows, Integer cols, 
                                Integer mat_rows, Integer mat_cols) const;
        void                check_partition(Integer size, const mr::integer_dense& part) const;
        KSP                 get_block_KSP(Integer block) const;
        KSP*                get_block_KSP_ref(Integer block) const;
        KSP                 get_block_KSP_fieldsplit(Integer block) const;
        KSP*                get_block_KSP_fieldsplit_ref(Integer block) const;

        // return old KSP which should be destoyed; Petsc reference couting is not perfect
        // and it depends on context whether object can be destroyed or cannot
        void                set_shell_ksp(KSP* ksp_ref, const linear_operator& solv, KSP& old);
        void                set_shell_ksp(KSP* ksp_ref, const linsolve_obj& solv, KSP& old);
        void                set_shell_ksp(KSP* ksp_ref, const ksp_solver& solv, KSP& old);

    private:
        //just to increase refcount; must be initialized before m_ksp_operator
        Matrix              m_mat_A;
        Matrix              m_mat_P;
        std::string         m_id;
        bool                m_has_right_null;
        bool                m_has_left_null;
        bool                m_ls_problem;
        mutable bool        m_setup_called;
        petsc_mat_null      m_null_right;
        petsc_mat_null      m_null_left;
        pcshell_vec         m_pc_shell_vec;
        block_data          m_block_data;
        matrix_vec          m_matrix_vec;
        ksp_vec             m_ksp_vec;

        //petsc_operator stores unique matrix
        petsc_matrix_ptr    m_petsc_mat_A;
        petsc_matrix_ptr    m_petsc_mat_P;

        notifier            m_notifier;
        petsc_lock          m_ksp_lock;
        ksp_solver_solver   m_ksp;        
        options             m_missing_opts;
        options             m_set_options;
        opt_handler_ptr     m_options_handler;
        mg_data             m_mg_data;

        ksp_solver_impl(const ksp_solver_impl&) = delete;
        ksp_solver_impl& operator=(const ksp_solver_impl&) = delete;
};
 
}};

#endif