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

#include "matcl-core/error/exception_classes.h"
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-matrep/matrix/matrix.h"

#include "matcl-core/error/safe_string_message.h"

#pragma warning(push)
#pragma warning(disable:4251)

//TODO: cleanup

namespace matcl { namespace error
{
    
class exception_message_linalg;
using exception_message_linalg_ptr  = std::shared_ptr<exception_message_linalg>;

MATCL_LINALG_EXPORT void                            set_global_messanger_linalg(exception_message_linalg_ptr msg);
MATCL_LINALG_EXPORT exception_message_linalg_ptr    get_global_messanger_linalg();

//TODO: warning filtering

class MATCL_LINALG_EXPORT exception_message_linalg
{
    public:
        virtual ~exception_message_linalg(){};
        
        virtual const char* error_general(const std::string& msg) = 0;		
        virtual const char* error_lu(const std::string& msg) = 0;	
        virtual const char* error_arpack(const std::string& msg) = 0;
        virtual const char* error_lsolve(Integer r1,Integer c1, Integer r2, Integer c2) = 0;
        virtual const char* error_lsolve_rev(Integer r1,Integer c1, Integer r2, Integer c2) = 0;
        virtual const char* error_fft(Integer status) = 0;
        virtual const char* error_singular() = 0;
        virtual const char* error_gen_sym_eigen_job() = 0;
        virtual const char* error_schur_incorrect_factors() = 0;
        virtual const char* error_schur_U_not_computed() = 0;        
        virtual const char* error_gschur_QZ_not_computed() = 0;
        virtual const char* error_gschur_incorrect_factors() = 0;
        virtual const char* error_size_qrs(Integer r1,Integer c1) = 0;
        virtual const char* error_qrs(const Integer code) = 0;
        virtual const char* error_size_eig(Integer r1,Integer c1) = 0;
        virtual const char* error_size_geig(Integer r1,Integer c1, Integer r2, Integer c2) = 0;
        virtual const char* error_size_sylv(Integer r1,Integer c1, Integer r2, Integer c2, 
                                Integer r3, Integer c3) = 0;
        virtual const char* error_size_gsylv(Integer rA1_, Integer cA1_, Integer rA2_, Integer cA2_,
                                Integer rB1_, Integer cB1_, Integer rB2_, Integer cB2_,
                                Integer rC1_, Integer cC1_, Integer rC2_, Integer cC2_) = 0;
        virtual const char* error_size_ricc() = 0;
        virtual const char* error_size_lyapunov() = 0;
        virtual const char* error_schur() = 0;
        virtual const char* error_ldl() = 0;
        virtual const char* error_size_ldl(Integer r1,Integer c1) = 0;
        virtual const char* error_sylv() = 0;
        virtual const char* error_gsylv() = 0;
        virtual const char* error_ricc(Integer info) = 0;
        virtual const char* error_lyapunov() = 0;
        virtual const char* error_gschur() = 0;
        virtual const char* error_norm_type() = 0;
        virtual const char* error_svd() = 0;
        
        virtual const char* error_hess_nonsq() = 0;
        virtual const char* error_gen_hess(Integer Ar, Integer Ac, Integer Br, Integer Bc) = 0;
        virtual const char* hess_eig_failed(const Matrix& VL, const Matrix& VR, const Matrix& fail_L, 
                                const Matrix& fail_R, bool eval_L, bool eval_R, Integer failed_eig) = 0;
        virtual const char* error_hess_eig_invalid_init(Integer I_rows, Integer I_cols, Integer H_size, 
                                Integer E_size, bool left_init) = 0;
        virtual const char* error_hess_eig_invalid_conj_pair(Integer pos, bool missing_conj) = 0;
        virtual const char* error_hess_eig_invalid_eig(Integer E_rows, Integer E_cols, Integer H_size) = 0;
        virtual const char* error_hess_eig_her_complex_eig(Integer pos) = 0;

        virtual const char* error_svd_inf_nan() = 0;
        virtual const char* error_schur_sel() = 0;
        virtual const char* error_schur_sel_comp() = 0;
        virtual const char* error_schur_cond_vec(Integer N, Integer M_sel, Integer N_vec, 
                                Integer M_vec, bool left_vector) = 0;
        virtual const char* error_schur_eig_cluster_size(Integer M, Integer N) = 0;
        virtual const char* error_gschur_sel_comp() = 0;
        virtual const char* error_nonsymh_chol() = 0;
        virtual const char* error_nonposdef(bool allow_semipos) = 0;
        virtual const char* error_size_tridiag_sym(Integer A_rows, Integer A_cols, 
                                Integer B_rows, Integer B_cols) = 0;
        virtual const char* error_size_tridiag_nsym(Integer Dm1_rows, Integer Dm1_cols,Integer D0_rows, 
                                Integer D0_cols, Integer Dp1_rows, Integer Dp1_cols) = 0;
        virtual const char* error_symeigen_complex_diagonal() = 0;

        virtual const char* invalid_permutation_vector() = 0;
        virtual const char* error_cholmod(const std::string& chol_msg) = 0;        

        virtual const char* error_complex_value_type_not_allowed() = 0;                
        virtual const char* symmetric_matrix_required() = 0;
        virtual const char* tridiagonal_matrix_required(Integer rows, Integer cols, Integer ld, Integer ud) = 0;
        virtual const char* error_upper_hess_matrix_required() = 0;
        virtual const char* invalid_unitary_matrix() = 0;
        virtual const char* error_invalid_chol_update_nonposdef() = 0;
        virtual const char* error_invalid_chol_factor(bool upper) = 0;
        virtual const char* error_invalid_qr_factor() = 0;
        virtual const char* error_invalid_rq_factor(Integer r, Integer c) = 0;
        virtual const char* error_invalid_chol_update(Integer T_rows, Integer w_rows, Integer w_cols) = 0;
        virtual const char* error_invalid_qr_update(Integer T_rows, Integer T_cols, Integer u_rows, Integer u_cols, 
                                Integer v_rows, Integer v_cols) = 0;
        virtual const char* error_invalid_update_sigma(Integer K_vec, Integer sigma_rows, 
                                Integer sigma_cols) = 0;

        virtual const char* error_invalid_givens_seq_size(Integer C_r, Integer C_c, Integer S_r, Integer S_c,
                                Integer Ind_r, Integer Ind_c, Integer mat_pos) = 0;
        virtual const char* error_invalid_givens_seq_value_code(value_code vc, Integer mat_pos) = 0;
        virtual const char* error_invalid_givens_seq_index(Integer row_pos, Integer col_pos, Integer ind, 
                                Integer mat_size) = 0;
        virtual const char* invalid_holseholder_col(Integer col, Integer max_col) = 0;

        virtual const char* uninitialized_object_used(const std::string& obj_name) = 0; 
        virtual const char* unable_to_continue_arnoldi() = 0;
        virtual const char* invalid_speigs_k(Integer k, Integer N) = 0;

        virtual const char* invalid_lu_factors() = 0;
        virtual const char* invalid_cholecky_factors() = 0;
        virtual const char* invalid_svd_factors() = 0;
        virtual const char* invalid_schur_factors() = 0;
        virtual const char* invalid_qr_factors() = 0;
        virtual const char* invalid_ldl_factors() = 0;
        virtual const char* invalid_hess_factors() = 0;
        virtual const char* invalid_bidiag_factors() = 0;
        virtual const char* invalid_diagonal_22() = 0;

        virtual const char* invalid_length_of_permutation(Integer p_length, Integer N_req) = 0;
        virtual const char* band_matrix_required(Integer ld, Integer ud, Integer req_ld, Integer req_ud) = 0;
        virtual const char* diagonal_matrix_required(Integer ld, Integer ud) = 0;
        virtual const char* triangular_matrix_required(Integer ld, Integer ud) = 0;

        virtual const char* band_matrix_with_main_diag_required(Integer fd, Integer ld) = 0;

        virtual const char* iterative_solver_not_converged() = 0;
        virtual const char* logdet_not_available() = 0;
        virtual const char* unable_convert_linear_operator_to_matrix() = 0;

        virtual const char* invalid_partition(Integer part_length, Integer mat_size) = 0;
        virtual const char* invalid_partition_block(Integer block_size) = 0;
        virtual const char* invalid_partition_size(Integer req_size, Integer part_size) = 0;
        virtual const char* invalid_aggregation_matrix(Integer rows, Integer nnz) = 0;

        virtual const char* scotch_error(const std::string& str) = 0;        
        virtual const char* metis_error(Integer err) = 0;
        virtual const char* invalid_dm_block(Integer block, Integer max_blocks, Integer type) = 0;
        virtual const char* gsvd_failed() = 0;
        virtual const char* inaccurate_solution(Integer k, Real fe, Integer req_prec) = 0;

        virtual void        warning_singular() = 0;
        virtual void        warning_cholmod(const std::string& chol_msg) = 0;
        virtual void        warning_arpack(const std::string& chol_msg) = 0;
        virtual void        warning_inaccurate_solution(Integer k, Real fe, Integer req_prec) = 0;
};

class MATCL_LINALG_EXPORT default_exception_message_linalg : public exception_message_linalg
{
    private:
        safe_string_message current_message;

    public:
        virtual ~default_exception_message_linalg(){};

        virtual const char* error_general(const std::string& msg) override;		
        virtual const char* error_lu(const std::string& msg) override;	
        virtual const char* error_arpack(const std::string& msg) override;
        virtual const char* error_lsolve(Integer r1,Integer c1, Integer r2, Integer c2) override;
        virtual const char* error_lsolve_rev(Integer r1,Integer c1, Integer r2, Integer c2) override;
        virtual const char* error_fft(Integer status) override;
        virtual const char* error_singular() override;
        virtual const char* error_gen_sym_eigen_job() override;
        virtual const char* error_schur_incorrect_factors() override;
        virtual const char* error_schur_U_not_computed() override;
        virtual const char* error_gschur_QZ_not_computed() override;
        virtual const char* error_gschur_incorrect_factors() override;
        virtual const char* error_size_qrs(Integer r1,Integer c1) override;
        virtual const char* error_qrs(const Integer code) override;
        virtual const char* error_size_eig(Integer r1,Integer c1) override;
        virtual const char* error_size_geig(Integer r1,Integer c1, Integer r2, Integer c2) override;
        virtual const char* error_size_sylv(Integer r1,Integer c1, Integer r2, Integer c2, 
                                    Integer r3, Integer c3) override;
        virtual const char* error_size_gsylv(Integer rA1_, Integer cA1_, Integer rA2_, Integer cA2_,
                                    Integer rB1_, Integer cB1_, Integer rB2_, Integer cB2_,
                                    Integer rC1_, Integer cC1_, Integer rC2_, Integer cC2_) override;
        virtual const char* error_size_ricc() override;
        virtual const char* error_size_lyapunov() override;
        virtual const char* error_sylv() override;
        virtual const char* error_gsylv() override;
        virtual const char* error_ricc(Integer info) override;
        virtual const char* error_lyapunov() override;
        virtual const char* error_ldl() override;
        virtual const char* error_size_ldl(Integer r1,Integer c1) override;
        virtual const char* error_schur() override;
        virtual const char* error_gschur() override;
        virtual const char* error_norm_type() override;
        virtual const char* error_svd() override;
        virtual const char* error_svd_inf_nan() override;
        
        virtual const char* error_hess_nonsq() override;
        virtual const char* error_gen_hess(Integer Ar, Integer Ac, Integer Br, Integer Bc) override;
        virtual const char* hess_eig_failed(const Matrix& VL, const Matrix& VR, const Matrix& fail_L, 
                                    const Matrix& fail_R, bool eval_L, bool eval_R, Integer failed_eig) override;
        virtual const char* error_hess_eig_invalid_init(Integer I_rows, Integer I_cols, Integer H_size, 
                                    Integer E_size, bool left_init) override;
        virtual const char* error_hess_eig_invalid_conj_pair(Integer pos, bool missing_conj) override;
        virtual const char* error_hess_eig_invalid_eig(Integer E_rows, Integer E_cols, Integer H_size) override;
        virtual const char* error_hess_eig_her_complex_eig(Integer pos) override;

        virtual const char* error_schur_sel() override;
        virtual const char* error_schur_sel_comp() override;
        virtual const char* error_schur_cond_vec(Integer N, Integer M_sel, Integer N_vec, 
                                    Integer M_vec, bool left_vector) override;
        virtual const char* error_schur_eig_cluster_size(Integer M, Integer N) override;
        virtual const char* error_gschur_sel_comp() override;
        virtual const char* error_nonsymh_chol() override;
        virtual const char* error_nonposdef(bool allow_semipos) override;
        virtual const char* error_size_tridiag_sym(Integer A_rows, Integer A_cols, Integer B_rows, 
                                    Integer B_cols) override;
        virtual const char* error_size_tridiag_nsym(Integer Dm1_rows, Integer Dm1_cols,Integer D0_rows, 
                                    Integer D0_cols, Integer Dp1_rows, Integer Dp1_cols) override;
        virtual const char* error_symeigen_complex_diagonal() override;
        virtual const char* error_cholmod(const std::string& chol_msg) override;

        virtual const char* invalid_permutation_vector() override;

        virtual const char* error_complex_value_type_not_allowed() override;        
        virtual const char* symmetric_matrix_required() override;
        virtual const char* tridiagonal_matrix_required(Integer rows, Integer cols, Integer ld, Integer ud) override;
        virtual const char* error_upper_hess_matrix_required() override;
        virtual const char* invalid_unitary_matrix() override;
        virtual const char* error_invalid_chol_update_nonposdef() override;
        virtual const char* error_invalid_chol_factor(bool upper) override;
        virtual const char* error_invalid_qr_factor() override;
        virtual const char* error_invalid_rq_factor(Integer r, Integer c) override;
        virtual const char* error_invalid_chol_update(Integer T_rows, Integer w_rows, 
                                    Integer w_cols) override;
        virtual const char* error_invalid_qr_update(Integer T_rows, Integer T_cols, Integer u_rows, Integer u_cols, 
                                    Integer v_rows, Integer v_cols) override;
        virtual const char* error_invalid_update_sigma(Integer K_vec, Integer sigma_rows, 
                                    Integer sigma_cols) override;

        virtual const char* error_invalid_givens_seq_size(Integer C_r, Integer C_c, Integer S_r, Integer S_c,
                                    Integer Ind_r, Integer Ind_c, Integer mat_pos) override;
        virtual const char* error_invalid_givens_seq_value_code(value_code vc, Integer mat_pos) override;
        virtual const char* error_invalid_givens_seq_index(Integer row_pos, Integer col_pos, Integer ind, 
                                    Integer mat_size) override;
        virtual const char* invalid_holseholder_col(Integer col, Integer max_col) override;
        virtual const char* uninitialized_object_used(const std::string& obj_name) override; 
        virtual const char* unable_to_continue_arnoldi() override;
        virtual const char* invalid_speigs_k(Integer k, Integer N) override;

        virtual const char* invalid_lu_factors() override;
        virtual const char* invalid_cholecky_factors() override;
        virtual const char* invalid_svd_factors() override;
        virtual const char* invalid_schur_factors() override;
        virtual const char* invalid_qr_factors() override;
        virtual const char* invalid_ldl_factors() override;
        virtual const char* invalid_hess_factors() override;
        virtual const char* invalid_bidiag_factors() override;
        virtual const char* invalid_diagonal_22() override;

        virtual const char* invalid_length_of_permutation(Integer p_length, Integer N_req) override;
        virtual const char* band_matrix_required(Integer ld, Integer ud, Integer req_ld, Integer req_ud) override;
        virtual const char* diagonal_matrix_required(Integer ld, Integer ud) override;
        virtual const char* triangular_matrix_required(Integer ld, Integer ud) override;
        virtual const char* band_matrix_with_main_diag_required(Integer fd, Integer ld) override;

        virtual const char* iterative_solver_not_converged() override;
        virtual const char* logdet_not_available() override;
        virtual const char* unable_convert_linear_operator_to_matrix() override;
        virtual const char* invalid_partition(Integer part_length, Integer mat_size) override;
        virtual const char* invalid_partition_block(Integer block_size) override;
        virtual const char* invalid_partition_size(Integer req_size, Integer part_size) override;
        virtual const char* invalid_aggregation_matrix(Integer rows, Integer nnz) override;

        virtual const char* scotch_error(const std::string& str) override;
        virtual const char* metis_error(Integer err) override;
        virtual const char* invalid_dm_block(Integer block, Integer max_blocks, Integer type) override;
        virtual const char* gsvd_failed() override;
        virtual const char* inaccurate_solution(Integer k, Real fe, Integer req_prec) override;

        virtual void        warning(const std::string& msg);
        virtual void        warning_singular() override;
        virtual void        warning_cholmod(const std::string& chol_msg) override;
        virtual void        warning_arpack(const std::string& chol_msg) override;
        virtual void        warning_inaccurate_solution(Integer k, Real fe, Integer req_prec) override;
};

MATCL_LINALG_EXPORT void                            set_global_messanger_linalg(exception_message_linalg_ptr msg);
MATCL_LINALG_EXPORT exception_message_linalg_ptr    get_global_messanger_linalg();

class MATCL_LINALG_EXPORT matcl_exception_linalg : public matcl_exception
{
    public:
        virtual const char* what() const throw() override;

    private:
        virtual const char* what(exception_message& em) const override;
        virtual const char* what(exception_message_linalg& em) const = 0;
};

class MATCL_LINALG_EXPORT error_lu : public matcl_exception_linalg
{
    public:
        std::string msg;

    public:
        error_lu(const std::string& msg) : msg(msg) {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_arpack : public matcl_exception_linalg
{
    public:
        std::string msg;

    public:
        error_arpack(const std::string& msg) : msg(msg) {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_lsolve : public matcl_exception_linalg
{
    private:
        Integer r1, r2, c1, c2;

    public:
        error_lsolve(Integer r1_,Integer c1_, Integer r2_, Integer c2_) 
            : r1(r1_),r2(r2_),c1(c1_),c2(c2_) {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_lsolve_rev : public matcl_exception_linalg
{
    private:
        Integer r1, r2, c1, c2;

    public:
        error_lsolve_rev(Integer r1_,Integer c1_, Integer r2_, Integer c2_) 
            : r1(r1_),r2(r2_),c1(c1_),c2(c2_) {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_fft : public matcl_exception_linalg
{
    private:
        Integer status;

    public:
        error_fft(Integer status_) 
            : status(status_){};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_singular : public matcl_exception_linalg
{
    public:
        error_singular() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_schur_incorrect_factors : public matcl_exception_linalg
{
    public:
        error_schur_incorrect_factors() {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_schur_U_not_computed : public matcl_exception_linalg
{
    public:
        error_schur_U_not_computed() {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_gschur_QZ_not_computed : public matcl_exception_linalg
{
    public:
        error_gschur_QZ_not_computed() {};

        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_gen_sym_eigen_job : public matcl_exception_linalg
{
    private:

    public:
        error_gen_sym_eigen_job() 
            {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_gschur_incorrect_factors : public matcl_exception_linalg
{
    private:

    public:
        error_gschur_incorrect_factors() 
            {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_size_eig : public matcl_exception_linalg
{
    private:
        Integer r1, c1;

    public:
        error_size_eig(Integer A_r, Integer A_c) 
            : r1(A_r),c1(A_c) {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_size_qrs : public matcl_exception_linalg
{
    private:
        Integer r1, c1;

    public:
        error_size_qrs(Integer A_r, Integer A_c) 
            : r1(A_r),c1(A_c) {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_qrs : public matcl_exception_linalg
{
    private:
        const Integer code;

    public:
        error_qrs(const Integer code) 
            : code(code) {};

        virtual const char* what(exception_message_linalg& em) const;

	private:
		error_qrs& operator=(const error_qrs&) = delete;
};
class MATCL_LINALG_EXPORT error_size_ldl : public matcl_exception_linalg
{
    private:
        Integer r1, c1;

    public:
        error_size_ldl(Integer A_r, Integer A_c) 
            : r1(A_r),c1(A_c) {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_ldl : public matcl_exception_linalg
{
    public:
        error_ldl() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_size_geig : public matcl_exception_linalg
{
    private:
        Integer r1, r2, c1, c2;

    public:
        error_size_geig(Integer A_r, Integer A_c, Integer B_r, Integer B_c) 
            : r1(A_r),r2(B_r),c1(A_c),c2(B_c) {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_size_sylv : public matcl_exception_linalg
{
    private:
        Integer r1, r2, c1, c2, r3, c3;

    public:
        error_size_sylv(Integer A_r, Integer A_c, Integer B_r, Integer B_c, Integer C_r, Integer C_c) 
            : r1(A_r),r2(B_r),c1(A_c),c2(B_c),r3(C_r),c3(C_c) {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_size_gsylv : public matcl_exception_linalg
{
    private:
        Integer rA1, rA2, rB1, rB2, rC1, rC2;
        Integer cA1, cA2, cB1, cB2, cC1, cC2;

    public:
        error_size_gsylv(Integer rA1_, Integer cA1_, Integer rA2_, Integer cA2_,
                         Integer rB1_, Integer cB1_, Integer rB2_, Integer cB2_,
                         Integer rC1_, Integer cC1_, Integer rC2_, Integer cC2_) 
            : rA1(rA1_), cA1(cA1_), rA2(rA2_), cA2(cA2_),
              rB1(rB1_), cB1(cB1_), rB2(rB2_), cB2(cB2_),
              rC1(rC1_), cC1(cC1_), rC2(rC2_), cC2(cC2_)
        {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_size_ricc : public matcl_exception_linalg
{
    public:
        error_size_ricc() {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_size_lyapunov : public matcl_exception_linalg
{
    public:
        error_size_lyapunov() {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_sylv : public matcl_exception_linalg
{
    public:
        error_sylv() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_gsylv : public matcl_exception_linalg
{
    public:
        error_gsylv() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_ricc : public matcl_exception_linalg
{
    private:
        Integer info;
    public:
        error_ricc(Integer info_in)
            : info(info_in) {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_lyapunov : public matcl_exception_linalg
{
    private:
    public:
        error_lyapunov() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_schur : public matcl_exception_linalg
{
    public:
        error_schur() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_gschur : public matcl_exception_linalg
{
    public:
        error_gschur() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_norm_type : public matcl_exception_linalg
{
    public:
        error_norm_type() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_svd : public matcl_exception_linalg
{
    public:
        error_svd() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_svd_inf_nan : public matcl_exception_linalg
{
    public:
        error_svd_inf_nan() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_hess_nonsq : public matcl_exception_linalg
{
    public:
        error_hess_nonsq() {};
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_gen_hess : public matcl_exception_linalg
{
    private:
        Integer m_Ar;
        Integer m_Ac;
        Integer m_Br;
        Integer m_Bc;

    public:
        error_gen_hess(Integer A_r, Integer A_c, Integer B_r, Integer B_c)
            :m_Ar(A_r), m_Ac(A_c), m_Br(B_r), m_Bc(B_c)
        {};
        
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_hess_eig_invalid_eig : public matcl_exception_linalg
{
    private:
        Integer m_E_rows;
        Integer m_E_cols;
        Integer m_H_size;

    public:
        error_hess_eig_invalid_eig(Integer E_r, Integer E_c, Integer N)
            :m_E_rows(E_r), m_E_cols(E_c), m_H_size(N)
        {};
        
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_hess_eig_invalid_init : public matcl_exception_linalg
{
    private:
        Integer m_I_rows;
        Integer m_I_cols;
        Integer m_H_size;
        Integer m_E_size;
        bool    m_left_init;

    public:
        error_hess_eig_invalid_init(Integer I_r, Integer I_c, Integer H_size, Integer E_size, 
                                    bool left)
            :m_I_rows(I_r), m_I_cols(I_c), m_H_size(H_size), m_E_size(E_size), m_left_init(left)
        {};
        
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_hess_eig_invalid_conj_pair : public matcl_exception_linalg
{
    private:
        Integer m_pos;
        bool    m_missing_conj;

    public:
        error_hess_eig_invalid_conj_pair(Integer pos, bool missing_conj)
            :m_pos(pos), m_missing_conj(missing_conj)
        {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_hess_eig_her_complex_eig : public matcl_exception_linalg
{
    private:
        Integer m_pos;

    public:
        error_hess_eig_her_complex_eig(Integer pos)
            :m_pos(pos)
        {};
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_schur_sel : public matcl_exception_linalg
{
    public:
        error_schur_sel() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_schur_sel_comp : public matcl_exception_linalg
{
    public:
        error_schur_sel_comp() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_schur_cond_vec : public matcl_exception_linalg
{
    private:
        Integer m_exp_N;
        Integer m_exp_M;
        Integer m_sup_N;
        Integer m_sup_M;
        bool    m_left_vector;

    public:
        error_schur_cond_vec(Integer N, Integer M_sel, Integer N_vec, Integer M_vec, bool left_vector)
            :m_exp_N(N), m_exp_M(M_sel), m_sup_N(N_vec), m_sup_M(M_vec), m_left_vector(left_vector)
        {};

        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_schur_eig_cluster_size : public matcl_exception_linalg
{
    private:
        Integer m_cluster;
        Integer m_size;

    public:
        error_schur_eig_cluster_size(Integer M, Integer N)
            :m_cluster(M), m_size(N)
        {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_gschur_sel_comp : public matcl_exception_linalg
{
    public:
        error_gschur_sel_comp() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_size_tridiag_sym : public matcl_exception_linalg
{
    private:
        Integer A_rows, A_cols, B_rows, B_cols;

    public:
        error_size_tridiag_sym(Integer A_rows_, Integer A_cols_, Integer B_rows_, 
                               Integer B_cols_)
            :A_rows(A_rows_), B_rows(B_rows_), A_cols(A_cols_), B_cols(B_cols_)
        {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_size_tridiag_nsym : public matcl_exception_linalg
{
    private:
        Integer D0_rows, D0_cols;
        Integer D_m1_rows, D_m1_cols;
        Integer D_p1_rows, D_p1_cols;

    public:
        error_size_tridiag_nsym(Integer Dm1_rows_, Integer Dm1_cols_, Integer D0_rows_, Integer D0_cols_, 
                            Integer Dp1_rows_, Integer Dp1_cols_)
            :D0_rows(D0_rows_),D0_cols(D0_cols_),D_m1_rows(Dm1_rows_),D_m1_cols(Dm1_cols_)
            ,D_p1_rows(Dp1_rows_),D_p1_cols(Dp1_cols_)
        {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_symeigen_complex_diagonal : public matcl_exception_linalg
{
    public:
        error_symeigen_complex_diagonal(){};
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_nonsymh_chol : public matcl_exception_linalg
{
    public:
        error_nonsymh_chol() {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_nonposdef : public matcl_exception_linalg
{
    private:
        bool    m_allow_semipos;

    public:
        error_nonposdef(bool allow_semipos) 
            :m_allow_semipos(allow_semipos)
        {};

        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_cholmod : public matcl_exception_linalg
{
    public:
        error_cholmod(const std::string& _chol_msg): chol_msg(_chol_msg) {};

        virtual const char* what(exception_message_linalg& em) const override;

    private:
        std::string chol_msg;
};

class MATCL_LINALG_EXPORT symmetric_matrix_required: public matcl_exception_linalg
{
    public:
        virtual const char* what(exception_message_linalg& em) const override;
};

class MATCL_LINALG_EXPORT tridiagonal_matrix_required: public matcl_exception_linalg
{
    private:
        Integer m_rows;
        Integer m_cols;
        Integer m_ld;
        Integer m_ud;

    public:
        tridiagonal_matrix_required(Integer r, Integer c, Integer ld, Integer ud)
            :m_rows(r), m_cols(c), m_ld(ld), m_ud(ud)
        {};
        virtual const char* what(exception_message_linalg& em) const override;
};

class MATCL_LINALG_EXPORT diagonal_matrix_required: public matcl_exception_linalg
{
    private:
        Integer m_ld;
        Integer m_ud;

    public:
        diagonal_matrix_required(Integer ld, Integer ud)
            :m_ld(ld), m_ud(ud)
        {};
        virtual const char* what(exception_message_linalg& em) const override;
};
class MATCL_LINALG_EXPORT triangular_matrix_required: public matcl_exception_linalg
{
    private:
        Integer m_ld;
        Integer m_ud;

    public:
        triangular_matrix_required(Integer ld, Integer ud)
            :m_ld(ld), m_ud(ud)
        {};
        virtual const char* what(exception_message_linalg& em) const override;
};

class MATCL_LINALG_EXPORT band_matrix_required: public matcl_exception_linalg
{
    private:
        Integer m_ld;
        Integer m_ud;
        Integer m_req_ld;   //-1 if any
        Integer m_req_ud;   //-1 if any

    public:
        band_matrix_required(Integer ld, Integer ud, Integer req_ld, Integer req_ud)
            :m_ld(ld), m_ud(ud), m_req_ld(req_ld), m_req_ud(req_ud)
        {};
        virtual const char* what(exception_message_linalg& em) const override;
};

class MATCL_LINALG_EXPORT invalid_length_of_permutation: public matcl_exception_linalg
{
    private:
        Integer m_p_length;
        Integer m_N_req;

    public:
        invalid_length_of_permutation(Integer pn, Integer N)
            :m_p_length(pn), m_N_req(N)
        {};
        virtual const char* what(exception_message_linalg& em) const override;
};

class MATCL_LINALG_EXPORT error_upper_hess_matrix_required: public matcl_exception_linalg
{
    public:
        error_upper_hess_matrix_required(){};
        virtual const char* what(exception_message_linalg& em) const override;
};

class MATCL_LINALG_EXPORT error_complex_value_type_not_allowed : public matcl_exception_linalg
{
    public:
        error_complex_value_type_not_allowed() {};

        virtual const char* what(exception_message_linalg& em) const override;
};

class MATCL_LINALG_EXPORT invalid_unitary_matrix : public matcl_exception_linalg
{
    public:
        invalid_unitary_matrix(){};

        virtual const char* what(exception_message_linalg& em) const;
};


class MATCL_LINALG_EXPORT error_invalid_chol_factor : public matcl_exception_linalg
{
    private:
        bool    m_upper;

    public:
        error_invalid_chol_factor(bool upper)
            :m_upper(upper)
        {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_invalid_qr_factor : public matcl_exception_linalg
{
    public:
        error_invalid_qr_factor()
        {};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_invalid_rq_factor : public matcl_exception_linalg
{
    private:
        Integer m_rows;
        Integer m_cols;

    public:
        error_invalid_rq_factor(Integer r, Integer c)
            :m_rows(r), m_cols(c)
        {};

        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_invalid_chol_update : public matcl_exception_linalg
{
    private:
        Integer     m_T_rows;
        Integer     m_w_rows;
        Integer     m_w_cols;

    public:
        error_invalid_chol_update(Integer T_rows, Integer w_rows, Integer w_cols)
            :m_T_rows(T_rows), m_w_rows(w_rows), m_w_cols(w_cols)
        {};
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT error_invalid_update_sigma : public matcl_exception_linalg
{
    private:
        Integer     m_K_vector;
        Integer     m_sigma_r;
        Integer     m_sigma_c;

    public:
        error_invalid_update_sigma(Integer K_vec, Integer sigma_r, Integer sigma_c)
            :m_K_vector(K_vec), m_sigma_r(sigma_r), m_sigma_c(sigma_c)
        {};

        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_invalid_chol_update_nonposdef : public matcl_exception_linalg
{
    public:
        error_invalid_chol_update_nonposdef(){};
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_invalid_qr_update : public matcl_exception_linalg
{
    private:
        Integer     m_T_rows;
        Integer     m_T_cols;
        Integer     m_u_rows;
        Integer     m_u_cols;
        Integer     m_v_rows;
        Integer     m_v_cols;

    public:
        error_invalid_qr_update(Integer T_rows, Integer T_cols, Integer u_rows, Integer u_cols, 
                                Integer v_rows, Integer v_cols)
            :m_T_rows(T_rows), m_T_cols(T_cols), m_u_rows(u_rows), m_u_cols(u_cols), m_v_rows(v_rows), 
            m_v_cols(v_cols)
        {};
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_invalid_givens_seq_size : public matcl_exception_linalg
{
    private:
        Integer m_C_r;
        Integer m_C_c;
        Integer m_S_r;
        Integer m_S_c;
        Integer m_Ind_r;
        Integer m_Ind_c;
        Integer m_invalid_mat_pos;

    public:
        error_invalid_givens_seq_size(Integer C_r, Integer C_c, Integer S_r, Integer S_c,
                                      Integer Ind_r, Integer Ind_c, Integer mat_pos)
            :m_C_r(C_r), m_C_c(C_c), m_S_r(S_r), m_S_c(S_c), m_Ind_r(Ind_r), m_Ind_c(Ind_c)
            ,m_invalid_mat_pos(mat_pos)
        {};

        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_invalid_givens_seq_value_code : public matcl_exception_linalg
{
    private:
        value_code  m_value_code;
        Integer     m_mat_pos;

    public:
        error_invalid_givens_seq_value_code(value_code vc, Integer mat_pos)
            :m_value_code(vc), m_mat_pos(mat_pos)
        {};

        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT error_invalid_givens_seq_index : public matcl_exception_linalg
{
    private:
        Integer m_pos_row;
        Integer m_pos_col;
        Integer m_index;
        Integer m_matrix_size;

    public:
        error_invalid_givens_seq_index(Integer row_pos, Integer col_pos, Integer ind, Integer mat_size)
            :m_pos_row(row_pos), m_pos_col(col_pos), m_index(ind), m_matrix_size(mat_size)
        {};

        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT invalid_holseholder_col : public matcl_exception_linalg
{
    private:
        Integer m_col;
        Integer m_max_col;

    public:
        invalid_holseholder_col(Integer col, Integer max_col)
            :m_col(col), m_max_col(max_col)
        {};

        virtual const char* what(exception_message_linalg& em) const;

};

class MATCL_LINALG_EXPORT uninitialized_object_used : public matcl_exception_linalg
{
    private:
        std::string     m_object;

    public:
        uninitialized_object_used(const std::string& obj)
            :m_object(obj)
        {};

        virtual const char* what(exception_message_linalg& em) const;

};

class MATCL_LINALG_EXPORT unable_to_continue_arnoldi : public matcl_exception_linalg
{
    public:
        unable_to_continue_arnoldi(){};

        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT invalid_speigs_k : public matcl_exception_linalg
{
    private:
        Integer     m_k;
        Integer     m_N;

    public:
        invalid_speigs_k(Integer k, Integer N)
            :m_k(k), m_N(N)
        {};

        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT invalid_lu_factors : public matcl_exception_linalg
{
    public:
        invalid_lu_factors(){};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT invalid_cholecky_factors : public matcl_exception_linalg
{
    public:
        invalid_cholecky_factors(){};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT invalid_svd_factors : public matcl_exception_linalg
{
    public:
        invalid_svd_factors(){};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT invalid_schur_factors : public matcl_exception_linalg
{
    public:
        invalid_schur_factors(){};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT invalid_qr_factors : public matcl_exception_linalg
{
    public:
        invalid_qr_factors(){};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT invalid_ldl_factors : public matcl_exception_linalg
{
    public:
        invalid_ldl_factors(){};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT invalid_hess_factors : public matcl_exception_linalg
{
    public:
        invalid_hess_factors(){};

        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT invalid_bidiag_factors : public matcl_exception_linalg
{
    public:
        invalid_bidiag_factors(){};

        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT invalid_diagonal_22 : public matcl_exception_linalg
{
    public:
        invalid_diagonal_22(){};

        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT band_matrix_with_main_diag_required
    : public matcl_exception_linalg
{
    private:
        Integer m_first_diag;
        Integer m_last_diag;

    public:
        band_matrix_with_main_diag_required(Integer fd, Integer ld)
            :m_first_diag(fd), m_last_diag(ld)
        {};

        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT iterative_solver_not_converged : public matcl_exception_linalg
{
    public:
        iterative_solver_not_converged(){};

    public:
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT logdet_not_available : public matcl_exception_linalg
{
    public:
        logdet_not_available(){};

    public:
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT unable_convert_linear_operator_to_matrix : public matcl_exception_linalg
{
    public:
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT invalid_partition : public matcl_exception_linalg
{
    private:
        Integer     m_part_length;
        Integer     m_mat_size;

    public:
        invalid_partition(Integer part_length, Integer mat_size)
            :m_part_length(part_length), m_mat_size(mat_size)
        {};
        
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT invalid_partition_block : public matcl_exception_linalg
{
    private:
        Integer     m_block;

    public:
        invalid_partition_block(Integer block)
            :m_block(block)
        {};
        
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT invalid_partition_size : public matcl_exception_linalg
{
    private:
        Integer     m_req_size;
        Integer     m_part_size;

    public:
        invalid_partition_size(Integer req_size, Integer part_size)
            :m_req_size(req_size), m_part_size(part_size)
        {};
        
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT scotch_error : public matcl_exception_linalg
{
    private:
        std::string m_msg;

    public:
        scotch_error(const std::string& msg)
            :m_msg(msg)
        {};
    public:
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT metis_error : public matcl_exception_linalg
{
    private:
        Integer m_err;

    public:
        metis_error(Integer err)
            :m_err(err)
        {};

    public:
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT invalid_aggregation_matrix : public matcl_exception_linalg
{
    private:
        Integer m_rows;
        Integer m_nnz;

    public:
        invalid_aggregation_matrix(Integer rows, Integer nnz)
            :m_rows(rows), m_nnz(nnz)
        {};

    public:
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT invalid_dm_block : public matcl_exception_linalg
{
    private:
        Integer m_block;
        Integer m_max_blocks;
        Integer m_type;

    public:
        invalid_dm_block(Integer block, Integer max_blocks, Integer type)
            :m_block(block), m_max_blocks(max_blocks), m_type(type)
        {};

    public:
        virtual const char* what(exception_message_linalg& em) const;
};
class MATCL_LINALG_EXPORT gsvd_failed : public matcl_exception_linalg
{
    public:
        gsvd_failed(){};

    public:
        virtual const char* what(exception_message_linalg& em) const;
};

class MATCL_LINALG_EXPORT inaccurate_solution : public matcl_exception_linalg
{
    private:
        Integer m_vec;
        Real    m_fe;
        Integer m_req_prec;

    public:
        inaccurate_solution(Integer k, Real fe, Integer req_prec)
            :m_vec(k), m_fe(fe), m_req_prec(req_prec)
        {};

    public:
        virtual const char* what(exception_message_linalg& em) const;
};

inline void check_lsolve(Integer r1, Integer c1, Integer r2, Integer c2)
{
    if (r1 != c1 || c1 != r2)
    {
        throw error_lsolve(r1, c1, r2, c2);
    };
};
inline void check_lsolve_rev(Integer r1, Integer c1, Integer r2, Integer c2)
{
    if ((r1 != c1) || (c1 != c2))
    {
        throw error_lsolve_rev(r1, c1, r2, c2);
    };
};

};};

#pragma warning(pop)