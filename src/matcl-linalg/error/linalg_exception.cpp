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

#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-core/IO/output_stream.h"
#include <iostream>

#include <algorithm>

namespace matcl { namespace error
{

static exception_message_linalg_ptr messanger(new default_exception_message_linalg());

const char* default_exception_message_linalg::error_general(const std::string& msg)
{
    current_message = msg;
    return current_message.c_str();
};

void default_exception_message_linalg::warning(const std::string& msg)
{
    std::stringstream of;

    of << "warning: " << msg << "\n";
    global_output_stream()->disp(of);
};

void default_exception_message_linalg::warning_cholmod(const std::string& chol_msg)
{
    warning("CHOLMOD: " + chol_msg);
};

void default_exception_message_linalg::warning_arpack(const std::string& chol_msg)
{
    warning("ARPACK: " + chol_msg);
};

const char* default_exception_message_linalg::error_cholmod(const std::string& chol_msg)
{
    current_message = "CHOLMOD: " + chol_msg;
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_lu(const std::string& msg)
{
    current_message = std::string() + "error during lu factorization: " + msg;
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_arpack(const std::string& msg)
{
    current_message = std::string() + "error during arpack excecution: " + msg;
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_lsolve(Integer r1,Integer c1, Integer r2, Integer c2)
{
    std::ostringstream os;

    os  << "nonsquare matrix or nonconformat sizes in `linsolve';\n"
        << "first operand: " << r1 << "x" << c1
        << ", second operand: "  << r2 << "x" << c2;

    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_lsolve_rev(Integer r1,Integer c1, Integer r2, Integer c2)
{
    std::ostringstream os;

    os  << "nonsquare matrix or nonconformat sizes in `linsolve_rev';\n"
        << "first operand: " << r1 << "x" << c1
        << ", second operand: "  << r2 << "x" << c2;

    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_fft(Integer status)
{
    std::ostringstream os;

    os  << "DFTI returned status" << status;

    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_singular()
{
    current_message = "matrix is singular";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_schur_incorrect_factors()
{
    std::ostringstream os;

    os  << "not a Schur factorization ";
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_schur_U_not_computed()
{
    std::ostringstream os;

    os  << "unitary matrix of schur decomposition is not available";
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_gschur_QZ_not_computed()
{
    std::ostringstream os;

    os  << "unitary matrices of generalized schur decomposition are not available";
    current_message = os.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_gen_sym_eigen_job()
{
    std::ostringstream os;

    os  << "Only job=1, job=2 or job=3 allowed";
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_gschur_incorrect_factors()
{
    std::ostringstream os;

    os  << " not a generalized Schur factorization ";
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_qrs(const Integer code)
{
    std::ostringstream os;

    os  << "Quern failed, returned code "<< code;
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_size_eig(Integer r1,Integer c1)
{
    std::ostringstream os;

    os  << "nonsquare matrix in eigenproblem;\n"
        << "first operand: " << r1 << "x" << c1;
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_size_qrs(Integer r1,Integer c1)
{
    std::ostringstream os;

    os  << "m < n required in qrs';\n"
        << "first operand: m = " << r1 << ", n = " << c1;
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_size_ldl(Integer r1,Integer c1)
{
    std::ostringstream os;

    os  << "nonsquare matrix in Bunch-Kaufman LDL';\n"
        << "first operand: " << r1 << "x" << c1;
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_ldl()
{
    current_message = "SYTRF/HETRF returne info < 0";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_size_geig(Integer r1,Integer c1, Integer r2, Integer c2)
{
    std::ostringstream os;

    os  << "nonsquare matrix or nonconformat sizes in generalized eigenproblem;\n"
        << "first operand: " << r1 << "x" << c1
        << ", second operand: "  << r2 << "x" << c2;
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_size_sylv(
                Integer r1,Integer c1, Integer r2, Integer c2, Integer r3, Integer c3)
{
    std::ostringstream os;

    os  << "nonsquare matrices or nonconformat sizes in Sylvester equation;\n"
        << "first operand: " << r1 << "x" << c1
        << ", second operand: "  << r2 << "x" << c2
        << ", RHS : "  << r3 << "x" << c3;
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_size_gsylv(Integer rA1, Integer cA1, Integer rA2, Integer cA2,
                                    Integer rB1, Integer cB1, Integer rB2, Integer cB2,
                                    Integer rC1, Integer cC1, Integer rC2, Integer cC2)
{
    std::ostringstream os;

    os  << "nonsquare matrices or nonconformat sizes in generalized Sylvester equation;\n"
        << "first equation: (" << rA1 << "x" << cA1 << ") * X + Y * (" << rA2 << "x" << cA2 << ") = ("  
                            << rC1 << "x" << cC1 << ")\n"
        << "second equation: (" << rB1 << "x" << cB1 << ") * X + Y * (" << rB2 << "x" << cB2 << ") = ("  
                            << rC2 << "x" << cC2 << ")\n";

    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_size_ricc()
{
    std::ostringstream os;

    os  << "nonsquare matrices or nonconformat sizes in algebraic Riccati equation";
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_size_lyapunov()
{
    std::ostringstream os;

    os  << "nonsquare matrices or nonconformat sizes in algebraic Lypunov equation";
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_schur()
{
    current_message = "Schur decomposition routine did not converge";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_sylv()
{
    current_message = "Sylvester equation routine did not converge (common eigenvalues?)";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_gsylv()
{
    current_message = "generalized Sylvester equation routine did not converge (common eigenvalues?)";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_ricc(Integer info)
{
    std::ostringstream os;
    os << "Riccati equation solver failed: ";
    switch(info)
    {
        case -1:
            os << "RCOND close to zero. RCOND is an estimate of the reciprocal of the condition number (in ";
            os << "the 1-norm) of the N-th order system of algebraic ";
            os << "equations from which the solution matrix X is obtained.";
            break;
        case -2:
            os << "Generalized spectrum too close to unit circle";
            break;
        // SLICOT messages
        case 1:
            os  << "the computed extended matrix pencil is singular, possibly due to rounding errors";
            break;
        case 2:
            os  << "the QZ (or QR) algorithm failed";
            break;
        case 3:
            os  << "reordering of the (generalized) eigenvalues failed";
            break;
        case 4:
            os  << "after reordering, roundoff changed values of " <<
                "some complex eigenvalues so that leading eigenvalues " <<
                "in the (generalized) Schur form no longer satisfy " <<
                "the stability condition; this could also be caused " <<
                "due to scaling";
            break;
        case 5:
            os  << "the computed dimension of the solution does not equal N;";
            break;
        case 6:
            os  << "a singular matrix was encountered during the computation of the solution matrix X";
            break;
        default:
            os  << "unknown error";
            break;
    }
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_lyapunov()
{
    std::ostringstream os;
    os << "Lyapunov equation solver failed";
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_gschur()
{
    current_message = "generalized Schur decomposition routine did not converge";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_norm_type()
{
    std::ostringstream os;
    os << "invalid matrix norm; admissible values are 1, 2, Inf, -1 (Frobenius), and -2 (max abs)";
    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_svd()
{
    current_message = "SVD did not converge";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_svd_inf_nan()
{
    current_message = "Input to SVD must not contain nan or Inf";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_hess_nonsq()
{
    current_message = "Input to hess must be square";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_gen_hess(Integer Ar, Integer Ac, Integer Br, Integer Bc)
{
    std::ostringstream msg;
    msg << "invalid matrices passed to generalized hessenberg routines; expecting square matrices of the same size;"
        << " first matrix size: " << Ar << "x" << Ac << ", second matrix size: " << Br << "x" << Bc;

    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_hess_eig_invalid_eig(Integer E_rows, Integer E_cols, 
                                                                         Integer H_size)
{
    std::ostringstream msg;
    msg << "invalid eigenvalue matrix passed to hess_eig; expecting vector of length at most " << H_size
        << "; supplied matrix has size " << E_rows << "x" << E_cols;

    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_hess_eig_her_complex_eig(Integer pos)
{
    std::ostringstream msg;
    msg << "invalid eigenvalue matrix passed to hess_eig for hermitian matrix; expecting real eigenvalues only;"
        << "; supplied matrix has complex eigenvalue at position " << pos;

    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::hess_eig_failed(const Matrix& VL, const Matrix& VR, 
                    const Matrix& fail_L, const Matrix& fail_R, bool eval_L, bool eval_R, Integer failed_eig)
{
    (void)VL;
    (void)VR;
    (void)fail_L;
    (void)fail_R;
    (void)eval_L;
    (void)eval_R;

    std::ostringstream msg;
    msg << "inverse iteration in hess_eig did not converge; there are " << failed_eig << " nonconverged eigenvectors;"
        << " additional information are stored in exception class 'hess_eig_failed'";

    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_hess_eig_invalid_init(Integer I_rows, Integer I_cols, 
                    Integer H_size, Integer E_size, bool left_init)
{
    std::ostringstream msg;
    msg << "invalid initial " << (left_init? "left" : "right") << " vectors are supplied to hess_eig;"
        <<" expecting matrix of size " << H_size << "x" << E_size
        <<" or empty matrix; supplied matrix has size " << I_rows << "x" << I_cols;

    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_hess_eig_invalid_conj_pair(Integer pos, bool missing_conj)
{
    std::ostringstream msg;
    msg << "invalid eigenvalue vector passed to hess_eig;";
    if (missing_conj)
        msg << " invalid complex conjugate pair at last position";
    else
        msg << " invalid complex conjugate pair at position " << pos;

    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_schur_sel()
{
    current_message = "invalid eigenvalue selector in (g)schur, size of the selector vector must be equal to ";
    current_message+= "size of matrices, and must contain 0 or 1 values";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_schur_sel_comp()
{
    current_message = "eigenvalues reordering in Schur decomposition routine failed";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_schur_cond_vec(Integer N, Integer M_sel, Integer N_vec, 
                                            Integer M_vec, bool left_vector)
{
    std::ostringstream msg;
    msg << "invalid " << (left_vector? "left" : "right") << " eigenvector matrix; expecting matrix of size: "
        << N << "x" << M_sel << "; supplied matrix has size: " << N_vec << "x" << M_vec;

    current_message = msg.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_schur_eig_cluster_size(Integer M, Integer N)
{
    std::ostringstream msg;
    msg << "invalid eigenvalue cluster size; expecting value in range [0," << N <<"]; supplied value is " << M; 

    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_gschur_sel_comp()
{
    current_message = "eigenvalues reordering in generalized Schur decomposition routine failed";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_nonsymh_chol()
{
    current_message = "input to chol must be real symmetric or complex hermitian";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_nonposdef(bool allow_semi)
{
    if (allow_semi)
        current_message = "matrix is not semi positive definite";
    else
        current_message = "matrix is not positive definite";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_size_tridiag_sym(Integer A_rows, Integer A_cols, 
                                            Integer B_rows, Integer B_cols)
{
    std::ostringstream os;
    os << "invalid diagonals of hermitian tridiagonal matrix, expecting two vectors of elements "
          "on diagonal and subdiagonal of a matrix, supplied matrices have size: "
       << A_rows << "x" << A_cols << ", and " << B_rows << "x" << B_cols;
    current_message = os.str();
    return current_message.c_str();
}
const char* default_exception_message_linalg::error_size_tridiag_nsym(Integer Dm1_rows, Integer Dm1_cols, 
                Integer D0_rows, Integer D0_cols, Integer Dp1_rows, Integer Dp1_cols)
{
    std::ostringstream os;
    os << "invalid diagonals of tridiagonal matrix, expecting three vectors of elements "
          "on diagonal and off-diagonals of a matrix, supplied matrices have size: "
       << Dm1_rows << "x" << Dm1_cols << ", " << D0_rows << "x" << D0_cols << ", and "
       << Dp1_rows << "x" << Dp1_cols;

    current_message = os.str();
    return current_message.c_str();
}

const char* default_exception_message_linalg::error_symeigen_complex_diagonal()
{
    current_message = "invalid matrix passed to symmetric eigenvalue problem, diagonal contains complex entries";
    return current_message.c_str();
};

const char* default_exception_message_linalg::invalid_permutation_vector()
{
    current_message = "invalid permutation vector";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_complex_value_type_not_allowed()
{
    current_message = "complex value type is not allowed";
    return current_message.c_str();
};
const char* default_exception_message_linalg::symmetric_matrix_required()
{
    std::ostringstream msg;
    msg << "symmetric matrix is required";
    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::tridiagonal_matrix_required(Integer rows, Integer cols,
                                                            Integer ld, Integer ud)
{
    std::ostringstream msg;
    msg << "tridiagonal matrix is required, matrix size is: " << rows << "x" << cols
        << ", with " << ld << " subdiagonals and " << ud << "superdiagonals";
    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_upper_hess_matrix_required()
{
    current_message = std::string() + "upper hessenberg matrix is required";
    return current_message.c_str();
};

const char* default_exception_message_linalg::invalid_unitary_matrix()
{
    current_message = "matrix is not a valid unitary matrix";
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_invalid_chol_update_nonposdef()
{
    current_message = "unable to make rank-k update of cholesky factor: resulting matrix is not"
                      " positive definite";
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_invalid_chol_factor(bool upper)
{
    std::ostringstream os;

    os << "input matrix must be square ";
    if (upper)
        os << "upper";
    else
        os << "lower";

    os << " triangular";

    current_message = os.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_invalid_qr_factor()
{
    std::ostringstream os;

    os << "input matrix must be upper triangular";

    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_invalid_rq_factor(Integer r, Integer c)
{
    std::ostringstream os;

    os << "input matrix must be square upper triangular, matrix size is: " << r << "x" << c;

    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::error_invalid_qr_update(Integer T_rows, Integer T_cols, 
        Integer u_rows, Integer u_cols, Integer v_rows, Integer v_cols)
{
    std::ostringstream os;

    os << "unable to perform rank-k update of qr factor: "
            "nonconformant size of qr factor and update vectors, "
        << "qr factor size: " << T_rows << "x" << T_cols << ", update vectors size are " << u_rows << "x" << u_cols
        << " and " << v_rows << "x" << v_cols;

    current_message = os.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_invalid_chol_update(Integer T_rows, Integer w_rows, 
                Integer w_cols)
{
    std::ostringstream os;

    os << "unable to perform rank-k update of cholesky factor: "
            "nonconformant size of cholesky factor and update vectors, "
        << "cholesky factor size: " << T_rows << ", update vectors size is " << w_rows << "x" << w_cols;

    current_message = os.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_invalid_update_sigma(Integer K_vec, 
            Integer sigma_r, Integer sigma_c)
{
    std::ostringstream os;

    os << "unable to perform rank-k update: "
            "invalid diagonal scaling, expecting scalar or vector of size " << K_vec
        << "; scaling size is: " << sigma_r << "x" << sigma_c;

    current_message = os.str();
    return current_message.c_str();
}

const char* default_exception_message_linalg::error_invalid_givens_seq_value_code(value_code vc, Integer mat_pos)
{
    std::ostringstream os;
    os << "invalid plane rotations sequence";
    
    if (mat_pos == 1)
    {
        os  << "; cosine array stores value of invalid type: expecting real scalars, matrix stores "
            << default_exception_message::val_to_string(vc) << " scalars";
    }
    else if (mat_pos == 2)
    {
        os  << "; sine array stores value of invalid type: expecting real or complex scalars, matrix stores "
            << default_exception_message::val_to_string(vc) << " scalars";
    }
    else
    {
        os  << "; index array stores value of invalid type: expecting integer scalars, matrix stores "
            << default_exception_message::val_to_string(vc) << " scalars";
    }

    current_message = os.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_invalid_givens_seq_index(Integer row_pos, Integer col_pos, 
                Integer ind, Integer mat_size)
{
    std::ostringstream os;
    os  << "invalid plane rotations sequence; index array contains invalid element at position [" << row_pos << ", " 
        << col_pos << "]" << "; expecting value in range [" << 1 << ", " << mat_size << "]; value is: "
        << ind;

    current_message = os.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::invalid_holseholder_col(Integer col, Integer max_col)
{
    std::ostringstream os;
    os  << "invalid number of unitary matrix column constructed from elementary reflectors; expecting value in range: "
        << 1 << "-" << max_col << "; supplied value : " << col;

    current_message = os.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::uninitialized_object_used(const std::string& obj_name)
{
    std::ostringstream os;
    os  << "uninitialized object used: " << obj_name;

    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message_linalg::unable_to_continue_arnoldi()
{
    std::ostringstream os;
    os  << "unable to continue Arnoldi iterations; Arnoldi iterations not started";

    current_message = os.str();
    return current_message.c_str();
};

const char* default_exception_message_linalg::error_invalid_givens_seq_size(Integer C_r, Integer C_c, 
                Integer S_r, Integer S_c, Integer Ind_r, Integer Ind_c, Integer mat_pos)
{
    std::ostringstream os;
    os  << "invalid plane rotations sequence";

    if (mat_pos == 1)
    {
        os << "; cosine array has invalid size; expecting vector, matrix size is: " << C_r << "x" << C_c;
    }
    else if (mat_pos == 2)
    {
        Integer K   = (C_r == 0 || C_c == 0) ? 0 : std::max(C_r, C_c);

        os  << "; sine array has invalid size; expecting vector of length " << K 
            << ", matrix size is: " << S_r << "x" << S_c;
    }
    else
    {
        Integer K   = (C_r == 0 || C_c == 0) ? 0 : std::max(C_r, C_c);

        os  << "; index array has invalid size; expecting matrix of size " << K << " x 2"
            << ", matrix size is: " << Ind_r << "x" << Ind_c;
    };

    current_message = os.str();
    return current_message.c_str(); 
};

const char* default_exception_message_linalg::invalid_speigs_k(Integer k, Integer N)
{
    std::ostringstream os;

    if (N <= 0)
        os  << "problem is too small for speigs, use direct method insted";
    else
    {
        os << "invalid number of eigenvalues requested; expecting number in range: " << 1 << "-" << N
           << "; requested number of eigenvalues is " << k;
    };

    current_message = os.str();
    return current_message.c_str(); 
};

const char* default_exception_message_linalg::invalid_lu_factors()
{
    current_message = "invalid lu factors";
    return current_message.c_str();
};
const char* default_exception_message_linalg::invalid_cholecky_factors()
{
    current_message = "invalid Cholesky factors";
    return current_message.c_str();
};
const char* default_exception_message_linalg::invalid_svd_factors()
{
    current_message = "invalid svd factors";
    return current_message.c_str();
};
const char* default_exception_message_linalg::invalid_schur_factors()
{
    current_message = "invalid Schur factors";
    return current_message.c_str();
};
const char* default_exception_message_linalg::invalid_qr_factors()
{
    current_message = "invalid qr factors";
    return current_message.c_str();
};
const char* default_exception_message_linalg::invalid_ldl_factors()
{
    current_message = "invalid ldl factors";
    return current_message.c_str();
};
const char* default_exception_message_linalg::invalid_hess_factors()
{
    current_message = "invalid hess factors";
    return current_message.c_str();
};
const char* default_exception_message_linalg::invalid_bidiag_factors()
{
    current_message = "invalid bidiagonal factors";
    return current_message.c_str();
};
const char* default_exception_message_linalg::invalid_diagonal_22()
{
    current_message = "invalid block diagonal matrix with blocks of size 1x1 or 2x2";
    return current_message.c_str();
};

const char* default_exception_message_linalg::diagonal_matrix_required(Integer ld, Integer ud)
{
    std::ostringstream os;
    
    os  << "diagonal matrix is required; supplied muatrix has " << ld << " subdiagonals" << " and " << ud
        << " superdiagonals";

    current_message = os.str();
    return current_message.c_str(); 
}
const char* default_exception_message_linalg::triangular_matrix_required(Integer ld, Integer ud)
{
    std::ostringstream os;
    
    os  << "upper or lower traingular matrix is required; supplied muatrix has " 
        << ld << " subdiagonals" << " and " << ud << " superdiagonals";

    current_message = os.str();
    return current_message.c_str(); 
};

const char* default_exception_message_linalg::band_matrix_with_main_diag_required(Integer fd, Integer ld)
{
    std::ostringstream os;
    
    os  << "band matrix with main diagonal is required; matrix first diagonal index: " << fd
        << ", last diagonal index: " << ld;

    current_message = os.str();
    return current_message.c_str(); 
};

const char* default_exception_message_linalg::iterative_solver_not_converged()
{
    current_message = "iterative solver did not converge";
    return current_message.c_str(); 
};
const char* default_exception_message_linalg::logdet_not_available()
{
    current_message = "log_det is not available";
    return current_message.c_str(); 
};
const char* default_exception_message_linalg::unable_convert_linear_operator_to_matrix()
{
    current_message = "linear operator cannot be converted to matrix";
    return current_message.c_str(); 
};

const char* default_exception_message_linalg::invalid_length_of_permutation(Integer p_length, Integer N_req)
{
    std::ostringstream os;
    
    os  << "invalid permutation, expecting permutation of length " << N_req << ", supplied permutation has length "
        << p_length;

    current_message = os.str();
    return current_message.c_str(); 
};
const char* default_exception_message_linalg::band_matrix_required(Integer ld, Integer ud, Integer req_ld, Integer req_ud)
{
    std::ostringstream os;

    bool cont = false;

    if (req_ud >= 0)    
    {
        os  << "expecting matrix with " << req_ud << " superdiagonals" << "; supplied matrix has " 
            << ud << "superdiagonals";
        cont = true;
    };
    if (req_ld >= 0)
    {
        if (cont)
            os << "; ";

        os  << "expecting matrix with " << req_ld << " subdiagonals" << "; supplied matrix has " 
            << ld << "subdiagonals";
    };

    current_message = os.str();
    return current_message.c_str(); 
};

const char* default_exception_message_linalg::invalid_partition(Integer part_length, Integer mat_size)
{
    std::ostringstream os;
    os << "invalid partition, partition size: " << part_length << ", matrix size: " << mat_size;

    current_message = os.str();
    return current_message.c_str(); 
}
const char* default_exception_message_linalg::invalid_partition_block(Integer block_size)
{
    std::ostringstream os;
    os << "invalid block size: " << block_size;

    current_message = os.str();
    return current_message.c_str(); 
}

const char* default_exception_message_linalg::invalid_partition_size(Integer req_size, Integer part_size)
{
    std::ostringstream os;
    os << "invalid particion, matrix size is: " << req_size << ", total size of blocks is: " << part_size;

    current_message = os.str();
    return current_message.c_str(); 
}

const char* default_exception_message_linalg::invalid_aggregation_matrix(Integer rows, Integer nnz)
{
    std::ostringstream os;
    os << "invalid aggregation matrix of size " << rows << " with " << nnz << " nonzeros"
       << "; expecting matrix with " << rows << " nonzeros";

    current_message = os.str();
    return current_message.c_str(); 
}

const char* default_exception_message_linalg::scotch_error(const std::string& str)
{
    std::ostringstream os;
    os << "scotch error: " << str;

    current_message = os.str();
    return current_message.c_str(); 
};

const char* default_exception_message_linalg::metis_error(Integer err)
{
    std::ostringstream os;
    os << "metis error: " << err;

    current_message = os.str();
    return current_message.c_str(); 
};

const char* default_exception_message_linalg::gsvd_failed()
{
    current_message = "gsvd failed";
    return current_message.c_str(); 
};

const char* default_exception_message_linalg::inaccurate_solution(Integer k, Real fe, Integer req_prec)
{
    std::ostringstream os;
    os << "inaccurate solution with " << k << "-th vector; estimated forward error is: " << fe
        <<", required precision is: " << req_prec;

    current_message = os.str();
    return current_message.c_str(); 
};

void default_exception_message_linalg::warning_inaccurate_solution(Integer k, Real fe, Integer req_prec)
{
    std::ostringstream os;
    os << "inaccurate solution with " << k << "-th vector; estimated forward error is: " << fe
        <<", required precision is: " << req_prec;

    warning(os.str());
};

const char* default_exception_message_linalg::invalid_dm_block(Integer block, Integer max_blocks, 
                                                               Integer type)
{
    std::ostringstream os;
    os << "invalid block number in dm_decomp: " << block << "; of type: ";
    if (type == 1)
        os << "overdetermined part";
    else if (type == 2)
        os << "exactly determined part";
    else
        os << "underdetermined part";

    os << "; maximum number of block index is " << max_blocks;

    current_message = os.str();
    return current_message.c_str(); 
};

const char* invalid_length_of_permutation::what(exception_message_linalg& em) const
{
    return em.invalid_length_of_permutation(m_p_length, m_N_req);
};

const char* band_matrix_required::what(exception_message_linalg& em) const
{
    return em.band_matrix_required(m_ld, m_ud, m_req_ld, m_req_ud);
};

const char* diagonal_matrix_required::what(exception_message_linalg& em) const
{
    return em.diagonal_matrix_required(m_ld, m_ud);
};

const char* triangular_matrix_required::what(exception_message_linalg& em) const
{
    return em.triangular_matrix_required(m_ld, m_ud);
};

void default_exception_message_linalg::warning_singular()
{
    std::cout << "warning: " << "matrix is singular" << "\n";
};

void set_global_messanger_linalg(exception_message_linalg_ptr msg)
{
    if (msg)
    {
        messanger = msg;
    };
};
exception_message_linalg_ptr get_global_messanger_linalg()
{
    return messanger;
};

const char* matcl_exception_linalg::what() const throw()
{
    m_message = what(*get_global_messanger_linalg());
    return m_message.c_str();
};

const char* matcl_exception_linalg::what(exception_message&) const
{
    m_message = "invalid exception_message";
    return m_message.c_str();
};

const char* error_lu::what(exception_message_linalg& em) const
{
    return em.error_lu(msg);
};

const char* error_arpack::what(exception_message_linalg& em) const
{
    return em.error_arpack(msg);
};

const char* error_lsolve::what(exception_message_linalg& em) const
{
    return em.error_lsolve(r1,c1,r2,c2);
};
const char* error_lsolve_rev::what(exception_message_linalg& em) const
{
    return em.error_lsolve_rev(r1,c1,r2,c2);
};

const char* error_fft::what(exception_message_linalg& em) const
{
    return em.error_fft(status);
};
const char* error_singular::what(exception_message_linalg& em) const
{
    return em.error_singular();
};
const char* error_schur_incorrect_factors::what(exception_message_linalg& em) const
{
    return em.error_schur_incorrect_factors();
};
const char* error_schur_U_not_computed::what(exception_message_linalg& em) const
{
    return em.error_schur_U_not_computed();
};
const char* error_gschur_QZ_not_computed::what(exception_message_linalg& em) const
{
    return em.error_gschur_QZ_not_computed();
};

const char* error_gen_sym_eigen_job::what(exception_message_linalg& em) const
{
    return em.error_gen_sym_eigen_job();
};
const char* error_gschur_incorrect_factors::what(exception_message_linalg& em) const
{
    return em.error_gschur_incorrect_factors();
};
const char* error_size_eig::what(exception_message_linalg& em) const
{
    return em.error_size_eig(r1,c1);
};
const char* error_size_qrs::what(exception_message_linalg& em) const
{
    return em.error_size_qrs(r1,c1);
};
const char* error_qrs::what(exception_message_linalg& em) const
{
    return em.error_qrs(code);
};
const char* error_size_ldl::what(exception_message_linalg& em) const
{
    return em.error_size_ldl(r1,c1);
};
const char* error_ldl::what(exception_message_linalg& em) const
{
    return em.error_ldl();
};
const char* error_size_geig::what(exception_message_linalg& em) const
{
    return em.error_size_geig(r1,c1,r2,c2);
};
const char* error_size_sylv::what(exception_message_linalg& em) const
{
    return em.error_size_sylv(r1,c1,r2,c2, r3,c3);
};
const char* error_size_gsylv::what(exception_message_linalg& em) const
{
    return em.error_size_gsylv(rA1,cA1, rA2, cA2, rB1, cB1, rB2, cB2, rC1, cC1, rC2, cC2);
};
const char* error_size_ricc::what(exception_message_linalg& em) const
{
    return em.error_size_ricc();
};
const char* error_size_lyapunov::what(exception_message_linalg& em) const
{
    return em.error_size_lyapunov();
};
const char* error_sylv::what(exception_message_linalg& em) const
{
    return em.error_sylv();
};
const char* error_gsylv::what(exception_message_linalg& em) const
{
    return em.error_gsylv();
};
const char* error_ricc::what(exception_message_linalg& em) const
{
    return em.error_ricc(info);
};
const char* error_lyapunov::what(exception_message_linalg& em) const
{
    return em.error_lyapunov();
};
const char* error_schur::what(exception_message_linalg& em) const
{
    return em.error_schur();
};
const char* error_gschur::what(exception_message_linalg& em) const
{
    return em.error_gschur();
};
const char* error_norm_type::what(exception_message_linalg& em) const
{
    return em.error_norm_type();
};
const char* error_svd::what(exception_message_linalg& em) const
{
    return em.error_svd();
};
const char* error_svd_inf_nan::what(exception_message_linalg& em) const
{
    return em.error_svd_inf_nan();
};
const char* error_hess_nonsq::what(exception_message_linalg& em) const
{
    return em.error_hess_nonsq();
};
const char* error_gen_hess::what(exception_message_linalg& em) const
{
    return em.error_gen_hess(m_Ar, m_Ac, m_Br, m_Bc);
};
const char* error_hess_eig_invalid_init::what(exception_message_linalg& em) const
{
    return em.error_hess_eig_invalid_init(m_I_rows, m_I_cols, m_H_size, m_E_size, m_left_init);
};
const char* error_hess_eig_invalid_conj_pair::what(exception_message_linalg& em) const
{
    return em.error_hess_eig_invalid_conj_pair(m_pos, m_missing_conj);
};

const char* error_hess_eig_invalid_eig::what(exception_message_linalg& em) const
{
    return em.error_hess_eig_invalid_eig(m_E_rows, m_E_cols, m_H_size);
};
const char* error_hess_eig_her_complex_eig::what(exception_message_linalg& em) const
{
    return em.error_hess_eig_her_complex_eig(m_pos);
};

const char* error_schur_sel::what(exception_message_linalg& em) const
{
    return em.error_schur_sel();
};
const char* error_schur_sel_comp::what(exception_message_linalg& em) const
{
    return em.error_schur_sel_comp();
};
const char* error_schur_cond_vec::what(exception_message_linalg& em) const
{
    return em.error_schur_cond_vec(m_exp_N,m_exp_M,m_sup_N,m_sup_M,m_left_vector);
};
const char* error_schur_eig_cluster_size::what(exception_message_linalg& em) const
{
    return em.error_schur_eig_cluster_size(m_cluster,m_size);
};

const char* error_gschur_sel_comp::what(exception_message_linalg& em) const
{
    return em.error_gschur_sel_comp();
};
const char* error_nonsymh_chol::what(exception_message_linalg& em) const
{
    return em.error_nonsymh_chol();
};
const char* error_nonposdef::what(exception_message_linalg& em) const
{
    return em.error_nonposdef(m_allow_semipos);
};

const char* error_size_tridiag_sym::what(exception_message_linalg& em) const
{
    return em.error_size_tridiag_sym(A_rows,A_cols,B_rows,B_cols);
};
const char* error_size_tridiag_nsym::what(exception_message_linalg& em) const
{
    return em.error_size_tridiag_nsym(D_m1_rows,D_m1_cols,D0_rows,D0_cols,D_p1_rows,D_p1_cols);
};

const char* error_symeigen_complex_diagonal::what(exception_message_linalg& em) const
{
    return em.error_symeigen_complex_diagonal();
};

const char* error_cholmod::what(exception_message_linalg& em) const
{
    return em.error_cholmod(chol_msg);
};
const char* error_complex_value_type_not_allowed::what(exception_message_linalg& em) const
{
    return em.error_complex_value_type_not_allowed();
};
const char* invalid_unitary_matrix::what(exception_message_linalg& em) const
{
    return em.invalid_unitary_matrix();
};

const char* symmetric_matrix_required::what(exception_message_linalg& em) const
{
    return em.symmetric_matrix_required();
};
const char* tridiagonal_matrix_required::what(exception_message_linalg& em) const
{
    return em.tridiagonal_matrix_required(m_rows, m_cols, m_ld, m_ud);
};

const char* error_upper_hess_matrix_required::what(exception_message_linalg& em) const
{
    return em.error_upper_hess_matrix_required();
};

const char* error_invalid_chol_update_nonposdef::what(exception_message_linalg& em) const
{
    return em.error_invalid_chol_update_nonposdef();
};
const char* error_invalid_chol_factor::what(exception_message_linalg& em) const
{
    return em.error_invalid_chol_factor(m_upper);
};
const char* error_invalid_qr_factor::what(exception_message_linalg& em) const
{
    return em.error_invalid_qr_factor();
};
const char* error_invalid_rq_factor::what(exception_message_linalg& em) const
{
    return em.error_invalid_rq_factor(m_rows, m_cols);
};

const char* error_invalid_chol_update::what(exception_message_linalg& em) const
{
    return em.error_invalid_chol_update(m_T_rows,m_w_rows,m_w_cols);
}
const char* error_invalid_update_sigma::what(exception_message_linalg& em) const
{
    return em.error_invalid_update_sigma(m_K_vector, m_sigma_r, m_sigma_c);
}
const char* error_invalid_qr_update::what(exception_message_linalg& em) const
{
    return em.error_invalid_qr_update(m_T_rows, m_T_cols, m_u_rows, m_u_cols, m_v_rows, m_v_cols);
}

const char* error_invalid_givens_seq_size::what(exception_message_linalg& em) const
{
    return em.error_invalid_givens_seq_size(m_C_r, m_C_c, m_S_r, m_S_c, m_Ind_r, m_Ind_c, m_invalid_mat_pos);
}

const char* error_invalid_givens_seq_value_code::what(exception_message_linalg& em) const
{
    return em.error_invalid_givens_seq_value_code(m_value_code, m_mat_pos);
}

const char* error_invalid_givens_seq_index::what(exception_message_linalg& em) const
{
    return em.error_invalid_givens_seq_index(m_pos_row, m_pos_col, m_index, m_matrix_size);
}
const char* invalid_holseholder_col::what(exception_message_linalg& em) const
{
    return em.invalid_holseholder_col(m_col, m_max_col);
}
const char* uninitialized_object_used::what(exception_message_linalg& em) const
{
    return em.uninitialized_object_used(m_object);
}
const char* unable_to_continue_arnoldi::what(exception_message_linalg& em) const
{
    return em.unable_to_continue_arnoldi();
}
const char* invalid_speigs_k::what(exception_message_linalg& em) const
{
    return em.invalid_speigs_k(m_k, m_N);
}

const char* invalid_lu_factors::what(exception_message_linalg& em) const
{
    return em.invalid_lu_factors();
};
const char* invalid_cholecky_factors::what(exception_message_linalg& em) const
{
    return em.invalid_cholecky_factors();
};
const char* invalid_svd_factors::what(exception_message_linalg& em) const
{
    return em.invalid_svd_factors();
};
const char* invalid_schur_factors::what(exception_message_linalg& em) const
{
    return em.invalid_schur_factors();
};
const char* invalid_qr_factors::what(exception_message_linalg& em) const
{
    return em.invalid_qr_factors();
};
const char* invalid_ldl_factors::what(exception_message_linalg& em) const
{
    return em.invalid_ldl_factors();
};
const char* invalid_hess_factors::what(exception_message_linalg& em) const
{
    return em.invalid_hess_factors();
};
const char* invalid_bidiag_factors::what(exception_message_linalg& em) const
{
    return em.invalid_bidiag_factors();
};
const char* invalid_diagonal_22::what(exception_message_linalg& em) const
{
    return em.invalid_diagonal_22();
};
const char* band_matrix_with_main_diag_required::what(exception_message_linalg& em) const
{
    return em.band_matrix_with_main_diag_required(m_first_diag, m_last_diag);
};
const char* iterative_solver_not_converged::what(exception_message_linalg& em) const
{
    return em.iterative_solver_not_converged();
};
const char* logdet_not_available::what(exception_message_linalg& em) const
{
    return em.iterative_solver_not_converged();
};
const char* unable_convert_linear_operator_to_matrix::what(exception_message_linalg& em) const
{
    return em.unable_convert_linear_operator_to_matrix();
};

const char* invalid_partition::what(exception_message_linalg& em) const
{
    return em.invalid_partition(m_part_length, m_mat_size);
};

const char* invalid_partition_block::what(exception_message_linalg& em) const
{
    return em.invalid_partition_block(m_block);
};

const char* invalid_partition_size::what(exception_message_linalg& em) const
{
    return em.invalid_partition_size(m_req_size, m_part_size);
};

const char* scotch_error::what(exception_message_linalg& em) const
{
    return em.scotch_error(m_msg);
};

const char* metis_error::what(exception_message_linalg& em) const
{
    return em.metis_error(m_err);
};

const char* invalid_aggregation_matrix::what(exception_message_linalg& em) const
{
    return em.invalid_aggregation_matrix(m_rows, m_nnz);
};

const char* invalid_dm_block::what(exception_message_linalg& em) const
{
    return em.invalid_dm_block(m_block, m_max_blocks, m_type);
};

const char* gsvd_failed::what(exception_message_linalg& em) const
{
    return em.gsvd_failed();
};

const char* inaccurate_solution::what(exception_message_linalg& em) const
{
    return em.inaccurate_solution(m_vec, m_fe, m_req_prec);
};


};};