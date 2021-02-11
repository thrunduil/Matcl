/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include "matcl-linalg/decompositions/eig_functions.h"
#include "matcl-linalg/special_matrices/matrix_functors.h"
#include "matcl-linalg/decompositions/speigs.h"

namespace matcl { namespace details
{

// Choice of ARPACK routines: dnaupd/dneupd, snaupd/sneupd, znaupd/zneupd respectively
enum class kernel_type
{
    nonsymmetric, symmetric, complex
};

enum class arpack_output
{
    eig,        // eigenvalues only
    schur,      // eigenvalues and schur vectors
    ritz,       // eigenvalues and ritz vectors
    schur_ritz, // eigenvalues, schur and ritz vectors
};

// Helper class to perform arpack calls
// This class locks an arpack mutex for as long as it is in scope
template<class T>
class arpack_wrapper
{
    private:
        using mutex_type    = matcl::default_mutex;
        using lock_type     = std::unique_lock<mutex_type>;
        using Mat_D         = raw::Matrix<T, struct_dense>;
        using TR            = typename details::real_type<T>::type;

    private:
        arpack_wrapper(const arpack_wrapper&) = delete;
        arpack_wrapper& operator=(const arpack_wrapper&) = delete;

    public:
        // create standard version for Av = lambda*v problem
        arpack_wrapper(const linear_operator& A, arpack_output out, Integer max_iter, Integer n_arnoldi,
                       Real tol);

        // create generalized version for Av = lambda*B*v problem
        arpack_wrapper(const linear_operator& A, const linear_operator& B, bool hermitian, 
                       arpack_output out, Integer max_iter, Integer n_arnoldi, Real tol);

        ~arpack_wrapper();

        // set initial vector
        void                set_x0(const Matrix& x0);

        // set select type
        void                set_cluster(cluster_type ec);

        // do all calculations; k is number of eigenvalues requested;
        // if return_nonconverged is true, then return also eigenvalues which have not converged
        void                calculate(Integer k, bool return_nonconverged);

        // size of the operators
        Integer             size() const;

        // schur vectors
        Matrix              get_U() const;

        // quasi-upper triangular matrix such that (A*B^-1) * U = U * T
        Matrix              get_T() const;

        // eigenvectors
        Matrix              get_V() const;

        // eigenvalues
        Matrix              get_D() const;

        // return true if all eigenvalues converged
        bool                get_converged() const;

        value_code          get_value_code() const;

        // return number of converged eigenvalues
        Integer             number_converged_eigenvalues() const;

        // profiling statistics
        Integer             number_Arnoldi_iterations() const;
        Integer             number_operations_Ax() const;
        Integer             number_operations_Bx() const;
        Integer             number_reorthogonalizations() const;

    private:
        static kernel_type  get_kernel_type(const linear_operator& A);
        static kernel_type  get_kernel_type(bool hermitian);
        void                prepare_dimmensions(const Integer k);
        void                prepare_integer_arrays();
        void                prepare_arrays();
        Integer             do_arpack_loop(Mat_D& v);        
        void                process_xyeupd(Integer xyaupd_info, bool return_nonconverged, Mat_D& v);
        void                call_xyaupd(Integer& ido, Integer& info, Mat_D& v);
        void                call_xyeupd(Integer& info, Mat_D& v, Mat_D& eig);
        void                extract_H(Integer results_to_return);

    private:
        lock_type           m_lock;
        Integer             m_N;
        TR                  m_tol;
        kernel_type         m_kernel_type;
        arpack_output       m_output;
        linear_operator     m_A;
        linear_operator     m_B;
        bool                m_has_x0;
        Integer             m_num_arnoldi;

        cluster_type        m_cluster;
        bool                m_calculate_called;

        Matrix              m_U_result;
        Matrix              m_V_result;
        Matrix              m_D;
        Matrix              m_T_result;
        bool                m_converged;    
        bool                m_generalized;

        // arpack arrays & workspaces
        Integer             m_num_eig;
        Integer             m_num_vec;
        Integer             m_max_iter;
        Matrix              m_v;
        Integer             m_ldv;
        Matrix              m_iparam;
        Integer*            m_iparam_ptr;
        Matrix              m_ipntr;
        Integer*            m_ipntr_ptr;
        Matrix              m_select;
        Integer*            m_select_ptr;
        Matrix              m_x0;
        T*                  m_x0_ptr;
        Matrix              m_workd;
        T*                  m_workd_ptr; 
        Integer             m_lworkl;
        Matrix              m_workl;
        T*                  m_workl_ptr;
        Matrix              m_rwork;
        TR*                 m_rwork_ptr;
        Matrix              m_workev;
        T*                  m_workev_ptr;

        Matrix              m_vec_schur;
        Matrix              m_vec_ritz;        
};

}};