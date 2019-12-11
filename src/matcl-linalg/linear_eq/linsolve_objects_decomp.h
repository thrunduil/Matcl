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

#include "linsolve_objects_struct.h"

namespace matcl { namespace details
{

class linsolve_obj_lu_factors : public linsolve_obj_seq_2
{
    public:
        linsolve_obj_lu_factors(const Matrix& A, const Matrix& L, const Matrix& U, 
                const permvec& P, const permvec& Q, const options& opts);

        ~linsolve_obj_lu_factors();
};

class linsolve_obj_chol : public linsolve_obj_seq_2
{
    public:
        linsolve_obj_chol(const Matrix& A, const Matrix& L, const permvec& p, bool upper,
                          const options& opts);

        ~linsolve_obj_chol();

        virtual Real            log_det() const override;
};

class linsolve_obj_svd : public linsolve_obj_seq_3
{
    public:
        linsolve_obj_svd(const Matrix& A, const unitary_matrix& U, const Matrix& S, const unitary_matrix& V,
                         const options& opts);

        ~linsolve_obj_svd();

        virtual Real            normest_2() const override;
        virtual Real            mat_normest_2() const override;
};

class linsolve_obj_schur : public linsolve_obj_seq_3
{
    public:
        linsolve_obj_schur(const Matrix& A, const unitary_matrix& U, const Matrix& T,
                           const options& opts);

        ~linsolve_obj_schur();

        virtual Real            normest_2() const override;
        virtual Real            mat_normest_2() const override;
};

class linsolve_obj_bidiag : public linsolve_obj_seq_3
{
    public:
        linsolve_obj_bidiag(const Matrix& A, const unitary_matrix& U, const Matrix& R, const unitary_matrix& V,
                            const options& opts);

        ~linsolve_obj_bidiag();

        virtual Real            normest_2() const override;
        virtual Real            mat_normest_2() const override;
};

class linsolve_obj_rq : public linsolve_obj_seq_2
{
    public:
        //RL is lower or upper triangular
        linsolve_obj_rq(const Matrix& A, const Matrix& RL, const unitary_matrix& Q,
                        const options& opts);

        ~linsolve_obj_rq();

        virtual Real            normest_2() const override;
        virtual Real            mat_normest_2() const override;
};

class linsolve_obj_qr : public linsolve_obj_seq_2
{
    public:
        //RL is lower or upper triangular
        linsolve_obj_qr(const Matrix& A, const unitary_matrix& Q, const Matrix& RL, const permvec& p,
                        const options& opts);

        ~linsolve_obj_qr();

        virtual Real            normest_2() const override;
        virtual Real            mat_normest_2() const override;
};

class linsolve_obj_hess : public linsolve_obj_seq_3
{
    public:
        linsolve_obj_hess(const Matrix& A, const unitary_matrix& U, const Matrix& H, const options& opts);

        ~linsolve_obj_hess();

        virtual Real            normest_2() const override;
        virtual Real            mat_normest_2() const override;
};

class linsolve_obj_ldl : public linsolve_obj_seq_3
{
    private:
        bool                    m_sym;

    public:
        linsolve_obj_ldl(const Matrix& A, const Matrix& L, const Matrix& D, const permvec& p, bool sym,
                         const options& opts);

        ~linsolve_obj_ldl();

        virtual Real            log_det() const override;

        virtual data_ptr        convert(value_code new_val_code) const override;

    private:
        linsolve_obj_ldl(const Matrix& A, const linsolve_obj& L1, const linsolve_obj& D, 
                         const linsolve_obj& L2, bool sym);
};

template<class V>
class linsolve_obj_lu_dense : public linsolve_obj_base
{
    private:
        using Mat       = raw::Matrix<V, struct_dense>;
        using Mat_I     = raw::Matrix<Integer, struct_dense>;        

    private:
        Mat                     m_A;
        Mat                     m_A_decomp;
        Mat_I                   m_piv;

    public:
        linsolve_obj_lu_dense(const Mat& A, const Mat& A_decomp, const Mat_I& ipiv,
                              const options& opts);
        linsolve_obj_lu_dense(const Mat& A, const Mat& A_decomp, const Mat_I& ipiv, bool modif, from_inv);

        ~linsolve_obj_lu_dense();

        virtual bool            is_hermitian() const override;
        virtual bool            is_posdef() const override;
        virtual bool            all_finite() const override     { return true; };
        virtual Real            log_det() const override;
        virtual Real            mat_normest_1() const override;
        virtual Real            mat_normest_inf() const override;

        virtual data_ptr        convert(value_code new_val_code) const override;

        // create inverse matrix
        virtual matcl::Matrix   inv() const override;
        virtual matcl::Matrix   base_matrix() const override;

        // solve op(A) * Y = X, where A is the matrix represented by this object
        virtual matcl::Matrix   solve(const Matrix& X, trans_type tA) const override;
        virtual matcl::Matrix   solve(Matrix&& X, trans_type tA) const override;

        // solve Y * op(A) = X, where A is the matrix represented by this object
        virtual matcl::Matrix   solve_rev(const Matrix& X, trans_type tA) const override;
        virtual matcl::Matrix   solve_rev(Matrix&& X, trans_type tA) const override;

        virtual Matrix          mmul_right(const Matrix& X, trans_type t) const override;
        virtual Matrix          mmul_right(Matrix&& X, trans_type t) const override;
        virtual Matrix          mmul_left(const Matrix& X, trans_type t) const override;
        virtual Matrix          mmul_left(Matrix&& X, trans_type t) const override;

        template<class T>
        data_ptr                convert_impl() const;

        template<class T>
        void                    solve_impl(Matrix& ret, Matrix& X, trans_type tA) const;
        template<class T>
        void                    solve_rev_impl(Matrix& ret, Matrix& X, trans_type tA) const;

    private:
        void                    correct_singular(Mat& m_A_decomp, const options& opts);
};

template<class V>
class linsolve_obj_lu_band : public linsolve_obj_base
{
    private:
        using Mat       = raw::Matrix<V, struct_banded>;
        using Mat_I     = raw::Matrix<Integer, struct_dense>;

    private:
        Mat                     m_A_decomp;
        Mat                     m_A;
        Mat_I                   m_piv;
        Integer                 m_N;
        Integer                 m_ldiags;
        Integer                 m_udiags;

    public:
        linsolve_obj_lu_band(const Mat& A, const Mat& A_decomp, const Mat_I& ipiv, 
                                Integer N, Integer ld, Integer ud, const options& opts);
        linsolve_obj_lu_band(const Mat& A, const Mat& A_decomp, const Mat_I& ipiv, 
                                Integer N, Integer ld, Integer ud, bool modif, from_inv);

        ~linsolve_obj_lu_band();

        virtual bool            is_hermitian() const override;
        virtual bool            is_posdef() const override;
        virtual bool            all_finite() const override     { return true; };
        virtual Real            log_det() const override;
        virtual Real            mat_normest_1() const override;
        virtual Real            mat_normest_inf() const override;

        virtual data_ptr        convert(value_code new_val_code) const override;

        virtual matcl::Matrix   inv() const override;
        virtual matcl::Matrix   base_matrix() const override;

        // solve op(A) * Y = X, where A is the matrix represented by this object
        virtual matcl::Matrix   solve(const Matrix& X, trans_type tA) const override;
        virtual matcl::Matrix   solve(Matrix&& X, trans_type tA) const override;

        // solve Y * op(A) = X, where A is the matrix represented by this object
        virtual matcl::Matrix   solve_rev(const Matrix& X, trans_type tA) const override;
        virtual matcl::Matrix   solve_rev(Matrix&& X, trans_type tA) const override;

        virtual Matrix          mmul_right(const Matrix& X, trans_type t) const override;
        virtual Matrix          mmul_right(Matrix&& X, trans_type t) const override;
        virtual Matrix          mmul_left(const Matrix& X, trans_type t) const override;
        virtual Matrix          mmul_left(Matrix&& X, trans_type t) const override;

        template<class T>
        data_ptr                convert_impl() const;

        template<class T>
        void                    solve_impl(Matrix& ret, Matrix& X, trans_type tA) const;
        template<class T>
        void                    solve_rev_impl(Matrix& ret, Matrix& X, trans_type tA) const;

    private:
        void                    correct_singular(Mat& m_A_decomp, const options& opts);
};

template<class V>
class linsolve_obj_chol_tridiag_fac : public linsolve_obj_base
{
    private:
        using VR                = typename real_type<V>::type;
        using Mat               = raw::Matrix<V,struct_dense>;
        using Mat_R             = raw::Matrix<VR,struct_dense>;
        using Mat_B             = raw::Matrix<V,struct_banded>;

    private:
        Mat_B                   m_A;
        Mat_R                   m_D0_i;
        Mat                     m_D1_i;

    public:
        linsolve_obj_chol_tridiag_fac(const Mat_B& A, const Mat_R& D0_i, const Mat& D1_i,
                                      const options& opts);
        linsolve_obj_chol_tridiag_fac(const Mat_B& A, const Mat_R& D0_i, const Mat& D1_i, 
                                      bool modif, from_inv);

        ~linsolve_obj_chol_tridiag_fac();

        virtual bool            is_hermitian() const override   { return true; };
        virtual bool            is_posdef() const override      { return true; };
        virtual bool            all_finite() const override     { return true; };
        virtual Real            log_det() const override;
        virtual Real            normest_1() const override;
        virtual Real            normest_inf() const override;

        virtual matcl::Matrix   inv() const override;

        virtual matcl::Matrix   solve(const Matrix& X, trans_type tA) const override;
        virtual matcl::Matrix   solve(Matrix&& X, trans_type tA) const override;

        virtual matcl::Matrix   solve_rev(const Matrix& X, trans_type tA) const override;
        virtual matcl::Matrix   solve_rev(Matrix&& X, trans_type tA) const override;

        template<class T>
        void                    solve_impl(Matrix& ret, Matrix& X, trans_type tA) const;
        template<class T>
        void                    solve_rev_impl(Matrix& ret, Matrix& X, trans_type tA) const;

        virtual Real            mat_normest_1() const override;
        virtual Real            mat_normest_inf() const override;

        virtual data_ptr        convert(value_code new_val_code) const override;
        virtual matcl::Matrix   base_matrix() const override;

        virtual Matrix          mmul_right(const Matrix& X, trans_type t) const override;
        virtual Matrix          mmul_right(Matrix&& X, trans_type t) const override;
        virtual Matrix          mmul_left(const Matrix& X, trans_type t) const override;
        virtual Matrix          mmul_left(Matrix&& X, trans_type t) const override;

        template<class T>
        data_ptr                convert_impl() const;

    private:        
        void                    test_singular() const;
};

template<class V>
class linsolve_obj_tridiag_fac : public linsolve_obj_base
{
    private:
        using VR                = typename real_type<V>::type;
        using Mat               = raw::Matrix<V,struct_dense>;
        using Mat_R             = raw::Matrix<VR,struct_dense>;
        using Mat_B             = raw::Matrix<V,struct_banded>;
        using Mat_I             = raw::Matrix<Integer,struct_dense>;

    private:
        Mat_B                   m_A;
        Mat                     m_Dm1_i;
        Mat                     m_D0_i;
        Mat                     m_Dp1_i;
        Mat                     m_Dp2_i;
        Mat_I                   m_piv;

    public:
        linsolve_obj_tridiag_fac(const Mat_B& A, const Mat& Dm1_i, const Mat& D0_i, const Mat& Dp1_i,
                                 const Mat& Dp2_i, const Mat_I& piv, const options& opts);
        linsolve_obj_tridiag_fac(const Mat_B& A, const Mat& Dm1_i, const Mat& D0_i, const Mat& Dp1_i,
                                 const Mat& Dp2_i, const Mat_I& piv, bool modif, from_inv);

        ~linsolve_obj_tridiag_fac();

        virtual bool            is_hermitian() const override   { return false; };
        virtual bool            is_posdef() const override      { return false; };
        virtual bool            all_finite() const override     { return true; };
        virtual Real            log_det() const override;

        virtual matcl::Matrix   inv() const override;

        virtual matcl::Matrix   solve(const Matrix& X, trans_type tA) const override;
        virtual matcl::Matrix   solve(Matrix&& X, trans_type tA) const override;

        virtual matcl::Matrix   solve_rev(const Matrix& X, trans_type tA) const override;
        virtual matcl::Matrix   solve_rev(Matrix&& X, trans_type tA) const override;

        template<class T>
        void                    solve_impl(Matrix& ret, Matrix& X, trans_type tA) const;
        template<class T>
        void                    solve_rev_impl(Matrix& ret, Matrix& X, trans_type tA) const;

        virtual Real            mat_normest_1() const override;
        virtual Real            mat_normest_inf() const override;

        virtual data_ptr        convert(value_code new_val_code) const override;
        virtual matcl::Matrix   base_matrix() const override;

        virtual Matrix          mmul_right(const Matrix& X, trans_type t) const override;
        virtual Matrix          mmul_right(Matrix&& X, trans_type t) const override;
        virtual Matrix          mmul_left(const Matrix& X, trans_type t) const override;
        virtual Matrix          mmul_left(Matrix&& X, trans_type t) const override;

        template<class T>
        data_ptr                convert_impl() const;
};

class linsolve_forwarding : public linsolve_obj_base
{
    private:
        linsolve_obj            m_base;

    public:
        linsolve_forwarding(const linsolve_obj& base);
        ~linsolve_forwarding(){};

        const linsolve_obj&     base() const                    { return m_base; };

        virtual bool            all_finite() const override     { return true; };
        virtual bool            is_hermitian() const override   { return m_base.is_hermitian(); };
        virtual bool            is_posdef() const override      { return m_base.is_posdef(); };
        virtual bool            is_modified() const override    { return m_base.is_modified(); };
        virtual bool            is_direct() const override      { return m_base.is_direct(); };
        virtual Matrix          inv() const override            { return m_base.inv(); };
        virtual Matrix          base_matrix() const override    { return m_base.base_matrix(); };
        virtual Real            log_det() const override        { return m_base.log_det(); };
        virtual Real            normest_1() const override      { return m_base.normest(basic_vector_norm::norm_1); };
        virtual Real            normest_2() const override      { return m_base.normest(basic_vector_norm::norm_2); };
        virtual Real            normest_inf() const override    { return m_base.normest(basic_vector_norm::norm_inf); };
        virtual Real            mat_normest_1() const override  { return m_base.mat_normest(basic_vector_norm::norm_1); };
        virtual Real            mat_normest_2() const override  { return m_base.mat_normest(basic_vector_norm::norm_2); };
        virtual Real            mat_normest_inf() const override{ return m_base.mat_normest(basic_vector_norm::norm_inf); };

        virtual data_ptr        convert(value_code new_val_code) const override;        
        virtual Matrix          solve(const Matrix& X, trans_type tA) const override;
        virtual Matrix          solve(Matrix&& X, trans_type tA) const override;
        virtual Matrix          solve_rev(const Matrix& X, trans_type tA) const override;
        virtual Matrix          solve_rev(Matrix&& X, trans_type tA) const override;        
        virtual Matrix          mmul_right(const Matrix& X, trans_type t) const override;
        virtual Matrix          mmul_right(Matrix&& X, trans_type t) const override;
        virtual Matrix          mmul_left(const Matrix& X, trans_type t) const override;
        virtual Matrix          mmul_left(Matrix&& X, trans_type t) const override;
};

class linsolve_perm : public linsolve_forwarding
{
    private:
        Matrix      m_A;
        permvec     m_p;
        permvec     m_q;
        bool        m_has_A;
        bool        m_symperm;

    public:
        // A can be complex; internally QR decomp is performed
        linsolve_perm(const Matrix& A, const linsolve_obj& Apq, const permvec& p, const permvec& q,
                      bool symperm);
        linsolve_perm(const linsolve_obj& Apq, const permvec& p, const permvec& q, bool symperm);
        ~linsolve_perm();

        virtual bool            is_hermitian() const override;
        virtual bool            is_posdef() const override;
        virtual Matrix          inv() const override;
        virtual Matrix          base_matrix() const override;

        virtual data_ptr        convert(value_code new_val_code) const override;        

        virtual Matrix          solve(const Matrix& X, trans_type tA) const override;
        virtual Matrix          solve(Matrix&& X, trans_type tA) const override;
        virtual Matrix          solve_rev(const Matrix& X, trans_type tA) const override;
        virtual Matrix          solve_rev(Matrix&& X, trans_type tA) const override;        
        virtual Matrix          mmul_right(const Matrix& X, trans_type t) const override;
        virtual Matrix          mmul_right(Matrix&& X, trans_type t) const override;
        virtual Matrix          mmul_left(const Matrix& X, trans_type t) const override;
        virtual Matrix          mmul_left(Matrix&& X, trans_type t) const override;
};

};};