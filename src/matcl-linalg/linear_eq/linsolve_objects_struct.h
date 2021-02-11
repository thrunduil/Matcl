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

#include "matcl-linalg/linear_eq/linsolve_object.h"
#include "matcl-linalg/special_matrices/unitary_matrix.h"
#include "matcl-core/lib_functions/constants.h"

namespace matcl { namespace details
{

struct from_inv{};

// All linsolve objects except linsolve_obj_nan cannot contain NAN values

class linsolve_obj_base : public linsolve_obj_data
{
    protected:
        Integer                 m_M;
        Integer                 m_N;
        value_code              m_vc;
        ti::ti_object           m_ti;
        bool                    m_modified;

    public:
        linsolve_obj_base(Integer M, Integer N, value_code vc, const ti::ti_object& ti,
                          bool modif);
        ~linsolve_obj_base();

        // get number of rows
        virtual Integer         rows() const override;

        // get number of columns
        virtual Integer         cols() const  override;

        // code of stored elements type
        virtual value_code      get_value_code() const override;

        // return type_info of stored elements
        virtual ti::ti_object   get_type() const override;

        // return true if small perturbations was added to make factors nonsingular
        virtual bool            is_modified() const override;

        virtual bool            is_direct() const override;
};

class linsolve_obj_empty : public linsolve_obj_base
{
    public:
        linsolve_obj_empty(value_code vc, const ti::ti_object& ti);
        ~linsolve_obj_empty();

        virtual bool            is_hermitian() const override   { return false; };
        virtual bool            is_posdef() const override      { return false; };
        virtual Real            log_det() const override        { return 0.0; };
        virtual Real            normest_1() const override      { return 0.0; };
        virtual Real            normest_2() const override      { return 0.0; };
        virtual Real            normest_inf() const override    { return 0.0; };
        virtual Real            mat_normest_1() const override  { return 0.0; };
        virtual Real            mat_normest_2() const override  { return 0.0; };
        virtual Real            mat_normest_inf() const override{ return 0.0; };
        virtual bool            all_finite() const override     { return true; };

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
};

class linsolve_obj_nan : public linsolve_obj_base
{
    public:
        linsolve_obj_nan(Integer N, value_code vc, const ti::ti_object& ti);
        ~linsolve_obj_nan();

        virtual bool            is_hermitian() const override   { return false; };
        virtual bool            is_posdef() const override      { return false; };
        virtual Real            log_det() const override        { return constants::nan(); };
        virtual Real            normest_1() const override      { return constants::nan(); };
        virtual Real            normest_2() const override      { return constants::nan(); };
        virtual Real            normest_inf() const override    { return constants::nan(); };
        virtual Real            mat_normest_1() const override  { return constants::nan(); };
        virtual Real            mat_normest_2() const override  { return constants::nan(); };
        virtual Real            mat_normest_inf() const override{ return constants::nan(); };
        virtual bool            all_finite() const override     { return false; };

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
};

class linsolve_obj_id : public linsolve_obj_base
{
    public:
        linsolve_obj_id(Integer N, value_code vc, const ti::ti_object& ti);
        ~linsolve_obj_id();

        virtual bool            is_hermitian() const override   { return true; };
        virtual bool            is_posdef() const override      { return true; };
        virtual Real            log_det() const override        { return 0.0; };
        virtual Real            normest_1() const override      { return 1.0; };
        virtual Real            normest_2() const override      { return 1.0; };
        virtual Real            normest_inf() const override    { return 1.0; };
        virtual Real            mat_normest_1() const override  { return 1.0; };
        virtual Real            mat_normest_2() const override  { return 1.0; };
        virtual Real            mat_normest_inf() const override{ return 1.0; };
        virtual bool            all_finite() const override     { return true; };

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
};

class linsolve_obj_scalar : public linsolve_obj_base
{
    private:
        struct from_inv{};

    private:
        Matrix                  m_A;

    public:
        linsolve_obj_scalar(const Matrix& A, const options& opts);

        ~linsolve_obj_scalar();

        virtual bool            all_finite() const override { return true; };

        virtual bool            is_hermitian() const override;
        virtual bool            is_posdef() const override;
        virtual Real            log_det() const override;
        virtual Real            normest_1() const override;
        virtual Real            normest_2() const override;
        virtual Real            normest_inf() const override;
        virtual Real            mat_normest_1() const override;
        virtual Real            mat_normest_2() const override;
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

    private:
        linsolve_obj_scalar(const Matrix& A, bool modif, from_inv);
};

class linsolve_obj_diag : public linsolve_obj_base
{
    private:
        struct from_inv{};

    public:
        struct from_vec_inv{};

    private:
        Matrix                  m_D;
        Matrix                  m_Di;
        bool                    m_is_posdef;

    public:
        linsolve_obj_diag(const Matrix& A, const options& opts);
        linsolve_obj_diag(const Matrix& D, from_vec_inv);
        ~linsolve_obj_diag();

        virtual bool            all_finite() const override { return true; };

        virtual bool            is_hermitian() const override;
        virtual bool            is_posdef() const override;
        virtual Real            log_det() const override;
        virtual Real            normest_1() const override;
        virtual Real            normest_2() const override;
        virtual Real            normest_inf() const override;
        virtual Real            mat_normest_1() const override;
        virtual Real            mat_normest_2() const override;
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

    private:
        linsolve_obj_diag(const Matrix& D, const Matrix& Di, bool is_posdef, bool modif, from_inv);

        void                    initialize(const Matrix& A, const options& opts);
        void                    initialize_vec_inv(const Matrix& A);
        void                    initialize_posdef(const Matrix& A);
};

class linsolve_obj_diag_22 : public linsolve_obj_base
{
    private:
        struct from_inv{};

    private:
        Matrix                  m_Ai;
        Matrix                  m_A;
        Real                    m_sig_min;
        Real                    m_sig_max;
        bool                    m_hermitian;
        bool                    m_posdef;
        Real                    m_logdet;

    public:
        linsolve_obj_diag_22(const Matrix& A, const options& opts);
        ~linsolve_obj_diag_22();

        virtual bool            is_hermitian() const override;
        virtual bool            is_posdef() const override;
        virtual Real            log_det() const override;
        virtual bool            all_finite() const override     { return true; };

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

        virtual Real            normest_2() const override;
        virtual Real            mat_normest_2() const override;

        //not reimplemented
        //virtual Real          normest_1() const override;
        //virtual Real          normest_inf() const override;
        //virtual Real          mat_normest_1() const override;
        //virtual Real          mat_normest_inf() const override;

    private:
        linsolve_obj_diag_22(const Matrix& A, const Matrix& Ai, bool hermitian, bool posdef,
                                Real logdet, Real sig_max, Real sig_min, bool modif, from_inv);

        void                    initialize(const Matrix& A, Real& logdet, Real& sig_max, Real& sig_min,
                                           const options& opts);
};

class linsolve_obj_unitary : public linsolve_obj_base
{
    private:
        struct from_inv{};

    private:
        Matrix                  m_A;

    public:
        linsolve_obj_unitary(const Matrix& A);
        ~linsolve_obj_unitary();

        virtual bool            is_hermitian() const override   { return false; };
        virtual bool            is_posdef() const override      { return false; };
        virtual Real            log_det() const override        { return 0.0;   };
        virtual Real            normest_2() const override      { return 1.0; };
        virtual Real            mat_normest_2() const override  { return 1.0; };
        virtual bool            all_finite() const override     { return true; };
        virtual Real            normest_1() const override;
        virtual Real            normest_inf() const override;
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

    private:
        linsolve_obj_unitary(const Matrix& A, from_inv);
};

class linsolve_obj_unitary_mat : public linsolve_obj_base
{
    private:
        struct from_inv{};

    private:
        unitary_matrix          m_A;

    public:
        linsolve_obj_unitary_mat(const unitary_matrix& A);
        ~linsolve_obj_unitary_mat();

        virtual bool            is_hermitian() const override   { return false; };
        virtual bool            is_posdef() const override      { return false; }
        virtual Real            log_det() const override        { return 0.0; };
        virtual Real            normest_2() const override      { return 1.0; };
        virtual Real            mat_normest_2() const override  { return 1.0; };
        virtual bool            all_finite() const override     { return true; };

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

        //not reimplemented
        //virtual Real          normest_1() const override;
        //virtual Real          normest_inf() const override;
        //virtual Real          mat_normest_1() const override;
        //virtual Real          mat_normest_inf() const override;
};

class linsolve_obj_triang : public linsolve_obj_base
{
    private:
        Matrix      m_LU;
        permvec     m_p;
        permvec     m_q;
        bool        m_tril;

    public:
        linsolve_obj_triang(const Matrix& LU, const options& opts);

        // solve op(LU(p^-1,q^-1)) * X = b
        linsolve_obj_triang(const Matrix& LU, const permvec& p, const permvec& q, const options& opts);
        linsolve_obj_triang(const Matrix& LU, const permvec& p, const permvec& q, bool modif, from_inv);
        ~linsolve_obj_triang();

        virtual bool            is_hermitian() const override   { return false; };
        virtual bool            is_posdef() const override      { return false; };
        virtual bool            all_finite() const override     { return true; };
        virtual Real            log_det() const override;

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

        // not reimplemented
        //virtual Real          normest_2() const override;
        //virtual Real          mat_normest_2() const override;
        //virtual Real          normest_1() const override;
        //virtual Real          normest_inf() const override;
        //virtual Real          mat_normest_1() const override;        
        //virtual Real          mat_normest_inf() const override;

    public:
        template<class Mat>
        void                    inv_impl(Matrix& ret) const;
        template<class V>
        void                    inv_scal_impl(Matrix& ret, const V& scal) const;

    private:
        void                    correct_singular(Matrix& A, const options& opts);
        void                    test_singular() const;
};

class linsolve_obj_seq_2 : public linsolve_obj_base
{
    protected:
        linsolve_obj            m_A1;
        linsolve_obj            m_A2;
        Matrix                  m_base_mat;
        bool                    m_has_base;

    public:
        linsolve_obj_seq_2();
        linsolve_obj_seq_2(const linsolve_obj& A1, const linsolve_obj& A2);
        linsolve_obj_seq_2(const Matrix& base, const linsolve_obj& A1, const linsolve_obj& A2);

        ~linsolve_obj_seq_2();

        virtual bool            is_direct() const override;

        void                    reset(const linsolve_obj& A1, const linsolve_obj& A2);
        void                    reset(const Matrix& A, const linsolve_obj& A1, const linsolve_obj& A2);

        virtual bool            all_finite() const override { return true; };

        virtual bool            is_hermitian() const override;
        virtual bool            is_posdef() const override;
        virtual Real            log_det() const override;

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

        virtual Real            mat_normest_1() const;
        virtual Real            mat_normest_inf() const;

        //not reimplemented
        //virtual Real          normest_2() const;
        //virtual Real          mat_normest_2() const;
        //virtual Real          normest_1() const;
        //virtual Real          normest_inf() const;
};

class linsolve_obj_seq_3 : public linsolve_obj_base
{
    protected:
        linsolve_obj            m_A1;
        linsolve_obj            m_A2;
        linsolve_obj            m_A3;
        Matrix                  m_base_mat;
        bool                    m_has_base;

    public:
        linsolve_obj_seq_3();
        linsolve_obj_seq_3(const linsolve_obj& A1, const linsolve_obj& A2, const linsolve_obj& A3);
        linsolve_obj_seq_3(const Matrix& base, const linsolve_obj& A1, const linsolve_obj& A2, 
                           const linsolve_obj& A3);

        ~linsolve_obj_seq_3();

        virtual bool            is_direct() const override;

        void                    reset(const linsolve_obj& A1, const linsolve_obj& A2, 
                                      const linsolve_obj& A3);
        void                    reset(const Matrix& base, const linsolve_obj& A1, const linsolve_obj& A2, 
                                      const linsolve_obj& A3);

        virtual bool            all_finite() const override     { return true; };
        virtual bool            is_hermitian() const override;
        virtual bool            is_posdef() const override;
        virtual Real            log_det() const override;
        virtual Real            mat_normest_1() const;
        virtual Real            mat_normest_inf() const;

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

        //not reimplemented
        //virtual Real          normest_2() const;
        //virtual Real          mat_normest_2() const;
        //virtual Real          normest_1() const;
        //virtual Real          normest_inf() const;
};

class linsolve_obj_uhess : public linsolve_obj_seq_2
{
    public:
        // A can be complex; internally QR decomp is performed
        linsolve_obj_uhess(const Matrix& A, const options& opts);
        ~linsolve_obj_uhess();
};

class linsolve_obj_lhess_dense : public linsolve_obj_seq_2
{
    public:
        // A can be complex; internally QR decomp is performed
        linsolve_obj_lhess_dense(const Matrix& A, const options& opts);
        ~linsolve_obj_lhess_dense();
};

};};