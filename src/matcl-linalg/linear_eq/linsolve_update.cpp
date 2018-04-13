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

#include "matcl-linalg/linear_eq/linsolve_object.h"
#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "linsolve_objects_decomp.h"
#include "matcl-linalg/norms_error/norm.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-linalg/special_matrices/struct_flag_linalg.h"

namespace matcl { namespace details
{

class SM_precomp
{
    public:
        struct ver_block{};

    private:
        Real            m_log_det;

        linsolve_obj    m_Q;

        //inv(A) * U * inv(Q)
        Matrix          m_AUQ;

        //V' * inv(A) = (inv(A)' * V)'
        Matrix          m_AV;

    public:
        bool all_finite() const
        {
            return m_Q.all_finite() && m_AUQ.all_finite() && m_AV.all_finite();
        };

        bool is_modified() const
        {
            return m_Q.is_modified();
        };

        bool is_direct() const
        {
            return m_Q.is_direct();
        };

        //I + V'*inv(A)*U
        SM_precomp(bool sym, const Matrix& U, const Matrix& V, const linsolve_obj& A, const options& opt)
        {
            value_code vc   = matrix_traits::unify_value_types(U.get_value_code(), A.get_value_code());
            vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

            if (sym == false)
                vc          = matrix_traits::unify_value_types(vc, V.get_value_code());

            //form: Q = I + V'*inv(A)*U

            Integer k   = U.cols();

            Matrix AU   = A.solve(U, trans_type::no_trans);

            if (k == 1)
            {
                Matrix Q;                
          
                if (sym == false)
                    Q = 1.0f + mmul(V, AU, trans_type::conj_trans, trans_type::no_trans);
                else
                    Q = 1.0f + mmul(U, AU, trans_type::conj_trans, trans_type::no_trans);

                m_Q         = make_linsolve_obj(Q,opt);
                m_AUQ       = AU * (1.0f / Q);
                m_log_det   = matcl::log(matcl::abs(Q)).get_scalar<Real>();
            }
            else
            {
                Matrix Q    = eye(k,k,vc);
          
                if (sym == false)
                    gemm(1.0f, V, AU, trans_type::conj_trans, trans_type::no_trans, 1.0f, Q);
                else
                    gemm(1.0f, U, AU, trans_type::conj_trans, trans_type::no_trans, 1.0f, Q);

                bool posdef = (sym == false)? false : A.is_posdef();

                if (posdef == true)
                {
                    Q.set_struct(predefined_struct_type::her);
                    Q.add_struct(posdef_flag());
                };

                linsolve_obj Q_mat   = make_linsolve_obj(Q,opt);

                m_Q         = Q_mat;
                m_AUQ       = Q_mat.solve_rev(std::move(AU), trans_type::no_trans);
                m_log_det   = Q_mat.log_det();
            };

            //V'*inv(A) = ((V'*inv(A))')' = (inv(A)'*V)'
            m_AV    = (sym == false) ? A.solve(V, trans_type::conj_trans) : A.solve(U, trans_type::conj_trans);
        };

        //inv(C) + V'*inv(A)*U
        SM_precomp(bool sym, const Matrix& U, const Matrix& V, const linsolve_obj& C, const linsolve_obj& A,
                   const options& opt)
        {
            value_code vc   = matrix_traits::unify_value_types(U.get_value_code(), A.get_value_code());
            vc              = matrix_traits::unify_value_types(vc, value_code::v_float);
            vc              = matrix_traits::unify_value_types(vc, C.get_value_code());

            if (sym == false)
                vc          = matrix_traits::unify_value_types(vc, V.get_value_code());

            //form: Q = inv(C) + V'*inv(A)*U

            Integer k   = U.cols();

            Matrix AU   = A.solve(U, trans_type::no_trans);

            if (k == 1)
            {
                Matrix Q;                
          
                if (sym == false)
                    Q = C.inv() + mmul(V, AU, trans_type::conj_trans, trans_type::no_trans);
                else
                    Q = C.inv() + mmul(U, AU, trans_type::conj_trans, trans_type::no_trans);

                m_Q         = make_linsolve_obj(Q,opt);
                m_AUQ       = AU * (1.0f / Q);
                m_log_det   = matcl::log(matcl::abs(Q)).get_scalar<Real>() + C.log_det();
            }
            else
            {
                Matrix Q    = C.inv();
          
                if (sym == false)
                    gemm(1.0f, V, AU, trans_type::conj_trans, trans_type::no_trans, 1.0f, Q);
                else
                    gemm(1.0f, U, AU, trans_type::conj_trans, trans_type::no_trans, 1.0f, Q);

                bool posdef = (sym == false)? false : A.is_posdef() && C.is_posdef();

                if (posdef == true)
                {
                    Q.set_struct(predefined_struct_type::her);
                    Q.add_struct(posdef_flag());
                };

                linsolve_obj Q_mat   = make_linsolve_obj(Q,opt);

                m_Q         = Q_mat;
                m_AUQ       = Q_mat.solve_rev(std::move(AU), trans_type::no_trans);
                m_log_det   = Q_mat.log_det() + C.log_det();
            };

            //V'*inv(A) = ((V'*inv(A))')' = (inv(A)'*V)'
            m_AV    = (sym == false) ? A.solve(V, trans_type::conj_trans) : A.solve(U, trans_type::conj_trans);
        };

        //inv(I + C * V' * inv(A) * U) * C
        SM_precomp(bool sym, const Matrix& U, const Matrix& V, const Matrix& C, const linsolve_obj& A,
                   const options& opt)
        {
            value_code vc   = matrix_traits::unify_value_types(U.get_value_code(), A.get_value_code());
            vc              = matrix_traits::unify_value_types(vc, value_code::v_float);
            vc              = matrix_traits::unify_value_types(vc, C.get_value_code());

            if (sym == false)
                vc          = matrix_traits::unify_value_types(vc, V.get_value_code());

            //form: Q = I + C * V'*inv(A)*U

            Matrix AU       = A.solve(U, trans_type::no_trans);
            Matrix AUV      = (sym == false)? mmul(V, AU, trans_type::conj_trans)
                                            : mmul(U, AU, trans_type::conj_trans);

            Integer N       = C.rows();
            Matrix Q        = eye(N,N,vc);
          
            gemm(1.0f, C, AUV, trans_type::no_trans, trans_type::no_trans, 1.0f, Q);

            linsolve_obj Q_mat   = make_linsolve_obj(Q,opt);

            m_Q         = Q_mat;
            m_AUQ       = Q_mat.solve_rev(std::move(AU), trans_type::no_trans) * C;
            m_log_det   = Q_mat.log_det();

            //V'*inv(A) = ((V'*inv(A))')' = (inv(A)'*V)'
            m_AV    = (sym == false) ? A.solve(V, trans_type::conj_trans) : A.solve(U, trans_type::conj_trans);
        };

        //inv(D - V' * inv(A) * U)
        SM_precomp(bool sym, const Matrix& U, const Matrix& V, const Matrix& D, const linsolve_obj& A, 
                   const options& opt, ver_block)
        {
            value_code vc   = matrix_traits::unify_value_types(U.get_value_code(), A.get_value_code());
            vc              = matrix_traits::unify_value_types(vc, value_code::v_float);
            vc              = matrix_traits::unify_value_types(vc, D.get_value_code());

            if (sym == false)
                vc          = matrix_traits::unify_value_types(vc, V.get_value_code());

            //form: Q = D - V'*inv(A)*U

            Matrix AU   = A.solve(U, trans_type::no_trans);
            Matrix VAU  = (sym == false)? mmul(V, AU, trans_type::conj_trans)
                                        : mmul(U, AU, trans_type::conj_trans);
            Matrix Q    = D - std::move(VAU);
          
            linsolve_obj Q_mat   = make_linsolve_obj(Q,opt);

            m_Q         = Q_mat;
            m_AUQ       = Q_mat.solve_rev(std::move(AU), trans_type::no_trans);
            m_log_det   = Q_mat.log_det();

            //V'*inv(A) = ((V'*inv(A))')' = (inv(A)'*V)'
            m_AV    = (sym == false) ? A.solve(V, trans_type::conj_trans) : A.solve(U, trans_type::conj_trans);
        };

        SM_precomp convert(value_code vc) const
        {
            Matrix AUQ          = details::convert_value(m_AUQ, vc);
            Matrix AV           = details::convert_value(m_AV, vc);
            linsolve_obj Qc  = m_Q.convert(vc);
            return SM_precomp(AUQ, AV, Qc, m_log_det);
        };

        Real log_det() const
        {
            return m_log_det;
        };

        const matcl::Matrix& get_AUQ() const
        {
            return m_AUQ;
        };
        const matcl::Matrix& get_AV() const
        {
            return m_AV;
        };
        const linsolve_obj& get_Q() const
        {
            return m_Q;
        };

        //form AUQ * X
        Matrix mult_AUQ(const Matrix& X, trans_type t) const
        {
            return mmul(m_AUQ, X, t);
        };

        //form X*AUQ
        Matrix mult_AUQ_rev(const Matrix& X, trans_type t) const
        {
            return mmul(X, m_AUQ, trans_type::no_trans, t);
        };

        //form AV' * X
        Matrix mult_AV(const Matrix& X, trans_type t) const
        {
            trans_type_ext te   = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(t), 
                                    trans_type_ext::conj_trans);

            return mmul(m_AV, X, te, trans_type_ext::no_trans);
        };

        //form X * op(AV')
        Matrix mult_AV_rev(const Matrix& X, trans_type t) const
        {
            trans_type_ext te   = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(t), 
                                    trans_type_ext::conj_trans);

            return mmul(X, m_AV, trans_type_ext::no_trans, te);
        };

        //form inv(Q) * X
        Matrix mult_Q(const Matrix& X, trans_type t) const
        {
            return m_Q.solve(X, t);
        };

        //form X*inv(Q)
        Matrix mult_Q_rev(const Matrix& X, trans_type t) const
        {
            return m_Q.solve_rev(X, t);
        };

        //form op(inv(A) * U * inv(Q) * V' * inv(A)) * X = op(AUQ * AV') * X

        Matrix mult(const Matrix& X, trans_type t) const
        {
            if (t == trans_type::no_trans)
                return m_AUQ * mmul(m_AV, X, trans_type::conj_trans);

            //op(AUQ * AV') * X = op(AV') * op(AUQ) * X
            Matrix y = mmul(m_AUQ, X, t);

            trans_type_ext te   = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(t), 
                                    trans_type_ext::conj_trans);

            return mmul(m_AUQ, std::move(y), te, trans_type_ext::no_trans);
        };

        //form X * op(inv(A) * U * inv(Q) * V' * inv(A)) = X * op(AUQ * AV')
        Matrix mult_rev(const Matrix& X, trans_type t) const
        {
            if (t == trans_type::no_trans)
            {
                Matrix y = mmul(X, m_AUQ, trans_type::no_trans, trans_type::no_trans);
                return  mmul(std::move(y), m_AV, trans_type::no_trans, trans_type::conj_trans);
            };

            //X * op(AUQ * AV') = X * op(AV') * op(AUQ)
            
            trans_type_ext te   = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(t), 
                                    trans_type_ext::conj_trans);

            Matrix y    = mmul(X, m_AV, trans_type_ext::no_trans, te);
            y           = mmul(std::move(y), m_AUQ, trans_type::no_trans, t);
            return y;
        };

    private:
        SM_precomp(const Matrix& AUQ, const Matrix& AV, const linsolve_obj& Qc, Real log_det)
            :m_AUQ(AUQ), m_AV(AV), m_log_det(log_det), m_Q(Qc)
        {};
};

// inv(A + U*V')  = inv(A) - inv(A) * U * inv( I + V'*inv(A)*U) ) * V' * inv(A)
class linsolve_update_SM : public linsolve_obj_data
{
    private:
        linsolve_obj    m_A;
        Matrix          m_U;
        Matrix          m_UV;
        value_code      m_vc;
        bool            m_sym;

        SM_precomp      m_precomp;

    public:
        linsolve_update_SM(const linsolve_obj& A, const Matrix& U, const Matrix& V, const options& opt,
                           bool& all_finite)
            :m_A(A), m_U(U), m_UV(V), m_sym(false), m_precomp(false, U, V, A, opt)
        {
            unify_check(false);
            all_finite = m_precomp.all_finite();
        };
        linsolve_update_SM(const linsolve_obj& A, const Matrix& U, const options& opt,bool& all_finite)
            :m_A(A), m_U(U), m_sym(true), m_precomp(true, U, U, A, opt)
        {
            unify_check(true);
            all_finite = m_precomp.all_finite();
        };

        virtual ~linsolve_update_SM(){};

        virtual bool all_finite() const override
        { 
            return true;
        };

        virtual bool is_modified() const override
        {
            return m_A.is_modified() || m_precomp.is_modified();
        };

        virtual bool is_direct() const override
        {
            return m_A.is_direct() && m_precomp.is_direct();
        };

        virtual Integer rows() const override
        {
            return m_A.rows();
        };

        virtual Integer cols() const override
        {
            return m_A.cols();
        };

        virtual value_code get_value_code() const override
        {
            return m_vc;
        };

        virtual ti::ti_object get_type() const override
        {
            return m_A.get_type();
        };

        virtual bool is_hermitian() const override
        { 
            return m_sym == false? false : m_A.is_hermitian();
        };

        virtual bool is_posdef() const override
        {
            return m_sym == false? false : m_A.is_posdef();
        }
        virtual data_ptr convert(value_code vc) const override
        {
            linsolve_obj Ac     = m_A.convert(vc);
            Matrix Uc           = details::convert_value(m_U, vc);
            Matrix UVc          = (m_sym == false) ? details::convert_value(m_UV, vc) : m_UV;
            SM_precomp pc       = m_precomp.convert(vc);

            return data_ptr(new linsolve_update_SM(Ac, Uc, UVc, pc, m_sym, vc) );
        };        

        virtual Real log_det() const override
        {
            Real det_A  = m_A.log_det();
            Real det_Q  = m_precomp.log_det();

            return det_A + det_Q;
        };

        virtual matcl::Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            // op(A + U*V') * Y = op(A)*Y + op(U*V') * Y
            Matrix y1       = m_A.mmul_right(X, t);

            if (t == trans_type::no_trans)
            {
                Matrix y2   = mmul(get_V(), X, trans_type::conj_trans);
                y2          = mmul(m_U, std::move(y2), trans_type::no_trans);

                return std::move(y1) + std::move(y2);
            };

            // op(U*V') * Y = op(V') * op(U) * Y
            Matrix y2       = mmul(m_U, X, t);

            trans_type_ext te   = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(t), 
                                    trans_type_ext::conj_trans);

            y2              = mmul(get_V(), std::move(y2), te, trans_type_ext::no_trans);

            return std::move(y1) + std::move(y2);
        };        

        virtual matcl::Matrix mmul_left(const Matrix& X, trans_type t) const override
        {
            // X * op(A + U*V') = X * op(A) + X * op(U*V')
            Matrix y1       = m_A.mmul_left(X, t);

            if (t == trans_type::no_trans)
            {
                Matrix y2   = mmul(X, m_U, trans_type::no_trans);
                y2          = mmul(std::move(y2), get_V(), trans_type::no_trans, trans_type::conj_trans);                

                return std::move(y1) + std::move(y2);
            };

            // X * op(U*V') = X * op(V') * op(U)            
            trans_type_ext te   = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(t), 
                                    trans_type_ext::conj_trans);

            Matrix y2       = mmul(X, get_V(), trans_type_ext::no_trans, te);
            y2              = mmul(std::move(y2), m_U, trans_type::no_trans, t);

            return std::move(y1) + std::move(y2);
        };     

        virtual matcl::Matrix solve(const Matrix& X, trans_type tA) const override
        {
            // op(A + U*V') * Y = X => Y = op(inv(A + U*V')) * X
            // op(inv(A + U*V')) = op(inv(A)) - op(AUQ * AV')

            Matrix y2       = m_precomp.mult(X, tA);
            Matrix y1       = m_A.solve(X, tA);            
            return std::move(y1) - std::move(y2);
        };
        virtual matcl::Matrix solve(Matrix&& X, trans_type tA) const override
        {
            // op(A + U*V') * Y = X => Y = op(inv(A + U*V')) * X
            // op(inv(A + U*V')) = op(inv(A)) - op(AUQ * AV')

            Matrix y2       = m_precomp.mult(X, tA);
            Matrix y1       = m_A.solve(std::move(X), tA);            
            return std::move(y1) - std::move(y2);
        };

        virtual matcl::Matrix solve_rev(const Matrix& X, trans_type tA) const override
        {
            // Y * op(A + U*V') = X => Y = X * op(inv(A + U*V'))
            // op(inv(A + U*V')) = op(inv(A)) - op(AUQ * AV')

            Matrix y2       = m_precomp.mult_rev(X, tA);
            Matrix y1       = m_A.solve_rev(X, tA);
            return std::move(y1) - std::move(y2);
        };

        virtual matcl::Matrix solve_rev(Matrix&& X, trans_type tA) const override
        {
            // Y * op(A + U*V') = X => Y = X * op(inv(A + U*V'))
            // op(inv(A + U*V')) = op(inv(A)) - op(AUQ * AV')

            Matrix y2       = m_precomp.mult_rev(X, tA);
            Matrix y1       = m_A.solve_rev(std::move(X), tA);
            return std::move(y1) - std::move(y2);
        };

        virtual matcl::Matrix inv() const override
        {
            //inv(A + U*V')    = inv(A) - AUQ * AV'
            Matrix y        = m_A.inv() - mmul(m_precomp.get_AUQ(), m_precomp.get_AV(), trans_type::no_trans,
                                               trans_type::conj_trans);

            if (is_hermitian() == true)
            {
                y.add_struct(predefined_struct_type::her);

                if (is_posdef() == true)
                    y.add_struct(posdef_flag());
            };

            return y;
        }

        virtual matcl::Matrix base_matrix() const
        {
            // A + U*V'
            Matrix A        = m_A.base_matrix() + m_U * get_V();
            return A;
        };        

        //not reimplemented
        //virtual Real normest_1() const override;
        //virtual Real normest_inf() const override
        //virtual matcl::Matrix mmul_right(Matrix&& X, trans_type t) const override;
        //virtual Real normest_2() const override;
        //virtual Real mat_normest_2() const override;

    private:
        //const from data
        linsolve_update_SM(const linsolve_obj& A, const Matrix& U, const Matrix& V, 
                           const SM_precomp& p, bool sym, value_code vc)
            : m_A(A), m_U(U), m_UV(V), m_precomp(p), m_sym(sym), m_vc(vc)
        {};

        const Matrix& get_V() const
        {
            return m_sym == false ? m_UV : m_U;
        };

        void unify_check(bool sym)
        {
            value_code vA   = m_A.get_value_code();
            value_code vU   = m_U.get_value_code();
            value_code vV   = m_UV.get_value_code();
            m_vc            = matrix_traits::unify_value_types(vA, vU);
            m_vc            = matrix_traits::unify_value_types(m_vc, value_code::v_float);

            if (sym == false)
                m_vc        = matrix_traits::unify_value_types(m_vc, vV);
            
            if (sym == false)
            {
                error::check_mul(m_U.rows(), m_U.cols(), m_UV.rows(), m_UV.cols(), trans_type::no_trans, 
                                 trans_type::conj_trans);

                error::check_eeop(m_A.rows(), m_A.cols(), m_U.rows(), m_UV.rows());
            }
            else
            {
                error::check_eeop(m_A.rows(), m_A.cols(), m_U.rows(), m_U.rows());
            };

            if (vA != m_vc)
                m_A         = m_A.convert(m_vc);

            if (vU != m_vc)
                m_U         = details::convert_value(m_U, m_vc);

            if (sym == false && vV != m_vc)
                m_UV        = details::convert_value(m_UV, m_vc);
        };
};

//------------------------------------------------------------------
//                      linsolve_update_W
//------------------------------------------------------------------
// inv(A + U*C*V')  = inv(A) - inv(A) * U * inv( inv(C) + V'*inv(A)*U) ) * V' * inv(A)
class linsolve_update_W : public linsolve_obj_data
{
    private:
        linsolve_obj    m_A;
        Matrix          m_U;
        linsolve_obj    m_C;
        Matrix          m_UV;
        value_code      m_vc;
        bool            m_sym;
        SM_precomp      m_precomp;

    public:
        linsolve_update_W(const linsolve_obj& A, const Matrix& U, const linsolve_obj& C, const Matrix& V,
                          const options& opt, bool& all_finite)
            :m_A(A), m_U(U), m_C(C), m_UV(V), m_sym(false), m_precomp(false, U, V, C, A, opt)
        {
            unify_check(false);
            all_finite = m_precomp.all_finite();
        };
        linsolve_update_W(const linsolve_obj& A, const Matrix& U, const linsolve_obj& C, const options& opt,
                          bool& all_finite)
            :m_A(A), m_U(U), m_C(C), m_sym(true), m_precomp(true, U, U, C, A, opt)
        {
            unify_check(true);
            all_finite = m_precomp.all_finite();
        };

        virtual ~linsolve_update_W(){};

        virtual bool all_finite() const override
        { 
            return true;
        };

        virtual bool is_modified() const override
        {
            return m_A.is_modified() || m_C.is_modified() || m_precomp.is_modified();
        };

        virtual bool is_direct() const override
        {
            return m_A.is_direct() && m_C.is_direct() && m_precomp.is_direct();
        };

        virtual Integer rows() const override
        {
            return m_A.rows();
        };

        virtual Integer cols() const override
        {
            return m_A.cols();
        };

        virtual value_code get_value_code() const override
        {
            return m_vc;
        };

        virtual ti::ti_object get_type() const override
        {
            return m_A.get_type();
        };

        virtual bool is_hermitian() const override
        { 
            return m_sym == false? false : m_A.is_hermitian() && m_C.is_hermitian();
        };

        virtual bool is_posdef() const override
        {
            return m_sym == false? false : m_A.is_posdef() && m_C.is_posdef();
        }

        virtual data_ptr convert(value_code vc) const override
        {
            linsolve_obj Ac     = m_A.convert(vc);
            linsolve_obj Cc     = m_C.convert(vc);
            Matrix Uc           = details::convert_value(m_U, vc);
            Matrix UVc          = (m_sym == false) ? details::convert_value(m_UV, vc) : m_UV;
            SM_precomp pc       = m_precomp.convert(vc);

            return data_ptr(new linsolve_update_W(Ac, Uc, UVc, Cc, pc, m_sym, vc) );
        };        

        virtual Real log_det() const override
        {
            Real det_A  = m_A.log_det();
            Real det_Q  = m_precomp.log_det();

            return det_A + det_Q;
        };

        virtual matcl::Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            // op(A + U*C*V') * Y = op(A)*Y + op(U*C*V') * Y
            Matrix y1       = m_A.mmul_right(X, t);

            if (t == trans_type::no_trans)
            {
                Matrix y2   = mmul(get_V(), X, trans_type::conj_trans);
                y2          = m_C.mmul_right(std::move(y2), trans_type::no_trans);
                y2          = mmul(m_U, std::move(y2), trans_type::no_trans);

                return std::move(y1) + std::move(y2);
            };

            // op(U*C*V') * Y = op(V') * op(C) * op(U) * Y
            Matrix y2       = mmul(m_U, X, t);
            y2              = m_C.mmul_right(std::move(y2), t);

            trans_type_ext te   = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(t), 
                                    trans_type_ext::conj_trans);

            y2              = mmul(get_V(), std::move(y2), te, trans_type_ext::no_trans);

            return std::move(y1) + std::move(y2);
        };        

        virtual matcl::Matrix mmul_left(const Matrix& X, trans_type t) const override
        {
            // X * op(A + U*C*V') = X * op(A) + X * op(U*C*V')
            Matrix y1       = m_A.mmul_left(X, t);

            if (t == trans_type::no_trans)
            {
                Matrix y2   = mmul(X, m_U, trans_type::no_trans, trans_type::no_trans);
                y2          = m_C.mmul_left(std::move(y2), trans_type::no_trans);
                y2          = mmul(std::move(y2), get_V(), trans_type::no_trans, trans_type::conj_trans);                

                return std::move(y1) + std::move(y2);
            };

            // X * op(U*C*V') = X * op(V') * op(C)*op(U)            
            trans_type_ext te   = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(t), 
                                    trans_type_ext::conj_trans);

            Matrix y2       = mmul(X, get_V(), trans_type_ext::no_trans, te);

            y2              = m_C.mmul_left(std::move(y2), t);
            y2              = mmul(std::move(y2), m_U, trans_type::no_trans, t);

            return std::move(y1) + std::move(y2);
        };     

        virtual matcl::Matrix solve(const Matrix& X, trans_type tA) const override
        {
            // op(A + U*C*V') * Y = X => Y = op(inv(A + U*C*V')) * X
            // op(inv(A + U*C*V')) = op(inv(A)) - op(AUQ * AV')

            Matrix y2       = m_precomp.mult(X, tA);
            Matrix y1       = m_A.solve(X, tA);            
            return std::move(y1) - std::move(y2);
        };

        virtual matcl::Matrix solve(Matrix&& X, trans_type tA) const override
        {
            // op(A + U*C*V') * Y = X => Y = op(inv(A + U*C*V')) * X
            // op(inv(A + U*C*V')) = op(inv(A)) - op(AUQ * AV')

            Matrix y2       = m_precomp.mult(X, tA);
            Matrix y1       = m_A.solve(std::move(X), tA);            
            return std::move(y1) - std::move(y2);
        };

        virtual matcl::Matrix solve_rev(const Matrix& X, trans_type tA) const override
        {
            // Y * op(A + U*C*V') = X => Y = X * op(inv(A + U*C*V'))
            // op(inv(A + U*C*V')) = op(inv(A)) - op(AUQ * AV')

            Matrix y2       = m_precomp.mult_rev(X, tA);
            Matrix y1       = m_A.solve_rev(X, tA);
            return std::move(y1) - std::move(y2);
        };

        virtual matcl::Matrix solve_rev(Matrix&& X, trans_type tA) const override
        {
            // Y * op(A + U*C*V') = X => Y = X * op(inv(A + U*C*V'))
            // op(inv(A + U*C*V')) = op(inv(A)) - op(AUQ * AV')

            Matrix y2       = m_precomp.mult_rev(X, tA);
            Matrix y1       = m_A.solve_rev(std::move(X), tA);
            return std::move(y1) - std::move(y2);
        };

        virtual matcl::Matrix inv() const override
        {
            //inv(A + U*C*V')   = inv(A) - AUQ * AV'
            Matrix y        = m_A.inv() - mmul(m_precomp.get_AUQ(), m_precomp.get_AV(), trans_type::no_trans,
                                               trans_type::conj_trans);

            if (is_hermitian() == true)
            {
                y.add_struct(predefined_struct_type::her);

                if (is_posdef() == true)
                    y.add_struct(posdef_flag());
            };

            return y;
        }

        virtual matcl::Matrix base_matrix() const
        {
            // A + U*V'
            Matrix A        = m_A.base_matrix() + chain_mult(m_U, m_C.base_matrix(), get_V());
            return A;
        };        

        //not reimplemented
        //virtual Real normest_1() const override;
        //virtual Real normest_inf() const override
        //virtual matcl::Matrix mmul_right(Matrix&& X, trans_type t) const override;
        //virtual Real normest_2() const override;
        //virtual Real mat_normest_2() const override;

    private:
        linsolve_update_W(const linsolve_obj& Ac, const Matrix& Uc, const Matrix& UVc, 
                          const linsolve_obj& Cc, const SM_precomp& pc, bool sym, value_code vc)
            : m_A(Ac), m_U(Uc), m_C(Cc), m_UV(UVc), m_precomp(pc), m_sym(sym), m_vc(vc)
        {};

        const Matrix& get_V() const
        {
            return m_sym == false ? m_UV : m_U;
        };

        void unify_check(bool sym)
        {
            value_code vA   = m_A.get_value_code();
            value_code vU   = m_U.get_value_code();
            value_code vC   = m_C.get_value_code();
            value_code vV   = m_UV.get_value_code();
            m_vc            = matrix_traits::unify_value_types(vA, vU);
            m_vc            = matrix_traits::unify_value_types(m_vc, vC);
            m_vc            = matrix_traits::unify_value_types(m_vc, value_code::v_float);

            if (sym == false)
                m_vc        = matrix_traits::unify_value_types(m_vc, vV);
            
            error::check_mul(m_U.rows(), m_U.cols(), m_C.rows(), m_C.cols(), trans_type::no_trans, 
                                trans_type::no_trans);

            if (sym == false)
            {
                error::check_mul(m_C.rows(), m_C.cols(), m_UV.rows(), m_UV.cols(), trans_type::no_trans, 
                                trans_type::conj_trans);
                error::check_eeop(m_A.rows(), m_A.cols(), m_U.rows(), m_UV.rows());
            }
            else
            {
                error::check_mul(m_C.rows(), m_C.cols(), m_U.rows(), m_U.cols(), trans_type::no_trans, 
                                trans_type::conj_trans);

                error::check_eeop(m_A.rows(), m_A.cols(), m_U.rows(), m_U.rows());
            };

            if (vA != m_vc)
                m_A         = m_A.convert(m_vc);

            if (vU != m_vc)
                m_U         = details::convert_value(m_U, m_vc);

            if (vC != m_vc)
                m_C         = m_C.convert(m_vc);

            if (sym == false && vV != m_vc)
                m_UV        = details::convert_value(m_UV, m_vc);
        };
};

//------------------------------------------------------------------
//                      linsolve_update_BI
//------------------------------------------------------------------
// inv(A + U*C*V')  = inv(A) - inv(A) * U * inv( I + C*V'*inv(A)*U) ) * C*V' * inv(A)
class linsolve_update_BI : public linsolve_obj_data
{
    private:
        linsolve_obj    m_A;
        Matrix          m_U;
        Matrix          m_C;
        Matrix          m_UV;
        value_code      m_vc;
        bool            m_sym;
        SM_precomp      m_precomp;

    public:
        linsolve_update_BI(const linsolve_obj& A, const Matrix& U, const Matrix& C, const Matrix& V,
                           const options& opt,bool& all_finite)
            :m_A(A), m_U(U), m_C(C), m_UV(V), m_sym(false), m_precomp(false, U, V, C, A, opt)
        {
            unify_check(false);
            all_finite = m_precomp.all_finite();
        };

        linsolve_update_BI(const linsolve_obj& A, const Matrix& U, const Matrix& C, const options& opt,
                           bool& all_finite)
            :m_A(A), m_U(U), m_C(C), m_sym(true), m_precomp(true, U, U, C, A, opt)
        {
            unify_check(true);
            all_finite = m_precomp.all_finite();
        };

        virtual ~linsolve_update_BI(){};

        virtual bool all_finite() const override
        { 
            return true;
        };

        virtual bool is_modified() const override
        {
            return m_A.is_modified() || m_precomp.is_modified();
        };

        virtual bool is_direct() const override
        {
            return m_A.is_direct() && m_precomp.is_direct();
        };

        virtual Integer rows() const override
        {
            return m_A.rows();
        };

        virtual Integer cols() const override
        {
            return m_A.cols();
        };

        virtual value_code get_value_code() const override
        {
            return m_vc;
        };

        virtual ti::ti_object get_type() const override
        {
            return m_A.get_type();
        };

        virtual bool is_hermitian() const override
        { 
            if (m_sym == false)
            {
                return false;
            }
            else
            {
                bool is_real    = matrix_traits::is_float_real(m_vc);
                bool is_herm_C  = m_C.get_struct().is_hermitian(m_C.is_square(), is_real);
                return m_A.is_hermitian() && is_herm_C;
            };
        };
        virtual bool is_posdef() const override
        { 
            if (m_sym == false)
            {
                return false;
            }
            else
            {
                bool is_real    = matrix_traits::is_float_real(m_vc);
                bool is_herm_C  = m_C.get_struct().is_hermitian(m_C.is_square(), is_real);
                bool is_posdef_C= matcl::is_posdef(m_C.get_struct());

                return m_A.is_posdef() && is_herm_C && is_posdef_C;
            };
        };

        virtual data_ptr convert(value_code vc) const override
        {
            linsolve_obj Ac     = m_A.convert(vc);
            Matrix Cc           = details::convert_value(m_C,vc);
            Matrix Uc           = details::convert_value(m_U, vc);
            Matrix UVc          = (m_sym == false) ? details::convert_value(m_UV, vc) : m_UV;
            SM_precomp pc       = m_precomp.convert(vc);

            return data_ptr(new linsolve_update_BI(Ac, Uc, UVc, Cc, pc, m_sym, vc) );
        };        

        virtual Real log_det() const override
        {
            Real det_A  = m_A.log_det();
            Real det_Q  = m_precomp.log_det();

            return det_A + det_Q;
        };

        virtual matcl::Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            // op(A + U*C*V') * Y = op(A)*Y + op(U*C*V') * Y
            Matrix y1       = m_A.mmul_right(X, t);

            if (t == trans_type::no_trans)
            {
                Matrix y2   = mmul(get_V(), X, trans_type::conj_trans);
                y2          = mmul(m_C, std::move(y2), trans_type::no_trans);
                y2          = mmul(m_U, std::move(y2), trans_type::no_trans);

                return std::move(y1) + std::move(y2);
            };

            // op(U*C*V') * Y = op(V') * op(C) * op(U) * Y
            Matrix y2       = mmul(m_U, X, t);
            y2              = mmul(m_C, std::move(y2), t);

            trans_type_ext te   = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(t), 
                                    trans_type_ext::conj_trans);

            y2              = mmul(get_V(), std::move(y2), te, trans_type_ext::no_trans);

            return std::move(y1) + std::move(y2);
        };        

        virtual matcl::Matrix mmul_left(const Matrix& X, trans_type t) const override
        {
            // X * op(A + U*C*V') = X * op(A) + X * op(U*C*V')
            Matrix y1       = m_A.mmul_left(X, t);

            if (t == trans_type::no_trans)
            {
                Matrix y2   = mmul(X, m_U, trans_type::no_trans, trans_type::no_trans);
                y2          = mmul(std::move(y2), m_C, trans_type::no_trans, trans_type::no_trans);
                y2          = mmul(std::move(y2), get_V(), trans_type::no_trans, trans_type::conj_trans);                

                return std::move(y1) + std::move(y2);
            };

            // X * op(U*C*V') = X * op(V') * op(C)*op(U)            
            trans_type_ext te   = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(t), 
                                    trans_type_ext::conj_trans);

            Matrix y2           = mmul(X, get_V(), trans_type_ext::no_trans, te);
            y2                  = mmul(std::move(y2), m_C, trans_type::no_trans, t);
            y2                  = mmul(std::move(y2), m_U, trans_type::no_trans, t);

            return std::move(y1) + std::move(y2);
        };   

        virtual matcl::Matrix solve(const Matrix& X, trans_type tA) const override
        {
            // op(A + U*C*V') * Y = X => Y = op(inv(A + U*C*V')) * X
            // op(inv(A + U*C*V')) = op(inv(A)) - op(AUQ * AV')

            Matrix y2       = m_precomp.mult(X, tA);
            Matrix y1       = m_A.solve(X, tA);            
            return std::move(y1) - std::move(y2);
        };

        virtual matcl::Matrix solve(Matrix&& X, trans_type tA) const override
        {
            // op(A + U*C*V') * Y = X => Y = op(inv(A + U*C*V')) * X
            // op(inv(A + U*C*V')) = op(inv(A)) - op(AUQ * AV')

            Matrix y2       = m_precomp.mult(X, tA);
            Matrix y1       = m_A.solve(std::move(X), tA);            
            return std::move(y1) - std::move(y2);
        };

        virtual matcl::Matrix solve_rev(const Matrix& X, trans_type tA) const override
        {
            // Y * op(A + U*C*V') = X => Y = X * op(inv(A + U*C*V'))
            // op(inv(A + U*C*V')) = op(inv(A)) - op(AUQ * AV')

            Matrix y2       = m_precomp.mult_rev(X, tA);
            Matrix y1       = m_A.solve_rev(X, tA);
            return std::move(y1) - std::move(y2);
        };

        virtual matcl::Matrix solve_rev(Matrix&& X, trans_type tA) const override
        {
            // Y * op(A + U*C*V') = X => Y = X * op(inv(A + U*C*V'))
            // op(inv(A + U*C*V')) = op(inv(A)) - op(AUQ * AV')

            Matrix y2       = m_precomp.mult_rev(X, tA);
            Matrix y1       = m_A.solve_rev(std::move(X), tA);
            return std::move(y1) - std::move(y2);
        };

        virtual matcl::Matrix inv() const override
        {
            //inv(A + U*C*V')   = inv(A) - AUQ * AV'
            Matrix y        = m_A.inv() - mmul(m_precomp.get_AUQ(), m_precomp.get_AV(), trans_type::no_trans,
                                               trans_type::conj_trans);

            if (is_hermitian() == true)
            {
                y.add_struct(predefined_struct_type::her);

                if (is_posdef() == true)
                    y.add_struct(posdef_flag());
            };

            return y;
        }

        virtual matcl::Matrix base_matrix() const
        {
            // A + U*C*V'
            Matrix A        = m_A.base_matrix() + chain_mult(m_U, m_C, get_V());
            return A;
        };        

        //not reimplemented
        //virtual Real normest_1() const override;
        //virtual Real normest_inf() const override
        //virtual matcl::Matrix mmul_right(Matrix&& X, trans_type t) const override;
        //virtual Real normest_2() const override;
        //virtual Real mat_normest_2() const override;

    private:
        linsolve_update_BI(const linsolve_obj& Ac, const Matrix& Uc, const Matrix& UVc, 
                          const Matrix& Cc, const SM_precomp& pc, bool sym, value_code vc)
            : m_A(Ac), m_U(Uc), m_C(Cc), m_UV(UVc), m_precomp(pc), m_sym(sym), m_vc(vc)
        {};

        const Matrix& get_V() const
        {
            return m_sym == false ? m_UV : m_U;
        };

        void unify_check(bool sym)
        {
            value_code vA   = m_A.get_value_code();
            value_code vU   = m_U.get_value_code();
            value_code vC   = m_C.get_value_code();
            value_code vV   = m_UV.get_value_code();
            m_vc            = matrix_traits::unify_value_types(vA, vU);
            m_vc            = matrix_traits::unify_value_types(m_vc, vC);
            m_vc            = matrix_traits::unify_value_types(m_vc, value_code::v_float);

            if (sym == false)
                m_vc        = matrix_traits::unify_value_types(m_vc, vV);
            
            error::check_mul(m_U.rows(), m_U.cols(), m_C.rows(), m_C.cols(), trans_type::no_trans, 
                                trans_type::no_trans);

            if (sym == false)
            {
                error::check_mul(m_C.rows(), m_C.cols(), m_UV.rows(), m_UV.cols(), trans_type::no_trans, 
                                trans_type::conj_trans);
                error::check_eeop(m_A.rows(), m_A.cols(), m_U.rows(), m_UV.rows());
            }
            else
            {
                error::check_mul(m_C.rows(), m_C.cols(), m_U.rows(), m_U.cols(), trans_type::no_trans, 
                                trans_type::conj_trans);

                error::check_eeop(m_A.rows(), m_A.cols(), m_U.rows(), m_U.rows());
            };

            if (vA != m_vc)
                m_A         = m_A.convert(m_vc);

            if (vU != m_vc)
                m_U         = details::convert_value(m_U, m_vc);

            if (vC != m_vc)
                m_C         = details::convert_value(m_C, m_vc);

            if (sym == false && vV != m_vc)
                m_UV        = details::convert_value(m_UV, m_vc);
        };
};

//------------------------------------------------------------------
//                      linsolve_insert_A
//------------------------------------------------------------------
class linsolve_insert_A : public linsolve_obj_data
{
    private:
        linsolve_obj    m_A;
        Matrix          m_B;
        Matrix          m_C;
        Matrix          m_D;
        value_code      m_vc;
        bool            m_sym;
        SM_precomp      m_precomp;

    public:
        linsolve_insert_A(const linsolve_obj& A, const Matrix& B, const Matrix& C, const Matrix& D,
                          const options& opt, bool& all_finite)
            :m_A(A), m_B(B), m_C(C), m_D(D), m_sym(false)
            , m_precomp(false, B, C, D, A, opt, SM_precomp::ver_block())
        {
            unify_check(false);
            all_finite = m_precomp.all_finite();
        };
        linsolve_insert_A(const linsolve_obj& A, const Matrix& B, const Matrix& D, const options& opt,
                          bool& all_finite)
            :m_A(A), m_B(B), m_D(D), m_sym(true)
            , m_precomp(false, B, B, D, A, opt, SM_precomp::ver_block())
        {
            unify_check(true);
            all_finite = m_precomp.all_finite();
        };

        virtual ~linsolve_insert_A(){};

        virtual bool all_finite() const override
        { 
            return true;
        };

        virtual bool is_modified() const override
        {
            return m_A.is_modified() || m_precomp.is_modified();
        };

        virtual bool is_direct() const override
        {
            return m_A.is_direct() && m_precomp.is_direct();
        };

        virtual Integer rows() const override
        {
            return m_A.rows() + m_D.rows();
        };

        virtual Integer cols() const override
        {
            return m_A.cols() + m_D.cols();
        };

        virtual value_code get_value_code() const override
        {
            return m_vc;
        };

        virtual ti::ti_object get_type() const override
        {
            return m_A.get_type();
        };

        virtual bool is_hermitian() const override
        { 
            if (m_sym == false)
            {
                return false;
            }
            else
            {
                bool is_real    = matrix_traits::is_float_real(m_vc);
                bool is_herm_D  = m_D.get_struct().is_hermitian(m_D.is_square(), is_real);
                return m_A.is_hermitian() && is_herm_D;
            };
        };

        virtual bool is_posdef() const override
        {
            return false;
        };

        virtual data_ptr convert(value_code vc) const override
        {
            linsolve_obj Ac     = m_A.convert(vc);
            Matrix Bc           = details::convert_value(m_B, vc);
            Matrix Cc           = (m_sym == false) ? details::convert_value(m_C, vc) : m_C;
            Matrix Dc           = details::convert_value(m_D, vc);
            SM_precomp pc       = m_precomp.convert(vc);

            return data_ptr(new linsolve_insert_A(Ac, Bc, Cc, Dc, pc, m_sym, vc) );
        };        

        virtual Real log_det() const override
        {
            Real det_A  = m_A.log_det();
            Real det_Q  = m_precomp.log_det();

            return det_A + det_Q;
        };

        virtual matcl::Matrix solve(const Matrix& X, trans_type tA) const override
        {
            Integer N1      = m_A.cols();
            Integer N2      = m_B.cols();

            Matrix X1       = X(colon(1, N1), colon());
            Matrix X2       = X(colon(N1+1, N1+N2), colon());

            if (tA == trans_type::no_trans)
            {
                // A_11 * Y = X => Y = inv(A_11) * X = inv(A)*X1 + AUQ * AV'*X1
                Matrix Z1   = m_precomp.mult_AV(X1, tA);
                Matrix Y11  = m_A.solve(X1, tA) + m_precomp.mult_AUQ(Z1, tA);

                // A_12 * Y = X => Y = inv(A_12) * X = -AUQ *X2
                Matrix Y12          = -m_precomp.mult_AUQ(X2, tA);

                // A_21 * X = X => Y = inv(A_21) * X = -inv(Q) * AV' *X1 = -inv(Q) * Z1
                Matrix Y21  = -m_precomp.mult_Q(Z1, tA);

                // A_22 * Y = X => Y = inv(A_22) * X = inv(Q) * X2
                Matrix Y22  = m_precomp.mult_Q(X2, tA);

                Matrix ret  = (mat_col(), Y11 + Y12, Y21 + Y22);
                return ret;
            }
            else
            {
                Matrix Z1   = m_precomp.mult_AUQ(X1,tA);
                Matrix Z2   = m_precomp.mult_Q(X2,tA);

                // op(A_11)*Y = X => Y = inv(op(A_11))*X = inv(op(A))*X1 + op(AV') * op(AUQ)*X1
                Matrix Y11  = m_A.solve_rev(X1, tA) + m_precomp.mult_AV(Z1, tA);

                // op(A_21)*Y = X => Y = inv(op(A_21))*X = -op(AV')*inv(op(Q))*X2
                Matrix Y12  = -m_precomp.mult_AV(Z2, tA);

                // op(A_12)*Y = X => Y = inv(op(A_12))*X = -op(AUQ)*X1
                Matrix Y21  = -Z1;

                // op(A_22)*Y = X => Y = inv(op(A_22))*X = inv(op(Q))*X
                Matrix Y22  = Z2;

                Matrix ret  = (mat_col(), Y11 + Y12, Y21 + Y22);
                return ret;
            };
        };

        virtual matcl::Matrix solve(Matrix&& X, trans_type tA) const override
        {
            return solve(X, tA);
        };

        virtual matcl::Matrix solve_rev(const Matrix& X, trans_type tA) const override
        {
            Integer N1      = m_A.rows();
            Integer N2      = m_D.rows();

            Matrix X1       = X(colon(), colon(1, N1));
            Matrix X2       = X(colon(), colon(N1+1, N1+N2));

            if (tA == trans_type::no_trans)
            {
                Matrix Z1   = m_precomp.mult_AUQ_rev(X1, tA);
                Matrix Z2   = m_precomp.mult_Q_rev(X2, tA);

                // Y * A_11 = X => Y = X*inv(A_11) = X1*inv(A) + X1*AUQ * AV                
                Matrix Y11  = m_A.solve(X1, tA) + m_precomp.mult_AV_rev(Z1, tA);

                // Y * A_12 = X => Y = X*inv(A_12) = -X1*AUQ = -Z1
                Matrix Y12  = -std::move(Z1);

                // Y * A_21 = X => Y = X*inv(A_21) = -X2 * inv(Q) * AV'               
                Matrix Y21  = -m_precomp.mult_AV_rev(Z2, tA);

                // Y * A_22 = X => Y = X * inv(A_22) = X2 * inv(Q) = Z2
                Matrix Y22  = Z2;

                Matrix ret  = (mat_row(), Y11 + Y21, Y12 + Y22);
                return ret;
            }
            else
            {
                Matrix Z1   = m_precomp.mult_AV_rev(X1, tA);

                // Y * op(A_11) = X => Y = X * inv(op(A_11)) = X1 * inv(op(A)) + X1*op(AV') * op(AUQ)
                Matrix Y11  = m_A.solve_rev(X1, tA) + m_precomp.mult_AUQ_rev(Z1, tA);

                // Y * op(A_12) = X => Y = X * inv(op(A_12)) = -X2 * op(AUQ)
                Matrix Y21  = -m_precomp.mult_AUQ_rev(X2, tA);

                // Y * op(A_21) = X => Y = X*inv(op(A_21)) = -X1*op(AV')*inv(op(Q))
                Matrix Y12  = -m_precomp.mult_Q_rev(std::move(Z1),tA);

                // Y * op(A_22) = X => Y = X*inv(op(A_22)) = X2 * inv(op(Q))
                Matrix Y22  = m_precomp.mult_Q_rev(X2, tA);

                Matrix ret  = (mat_row(), Y11 + Y21, Y12 + Y22);
                return ret;
            };
        };

        virtual matcl::Matrix solve_rev(Matrix&& X, trans_type tA) const override
        {
            return solve_rev(X,tA);
        };

        virtual matcl::Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            Integer N1      = m_A.cols();
            Integer N2      = m_B.cols();

            Matrix X1       = X(colon(1, N1), colon());
            Matrix X2       = X(colon(N1+1, N1+N2), colon());

            if (t == trans_type::no_trans)
            {
                Matrix Y1   = m_A.mmul_right(X1, t) + mmul(m_B, X2);
                Matrix Y2   = (m_sym == false)  ? mmul(m_C, X1) + mmul(m_D, X2)
                                                : mmul(m_B, X1, trans_type::conj_trans) + mmul(m_D, X2);

                return (mat_col(), Y1, Y2);
            }
            else
            {
                Matrix Y1;

                if (m_sym == false)
                {
                    Y1      = m_A.mmul_right(X1, t) + mmul(m_C, X2, t);
                }
                else
                {
                    trans_type_ext te   = details::trans_manip::link_trans(
                                            details::trans_manip::convert_trans(t), 
                                            trans_type_ext::conj_trans);

                    Matrix Y11  = m_A.mmul_right(X1, t);
                    Matrix Y12  =  mmul(m_B, X2, te, trans_type_ext::no_trans);
                    Y1          = std::move(Y11) + std::move(Y12);
                };

                Matrix Y2   = mmul(m_B, X1, t) + mmul(m_D, X2, t);

                return (mat_col(), Y1, Y2);
            };
        };

        virtual matcl::Matrix mmul_left(const Matrix& X, trans_type t) const override
        {
            Integer N1      = m_A.rows();
            Integer N2      = m_D.rows();

            Matrix X1       = X(colon(), colon(1, N1));
            Matrix X2       = X(colon(), colon(N1+1, N1+N2));

            if (t == trans_type::no_trans)
            {
                Matrix Y1   = (m_sym == false)  ? mmul(X2, m_C, trans_type::no_trans, trans_type::no_trans)
                                                : mmul(X2, m_B, trans_type::no_trans, trans_type::conj_trans);
                Y1          = std::move(Y1) + m_A.mmul_left(X1, t);
                Matrix Y2   = mmul(X1, m_C) + mmul(X2, m_D);

                return (mat_row(), Y1, Y2);
            }
            else
            {
                Matrix Y1   = m_A.mmul_left(X1, t) + mmul(X2, m_B, trans_type::no_trans, t);
                Matrix Y2;

                if (m_sym == false)
                {
                    Y2      = mmul(X1, m_C, trans_type::no_trans, t) + mmul(X2, m_D, trans_type::no_trans, t);
                }
                else
                {
                    trans_type_ext te   = details::trans_manip::link_trans(
                                            details::trans_manip::convert_trans(t), 
                                            trans_type_ext::conj_trans);

                    Matrix Y21  = mmul(X2, m_D, trans_type::no_trans, t);
                    Matrix Y22  = mmul(X1, m_B, trans_type_ext::no_trans, te);
                    Y2          = std::move(Y21) + std::move(Y22);
                };

                return (mat_row(), Y1, Y2);
            };
        };

        virtual matcl::Matrix mmul_right(Matrix&& X, trans_type t) const override
        {
            return mmul_right(X, t);
        };

        virtual matcl::Matrix inv() const override
        {
            const Matrix& AUQ   = m_precomp.get_AUQ();
            const Matrix& AV    = m_precomp.get_AV();

            // A_11 * Y = X => Y = inv(A_11) * X = inv(A)*X1 + AUQ * AV'*X1
            Matrix Y11      = m_A.inv() + mmul(AUQ, AV, trans_type::no_trans, trans_type::conj_trans);

            // A_12 * Y = X => Y = inv(A_12) * X = -AUQ *X2
            Matrix Y12      = -AUQ;

            // A_22 * Y = X => Y = inv(A_22) * X = inv(Q) * X2
            Matrix Y22      = m_precomp.get_Q().inv();

            // A_21 * X = X => Y = inv(A_21) * X = -inv(Q) * AV' *X1
            Matrix Y21      = -mmul(Y22, AV, trans_type::no_trans, trans_type::conj_trans);

            Matrix ret      = (mat_row(), (mat_col(), Y11, Y21), (mat_col(), Y12, Y22) );

            if (is_hermitian() == true)
            {
                ret.add_struct(predefined_struct_type::her);
            };

            return ret;
        };

        virtual matcl::Matrix base_matrix() const
        {
            Matrix A;
            if (m_sym == false)
                A       = (mat_col(), (mat_row(), m_A.base_matrix(), m_B), (mat_row(), m_C, m_D));
            else
                A       = (mat_col(), (mat_row(), m_A.base_matrix(), m_B), (mat_row(), ctrans(m_B), m_D));
            return A;
        };

        //not reimplemented
        //virtual Real normest_1() const override;
        //virtual Real normest_inf() const override
        //virtual Real normest_2() const override;
        //virtual Real mat_normest_2() const override;

    private:
        linsolve_insert_A(const linsolve_obj& Ac, const Matrix& Bc, const Matrix& Cc, 
                          const Matrix& Dc, const SM_precomp& pc, bool sym, value_code vc)
            : m_A(Ac), m_B(Bc), m_C(Cc), m_D(Dc), m_precomp(pc), m_sym(sym), m_vc(vc)
        {};

        void unify_check(bool sym)
        {
            value_code vA   = m_A.get_value_code();
            value_code vB   = m_B.get_value_code();
            value_code vC   = m_C.get_value_code();
            value_code vD   = m_D.get_value_code();
            m_vc            = matrix_traits::unify_value_types(vA, vB);
            m_vc            = matrix_traits::unify_value_types(m_vc, vD);
            m_vc            = matrix_traits::unify_value_types(m_vc, value_code::v_float);

            if (sym == false)
                m_vc        = matrix_traits::unify_value_types(m_vc, vC);

            if (m_D.rows() != m_D.cols())
                throw error::square_matrix_required(m_D.rows(), m_D.cols());

            error::check_horzcat(m_A.rows(), m_A.cols(), m_B.rows(), m_B.cols());
            error::check_vertcat(m_B.rows(), m_B.cols(), m_D.rows(), m_D.cols());

            if (sym == false)
            {
                error::check_horzcat(m_C.rows(), m_C.cols(), m_D.rows(), m_D.cols());
                error::check_vertcat(m_A.rows(), m_A.cols(), m_C.rows(), m_C.cols());            
            };

            if (vA != m_vc)
                m_A         = m_A.convert(m_vc);

            if (vB != m_vc)
                m_B         = details::convert_value(m_B, m_vc);

            if (vD != m_vc)
                m_D         = details::convert_value(m_D, m_vc);

            if (sym == false && vC != m_vc)
                m_C         = details::convert_value(m_C, m_vc);
        };
};

//------------------------------------------------------------------
//                      linsolve_insert_C
//------------------------------------------------------------------
class linsolve_insert_C : public linsolve_obj_data
{
    private:
        Matrix          m_A;
        Matrix          m_B;
        Matrix          m_C;
        linsolve_obj    m_D;
        value_code      m_vc;
        bool            m_sym;
        SM_precomp      m_precomp;

    public:
        linsolve_insert_C(const Matrix& A, const Matrix& B, const Matrix& C, const linsolve_obj& D,
                          const options& opt, bool& all_finite)
            :m_A(A), m_B(B), m_C(C), m_D(D), m_sym(false)
            , m_precomp(false, C, B, A, D, opt, SM_precomp::ver_block())
        {
            unify_check(false);
            all_finite = m_precomp.all_finite();
        };
        linsolve_insert_C(const Matrix& A, const Matrix& B, const linsolve_obj& D, const options& opt,
                          bool& all_finite)
            :m_A(A), m_B(B), m_D(D), m_sym(true)
            , m_precomp(false, B, B, A, D, opt, SM_precomp::ver_block())
        {
            unify_check(true);
            all_finite = m_precomp.all_finite();
        };

        virtual ~linsolve_insert_C(){};

        virtual bool all_finite() const override
        { 
            return true;
        };

        virtual bool is_modified() const override
        {
            return m_D.is_modified() || m_precomp.is_modified();
        };

        virtual bool is_direct() const override
        {
            return m_D.is_direct() && m_precomp.is_direct();
        };

        virtual Integer rows() const override
        {
            return m_A.rows() + m_D.rows();
        };

        virtual Integer cols() const override
        {
            return m_A.cols() + m_D.cols();
        };

        virtual value_code get_value_code() const override
        {
            return m_vc;
        };

        virtual ti::ti_object get_type() const override
        {
            return m_A.get_type();
        };

        virtual bool is_hermitian() const override
        { 
            if (m_sym == false)
            {
                return false;
            }
            else
            {
                bool is_real    = matrix_traits::is_float_real(m_vc);
                bool is_herm_A  = m_A.get_struct().is_hermitian(m_A.is_square(), is_real);
                return m_D.is_hermitian() && is_herm_A;
            };
        };

        virtual bool is_posdef() const override
        {
            return false;
        };

        virtual data_ptr convert(value_code vc) const override
        {
            linsolve_obj Dc     = m_D.convert(vc);
            Matrix Bc           = details::convert_value(m_B, vc);
            Matrix Cc           = (m_sym == false) ? details::convert_value(m_C, vc) : m_C;
            Matrix Ac           = details::convert_value(m_A, vc);
            SM_precomp pc       = m_precomp.convert(vc);

            return data_ptr(new linsolve_insert_C(Ac, Bc, Cc, Dc, pc, m_sym, vc) );
        };        

        virtual Real log_det() const override
        {
            Real det_D  = m_D.log_det();
            Real det_Q  = m_precomp.log_det();

            return det_D + det_Q;
        };

        virtual matcl::Matrix solve(const Matrix& X, trans_type tA) const override
        {
            Integer N1      = m_A.cols();
            Integer N2      = m_B.cols();

            Matrix X1       = X(colon(1, N1), colon());
            Matrix X2       = X(colon(N1+1, N1+N2), colon());

            if (tA == trans_type::no_trans)
            {
                Matrix Z1   = m_precomp.mult_AV(X2, tA);

                // A_11 * Y = X => Y = inv(A_11) * X = inv(Q) * X1
                Matrix Y11  = m_precomp.mult_Q(X1, tA);

                // A_12 * X = X => Y = inv(A_12) * X = -inv(Q) * AV' *X2 = -inv(Q) * Z1
                Matrix Y21  = -m_precomp.mult_Q(Z1, tA);

                // A_21 * Y = X => Y = inv(A_21) * X = -AUQ *X1
                Matrix Y12          = -m_precomp.mult_AUQ(X1, tA);

                // A_22 * Y = X => Y = inv(A_22) * X = inv(D)*X2 + AUQ * AV'*X2                
                Matrix Y22  = m_D.solve(X1, tA) + m_precomp.mult_AUQ(Z1, tA);

                Matrix ret  = (mat_col(), Y11 + Y12, Y21 + Y22);
                return ret;
            }
            else
            {
                Matrix Z1   = m_precomp.mult_AUQ(X2,tA);
                Matrix Z2   = m_precomp.mult_Q(X1,tA);

                // op(A_11)*Y = X => Y = inv(op(A_11))*X = inv(op(Q))*X1
                Matrix Y11  = Z2;

                // op(A_21)*Y = X => Y = inv(op(A_21))*X = -op(AUQ)*X2
                Matrix Y12  = -Z1;

                // op(A_12)*Y = X => Y = inv(op(A_12))*X = -op(AV')*inv(op(Q))*X1
                Matrix Y21  = -m_precomp.mult_AV(Z2, tA);

                // op(A_22)*Y = X => Y = inv(op(A_22))*X = inv(op(D))*X2 + op(AV') * op(AUQ)*X2
                Matrix Y22  = m_D.solve_rev(X2, tA) + m_precomp.mult_AV(Z1, tA);

                Matrix ret  = (mat_col(), Y11 + Y12, Y21 + Y22);
                return ret;
            };
        };

        virtual matcl::Matrix solve(Matrix&& X, trans_type tA) const override
        {
            return solve(X, tA);
        };

        virtual matcl::Matrix solve_rev(const Matrix& X, trans_type tA) const override
        {
            Integer N1      = m_A.rows();
            Integer N2      = m_D.rows();

            Matrix X1       = X(colon(), colon(1, N1));
            Matrix X2       = X(colon(), colon(N1+1, N1+N2));

            if (tA == trans_type::no_trans)
            {
                Matrix Z1   = m_precomp.mult_AUQ_rev(X2, tA);
                Matrix Z2   = m_precomp.mult_Q_rev(X1, tA);

                // Y * A_11 = X => Y = X * inv(A_11) = X1 * inv(Q) = Z2
                Matrix Y11  = Z2;

                // Y * A_12 = X => Y = X*inv(A_12) = -X1 * inv(Q) * AV'               
                Matrix Y12  = -m_precomp.mult_AV_rev(Z2, tA);

                // Y * A_21 = X => Y = X*inv(A_21) = -X2*AUQ = -Z1
                Matrix Y21  = -std::move(Z1);

                // Y * A_22 = X => Y = X*inv(A_22) = X2*inv(D) + X2*AUQ * AV                
                Matrix Y22  = m_D.solve(X2, tA) + m_precomp.mult_AV_rev(Z1, tA);

                Matrix ret  = (mat_row(), Y11 + Y21, Y12 + Y22);
                return ret;
            }
            else
            {
                Matrix Z1   = m_precomp.mult_AV_rev(X2, tA);

                // Y * op(A_11) = X => Y = X*inv(op(A_11)) = X1 * inv(op(Q))
                Matrix Y11  = m_precomp.mult_Q_rev(X1, tA);

                // Y * op(A_12) = X => Y = X*inv(op(A_12)) = -X2*op(AV')*inv(op(Q))
                Matrix Y21  = -m_precomp.mult_Q_rev(std::move(Z1),tA);

                // Y * op(A_21) = X => Y = X * inv(op(A_21)) = -X1 * op(AUQ)
                Matrix Y12  = -m_precomp.mult_AUQ_rev(X1, tA);

                // Y * op(A_22) = X => Y = X * inv(op(A_22)) = X2 * inv(op(A)) + X2*op(AV') * op(AUQ)
                Matrix Y22  = m_D.solve_rev(X2, tA) + m_precomp.mult_AUQ_rev(Z1, tA);

                Matrix ret  = (mat_row(), Y11 + Y21, Y12 + Y22);
                return ret;
            };
        };

        virtual matcl::Matrix solve_rev(Matrix&& X, trans_type tA) const override
        {
            return solve_rev(X,tA);
        };

        virtual matcl::Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            Integer N1      = m_A.cols();
            Integer N2      = m_B.cols();

            Matrix X1       = X(colon(1, N1), colon());
            Matrix X2       = X(colon(N1+1, N1+N2), colon());

            if (t == trans_type::no_trans)
            {
                Matrix Y1   = mmul(m_A, X1) + mmul(m_B, X2);
                Matrix Y2   = (m_sym == false)  ? m_D.mmul_right(X2, t) + mmul(m_C, X1)
                                                : m_D.mmul_right(X2, t) + mmul(m_B, X1, trans_type::conj_trans);

                return (mat_col(), Y1, Y2);
            }
            else
            {
                Matrix Y2   = mmul(m_B, X1, t) + m_D.mmul_right(X2, t);
                Matrix Y1;

                if (m_sym == false)
                {
                    Y1      = mmul(m_A, X1, t)+ mmul(m_C, X2, t);;
                }
                else
                {
                    trans_type_ext te   = details::trans_manip::link_trans(
                                            details::trans_manip::convert_trans(t), 
                                            trans_type_ext::conj_trans);

                    Matrix Y11  = mmul(m_A, X1, t);
                    Matrix Y12  = mmul(m_B, X2, te, trans_type_ext::no_trans);
                    Y1          = std::move(Y11) + std::move(Y12);
                };                

                return (mat_col(), Y1, Y2);
            };
        };

        virtual matcl::Matrix mmul_right(Matrix&& X, trans_type t) const override
        {
            return mmul_right(X, t);
        };

        virtual matcl::Matrix mmul_left(const Matrix& X, trans_type t) const override
        {
            Integer N1      = m_A.rows();
            Integer N2      = m_D.rows();

            Matrix X1       = X(colon(), colon(1, N1));
            Matrix X2       = X(colon(), colon(N1+1, N1+N2));

            if (t == trans_type::no_trans)
            {
                Matrix Y1   = (m_sym == false)  ? mmul(X2, m_C, trans_type::no_trans, trans_type::no_trans)
                                                : mmul(X2, m_B, trans_type::no_trans, trans_type::conj_trans);
                Y1          = std::move(Y1) + mmul(X1, m_A, trans_type::no_trans, t);
                Matrix Y2   = mmul(X1, m_C) + m_D.mmul_left(X2, t);

                return (mat_row(), Y1, Y2);
            }
            else
            {
                Matrix Y1   = mmul(X1, m_A, trans_type::no_trans, t) + mmul(X2, m_B, trans_type::no_trans, t);
                Matrix Y2;

                if (m_sym == false)
                {
                    Y2      = mmul(X1, m_C, trans_type::no_trans, t) + m_D.mmul_left(X2, t);
                }
                else
                {
                    trans_type_ext te   = details::trans_manip::link_trans(
                                            details::trans_manip::convert_trans(t), 
                                            trans_type_ext::conj_trans);

                    Matrix Y21  = m_D.mmul_left(X2, t);
                    Matrix Y22  = mmul(X1, m_B, trans_type_ext::no_trans, te);
                    Y2          = std::move(Y21) + std::move(Y22);
                };

                return (mat_row(), Y1, Y2);
            };
        };

        virtual matcl::Matrix inv() const override
        {
            const Matrix& AUQ   = m_precomp.get_AUQ();
            const Matrix& AV    = m_precomp.get_AV();

            // A_11 * Y = X => Y = inv(A_11) * X = inv(Q) * X1
            Matrix Y11      = m_precomp.get_Q().inv();

            // A_12 * X = X => Y = inv(A_12) * X = -inv(Q) * AV' *X2
            Matrix Y12      = -mmul(Y11, AV, trans_type::no_trans, trans_type::conj_trans);

            // A_21 * Y = X => Y = inv(A_21) * X = -AUQ *X2
            Matrix Y21      = -AUQ;

            // A_22 * Y = X => Y = inv(A_22) * X = inv(A)*X2 + AUQ * AV'*X2
            Matrix Y22      = m_D.inv() + mmul(AUQ, AV, trans_type::no_trans, trans_type::conj_trans);

            Matrix ret      = (mat_row(), (mat_col(), Y11, Y21), (mat_col(), Y12, Y22) );

            if (is_hermitian() == true)
            {
                ret.add_struct(predefined_struct_type::her);
            };

            return ret;
        };

        virtual matcl::Matrix base_matrix() const
        {
            Matrix A;
            if (m_sym == false)
                A       = (mat_col(), (mat_row(), m_A, m_B), (mat_row(), m_C, m_D.base_matrix()));
            else
                A       = (mat_col(), (mat_row(), m_A, m_B), (mat_row(), ctrans(m_B), m_D.base_matrix()));
            return A;
        };

        //not reimplemented
        //virtual Real normest_1() const override;
        //virtual Real normest_inf() const override
        //virtual Real normest_2() const override;
        //virtual Real mat_normest_2() const override;

    private:
        linsolve_insert_C(const Matrix& Ac, const Matrix& Bc, const Matrix& Cc, 
                          const linsolve_obj& Dc, const SM_precomp& pc, bool sym, value_code vc)
            : m_A(Ac), m_B(Bc), m_C(Cc), m_D(Dc), m_precomp(pc), m_sym(sym), m_vc(vc)
        {};

        void unify_check(bool sym)
        {
            value_code vA   = m_A.get_value_code();
            value_code vB   = m_B.get_value_code();
            value_code vC   = m_C.get_value_code();
            value_code vD   = m_D.get_value_code();
            m_vc            = matrix_traits::unify_value_types(vA, vB);
            m_vc            = matrix_traits::unify_value_types(m_vc, vD);
            m_vc            = matrix_traits::unify_value_types(m_vc, value_code::v_float);

            if (sym == false)
                m_vc        = matrix_traits::unify_value_types(m_vc, vC);

            if (m_A.rows() != m_A.cols())
                throw error::square_matrix_required(m_A.rows(), m_A.cols());

            error::check_horzcat(m_A.rows(), m_A.cols(), m_B.rows(), m_B.cols());
            error::check_vertcat(m_B.rows(), m_B.cols(), m_D.rows(), m_D.cols());

            if (sym == false)
            {
                error::check_horzcat(m_C.rows(), m_C.cols(), m_D.rows(), m_D.cols());
                error::check_vertcat(m_A.rows(), m_A.cols(), m_C.rows(), m_C.cols());            
            };

            if (vA != m_vc)
                m_A         = details::convert_value(m_A, m_vc);
                
            if (vB != m_vc)
                m_B         = details::convert_value(m_B, m_vc);

            if (vD != m_vc)
                m_D         = m_D.convert(m_vc);

            if (sym == false && vC != m_vc)
                m_C         = details::convert_value(m_C, m_vc);
        };
};

}};

namespace matcl
{

linsolve_obj matcl::update(const linsolve_obj& A, const Matrix& U, const Matrix& V,
                              const options& opt)
{
    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    value_code vA   = A.get_value_code();
    value_code vU   = U.get_value_code();
    value_code vV   = V.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vA, vU);
    vc              = matrix_traits::unify_value_types(vc, vV);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

    Integer N       = A.rows();

    if (N == 0)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A.get_type())));
    };

    if (A.all_finite() == false || U.all_finite() == false || V.all_finite() == false)
    {
        using data_ptr  = linsolve_obj::linsolve_data_ptr;

        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    };

    bool all_finite = false;
    linsolve_obj ret(impl_ptr(new details::linsolve_update_SM(A, U, V, opt,all_finite)));

    if (all_finite == false)
        return linsolve_obj(impl_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    else
        return ret;
};
linsolve_obj matcl::update(const linsolve_obj& A, const Matrix& U, const options& opt)
{
    value_code vA   = A.get_value_code();
    value_code vU   = U.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vA, vU);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

    Integer N       = A.rows();

    if (N == 0)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A.get_type())));
    };

    if (A.all_finite() == false || U.all_finite() == false)
    {
        using data_ptr  = linsolve_obj::linsolve_data_ptr;

        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    };

    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    bool all_finite = false;
    linsolve_obj ret(impl_ptr(new details::linsolve_update_SM(A, U, opt,all_finite)));

    if (all_finite == false)
        return linsolve_obj(impl_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    else
        return ret;
};

linsolve_obj matcl::update(const linsolve_obj& A, const Matrix& U, 
                            const linsolve_obj& C, const Matrix& V, const options& opt)
{
    value_code vA   = A.get_value_code();
    value_code vU   = U.get_value_code();
    value_code vC   = C.get_value_code();
    value_code vV   = V.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vA, vU);
    vc              = matrix_traits::unify_value_types(vc, vV);
    vc              = matrix_traits::unify_value_types(vc, vC);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

    Integer N       = A.rows();

    if (N == 0)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A.get_type())));
    };

    if (A.all_finite() == false || U.all_finite() == false || C.all_finite() == false
        || V.all_finite() == false)
    {
        using data_ptr  = linsolve_obj::linsolve_data_ptr;

        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    };

    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    bool all_finite = false;
    linsolve_obj ret(impl_ptr(new details::linsolve_update_W(A, U, C, V, opt,all_finite)));

    if (all_finite == false)
        return linsolve_obj(impl_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    else
        return ret;
};

linsolve_obj matcl::update(const linsolve_obj& A, const Matrix& U, 
                            const linsolve_obj& C, const options& opt)
{
    value_code vA   = A.get_value_code();
    value_code vU   = U.get_value_code();
    value_code vC   = C.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vA, vU);
    vc              = matrix_traits::unify_value_types(vc, vC);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

    Integer N       = A.rows();

    if (N == 0)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A.get_type())));
    };

    if (A.all_finite() == false || U.all_finite() == false || C.all_finite() == false)
    {
        using data_ptr  = linsolve_obj::linsolve_data_ptr;

        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    };

    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    bool all_finite = false;
    linsolve_obj ret(impl_ptr(new details::linsolve_update_W(A, U, C, opt, all_finite)));

    if (all_finite == false)
        return linsolve_obj(impl_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    else
        return ret;
};

linsolve_obj matcl::update_BI(const linsolve_obj& A, const Matrix& U, 
                              const Matrix& C, const Matrix& V, const options& opt)
{
    value_code vA   = A.get_value_code();
    value_code vU   = U.get_value_code();
    value_code vC   = C.get_value_code();
    value_code vV   = V.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vA, vU);
    vc              = matrix_traits::unify_value_types(vc, vV);
    vc              = matrix_traits::unify_value_types(vc, vC);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

    Integer N       = A.rows();

    if (N == 0)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A.get_type())));
    };

    if (A.all_finite() == false || U.all_finite() == false || C.all_finite() == false
        || V.all_finite() == false)
    {
        using data_ptr  = linsolve_obj::linsolve_data_ptr;

        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    };

    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    bool all_finite = false;
    linsolve_obj ret(impl_ptr(new details::linsolve_update_BI(A, U, C, V, opt, all_finite)));

    if (all_finite == false)
        return linsolve_obj(impl_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    else
        return ret;
};
linsolve_obj matcl::update_BI(const linsolve_obj& A, const Matrix& U, 
                              const Matrix& C,const options& opt)
{
    value_code vA   = A.get_value_code();
    value_code vU   = U.get_value_code();
    value_code vC   = C.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vA, vU);
    vc              = matrix_traits::unify_value_types(vc, vC);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

    Integer N       = A.rows();

    if (N == 0)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A.get_type())));
    };

    if (A.all_finite() == false || U.all_finite() == false || C.all_finite() == false)
    {
        using data_ptr  = linsolve_obj::linsolve_data_ptr;

        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    };

    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    bool all_finite = false;
    linsolve_obj ret(impl_ptr(new details::linsolve_update_BI(A, U, C, opt, all_finite)));

    if (all_finite == false)
        return linsolve_obj(impl_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    else
        return ret;
};

linsolve_obj matcl::block_invert(const linsolve_obj& A, const Matrix& B, 
                                    const Matrix& C, const Matrix& D, const options& opt)
{
    value_code vA   = A.get_value_code();
    value_code vB   = B.get_value_code();
    value_code vC   = C.get_value_code();
    value_code vD   = D.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vA, vB);
    vc              = matrix_traits::unify_value_types(vc, vC);
    vc              = matrix_traits::unify_value_types(vc, vD);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

    Integer N       = A.rows() + D.rows();

    if (N == 0)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A.get_type())));
    };

    if (A.all_finite() == false || B.all_finite() == false || C.all_finite() == false
        || D.all_finite() == false)
    {
        using data_ptr  = linsolve_obj::linsolve_data_ptr;

        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    };

    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    bool all_finite = false;
    linsolve_obj ret(impl_ptr(new details::linsolve_insert_A(A, B, C, D, opt, all_finite)));

    if (all_finite == false)
        return linsolve_obj(impl_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    else
        return ret;
};
linsolve_obj matcl::block_invert(const linsolve_obj& A, const Matrix& B, 
                                    const Matrix& D, const options& opt)
{
    value_code vA   = A.get_value_code();
    value_code vB   = B.get_value_code();
    value_code vD   = D.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vA, vB);
    vc              = matrix_traits::unify_value_types(vc, vD);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

    Integer N       = A.rows() + D.rows();

    if (N == 0)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A.get_type())));
    };

    if (A.all_finite() == false || B.all_finite() == false || D.all_finite() == false)
    {
        using data_ptr  = linsolve_obj::linsolve_data_ptr;

        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    };

    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    bool all_finite = false;
    linsolve_obj ret(impl_ptr(new details::linsolve_insert_A(A, B, D, opt, all_finite)));

    if (all_finite == false)
        return linsolve_obj(impl_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    else
        return ret;
};

linsolve_obj matcl::block_invert(const Matrix& A, const Matrix& B, 
                                    const Matrix& C, const linsolve_obj& D, const options& opt)
{
    value_code vA   = A.get_value_code();
    value_code vB   = B.get_value_code();
    value_code vC   = C.get_value_code();
    value_code vD   = D.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vA, vB);
    vc              = matrix_traits::unify_value_types(vc, vC);
    vc              = matrix_traits::unify_value_types(vc, vD);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

    Integer N       = A.rows() + D.rows();

    if (N == 0)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A.get_type())));
    };

    if (A.all_finite() == false || B.all_finite() == false || C.all_finite() == false
        || D.all_finite() == false)
    {
        using data_ptr  = linsolve_obj::linsolve_data_ptr;

        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    };

    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    bool all_finite = false;
    linsolve_obj ret(impl_ptr(new details::linsolve_insert_C(A, B, C, D, opt, all_finite)));
    
    if (all_finite == false)
        return linsolve_obj(impl_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    else
        return ret;
};

linsolve_obj matcl::block_invert(const Matrix& A, const Matrix& B, const linsolve_obj& D, const options& opt)
{
    value_code vA   = A.get_value_code();
    value_code vB   = B.get_value_code();
    value_code vD   = D.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vA, vB);
    vc              = matrix_traits::unify_value_types(vc, vD);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

    Integer N       = A.rows() + D.rows();

    if (N == 0)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A.get_type())));
    };

    if (A.all_finite() == false || B.all_finite() == false || D.all_finite() == false)
    {
        using data_ptr  = linsolve_obj::linsolve_data_ptr;

        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    };

    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    bool all_finite = false;
    linsolve_obj ret(impl_ptr(new details::linsolve_insert_C(A, B, D, opt, all_finite)));

    if (all_finite == false)
        return linsolve_obj(impl_ptr(new details::linsolve_obj_nan(N,vc, A.get_type())));
    else
        return ret;
};

};