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
#include "matcl-linalg/general/linalg_exception.h"
#include "linsolve_objects_decomp.h"
#include "matcl-linalg/norms_error/norm.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-linalg/special_matrices/struct_flag_linalg.h"
#include "matcl-linalg/linear_eq/linsolve.h"
#include "iterative_refinement.h"

namespace matcl { namespace details
{

// nan values are not allowed, return linsolve_nan in this case

//------------------------------------------------------------------
//                      linop_linsolve_mat
//------------------------------------------------------------------
class linop_linsolve_mat : public linear_operator_data
{
    private:
        linsolve_obj m_mat;

    public:
        linop_linsolve_mat(const linsolve_obj& mat)
            :m_mat(mat)
        {};

        virtual ~linop_linsolve_mat() 
        {};

        virtual value_code get_value_code() const override
        {
            return m_mat.get_value_code();
        }

        virtual data_ptr convert(value_code new_val_code) const override
        {
            return data_ptr(new linop_linsolve_mat(m_mat.convert(new_val_code)));
        };

        virtual Integer rows() const override
        {
            return m_mat.rows();
        };

        virtual Integer cols() const override
        {
            return m_mat.cols();
        };

        virtual bool is_hermitian() const override
        {
            return m_mat.is_hermitian();
        };

        virtual Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            return m_mat.mmul_right(X, t);
        };

        virtual Matrix mmul_right(Matrix&& X, trans_type t) const override
        {
            return m_mat.mmul_right(std::move(X), t);
        }
};

//------------------------------------------------------------------
//                      un_linsolve_obj
//------------------------------------------------------------------
// implement inv(a*A) or -inv(A)
class un_linsolve_obj : public linsolve_obj_data
{
    public:
        enum op_type
        {
            op_uminus, op_scale
        };

    private:
        linsolve_obj    m_A;
        op_type         m_op;
        value_code      m_vc;
        Matrix          m_alpha;    //used only if op_scale
        bool            m_modified;

    public:
        un_linsolve_obj(const linsolve_obj& A, op_type op, bool modif)
            :m_A(A), m_op(op), m_vc(A.get_value_code())
        {
            m_modified  = A.is_modified() || modif;
        };

        // op = op_scale
        un_linsolve_obj(const matcl::Matrix& alpha, const linsolve_obj& A, bool modif)
            :m_A(A), m_op(op_scale), m_vc(A.get_value_code()), m_alpha(alpha)
        {
            m_modified  = A.is_modified() || modif;
            unify_check();
        };

        virtual ~un_linsolve_obj(){};

        virtual bool all_finite() const override
        { 
            return true;
        };

        virtual bool is_direct() const override
        {
            return m_A.is_direct();
        };
        virtual Integer rows() const override
        {
            return m_A.rows();
        };

        virtual Integer cols() const override
        {
            return m_A.cols();
        };

        virtual bool is_modified() const override
        {
            return m_modified;
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
            return m_A.is_hermitian();
        };
        virtual bool is_posdef() const override
        {
            if (m_op == op_uminus)
                return false;
            else if (m_A.is_posdef() == false)
                return false;

            if (matrix_traits::is_float_real(m_alpha.get_value_code()) == false)
                return false;
            
            return bool(m_alpha > 0.0f);
        }

        virtual data_ptr convert(value_code vc) const override
        {
            return data_ptr(new un_linsolve_obj(m_A.convert(vc), m_op, vc, 
                                    details::convert_value(m_alpha, vc), m_modified));
        };

        virtual Real log_det() const
        {
            Real val    = m_A.log_det();

            if (m_op == op_uminus)
                return val;

            Real val_a  = matcl::log(abs(m_alpha)).get_scalar<Real>();

            return val + val_a;
        };

        virtual matcl::Matrix solve(const Matrix& X, trans_type tA) const
        {
            if (m_op == op_uminus)
                return -m_A.solve(X, tA);
            else
                return (1.0f / m_alpha) * m_A.solve(X, tA);
        };
        virtual matcl::Matrix solve(Matrix&& X, trans_type tA) const
        {
            if (m_op == op_uminus)
                return -m_A.solve(std::move(X), tA);
            else
                return (1.0f/m_alpha) * m_A.solve(std::move(X), tA);
        };
        virtual matcl::Matrix solve_rev(const Matrix& X, trans_type tA) const
        {
            if (m_op == op_uminus)
                return -m_A.solve_rev(X, tA);
            else
                return (1.0f/m_alpha) * m_A.solve_rev(X, tA);
        };
        virtual matcl::Matrix solve_rev(Matrix&& X, trans_type tA) const
        {
            if (m_op == op_uminus)
                return -m_A.solve_rev(std::move(X), tA);
            else
                return (1.0f/m_alpha) * m_A.solve_rev(std::move(X), tA);
        };
        virtual Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            if (m_op == op_uminus)
                return -m_A.mmul_right(X, t);
            else
                return m_alpha * m_A.mmul_right(X, t);
        };
        virtual Matrix mmul_right(Matrix&& X, trans_type t) const override
        {
            if (m_op == op_uminus)
                return -m_A.mmul_right(std::move(X), t);
            else
                return m_alpha * m_A.mmul_right(std::move(X), t);
        };

        virtual Matrix mmul_left(const Matrix& X, trans_type t) const override
        {
            if (m_op == op_uminus)
                return -m_A.mmul_left(X, t);
            else
                return m_alpha * m_A.mmul_left(X, t);
        };
        virtual Matrix mmul_left(Matrix&& X, trans_type t) const override
        {
            if (m_op == op_uminus)
                return -m_A.mmul_left(std::move(X), t);
            else
                return m_alpha * m_A.mmul_left(std::move(X), t);
        };

        virtual matcl::Matrix base_matrix() const override
        {
            Matrix ret  = m_A.base_matrix();

            if (m_op == op_uminus)
                ret = -ret;
            else
                ret = m_alpha * ret;

            if (is_hermitian() == true)
            {
                ret.add_struct(predefined_struct_type::her);

                if (is_posdef() == true)
                    ret.add_struct(posdef_flag());
            };

            return ret;
        };
        virtual matcl::Matrix inv() const override
        {
            Matrix ret;
            if (m_op == op_uminus)
                ret = -m_A.inv();
            else
                ret = (1.0f/m_alpha) * m_A.inv();

            if (is_hermitian() == true)
            {
                ret.add_struct(predefined_struct_type::her);

                if (is_posdef() == true)
                    ret.add_struct(posdef_flag());
            };

            return ret;
        };

        virtual Real mat_normest_1() const override
        {
            Real norm   = m_A.mat_normest(basic_vector_norm::norm_1);
            if (m_op == op_scale)
                return matcl::abs(m_alpha).get_scalar<Real>() * norm;
            else
                return norm;
        };
        virtual Real mat_normest_2() const override
        {
            Real norm   = m_A.mat_normest(basic_vector_norm::norm_2);
            if (m_op == op_scale)
                return matcl::abs(m_alpha).get_scalar<Real>() * norm;
            else
                return norm;
        };
        virtual Real mat_normest_inf() const override
        {
            Real norm   = m_A.mat_normest(basic_vector_norm::norm_inf);
            if (m_op == op_scale)
                return matcl::abs(m_alpha).get_scalar<Real>() * norm;
            else
                return norm;
        };

        virtual Real normest_1() const override
        {
            Real norm   = m_A.normest(basic_vector_norm::norm_1);
            if (m_op == op_scale)
                return (1.0/matcl::abs(m_alpha).get_scalar<Real>()) * norm;
            else
                return norm;
        };
        virtual Real normest_2() const override
        {
            Real norm   = m_A.normest(basic_vector_norm::norm_2);
            if (m_op == op_scale)
                return (1.0/matcl::abs(m_alpha).get_scalar<Real>()) * norm;
            else
                return norm;
        };
        virtual Real normest_inf() const override
        {
            Real norm   = m_A.normest(basic_vector_norm::norm_inf);
            if (m_op == op_scale)
                return (1.0/matcl::abs(m_alpha).get_scalar<Real>()) * norm;
            else
                return norm;
        };

    private:
        un_linsolve_obj(const linsolve_obj& A, op_type op, value_code vc, const Matrix& alpha, bool modif)
            :m_A(A), m_op(op), m_vc(vc), m_alpha(alpha), m_modified(modif)
        {};

        void unify_check()
        {
            value_code vA   = m_A.get_value_code();
            value_code vB   = m_alpha.get_value_code();
            m_vc            = matrix_traits::unify_value_types(vA, vB);
            m_vc            = matrix_traits::unify_value_types(m_vc, value_code::v_float);
            
            if (m_alpha.is_scalar() == false)
                throw error::scalar_required(m_alpha.rows(), m_alpha.cols());

            if (vA != m_vc)
                m_A         = m_A.convert(m_vc);

            if (vB != m_vc)
                m_alpha     = details::convert_value(m_alpha, m_vc);

            if (m_alpha == 0)
                throw error::error_singular();
        };
};

//------------------------------------------------------------------
//                      trans_linsolve_obj
//------------------------------------------------------------------
class trans_linsolve_obj : public linsolve_obj_data
{
    private:
        linsolve_obj    m_A;
        trans_type_ext  m_trans;
        bool            m_modified;

    public:
        trans_linsolve_obj(const linsolve_obj& A, trans_type_ext t, bool modif)
            :m_A(A), m_trans(t)
        {
            m_modified  = m_A.is_modified() || modif;
        };

        virtual ~trans_linsolve_obj(){};

        virtual bool all_finite() const override
        { 
            return true;
        };

        virtual bool is_modified() const override
        {
            return m_modified;
        };
        virtual bool is_direct() const override
        {
            return m_A.is_direct();
        };

        virtual Integer rows() const override
        {
            if (m_trans == trans_type_ext::no_trans || m_trans == trans_type_ext::conj)
                return m_A.rows();
            else
                return m_A.cols();
        };

        virtual Integer cols() const override
        {
            if (m_trans == trans_type_ext::no_trans || m_trans == trans_type_ext::conj)
                return m_A.cols();
            else
                return m_A.rows();
        };

        virtual value_code get_value_code() const override
        {
            return m_A.get_value_code();
        };

        virtual ti::ti_object get_type() const override
        {
            return m_A.get_type();
        };

        virtual bool is_hermitian() const override
        { 
            return m_A.is_hermitian();
        };

        virtual bool is_posdef() const override
        {
            return m_A.is_posdef();
        }

        virtual data_ptr convert(value_code vc) const override
        {
            return data_ptr(new trans_linsolve_obj(m_A.convert(vc), m_trans, m_modified));
        };

        virtual Real log_det() const override
        {
            return m_A.log_det();
        };

        virtual matcl::Matrix solve(const Matrix& X, trans_type tA) const override
        {
            trans_type_ext te   = details::trans_manip::link_trans(m_trans, 
                                    details::trans_manip::convert_trans(tA));

            bool conj;
            trans_type t2       = details::trans_manip::convert_trans(te, conj);
            bool is_real        = (conj == false) && matrix_traits::is_float_real(get_value_code());

            if (conj == false || is_real == true)
                return m_A.solve(X,t2);

            Matrix y = m_A.solve(matcl::conj(X), trans_type::no_trans);
            return matcl::conj(std::move(y));
        };

        virtual matcl::Matrix solve(Matrix&& X, trans_type tA) const override
        {
            trans_type_ext te   = details::trans_manip::link_trans(m_trans, 
                                    details::trans_manip::convert_trans(tA));

            bool conj;
            trans_type t2       = details::trans_manip::convert_trans(te, conj);
            bool is_real        = (conj == false) && matrix_traits::is_float_real(get_value_code());

            if (conj == false || is_real == true)
                return m_A.solve(std::move(X),t2);

            Matrix y = m_A.solve(matcl::conj(std::move(X)), trans_type::no_trans);
            return matcl::conj(std::move(y));
        };

        virtual matcl::Matrix solve_rev(const Matrix& X, trans_type tA) const override
        {
            trans_type_ext te   = details::trans_manip::link_trans(m_trans, 
                                    details::trans_manip::convert_trans(tA));

            bool conj;
            trans_type t2       = details::trans_manip::convert_trans(te, conj);
            bool is_real        = (conj == false) && matrix_traits::is_float_real(get_value_code());

            if (conj == false || is_real == true)
                return m_A.solve_rev(X,t2);

            Matrix y = m_A.solve_rev(matcl::conj(X), trans_type::no_trans);
            return matcl::conj(std::move(y));
        };

        virtual matcl::Matrix solve_rev(Matrix&& X, trans_type tA) const override
        {
            trans_type_ext te   = details::trans_manip::link_trans(m_trans, 
                                    details::trans_manip::convert_trans(tA));

            bool conj;
            trans_type t2       = details::trans_manip::convert_trans(te, conj);
            bool is_real        = (conj == false) && matrix_traits::is_float_real(get_value_code());

            if (conj == false || is_real == true)
                return m_A.solve_rev(std::move(X),t2);

            Matrix y = m_A.solve_rev(matcl::conj(std::move(X)), trans_type::no_trans);
            return matcl::conj(std::move(y));
        };

        virtual matcl::Matrix mmul_right(const Matrix& X, trans_type tA) const override
        {
            trans_type_ext te   = details::trans_manip::link_trans(m_trans, 
                                    details::trans_manip::convert_trans(tA));

            bool conj;
            trans_type t2       = details::trans_manip::convert_trans(te, conj);
            bool is_real        = (conj == false) && matrix_traits::is_float_real(get_value_code());

            if (conj == false || is_real == true)
                return m_A.mmul_right(X,t2);

            Matrix y = m_A.mmul_right(matcl::conj(X), trans_type::no_trans);
            return matcl::conj(std::move(y));
        };

        virtual matcl::Matrix mmul_right(Matrix&& X, trans_type tA) const override
        {
            trans_type_ext te   = details::trans_manip::link_trans(m_trans, 
                                    details::trans_manip::convert_trans(tA));

            bool conj;
            trans_type t2       = details::trans_manip::convert_trans(te, conj);
            bool is_real        = (conj == false) && matrix_traits::is_float_real(get_value_code());

            if (conj == false || is_real == true)
                return m_A.mmul_right(std::move(X),t2);

            Matrix y = m_A.mmul_right(matcl::conj(std::move(X)), trans_type::no_trans);
            return matcl::conj(std::move(y));
        };

        virtual matcl::Matrix mmul_left(const Matrix& X, trans_type tA) const override
        {
            trans_type_ext te   = details::trans_manip::link_trans(m_trans, 
                                    details::trans_manip::convert_trans(tA));

            bool conj;
            trans_type t2       = details::trans_manip::convert_trans(te, conj);
            bool is_real        = (conj == false) && matrix_traits::is_float_real(get_value_code());

            if (conj == false || is_real == true)
                return m_A.mmul_left(X,t2);

            Matrix y = m_A.mmul_left(matcl::conj(X), trans_type::no_trans);
            return matcl::conj(std::move(y));
        };
        virtual matcl::Matrix mmul_left(Matrix&& X, trans_type tA) const override
        {
            trans_type_ext te   = details::trans_manip::link_trans(m_trans, 
                                    details::trans_manip::convert_trans(tA));

            bool conj;
            trans_type t2       = details::trans_manip::convert_trans(te, conj);
            bool is_real        = (conj == false) && matrix_traits::is_float_real(get_value_code());

            if (conj == false || is_real == true)
                return m_A.mmul_left(std::move(X),t2);

            Matrix y = m_A.mmul_left(matcl::conj(std::move(X)), trans_type::no_trans);
            return matcl::conj(std::move(y));
        };

        virtual matcl::Matrix base_matrix() const override
        {
            Matrix ret      = m_A.base_matrix();

            bool conj;
            trans_type t    = details::trans_manip::convert_trans(m_trans, conj);

            if (conj == false)
                ret = matcl::trans(std::move(ret), t);
            else
                ret = matcl::conj(std::move(ret));

            if (is_hermitian() == true)
            {
                ret.add_struct(predefined_struct_type::her);

                if (is_posdef() == true)
                    ret.add_struct(posdef_flag());
            };

            return ret;
        };

        virtual matcl::Matrix inv() const override
        {
            Matrix ret      = m_A.inv();

            bool conj;
            trans_type t    = details::trans_manip::convert_trans(m_trans, conj);

            if (conj == false)
                ret = matcl::trans(std::move(ret), t);
            else
                ret = matcl::conj(std::move(ret));

            if (is_hermitian() == true)
            {
                ret.add_struct(predefined_struct_type::her);

                if (is_posdef() == true)
                    ret.add_struct(posdef_flag());
            };

            return ret;
        };

        virtual Real mat_normest_1() const override
        {
            if (m_trans == trans_type_ext::no_trans || m_trans == trans_type_ext::conj)
                return m_A.mat_normest(basic_vector_norm::norm_1);
            else
                return m_A.mat_normest(basic_vector_norm::norm_inf);
        };
        virtual Real mat_normest_2() const override
        {
            return m_A.mat_normest(basic_vector_norm::norm_2);
        };
        virtual Real mat_normest_inf() const override
        {
            if (m_trans == trans_type_ext::no_trans || m_trans == trans_type_ext::conj)
                return m_A.mat_normest(basic_vector_norm::norm_inf);
            else
                return m_A.mat_normest(basic_vector_norm::norm_1);
        };

        virtual Real normest_1() const override
        {
            if (m_trans == trans_type_ext::no_trans || m_trans == trans_type_ext::conj)
                return m_A.normest(basic_vector_norm::norm_1);
            else
                return m_A.normest(basic_vector_norm::norm_inf);
        };
        virtual Real normest_2() const override
        {
            return m_A.normest(basic_vector_norm::norm_2);
        };
        virtual Real normest_inf() const override
        {
            if (m_trans == trans_type_ext::no_trans || m_trans == trans_type_ext::conj)
                return m_A.normest(basic_vector_norm::norm_inf);
            else
                return m_A.normest(basic_vector_norm::norm_1);
        };
};

//------------------------------------------------------------------
//                      bin_linsolve_obj
//------------------------------------------------------------------
class bin_linsolve_obj : public linsolve_obj_data
{
    private:
        linsolve_obj    m_A;
        linsolve_obj    m_B;
        value_code      m_vc;
        trans_type      m_tA;
        trans_type      m_tB;
        bool            m_modified;

    public:
        bin_linsolve_obj(const linsolve_obj& A, const linsolve_obj& B, trans_type tA, 
                         trans_type tB, bool modified)
            :m_A(A), m_B(B), m_tA(tA), m_tB(tB)
        {
            m_modified  = A.is_modified() || B.is_modified() || modified;
            unify_check();
        };        

        virtual ~bin_linsolve_obj(){};

        virtual bool all_finite() const override
        { 
            return true;
        };

        virtual bool is_modified() const override
        {
            return m_modified;
        };
        virtual bool is_direct() const override
        {
            return m_A.is_direct() && m_B.is_direct();
        };

        virtual Integer rows() const override
        {
            if (m_tA == trans_type::no_trans)
                return m_A.rows();
            else
                return m_A.cols();
        };

        virtual Integer cols() const override
        {
            if (m_tB == trans_type::no_trans)
                return m_B.cols();
            else
                return m_B.rows();
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
            return false;
        };

        virtual bool is_posdef() const override
        {
            return false;
        }

        virtual data_ptr convert(value_code vc) const override
        {
            return data_ptr(new bin_linsolve_obj(m_A.convert(vc), m_B.convert(vc), m_tA, 
                                                 m_tB, m_modified));
        };

        virtual Real log_det() const override
        {
            return m_A.log_det() + m_B.log_det();
        };

        virtual matcl::Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            if (t == trans_type::no_trans)
            {
                //op(A) * op(B) * X
                Matrix y        = m_B.mmul_right(X, m_tB);
                y               = m_A.mmul_right(std::move(y), m_tA);
                return y;
            };

            //op(op(A) * op(B)) * X = op(op(B)) * op(op(A)) * X
            trans_type_ext tAe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tA), 
                                    details::trans_manip::convert_trans(t));
            trans_type_ext tBe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tB), 
                                    details::trans_manip::convert_trans(t));

            bool is_real        = matrix_traits::is_float_real(m_vc);

            bool conj_A, conj_B;
            trans_type tA2      = details::trans_manip::convert_trans(tAe, conj_A);
            trans_type tB2      = details::trans_manip::convert_trans(tBe, conj_B);
            bool need_conj      = false;

            Matrix y;

            if (conj_A == false || is_real)
            {
                y               = m_A.mmul_right(X, tA2);
            }
            else
            {
                y               = m_A.mmul_right(conj(X), tA2);
                need_conj       = true;
            }

            if (conj_B == false || is_real)
            {
                if (need_conj == true)
                    y           = m_B.mmul_right(conj(std::move(y)), tB2);
                else
                    y           = m_B.mmul_right(std::move(y), tB2);

                return y;
            }
            else
            {
                if (need_conj == false)
                    y           = m_B.mmul_right(conj(std::move(y)), tB2);
                else
                    y           = m_B.mmul_right(std::move(y), tB2);

                return matcl::conj(std::move(y));
            }
        };

        virtual matcl::Matrix mmul_right(Matrix&& X, trans_type t) const override
        {
            if (t == trans_type::no_trans)
            {
                //op(A) * op(B) * X
                Matrix y        = m_B.mmul_right(std::move(X), m_tB);
                y               = m_A.mmul_right(std::move(y), m_tA);
                return y;
            };

            //op(op(A) * op(B)) * X = op(op(B)) * op(op(A)) * X
            trans_type_ext tAe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tA), 
                                    details::trans_manip::convert_trans(t));
            trans_type_ext tBe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tB), 
                                    details::trans_manip::convert_trans(t));

            bool is_real        = matrix_traits::is_float_real(m_vc);

            bool conj_A, conj_B;
            trans_type tA2      = details::trans_manip::convert_trans(tAe, conj_A);
            trans_type tB2      = details::trans_manip::convert_trans(tBe, conj_B);
            bool need_conj      = false;

            Matrix y;

            if (conj_A == false || is_real)
            {
                y               = m_A.mmul_right(std::move(X), tA2);
            }
            else
            {
                y               = m_A.mmul_right(conj(std::move(X)), tA2);
                need_conj       = true;
            }

            if (conj_B == false || is_real)
            {
                if (need_conj == true)
                    y           = m_B.mmul_right(conj(std::move(y)), tB2);
                else
                    y           = m_B.mmul_right(std::move(y), tB2);

                return y;
            }
            else
            {
                if (need_conj == false)
                    y           = m_B.mmul_right(conj(std::move(y)), tB2);
                else
                    y           = m_B.mmul_right(std::move(y), tB2);

                return matcl::conj(std::move(y));
            }
        };

        virtual matcl::Matrix mmul_left(const Matrix& X, trans_type t) const override
        {
            if (t == trans_type::no_trans)
            {
                //X * op(A) * op(B)
                Matrix y        = m_A.mmul_left(X, m_tA);
                y               = m_B.mmul_left(std::move(y), m_tB);
                return y;
            };

            //X * op(op(A) * op(B)) = X * op(op(B)) * op(op(A))
            trans_type_ext tAe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tA), 
                                    details::trans_manip::convert_trans(t));
            trans_type_ext tBe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tB), 
                                    details::trans_manip::convert_trans(t));

            bool is_real        = matrix_traits::is_float_real(m_vc);

            bool conj_A, conj_B;
            trans_type tA2      = details::trans_manip::convert_trans(tAe, conj_A);
            trans_type tB2      = details::trans_manip::convert_trans(tBe, conj_B);
            bool need_conj      = false;

            Matrix y;

            if (conj_B == false || is_real)
            {
                y               = m_B.mmul_right(X, tB2);
            }
            else
            {
                y               = m_B.mmul_right(conj(X), tB2);
                need_conj       = true;
            }

            if (conj_A == false || is_real)
            {
                if (need_conj == true)
                    y           = m_A.mmul_right(conj(std::move(y)), tA2);
                else
                    y           = m_A.mmul_right(std::move(y), tA2);

                return y;
            }
            else
            {
                if (need_conj == false)
                    y           = m_A.mmul_right(conj(std::move(y)), tA2);
                else
                    y           = m_A.mmul_right(std::move(y), tA2);

                return matcl::conj(std::move(y));
            }
        };

        virtual matcl::Matrix mmul_left(Matrix&& X, trans_type t) const override
        {
            if (t == trans_type::no_trans)
            {
                //X * op(A) * op(B)
                Matrix y        = m_A.mmul_left(std::move(X), m_tA);
                y               = m_B.mmul_left(std::move(y), m_tB);
                return y;
            };

            //X * op(op(A) * op(B)) = X * op(op(B)) * op(op(A))
            trans_type_ext tAe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tA), 
                                    details::trans_manip::convert_trans(t));
            trans_type_ext tBe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tB), 
                                    details::trans_manip::convert_trans(t));

            bool is_real        = matrix_traits::is_float_real(m_vc);

            bool conj_A, conj_B;
            trans_type tA2      = details::trans_manip::convert_trans(tAe, conj_A);
            trans_type tB2      = details::trans_manip::convert_trans(tBe, conj_B);
            bool need_conj      = false;

            Matrix y;

            if (conj_B == false || is_real)
            {
                y               = m_B.mmul_right(std::move(X), tB2);
            }
            else
            {
                y               = m_B.mmul_right(conj(std::move(X)), tB2);
                need_conj       = true;
            }

            if (conj_A == false || is_real)
            {
                if (need_conj == true)
                    y           = m_A.mmul_right(conj(std::move(y)), tA2);
                else
                    y           = m_A.mmul_right(std::move(y), tA2);

                return y;
            }
            else
            {
                if (need_conj == false)
                    y           = m_A.mmul_right(conj(std::move(y)), tA2);
                else
                    y           = m_A.mmul_right(std::move(y), tA2);

                return matcl::conj(std::move(y));
            }
        };

        virtual matcl::Matrix solve(const Matrix& X, trans_type t) const override
        {
            if (t == trans_type::no_trans)
            {
                //op(A) * op(B) * Y = X => Y = inv(op(B)) * inv(op(A)) * X
                Matrix y        = m_A.solve(X, m_tA);
                y               = m_B.solve(std::move(y), m_tB);
                return y;
            };

            //op(op(A) * op(B)) * Y = X => op(op(B)) * op(op(A)) * Y = X
            //Y = inv(op(op(A))) * inv(op(op(B))) * X
            trans_type_ext tAe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tA), 
                                    details::trans_manip::convert_trans(t));
            trans_type_ext tBe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tB), 
                                    details::trans_manip::convert_trans(t));

            bool is_real        = matrix_traits::is_float_real(m_vc);

            bool conj_A, conj_B;
            trans_type tA2      = details::trans_manip::convert_trans(tAe, conj_A);
            trans_type tB2      = details::trans_manip::convert_trans(tBe, conj_B);
            bool need_conj      = false;

            Matrix y;

            if (conj_B == false || is_real)
            {
                y               = m_B.solve(X, tB2);
            }
            else
            {
                y               = m_B.solve(conj(X), tB2);
                need_conj       = true;
            }

            if (conj_A == false || is_real)
            {
                if (need_conj == true)
                    y           = m_A.solve(conj(std::move(y)), tA2);
                else
                    y           = m_A.solve(std::move(y), tA2);

                return y;
            }
            else
            {
                if (need_conj == false)
                    y           = m_A.mmul_right(conj(std::move(y)), tB2);
                else
                    y           = m_A.mmul_right(std::move(y), tB2);

                return matcl::conj(std::move(y));
            };
        };

        virtual matcl::Matrix solve(Matrix&& X, trans_type t) const override
        {
            if (t == trans_type::no_trans)
            {
                //op(A) * op(B) * Y = X => Y = inv(op(B)) * inv(op(A)) * X
                Matrix y        = m_A.solve(std::move(X), m_tA);
                y               = m_B.solve(std::move(y), m_tB);
                return y;
            };

            //op(op(A) * op(B)) * Y = X => op(op(B)) * op(op(A)) * Y = X
            //Y = inv(op(op(A))) * inv(op(op(B))) * X
            trans_type_ext tAe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tA), 
                                    details::trans_manip::convert_trans(t));
            trans_type_ext tBe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tB), 
                                    details::trans_manip::convert_trans(t));

            bool is_real        = matrix_traits::is_float_real(m_vc);

            bool conj_A, conj_B;
            trans_type tA2      = details::trans_manip::convert_trans(tAe, conj_A);
            trans_type tB2      = details::trans_manip::convert_trans(tBe, conj_B);
            bool need_conj      = false;

            Matrix y;

            if (conj_B == false || is_real)
            {
                y               = m_B.solve(std::move(X), tB2);
            }
            else
            {
                y               = m_B.solve(conj(std::move(X)), tB2);
                need_conj       = true;
            }

            if (conj_A == false || is_real)
            {
                if (need_conj == true)
                    y           = m_A.solve(conj(std::move(y)), tA2);
                else
                    y           = m_A.solve(std::move(y), tA2);

                return y;
            }
            else
            {
                if (need_conj == false)
                    y           = m_A.mmul_right(conj(std::move(y)), tB2);
                else
                    y           = m_A.mmul_right(std::move(y), tB2);

                return matcl::conj(std::move(y));
            };
        };

        virtual matcl::Matrix solve_rev(const Matrix& X, trans_type t) const override
        {
            if (t == trans_type::no_trans)
            {
                //Y * op(A) * op(B) = X => Y = X * inv(op(B)) * inv(op(A))
                Matrix y        = m_B.solve_rev(X, m_tB);
                y               = m_A.solve_rev(std::move(y), m_tA);
                return y;
            };

            //Y * op(op(A) * op(B)) = X => Y * op(op(B)) * op(op(A)) = X
            //Y = X * inv(op(op(A))) * inv(op(op(B)))
            trans_type_ext tAe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tA), 
                                    details::trans_manip::convert_trans(t));
            trans_type_ext tBe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tB), 
                                    details::trans_manip::convert_trans(t));

            bool is_real        = matrix_traits::is_float_real(m_vc);

            bool conj_A, conj_B;
            trans_type tA2      = details::trans_manip::convert_trans(tAe, conj_A);
            trans_type tB2      = details::trans_manip::convert_trans(tBe, conj_B);
            bool need_conj      = false;

            Matrix y;

            if (conj_A == false || is_real)
            {
                y               = m_A.solve_rev(X, tA2);
            }
            else
            {
                y               = m_A.solve_rev(conj(X), tA2);
                need_conj       = true;
            }

            if (conj_B == false || is_real)
            {
                if (need_conj == true)
                    y           = m_B.solve_rev(conj(std::move(y)), tB2);
                else
                    y           = m_B.solve_rev(std::move(y), tB2);

                return y;
            }
            else
            {
                if (need_conj == false)
                    y           = m_B.solve_rev(conj(std::move(y)), tB2);
                else
                    y           = m_B.solve_rev(std::move(y), tB2);

                return matcl::conj(std::move(y));
            }
        };

        virtual matcl::Matrix solve_rev(Matrix&& X, trans_type t) const override
        {
            if (t == trans_type::no_trans)
            {
                //Y * op(A) * op(B) = X => Y = X * inv(op(B)) * inv(op(A))
                Matrix y        = m_B.solve_rev(std::move(X), m_tB);
                y               = m_A.solve_rev(std::move(y), m_tA);
                return y;
            };

            //Y * op(op(A) * op(B)) = X => Y * op(op(B)) * op(op(A)) = X
            //Y = X * inv(op(op(A))) * inv(op(op(B)))
            trans_type_ext tAe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tA), 
                                    details::trans_manip::convert_trans(t));
            trans_type_ext tBe  = details::trans_manip::link_trans(
                                    details::trans_manip::convert_trans(m_tB), 
                                    details::trans_manip::convert_trans(t));

            bool is_real        = matrix_traits::is_float_real(m_vc);

            bool conj_A, conj_B;
            trans_type tA2      = details::trans_manip::convert_trans(tAe, conj_A);
            trans_type tB2      = details::trans_manip::convert_trans(tBe, conj_B);
            bool need_conj      = false;

            Matrix y;

            if (conj_A == false || is_real)
            {
                y               = m_A.solve_rev(std::move(X), tA2);
            }
            else
            {
                y               = m_A.solve_rev(conj(std::move(X)), tA2);
                need_conj       = true;
            }

            if (conj_B == false || is_real)
            {
                if (need_conj == true)
                    y           = m_B.solve_rev(conj(std::move(y)), tB2);
                else
                    y           = m_B.solve_rev(std::move(y), tB2);

                return y;
            }
            else
            {
                if (need_conj == false)
                    y           = m_B.solve_rev(conj(std::move(y)), tB2);
                else
                    y           = m_B.solve_rev(std::move(y), tB2);

                return matcl::conj(std::move(y));
            }
        };

        virtual matcl::Matrix base_matrix() const override
        {
            Matrix A    = m_A.base_matrix();
            Matrix B    = m_B.base_matrix();

            return trans(std::move(A), m_tA) * trans(std::move(B), m_tB);
        };

        virtual matcl::Matrix inv() const override
        {
            //(op(A) * op(B))^-1 = op(B)^-1 * op(A^-1)
            Matrix X    = trans(m_A.inv(), m_tA);
            X           = m_B.solve(std::move(X),m_tB);
            return X;
        };

        //not reimplemented
        //virtual Real normest_1() const override;
        //virtual Real normest_2() const override;
        //virtual Real normest_inf() const override
        //virtual Real mat_normest_1() const override;
        //virtual Real mat_normest_2() const override;
        //virtual Real mat_normest_inf() const override

    private:
        void unify_check()
        {
            value_code vA   = m_A.get_value_code();
            value_code vB   = m_B.get_value_code();
            m_vc            = matrix_traits::unify_value_types(vA, vB);
            m_vc            = matrix_traits::unify_value_types(m_vc, value_code::v_float);
            
            error::check_mul(m_A.rows(), m_A.cols(), m_B.rows(), m_B.cols(), m_tA, m_tB);

            if (vA != m_vc)
                m_A         = m_A.convert(m_vc);

            if (vB != m_vc)
                m_B         = m_B.convert(m_vc);
        };
};

}};

namespace matcl
{

//------------------------------------------------------------------
//                      linsolve_obj_data
//------------------------------------------------------------------
Real linsolve_obj_data::normest_1() const
{
    data_ptr ptr(const_cast<linsolve_obj_data*>(this)->shared_from_this());
    linsolve_obj obj(ptr);

    return matcl::normest_1(obj);
};
Real linsolve_obj_data::normest_2() const
{
    data_ptr ptr(const_cast<linsolve_obj_data*>(this)->shared_from_this());
    linsolve_obj obj(ptr);

    return matcl::normest_2(obj);
};
Real linsolve_obj_data::normest_inf() const
{
    data_ptr ptr(const_cast<linsolve_obj_data*>(this)->shared_from_this());
    linsolve_obj obj(ptr);

    return matcl::normest_inf(obj);
};

Real linsolve_obj_data::mat_normest_1() const
{
    data_ptr ptr(const_cast<linsolve_obj_data*>(this)->shared_from_this());
    linsolve_obj obj(ptr);
    linear_operator op(obj.linear_operator_mat());

    return matcl::normest_inf(op);
};
Real linsolve_obj_data::mat_normest_2() const
{
    data_ptr ptr(const_cast<linsolve_obj_data*>(this)->shared_from_this());
    linsolve_obj obj(ptr);
    linear_operator op(obj.linear_operator_mat());

    return matcl::normest_2(op);
};
Real linsolve_obj_data::mat_normest_inf() const
{
    data_ptr ptr(const_cast<linsolve_obj_data*>(this)->shared_from_this());
    linsolve_obj obj(ptr);
    linear_operator op(obj.linear_operator_mat());

    return matcl::normest_inf(obj);
};

Matrix linsolve_obj_data::mmul_right(Matrix&& X, trans_type t) const
{
    return this->mmul_right(X,t);
};

Matrix linsolve_obj_data::mmul_left(Matrix&& X, trans_type t) const
{
    return this->mmul_left(X,t);
};
matcl::Matrix linsolve_obj_data::base_matrix() const
{
    Integer N       = this->rows();
    value_code vc   = this->get_value_code();
    Matrix I        = speye(N, N, vc);
    return this->mmul_right(I, trans_type::no_trans);
};

matcl::Matrix linsolve_obj_data::solve_rev(const Matrix& X, trans_type tA) const
{
    if (tA == trans_type::no_trans)
    {
        //Y * A = X => A' * Y' = X'
        Matrix Y    = this->solve(ctrans(X), trans_type::conj_trans);
        return ctrans(std::move(Y));
    }
    else
    {
        //Y * op(A) = X => A * op(Y) = op(X)
        Matrix Y    = this->solve(trans(X, tA), trans_type::no_trans);
        return trans(std::move(Y), tA);
    };
}
matcl::Matrix linsolve_obj_data::solve_rev(Matrix&& X, trans_type tA) const
{
    Matrix X2(std::move(X));
    return solve_rev(X2,tA);
}

matcl::Matrix linsolve_obj_data::inv() const
{
    Matrix I = speye(this->cols(), this->cols(), this->get_value_code());
    return this->solve(I, trans_type::no_trans);
};

ti::ti_object linsolve_obj_data::get_type() const
{
    return ti::ti_object_type(this->get_value_code());
};

//------------------------------------------------------------------
//                      linsolve_obj
//------------------------------------------------------------------
linsolve_obj::linsolve_obj()
    :m_impl(new details::linsolve_obj_scalar(1.0, options()))
{};

linsolve_obj::linsolve_obj(const linsolve_data_ptr& rep)
    :m_impl(rep)
{
    if (!m_impl)
        throw error::uninitialized_object_used("linsolve_obj");
};

linsolve_obj::linsolve_obj(linsolve_data_ptr&& rep)
    :m_impl(std::move(rep))
{
    if (!m_impl)
        throw error::uninitialized_object_used("linsolve_obj");
};

linsolve_obj::linsolve_obj(const linsolve_obj& mat)
    :m_impl(mat.m_impl)
{};

linsolve_obj::linsolve_obj(linsolve_obj&& mat)
    :m_impl(std::move(mat.m_impl))
{};

linsolve_obj& linsolve_obj::operator=(const linsolve_obj& other) &
{
    m_impl  = other.m_impl;
    return *this;
}
linsolve_obj& linsolve_obj::operator=(linsolve_obj&& other) &
{
    m_impl  = std::move(other.m_impl);
    return *this;
};

linsolve_obj::~linsolve_obj()
{};

Integer linsolve_obj::rows() const
{
    return m_impl->rows();
};
Integer linsolve_obj::cols() const
{
    return m_impl->cols();
};
Integer linsolve_obj::length() const
{
    Integer r = this->rows();
    Integer c = this->cols();

    if (r == 0 || c == 0)
        return 0;

    return std::max(r,c);
};
Real linsolve_obj::numel() const
{
    return Real(this->rows()) * Real(this->cols());
};

bool linsolve_obj::all_finite() const
{
    return m_impl->all_finite();
};

bool linsolve_obj::is_empty() const
{
    Integer r = this->rows();
    Integer c = this->cols();

    if (r == 0 || c == 0)
        return true;
    else
        return false;
};

bool linsolve_obj::is_scalar() const
{
    Integer r = this->rows();
    Integer c = this->cols();

    if (r == 1 && c == 1)
        return true;
    else
        return false;
};

bool linsolve_obj::is_square() const
{
    Integer r = this->rows();
    Integer c = this->cols();

    if (r == c)
        return true;
    else
        return false;
};

bool linsolve_obj::is_vector() const
{
    Integer r = this->rows();
    Integer c = this->cols();

    if (r == 1 || c == 1 || r == 0 || c == 0)
        return true;
    else
        return false;
};

value_code linsolve_obj::get_value_code() const
{
    return m_impl->get_value_code();
};

ti::ti_object linsolve_obj::get_type() const
{
    return m_impl->get_type();
};

matcl::Matrix linsolve_obj::base_matrix() const
{
    return m_impl->base_matrix();
};

matcl::Matrix linsolve_obj::inv() const
{
    return m_impl->inv();
}

matcl::Matrix linsolve_obj::solve(const Matrix& X, trans_type tA) const
{
    return m_impl->solve(X, tA);
};
matcl::Matrix linsolve_obj::solve(Matrix&& X, trans_type tA) const
{
    return m_impl->solve(std::move(X), tA);
};
matcl::Matrix linsolve_obj::solve_rev(const Matrix& X, trans_type tA) const
{
    return m_impl->solve_rev(X, tA);
};
matcl::Matrix linsolve_obj::solve_rev(Matrix&& X, trans_type tA) const
{
    return m_impl->solve_rev(std::move(X), tA);
};

linsolve_obj linsolve_obj::convert(value_code new_val_code) const
{
    if (new_val_code == this->get_value_code())
        return *this;

    linsolve_data_ptr dp = m_impl->convert(new_val_code);
    return linsolve_obj(dp);
};

bool linsolve_obj::is_hermitian() const
{
    return m_impl->is_hermitian();
};
bool linsolve_obj::is_posdef() const
{
    return m_impl->is_posdef();
};
bool linsolve_obj::is_modified() const
{
    return m_impl->is_modified();
};
bool linsolve_obj::is_direct() const
{
    return m_impl->is_direct();
};
static basic_vector_norm trans_norm(basic_vector_norm p)
{
    if (p == basic_vector_norm::norm_2)
        return p;

    if (p == basic_vector_norm::norm_1)
        return basic_vector_norm::norm_inf;
    else
        return basic_vector_norm::norm_1;
};

Real linsolve_obj::normest(basic_vector_norm p, trans_type t) const
{
    if (t != trans_type::no_trans)
        p  = trans_norm(p);

    switch (p)
    {
        case basic_vector_norm::norm_1:     return m_impl->normest_1();
        case basic_vector_norm::norm_2:     return m_impl->normest_2();
        case basic_vector_norm::norm_inf:   return m_impl->normest_inf();

        default:
        {
            matcl_assert(0,"invalid case");
            return m_impl->normest_1();
        }
    }
};

Real linsolve_obj::mat_normest(basic_vector_norm p, trans_type t) const
{
    if (t != trans_type::no_trans)
        p  = trans_norm(p);

    switch (p)
    {
        case basic_vector_norm::norm_1:     return m_impl->mat_normest_1();
        case basic_vector_norm::norm_2:     return m_impl->mat_normest_2();
        case basic_vector_norm::norm_inf:   return m_impl->mat_normest_inf();

        default:
        {
            matcl_assert(0,"invalid case");
            return m_impl->mat_normest_1();
        }
    }

};

Real linsolve_obj::rcond_est(basic_vector_norm p) const
{
    Real norm_A = this->mat_normest(p);

    if (norm_A == 0)
        return 0.0;

    Real norm_Ai    = this->normest(p);

    if (norm_Ai == 0)
        return 0.0;

    return 1.0/(norm_A * norm_Ai);
};

Real linsolve_obj::log_det() const
{
    return m_impl->log_det();
};

Matrix linsolve_obj::mmul_right(const Matrix& X, trans_type t) const
{
    return m_impl->mmul_right(X, t);
};
Matrix linsolve_obj::mmul_right(Matrix&& X, trans_type t) const
{
    return m_impl->mmul_right(std::move(X), t);
};

Matrix linsolve_obj::mmul_left(const Matrix& X, trans_type t) const
{
    return m_impl->mmul_left(X, t);
};
Matrix linsolve_obj::mmul_left(Matrix&& X, trans_type t) const
{
    return m_impl->mmul_left(std::move(X), t);
};

linear_operator linsolve_obj::linear_operator_inv() const
{
    return linear_operator(*this);
};

linear_operator linsolve_obj::linear_operator_mat() const
{
    using data_ptr = linear_operator::data_ptr;
    data_ptr impl = data_ptr(new details::linop_linsolve_mat(*this));

    return linear_operator(impl);
};

Matrix linsolve_obj::resid(const Matrix& Y, const Matrix& X, trans_type t) const
{
    Matrix r = this->mmul_right(Y, t) - X;
    return r;
}
Matrix linsolve_obj::resid_rev(const Matrix& Y, const Matrix& X, trans_type t) const
{
    Matrix r = this->mmul_left(Y, t) - X;
    return r;
};

Matrix linsolve_obj::resid_norm(const Matrix& Y, const Matrix& X, basic_vector_norm p, 
                                   trans_type t) const
{
    Matrix r = this->resid(Y,X, t);
    return norm_vec(r, p, 1);
};

Matrix linsolve_obj::resid_rev_norm(const Matrix& Y, const Matrix& X, basic_vector_norm p, 
                                       trans_type t) const
{
    Matrix r = this->resid_rev(Y,X,t);
    return norm_vec(r, p, 2);
};

mat_tup_2 linsolve_obj::norm_error(const Matrix& Y, const Matrix& X, basic_vector_norm p, 
                                           trans_type t) const
{
    Matrix rn       = this->resid_norm(Y, X, p, t);
    Matrix ry       = norm_vec(Y, p, 1);
    Matrix rx       = norm_vec(X, p, 1);

    Real rA         = mat_normest(p, t);
    
    Matrix bN       = div_0(rn, rA * ry + rx);

    value_code vc   = Y.get_value_code();
    bool is_float   = (vc == value_code::v_float || vc == value_code::v_float_complex);
    Real gam        = is_float ? constants::eps<Float>() : constants::eps<Real>();

    // we cannot estimate number of nonzeros of A, taking maximum value (i.e. size of A)
    // can give much too large estimation of forward error, on the other hand not adding
    // error resulting from computation of residuals (as done usually, see LAPACK) can 
    // underestimate forward error substantially; we take the smallest possible estimation
    // of nz per row (i.e. 1)
    gam             = gam * 2;

    // take into account, that residuals are calculated inexactly
    bN              = bN + gam;

    if (t != trans_type::no_trans)
        p           = trans_norm(p);

    Real rc         = rcond_est(p);
    Matrix kA       = div_0(bN, rc);

    Matrix fN       = div_0(2 * kA, 1 - kA);

    return mat_tup_2(bN, fN);
};

mat_tup_2 linsolve_obj::norm_error_rev(const Matrix& Y, const Matrix& X, basic_vector_norm p, 
                                               trans_type t) const
{
    Matrix rn   = this->resid_rev_norm(Y, X, p, t);
    Matrix ry   = norm_vec(Y, p, 2);
    Matrix rx   = norm_vec(X, p, 2);

    trans_type_ext te   = details::trans_manip::convert_trans(t);
    trans_type_ext tt   = details::trans_manip::link_trans(te, trans_type_ext::conj_trans);

    bool conj_A;
    trans_type t2       = details::trans_manip::convert_trans(tt, conj_A);

    Real rA     = mat_normest(p, t2);
    
    Matrix bN   = div_0(rn, rA * ry + rx);

    value_code vc   = Y.get_value_code();
    bool is_float   = (vc == value_code::v_float || vc == value_code::v_float_complex);
    Real gam        = is_float ? constants::eps<Float>() : constants::eps<Real>();

    // we cannot estimate number of nonzeros of A, taking maximum value (i.e. size of A)
    // can give much too large estimation of forward error, on the other hand not adding
    // error resulting from computation of residuals (as done usually, see LAPACK) can 
    // underestimate forward error substantially; we take the smallest possible estimation
    // of nz per row (i.e. 1)
    gam             = gam * 2;

    // take into account, that residuals are calculated inexactly
    bN              = bN + gam;

    if (t2 != trans_type::no_trans)
        p       = trans_norm(p);

    Real rc     = rcond_est(p);
    Matrix kA   = div_0(bN, rc);

    Matrix fN   = div_0(2 * kA, 1 - kA);

    return mat_tup_2(bN, fN);
};

//TODO: move to matcl

/// form abs(X) * abs(op(A)) + abs(Y)
Matrix mmul_abs_rev(const Matrix& A, const Matrix& X, trans_type tA, const Matrix& Y)
{
    return mmul(abs(X), abs(A), trans_type::no_trans, tA) + abs(Y);
};

Matrix linsolve_obj::comp_bacward_error(const Matrix& Y, const Matrix& X, const Matrix& E, const Matrix& F,
                                    trans_type t) const
{
    value_code vc;
    Matrix rt;
    {
        Matrix r    = abs(this->resid(Y,X,t));
        vc          = r.get_value_code();

        Matrix d    = mmul_abs(E, Y, t, F);
        rt          = div_0(std::move(r), std::move(d));
    };

    Matrix bc       = max_d(rt, 1);

    Integer M       = std::max(1, E.rows());
    Real nz_est     = Real(E.structural_nnz())/ real(M);    
    bool is_float   = (vc == value_code::v_float || vc == value_code::v_float_complex);
    Real gam        = is_float ? constants::eps<Float>() : constants::eps<Real>();
    gam             = gam * (nz_est + 1.0);

    // take into account, that residuals are calculated inexactly
    return bc + gam;
};

Matrix linsolve_obj::comp_bacward_error_rev(const Matrix& Y, const Matrix& X,
                                    const Matrix& E, const Matrix& F, trans_type t) const
{
    value_code vc;
    Matrix rt;
    {
        Matrix r    = abs(this->resid_rev(Y,X,t));
        vc          = r.get_value_code();
        Matrix d    = mmul_abs_rev(E, Y, t, F);
        //Matrix d  = mmul(abs(Y), E, trans_type::no_trans, t) + F;
        rt          = div_0(std::move(r), std::move(d));        
    };

    Matrix bc       = max_d(rt, 2);
    
    Integer M       = std::max(1, E.rows());
    Real nz_est     = Real(E.structural_nnz())/ real(M);    
    bool is_float   = (vc == value_code::v_float || vc == value_code::v_float_complex);
    Real gam        = is_float ? constants::eps<Float>() : constants::eps<Real>();
    gam             = gam * (nz_est + 1.0);

    // take into account, that residuals are calculated inexactly
    return bc + gam;
};

Real linsolve_obj::skeel_gen_cond(const Matrix& E, basic_vector_norm p, trans_type tA) const
{
    Matrix inv  = this->inv();
    Matrix AE   = abs(trans(inv,tA)) * abs(trans(E,tA));

    if (p == basic_vector_norm::norm_2)
        return matcl::normest_2(AE);
    else
        return matcl::norm(AE, p);
};

Real linsolve_obj::skeel_cond(basic_vector_norm p, trans_type tA) const
{
    return skeel_gen_cond(this->base_matrix(), p, tA);
};
Matrix linsolve_obj::skeel_vec_cond(const Matrix& X, trans_type tA) const
{
    Matrix Y            = mmul(abs(this->base_matrix()), abs(X), tA);
    linear_operator op  = linop_trans(this->linear_operator_inv(), tA);
    Matrix n_Ar         = abs_normest_r_vec_inf(op, Y);
    Matrix n_X          = norm_vec(X,basic_vector_norm::norm_inf);
    return div_0(n_Ar, n_X);
};
Matrix linsolve_obj::skeel_gen_vec_cond(const Matrix& X, trans_type tA) const
{
    linear_operator op  = linop_trans(this->linear_operator_inv(), tA);
    Matrix n_Ar         = abs_normest_r_vec_inf(op, X);
    Matrix n_X          = norm_vec(X,basic_vector_norm::norm_inf);
    return div_0(n_Ar, n_X);
};

mat_tup_2 linsolve_obj::comp_error(const Matrix& Y, const Matrix& X, trans_type t) const
{    
    Matrix A        = this->base_matrix();

    Matrix r        = abs(this->resid(Y,X,t));
    Matrix d        = mmul_abs(A, Y, t, X);
    Matrix rt       = div_0(r, d);

    Matrix bC       = max_d(rt, 1);

    Integer M       = std::max(1, A.rows());
    Real nz_est     = Real(A.structural_nnz())/ real(M);
    value_code vc   = r.get_value_code();
    bool is_float   = (vc == value_code::v_float || vc == value_code::v_float_complex);
    Real gam        = is_float ? constants::eps<Float>() : constants::eps<Real>();
    gam             = gam * (nz_est + 1.0);

    linear_operator op  = linop_trans(this->linear_operator_inv(), t);
    Matrix n_Ar     = abs_normest_r_vec_inf(op, std::move(r) + gam * std::move(d));
    Matrix n_x      = norm(Y, basic_vector_norm::norm_inf);
    Matrix fC       = div_0(n_Ar, n_x);

    return mat_tup_2(bC, fC);
};

mat_tup_2 linsolve_obj::comp_error_rev(const Matrix& Y, const Matrix& X, trans_type t) const
{
    Matrix A        = this->base_matrix();
    Matrix r        = abs(this->resid_rev(Y,X,t));
    Matrix d        = mmul_abs_rev(A, Y, t, X);
    Matrix rt       = div_0(r, d);

    Matrix bC       = max_d(rt, 2);

    Integer M       = std::max(1, A.rows());
    Real nz_est     = Real(A.structural_nnz())/ real(M);
    value_code vc   = r.get_value_code();
    bool is_float   = (vc == value_code::v_float || vc == value_code::v_float_complex);
    Real gam        = is_float ? constants::eps<Float>() : constants::eps<Real>();
    gam             = gam * (nz_est + 1.0);

    linear_operator op  = linop_trans(this->linear_operator_inv(), t);

    Matrix n_Ar     = abs_normest_l_vec_inf(op, std::move(r) + gam * std::move(d));
    Matrix n_x      = norm(Y, basic_vector_norm::norm_inf);
    Matrix fC       = div_0(n_Ar, n_x);

    return mat_tup_2(bC, fC);
};

matcl::Matrix linsolve_obj::iterative_refinement(const Matrix& Y, const Matrix& X, trans_type tA,
                                                 const options& opts) const
{
    Matrix Y_sol(Y);

    details::iterative_refinement(X, tA, *this, opts).eval(Y_sol, false);
    return Y_sol;
}

matcl::Matrix linsolve_obj::iterative_refinement(Matrix&& Y, const Matrix& X, trans_type tA, 
                                                 const options& opts) const
{
    Matrix Y_sol(std::move(Y));

    details::iterative_refinement(X, tA, *this, opts).eval(Y_sol, false);
    return Y_sol;
}

matcl::Matrix linsolve_obj::iterative_refinement_rev(const Matrix& Y, const Matrix& X, trans_type tA,
                                                     const options& opts) const
{
    Matrix Y_sol(Y);

    details::iterative_refinement(X, tA, *this, opts).eval(Y_sol, true);
    return Y_sol;
};
matcl::Matrix linsolve_obj::iterative_refinement_rev(Matrix&& Y, const Matrix& X, trans_type tA,
                                                     const options& opts) const
{
    Matrix Y_sol(std::move(Y));

    details::iterative_refinement(X, tA, *this, opts).eval(Y_sol, true);
    return Y_sol;
};

matcl::Matrix linsolve_obj::improve_inv(const Matrix& I) const
{
    Integer N       = this->rows();
    value_code vc   = this->get_value_code();
    Matrix I2       = mmul(I, (2.0f * speye(N, N, vc) - this->mmul_right(I)));
    return I2;
};

//------------------------------------------------------------------
//                      operations
//------------------------------------------------------------------
linsolve_obj matcl::trans(const linsolve_obj& mat)
{
    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    if (mat.rows() == 0)
        return mat;

    if (mat.all_finite() == false)
        return linsolve_nan(mat.rows(), mat.get_value_code());

    linsolve_obj ret(impl_ptr(new details::trans_linsolve_obj(mat, trans_type_ext::trans, false)));
    return ret;
};

linsolve_obj matcl::ctrans(const linsolve_obj& mat)
{
    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    if (mat.rows() == 0)
        return mat;
    if (mat.all_finite() == false)
        return linsolve_nan(mat.rows(), mat.get_value_code());

    linsolve_obj ret(impl_ptr(new details::trans_linsolve_obj(mat, trans_type_ext::conj_trans, false)));
    return ret;
};

linsolve_obj matcl::trans(const linsolve_obj& mat, trans_type t)
{
    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    if (mat.rows() == 0)
        return mat;
    if (mat.all_finite() == false)
        return linsolve_nan(mat.rows(), mat.get_value_code());

    linsolve_obj ret(impl_ptr(new details::trans_linsolve_obj(mat, details::trans_manip::convert_trans(t),false)));
    return ret;
};

linsolve_obj matcl::conj(const linsolve_obj& mat)
{
    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    if (mat.rows() == 0)
        return mat;
    if (mat.all_finite() == false)
        return linsolve_nan(mat.rows(), mat.get_value_code());

    linsolve_obj ret(impl_ptr(new details::trans_linsolve_obj(mat, trans_type_ext::conj, false)));
    return ret;
};

linsolve_obj matcl::operator-(const linsolve_obj& A)
{
    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    if (A.rows() == 0)
        return A;
    if (A.all_finite() == false)
        return linsolve_nan(A.rows(), A.get_value_code());

    linsolve_obj ret(impl_ptr(new details::un_linsolve_obj(A, details::un_linsolve_obj::op_uminus, false)));
    return ret;
};

linsolve_obj matcl::uminus(const linsolve_obj& A)
{
    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    if (A.rows() == 0)
        return A;
    if (A.all_finite() == false)
        return linsolve_nan(A.rows(), A.get_value_code());

    linsolve_obj ret(impl_ptr(new details::un_linsolve_obj(A, details::un_linsolve_obj::op_uminus, false)));
    return ret;
};

linsolve_obj matcl::scale(const Matrix& alpha, const linsolve_obj& A)
{
    using impl_ptr = linsolve_obj::linsolve_data_ptr;

    if (A.rows() == 0)
        return A;

    if (alpha.all_finite() == false || A.all_finite() == false)
    {
        using data_ptr  = linsolve_obj::linsolve_data_ptr;
        value_code vc   = matrix_traits::unify_value_types(alpha.get_value_code(), A.get_value_code());

        return linsolve_nan(A.rows(), vc);
    };

    linsolve_obj ret(impl_ptr(new details::un_linsolve_obj(alpha, A, false)));
    return ret;
};

linsolve_obj matcl::operator*(const linsolve_obj& A, const linsolve_obj& B)
{
    return mmul(A, B, trans_type::no_trans, trans_type::no_trans);
};

linsolve_obj matcl::mmul(const linsolve_obj& A, const linsolve_obj& B,
                            trans_type tA, trans_type tB)
{
    using impl_ptr = linsolve_obj::linsolve_data_ptr;
            
    error::check_mul(A.rows(), A.cols(), B.rows(), B.cols(), tA, tB);

    if (A.rows() == 0)
        return A;

    if (A.all_finite() == false || B.all_finite() == false)
    {
        using data_ptr  = linsolve_obj::linsolve_data_ptr;
        value_code vc   = matrix_traits::unify_value_types(A.get_value_code(), B.get_value_code());

        return linsolve_nan(A.rows(), vc);
    };

    linsolve_obj ret(impl_ptr(new details::bin_linsolve_obj(A, B, tA, tB, false)));
    return ret;
};

};