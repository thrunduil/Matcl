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

#include "matcl-linalg/special_matrices/matrix_functors.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-linalg/linear_eq/linsolve_object.h"

namespace matcl { namespace details
{

//------------------------------------------------------------------
//                      linear_operator_mat
//------------------------------------------------------------------
class linear_operator_mat : public linear_operator_data
{
    private:
        Matrix              m_mat;

    public:
        linear_operator_mat(const Matrix& mat)
            :m_mat(mat)
        {};

        virtual ~linear_operator_mat() 
        {};

        virtual value_code get_value_code() const override
        {
            return m_mat.get_value_code();
        }

        virtual data_ptr convert(value_code new_val_code) const override
        {
            mat_code mc     = matrix_traits::get_matrix_type(new_val_code, m_mat.get_struct_code());
            Matrix mat_conv = matcl::convert(m_mat, mc);
            return data_ptr(new linear_operator_mat(mat_conv));
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
            bool is_real = matrix_traits::is_float_real(m_mat.get_value_code());
            return m_mat.get_struct().is_hermitian(m_mat.is_square(), is_real);
        };

        virtual Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            return mmul(m_mat, X, t);
        };

        virtual Matrix mmul_right(Matrix&& X, trans_type t) const override
        {
            return mmul(m_mat, std::move(X), t);
        }
        virtual void mmul_right(const Matrix& X, trans_type t, Matrix& y) const override
        {
            gemm(1.0f,m_mat,X,t,trans_type::no_trans, 0.0f, y);
        }
};

//------------------------------------------------------------------
//                      linear_operator_umat
//------------------------------------------------------------------
class linear_operator_umat : public linear_operator_data
{
    private:
        unitary_matrix  m_mat;

    public:
        linear_operator_umat(const unitary_matrix& mat)
            :m_mat(mat)
        {};

        virtual ~linear_operator_umat() 
        {};

        virtual value_code get_value_code() const override
        {
            return m_mat.get_value_code();
        }

        virtual data_ptr convert(value_code new_val_code) const override
        {
            return data_ptr(new linear_operator_umat(m_mat.convert(new_val_code)));
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
            return false;
        };

        virtual Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            return mmul(m_mat, X, t);
        };

        virtual Matrix mmul_right(Matrix&& X, trans_type t) const override
        {
            return mmul(m_mat, std::move(X), t);
        }
};

//------------------------------------------------------------------
//                      linear_operator_linsolve
//------------------------------------------------------------------
class linear_operator_linsolve : public linear_operator_data
{
    private:
        linsolve_obj m_mat;

    public:
        linear_operator_linsolve(const linsolve_obj& mat)
            :m_mat(mat)
        {};

        virtual ~linear_operator_linsolve() 
        {};

        virtual value_code get_value_code() const override
        {
            return m_mat.get_value_code();
        }

        virtual data_ptr convert(value_code new_val_code) const override
        {
            return data_ptr(new linear_operator_linsolve(m_mat.convert(new_val_code)));
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
            return m_mat.solve(X, t);
        };

        virtual Matrix mmul_right(Matrix&& X, trans_type t) const override
        {
            return m_mat.solve(std::move(X), t);
        }
};

//------------------------------------------------------------------
//                      bin_linop
//------------------------------------------------------------------
class bin_linop : public linear_operator_data
{
    public:
        enum op_type
        {
            op_plus, op_minus, op_mult
        };

    private:
        linear_operator m_A;
        linear_operator m_B;
        op_type         m_op;
        value_code      m_vc;

        //used only by op_mult
        trans_type      m_tA;
        trans_type      m_tB;

    public:
        bin_linop(const linear_operator& A, const linear_operator& B, op_type op)
            :m_A(A), m_B(B), m_op(op)
        {
            unify_check();
        };

        //represent mult
        bin_linop(const linear_operator& A, const linear_operator& B, trans_type tA, trans_type tB)
            :m_A(A), m_B(B), m_op(op_mult), m_tA(tA), m_tB(tB)
        {
            unify_check();
        };

        virtual ~bin_linop() 
        {};

        virtual value_code get_value_code() const override
        {
            return m_vc;
        }

        virtual data_ptr convert(value_code vc) const override
        {
            return data_ptr(new bin_linop(m_A.convert(vc), m_B.convert(vc), vc, m_op, m_tA, m_tB));
        };

        virtual Integer rows() const override
        {
            if (m_op == op_mult)
            {
                if (m_tA == trans_type::no_trans)
                    return m_A.rows();
                else
                    return m_A.cols();
            }
            else
            {
                return m_A.rows();
            }
        };

        virtual Integer cols() const override
        {
            if (m_op == op_mult)
            {
                if (m_tB == trans_type::no_trans)
                    return m_B.cols();
                else
                    return m_B.rows();
            }
            else
            {
                return m_B.cols();
            }
        };

        virtual bool is_hermitian() const override
        {
            if (m_op == op_mult)
                return false;

            return m_A.is_hermitian() && m_B.is_hermitian();
        };

        virtual Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            switch (m_op)
            {
                case op_plus:
                {
                    return m_A.mmul_right(X, t) + m_B.mmul_right(X, t);
                }
                case op_minus:
                {
                    return m_A.mmul_right(X, t) - m_B.mmul_right(X, t);
                }
                case op_mult:
                {
                    return eval_mult(X, t);
                }
                default:
                {
                    matcl_assert(0,"unknown case");
                    throw error::error_general("unknown case");
                }
            };
        };

        virtual Matrix mmul_right(Matrix&& X, trans_type t) const override
        {
            switch (m_op)
            {
                case op_plus:
                {
                    Matrix y1 = m_A.mmul_right(X, t);
                    Matrix y2 = m_B.mmul_right(std::move(X), t);
                    return  std::move(y1) + std::move(y2);
                }
                case op_minus:
                {
                    Matrix y1 = m_A.mmul_right(X, t);
                    Matrix y2 = m_B.mmul_right(std::move(X), t);
                    return  std::move(y1) - std::move(y2);
                }
                case op_mult:
                {
                    return eval_mult(std::move(X), t);
                }
                default:
                {
                    matcl_assert(0,"unknown case");
                    throw error::error_general("unknown case");
                }
            };
        }

        virtual void mmul_right(const Matrix& X, trans_type t, Matrix& y) const override
        {
            y(colon())  = mmul_right(X, t);
        }

    private:
        //create from data from other bin_linop
        bin_linop(const linear_operator& A, const linear_operator& B, value_code vc, op_type op, 
                  trans_type tA, trans_type tB)
            :m_A(A), m_B(B), m_vc(vc), m_op(op), m_tA(tA), m_tB(tB)
        {};

        void unify_check()
        {
            value_code vA   = m_A.get_value_code();
            value_code vB   = m_B.get_value_code();
            m_vc            = matrix_traits::unify_value_types(vA, vB);
            
            if (m_op == op_mult)
                error::check_mul(m_A.rows(), m_A.cols(), m_B.rows(), m_B.cols(), m_tA, m_tB);
            else
                error::check_eeop(m_A.rows(), m_A.cols(), m_B.rows(), m_B.cols());

            if (vA != m_vc)
                m_A         = m_A.convert(m_vc);

            if (vB != m_vc)
                m_B         = m_B.convert(m_vc);
        };

        matcl::Matrix eval_mult(const matcl::Matrix& X, trans_type t) const
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
        }

        matcl::Matrix eval_mult(matcl::Matrix&& X, trans_type t) const
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
        }
};

//------------------------------------------------------------------
//                      un_linop
//------------------------------------------------------------------
class un_linop : public linear_operator_data
{
    public:
        enum op_type
        {
            op_uminus, op_scale
        };

    private:
        linear_operator m_A;
        op_type         m_op;
        value_code      m_vc;
        Matrix          m_alpha;    //used only if m_op == op_scale

    public:
        un_linop(const linear_operator& A, op_type op)
            :m_A(A), m_op(op), m_vc(A.get_value_code())
        {};

        un_linop(const linear_operator& A, const Matrix& alpha)
            :m_A(A), m_op(op_scale), m_alpha(alpha)
        {
            unify_check();
        };

        virtual ~un_linop() 
        {};

        virtual value_code get_value_code() const override
        {
            return m_vc;
        }

        virtual data_ptr convert(value_code vc) const override
        {
            if (m_op != op_scale)
                return data_ptr(new un_linop(m_A.convert(vc), m_op));

            Matrix alpha    = details::convert_value(m_alpha, vc);
            return data_ptr(new un_linop(m_A.convert(vc), alpha, vc, m_op));
        };

        virtual Integer rows() const override
        {
            return m_A.rows();
        };

        virtual Integer cols() const override
        {
            return m_A.cols();
        };

        virtual bool is_hermitian() const override
        {
            return m_A.is_hermitian();
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
        }

    private:
        //create from data from other un_linop
        un_linop(const linear_operator& A, const Matrix& alpha, value_code vc, op_type op)
            :m_A(A), m_alpha(alpha), m_vc(vc), m_op(op)
        {};

        void unify_check()
        {
            value_code vA   = m_A.get_value_code();
            value_code vB   = m_alpha.get_value_code();
            m_vc            = matrix_traits::unify_value_types(vA, vB);
            
            if (m_alpha.is_scalar() == false)
                throw error::scalar_required(m_alpha.rows(), m_alpha.cols());

            if (vA != m_vc)
                m_A         = m_A.convert(m_vc);

            if (vB != m_vc)
                m_alpha     = details::convert_value(m_alpha, m_vc);
        };
};

//------------------------------------------------------------------
//                      symsum_linop
//------------------------------------------------------------------
class symsum_linop : public linear_operator_data
{
    private:
        linear_operator m_A;
        linear_operator m_At;
        value_code      m_vc;

    public:
        symsum_linop(const linear_operator& A)
            :m_A(A), m_vc(A.get_value_code()), m_At(linop_ctrans(A))
        {};

        virtual ~symsum_linop() 
        {};

        virtual value_code get_value_code() const override
        {
            return m_vc;
        }

        virtual data_ptr convert(value_code vc) const override
        {
            return data_ptr(new symsum_linop(m_A.convert(vc)));
        };

        virtual Integer rows() const override
        {
            return m_A.rows();
        };

        virtual Integer cols() const override
        {
            return m_A.rows();
        };

        virtual bool is_hermitian() const override
        {
            return true;
        };

        virtual Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            //op(A)*X + op(A')*X
            Matrix S1   = m_A.mmul_right(X, t);
            Matrix S2   = m_At.mmul_right(X, t);
            return std::move(S1) + std::move(S2);
        }

        virtual Matrix mmul_right(Matrix&& X, trans_type t) const override
        {
            //op(A)*X + op(A')*X
            Matrix S2   = m_At.mmul_right(X, t);
            Matrix S1   = m_A.mmul_right(std::move(X), t);            
            return std::move(S1) + std::move(S2);
        };
};

//------------------------------------------------------------------
//                      symprod_linop
//------------------------------------------------------------------
// A * A' if trans = false 
// A' * A if trans = true
class symprod_linop : public linear_operator_data
{
    private:
        linear_operator m_AAt;

    public:
        symprod_linop(const linear_operator& A, bool trans)
        {
            if (trans == false)
                m_AAt   = linop_mmul(A, A, trans_type::no_trans, trans_type::conj_trans);
            else
                m_AAt   = linop_mmul(A, A, trans_type::conj_trans, trans_type::no_trans);
        };

        virtual ~symprod_linop() 
        {};

        virtual value_code get_value_code() const override
        {
            return m_AAt.get_value_code();
        }

        virtual data_ptr convert(value_code vc) const override
        {
            return data_ptr(new symprod_linop(m_AAt.convert(vc)));
        };

        virtual Integer rows() const override
        {
            return m_AAt.rows();
        };

        virtual Integer cols() const override
        {
            return m_AAt.cols();
        };

        virtual bool is_hermitian() const override
        {
            return true;
        };

        virtual Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            return m_AAt.mmul_right(X,t);
        }

        virtual Matrix mmul_right(Matrix&& X, trans_type t) const override
        {
            return m_AAt.mmul_right(std::move(X),t);
        };

    private:
        symprod_linop(const linear_operator& AAt)
            :m_AAt(AAt)
        {};
};

//------------------------------------------------------------------
//                      symcat_linop
//------------------------------------------------------------------
// [0 A'; A 0] if trans = false
// [0 A; A' 0] if trans = true
class symcat_linop : public linear_operator_data
{
    private:
        linear_operator m_A;
        linear_operator m_At;
        bool            m_trans;
        value_code      m_vc;

    public:
        symcat_linop(const linear_operator& A, bool trans)
            :m_A(A), m_trans(trans), m_vc(A.get_value_code())
            ,m_At(linop_ctrans(A))
        {};

        virtual ~symcat_linop() 
        {};

        virtual value_code get_value_code() const override
        {
            return m_vc;
        }

        virtual data_ptr convert(value_code vc) const override
        {
            return data_ptr(new symcat_linop(m_A.convert(vc), m_trans));
        };

        virtual Integer rows() const override
        {
            return m_A.rows() + m_A.cols();
        };

        virtual Integer cols() const override
        {
            return m_A.rows() + m_A.cols();
        };

        virtual bool is_hermitian() const override
        {
            return true;
        };

        virtual Matrix mmul_right(const Matrix& X, trans_type t) const override
        {
            Integer M   = m_A.rows();
            Integer N   = m_A.cols();

            if (m_trans == false)
            {
                if (t == trans_type::no_trans)
                {
                    //[0 A'] * [X1] = [A'*X2]
                    //[A 0 ]   [X2]   [A*X1]

                    Matrix X1   = X(colon(1,N), colon());
                    Matrix X2   = X(colon(N+1,N+M), colon());
                    Matrix P1   = m_At.mmul_right(std::move(X2), trans_type::no_trans);
                    Matrix P2   = m_A.mmul_right(std::move(X1), trans_type::no_trans);
                    return mat_col().add(P1).add(P2);
                }
                else
                {
                    //[0      op(A)] * [X1] = [op(A)*X2]
                    //[op(A') 0    ]   [X2]   [op(A')*X1]

                    Matrix X1   = X(colon(1,N), colon());
                    Matrix X2   = X(colon(N+1,N+M), colon());
                    Matrix P1   = m_A.mmul_right(std::move(X2), t);
                    Matrix P2   = m_At.mmul_right(std::move(X1), t);
                    return mat_col().add(P1).add(P2);
                }
            }
            else
            {
                if (t == trans_type::no_trans)
                {
                    //[0 A]   * [X1] = [A*X2]
                    //[A' 0 ]   [X2]   [A'*X1]

                    Matrix X1   = X(colon(1,M), colon());
                    Matrix X2   = X(colon(M+1,N+M), colon());
                    Matrix P1   = m_A.mmul_right(std::move(X2), trans_type::no_trans);
                    Matrix P2   = m_At.mmul_right(std::move(X1), trans_type::no_trans);
                    return mat_col().add(P1).add(P2);
                }
                else
                {
                    //[0 op(A')] * [X1] = [op(A')*X2]
                    //[op(A) 0 ]   [X2]   [op(A)*X1]

                    Matrix X1   = X(colon(1,M), colon());
                    Matrix X2   = X(colon(M+1,N+M), colon());
                    Matrix P1   = m_At.mmul_right(std::move(X2), t);
                    Matrix P2   = m_A.mmul_right(std::move(X1), t);
                    return mat_col().add(P1).add(P2);
                };
            };
        };
        virtual Matrix mmul_right(Matrix&& X, trans_type t) const override
        {
            return mmul_right(X,t);
        };
};

}};

namespace matcl
{

//------------------------------------------------------------------
//                      linear_operator_data
//------------------------------------------------------------------
linear_operator_data::data_ptr 
linear_operator_data::convert(value_code new_val_code) const
{
    (void)new_val_code;
    return data_ptr(const_cast<linear_operator_data*>(this)->shared_from_this());
};

void linear_operator_data::mmul_right(const Matrix& X, trans_type t, Matrix& y) const
{
    Matrix z    = this->mmul_right(X,t);
    y(colon())  = z;
};

//------------------------------------------------------------------
//                      linear_operator
//------------------------------------------------------------------
linear_operator::linear_operator()
    :m_trans(trans_type_ext::no_trans)
{
    from_matrix(0.0);
};

linear_operator::linear_operator(const data_ptr& rep, trans_type_ext t)
    :m_trans(t), m_impl(rep)
{
    if (!m_impl)
        throw error::uninitialized_object_used("linear_operator");

};
linear_operator::linear_operator(data_ptr&& rep, trans_type_ext t)
    :m_impl(std::move(rep)), m_trans(t)
{
    if (!m_impl)
        throw error::uninitialized_object_used("linear_operator");
};

linear_operator::linear_operator(const linear_operator& mat)
    :m_impl(mat.m_impl), m_trans(mat.m_trans)
{};
linear_operator::linear_operator(linear_operator&& mat)
    :m_impl(std::move(mat.m_impl)), m_trans(mat.m_trans)
{};

linear_operator& linear_operator::operator=(const linear_operator& other) &
{
    m_impl  = other.m_impl;
    m_trans = other.m_trans;
    return *this;
}
linear_operator& linear_operator::operator=(linear_operator&& other) &
{
    m_impl  = std::move(other.m_impl);
    m_trans = other.m_trans;
    return *this;
}

linear_operator::~linear_operator()
{};

value_code linear_operator::get_value_code() const
{
    return m_impl->get_value_code();
};

linear_operator linear_operator::convert(value_code new_val_code) const
{
    data_ptr dp = m_impl->convert(new_val_code);
    return linear_operator(dp, m_trans);
};

Integer linear_operator::rows() const
{
    return (m_trans == trans_type_ext::no_trans || m_trans == trans_type_ext::conj)? 
                m_impl->rows() : m_impl->cols();
};

Integer linear_operator::cols() const
{
    return (m_trans == trans_type_ext::no_trans || m_trans == trans_type_ext::conj)? 
                m_impl->cols() : m_impl->rows();
};

bool linear_operator::is_hermitian() const
{
    return m_impl->is_hermitian();
};

Matrix linear_operator::mmul_right(const Matrix& X, trans_type t) const
{
    trans_type_ext te   = details::trans_manip::link_trans(m_trans, 
                            details::trans_manip::convert_trans(t));

    bool conj;
    trans_type t2       = details::trans_manip::convert_trans(te, conj);
    bool is_real        = (conj == false) && matrix_traits::is_float_real(get_value_code());

    if (conj == false || is_real == true)
        return m_impl->mmul_right(X,t2);

    Matrix y = m_impl->mmul_right(matcl::conj(X), trans_type::no_trans);
    return matcl::conj(std::move(y));
};

Matrix linear_operator::mmul_right(Matrix&& X, trans_type t) const
{
    trans_type_ext te   = details::trans_manip::link_trans(m_trans, 
                            details::trans_manip::convert_trans(t));

    bool conj;
    trans_type t2       = details::trans_manip::convert_trans(te, conj);
    bool is_real        = (conj == false) && matrix_traits::is_float_real(get_value_code());

    if (conj == false || is_real == true)
        return m_impl->mmul_right(std::move(X),t2);

    Matrix y = m_impl->mmul_right(matcl::conj(std::move(X)), trans_type::no_trans);
    return matcl::conj(std::move(y));
};
void linear_operator::mmul_right(const Matrix& X, trans_type t, Matrix& y) const
{
    trans_type_ext te   = details::trans_manip::link_trans(m_trans, 
                            details::trans_manip::convert_trans(t));

    bool conj;
    trans_type t2       = details::trans_manip::convert_trans(te, conj);
    bool is_real        = (conj == false) && matrix_traits::is_float_real(get_value_code());

    if (conj == false || is_real == true)
        return m_impl->mmul_right(X,t2,y);

    Matrix z = m_impl->mmul_right(matcl::conj(X), trans_type::no_trans);
    y(colon()) = matcl::conj(std::move(z));
};

void linear_operator::from_matrix(const Matrix& mat)
{
    m_impl = data_ptr(new details::linear_operator_mat(mat));
    m_trans = trans_type_ext::no_trans;
};

void linear_operator::from_unitary(const unitary_matrix& mat)
{
    if (mat.get_impl()->is_matrix())
        m_impl = data_ptr(new details::linear_operator_mat(mat.to_matrix()));
    else
        m_impl = data_ptr(new details::linear_operator_umat(mat));

    m_trans = trans_type_ext::no_trans;
};
void linear_operator::from_linsolve(const linsolve_obj& mat)
{
    m_impl  = data_ptr(new details::linear_operator_linsolve(mat));
    m_trans = trans_type_ext::no_trans;
};

linear_operator matcl::linop_trans(const linear_operator& mat)
{
    linear_operator ret(mat);
    ret.set_trans(details::trans_manip::link_trans(mat.get_trans(), trans_type_ext::trans));
    return ret;
}

linear_operator matcl::linop_ctrans(const linear_operator& mat)
{
    linear_operator ret(mat);
    ret.set_trans(details::trans_manip::link_trans(mat.get_trans(), trans_type_ext::conj_trans));
    return ret;
}

linear_operator matcl::linop_conj(const linear_operator& mat)
{
    linear_operator ret(mat);
    ret.set_trans(details::trans_manip::link_trans(mat.get_trans(), trans_type_ext::conj));
    return ret;
};

linear_operator matcl::linop_trans(const linear_operator& mat, trans_type t)
{
    linear_operator ret(mat);
    ret.set_trans(details::trans_manip::link_trans(mat.get_trans(), details::trans_manip::convert_trans(t)));
    return ret;
};

linear_operator matcl::linop_plus(const linear_operator& A, const linear_operator& B)
{
    using impl_ptr = linear_operator::data_ptr;

    linear_operator ret(impl_ptr(new details::bin_linop(A, B, details::bin_linop::op_plus)));
    return ret;
};

linear_operator matcl::linop_minus(const linear_operator& A, const linear_operator& B)
{
    using impl_ptr = linear_operator::data_ptr;

    linear_operator ret(impl_ptr(new details::bin_linop(A, B, details::bin_linop::op_minus)));
    return ret;
};

linear_operator matcl::linop_uminus(const linear_operator& A)
{
    using impl_ptr = linear_operator::data_ptr;

    linear_operator ret(impl_ptr(new details::un_linop(A, details::un_linop::op_uminus)));
    return ret;
};

linear_operator matcl::linop_mmul(const linear_operator& A, const linear_operator& B,trans_type tA,trans_type tB)
{
    using impl_ptr = linear_operator::data_ptr;

    linear_operator ret(impl_ptr(new details::bin_linop(A, B, tA, tB)));
    return ret;
};

linear_operator matcl::linop_scale(const Matrix& alpha, const linear_operator& A)
{
    using impl_ptr = linear_operator::data_ptr;

    linear_operator ret(impl_ptr(new details::un_linop(A, alpha)));
    return ret;
};

linear_operator matcl::linop_symprod(const linear_operator& A, bool trans)
{
    using impl_ptr = linear_operator::data_ptr;

    linear_operator ret(impl_ptr(new details::symprod_linop(A, trans)));
    return ret;
};

linear_operator matcl::linop_symsum(const linear_operator& A)
{
    if (A.rows() != A.cols())
        throw error::square_matrix_required(A.rows(), A.cols());

    using impl_ptr = linear_operator::data_ptr;

    linear_operator ret(impl_ptr(new details::symsum_linop(A)));
    return ret;
};

linear_operator matcl::linop_symcat(const linear_operator& A, bool trans)
{
    using impl_ptr = linear_operator::data_ptr;

    linear_operator ret(impl_ptr(new details::symcat_linop(A, trans)));
    return ret;
};

};
