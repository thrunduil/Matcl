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

#include "matcl-linalg/special_matrices/unitary_matrix.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-matrep/IO/matrix_io.h"
#include "matcl-linalg/decompositions/householder_q.h"
#include "matcl-linalg/utils/linalg_utils.h"

#include <algorithm>

namespace matcl { namespace details
{

namespace mdyp = matcl::dynamic::functions;

static bool is_real_code(value_code vc)
{
    switch (vc)
    {
        case value_code::v_integer:
        case value_code::v_float:
        case value_code::v_real:
            return true;
        default:
            return false;
    };
};
class unitary_matrix_from_mat : public unitary_matrix_data
{
    private:
        matcl::Matrix   m_mat;

    public:
        unitary_matrix_from_mat(const matcl::Matrix& mat)   :m_mat(mat) {};
        unitary_matrix_from_mat(matcl::Matrix&& mat)        :m_mat(std::move(mat)) {};

        virtual ~unitary_matrix_from_mat() override {};

    public:
        virtual Integer rows() const override
        {
            return m_mat.rows();
        };

        virtual Integer cols() const override
        {
            return m_mat.cols();
        };
        virtual value_code get_value_code() const override
        {
            return m_mat.get_value_code();
        };
        virtual ti::ti_object get_type() const override
        {
            return m_mat.get_type();
        };

        virtual bool all_finite() const override
        {
            return true;
        };

        virtual bool is_matrix() const override
        { 
            return true;
        };

        virtual bool is_id() const override
        { 
            return m_mat.get_struct().is_id(); 
        };

        virtual void to_matrix(matcl::Matrix& ret) const override
        {
            ret = m_mat;
        };

        virtual data_ptr convert(value_code new_val_code) const override
        {
            data_ptr ret(new unitary_matrix_from_mat(details::convert_value(m_mat, new_val_code)));
            return ret;
        };

        virtual void mult_right(matcl::Matrix& ret, const matcl::Matrix& mat, 
                                trans_type t_u) const override
        {
            ret = mmul(m_mat, mat, t_u);
        };
        
        virtual void mult_left(matcl::Matrix& ret, const matcl::Matrix& mat, 
                                trans_type t_u) const override
        {
            ret = mmul(mat, m_mat, trans_type::no_trans, t_u);
        };

        virtual serialization_helper<unitary_matrix_data>*
                                get_serialization_helper() const override
        {
            return serialization_helper<unitary_matrix_data>
                    ::get<unitary_matrix_from_mat>("unitary_matrix_from_mat");
        };

        virtual void save(oarchive& os) const override
        {
            matcl::save(os, m_mat);
        };
        virtual void save(std::ostream& os) const override
        {
            os << m_mat;
        }

        static unitary_matrix_data* load(std::istream& is)
        {
            Matrix ret;
            is >> ret;
            return new unitary_matrix_from_mat(ret);
        };

        static unitary_matrix_data* load(iarchive& ar)
        {
            Matrix ret;
            matcl::load(ar, ret);
            return new unitary_matrix_from_mat(ret);
        };
};

class unitary_matrix_mult : public unitary_matrix_data
{
    private:
        unitary_matrix  m_A;
        unitary_matrix  m_B;
        trans_type      m_tA;
        trans_type      m_tB;

    public:
        unitary_matrix_mult(const unitary_matrix& A, const unitary_matrix& B, trans_type t_A, trans_type t_B)
            : m_A(A), m_B(B), m_tA(t_A), m_tB(t_B)
        {};

        virtual ~unitary_matrix_mult() override {};

    public:
        virtual Integer rows() const override
        {
            return (m_tA == trans_type::no_trans)? m_A.rows() : m_A.cols();
        };
        virtual Integer cols() const override
        {
            return (m_tB == trans_type::no_trans)? m_B.cols() : m_B.rows();
        };
        virtual bool all_finite() const override
        {
            return true;
        };

        virtual value_code get_value_code() const override
        {
            value_code v1 = m_A.get_value_code();
            value_code v2 = m_B.get_value_code();
            return matrix_traits::unify_value_types(v1,v2);
        };
        virtual ti::ti_object get_type() const override
        {
            ti::ti_object t1 = m_A.get_type();
            ti::ti_object t2 = m_B.get_type();

            ti::ti_object ret_ti = ti::get_return_ti<ti::ti_object>(mdyp::op_mul::eval(), t1, t2);
            return ret_ti;
        };

        virtual void to_matrix(matcl::Matrix& ret) const override
        {
            ret = mmul(m_A.to_matrix(),m_B.to_matrix(),m_tA, m_tB);
        };

        virtual data_ptr convert(value_code new_val_code) const override
        {           
            unitary_matrix Ac = m_A.convert(new_val_code);
            unitary_matrix Bc = m_B.convert(new_val_code);

            data_ptr ret(new unitary_matrix_mult(Ac, Bc, m_tA, m_tB));
            return ret;
        };

        //m_A * m_B * mat
        virtual void mult_right(matcl::Matrix& ret, const matcl::Matrix& mat, 
                                trans_type t_u) const override
        {
            switch(t_u)
            {
                case trans_type::no_trans:
                {
                    Matrix tmp = mmul(m_B, mat, m_tB, trans_type::no_trans);
                    ret = mmul(m_A, tmp, m_tA, trans_type::no_trans);
                    return;
                }
                case trans_type::conj_trans:
                {
                    //(op(A)*op(B))^CT*X = op(B)^CT * op(A)^CT * X
                    Matrix tmp = mmul(ctrans(m_A), mat, m_tA, trans_type::no_trans);
                    ret = mmul(ctrans(m_B), tmp, m_tB, trans_type::no_trans);
                    return;
                }
                case trans_type::trans:
                {
                    //(op(A)*op(B))^T*X = op(B)^T * op(A)^T * X
                    Matrix tmp = mmul(trans(m_A), mat, m_tA, trans_type::no_trans);
                    ret = mmul(trans(m_B), tmp, m_tB, trans_type::no_trans);
                    return;
                }  
                default:
                    matcl_assert(0,"unknown case");
                    throw error::error_general("invalid case");
            }
        };
        
        virtual void mult_left(matcl::Matrix& ret, const matcl::Matrix& mat, 
                                trans_type t_u) const override
        {
            //mat * A * B
            switch(t_u)
            {
                case trans_type::no_trans:
                {
                    Matrix tmp = mmul(mat, m_A, trans_type::no_trans, m_tA);
                    ret = mmul(tmp, m_B, trans_type::no_trans, m_tB);
                    return;
                }
                case trans_type::conj_trans:
                {
                    //X * (op(A)*op(B))^CT = X * op(B)^CT * op(A)^CT
                    Matrix tmp = mmul(mat, ctrans(m_B), trans_type::no_trans, m_tB);
                    ret = mmul(tmp, ctrans(m_A), trans_type::no_trans, m_tA);
                    return;
                }
                case trans_type::trans:
                {
                    //X * (op(A)*op(B))^T = X * op(B)^T * op(A)^T
                    Matrix tmp = mmul(mat, trans(m_B), trans_type::no_trans, m_tB);
                    ret = mmul(tmp, trans(m_A), trans_type::no_trans, m_tA);
                    return;
                }  
                default:
                    matcl_assert(0,"unknown case");
                    throw error::error_general("invalid case");
            }
        };

        virtual serialization_helper<unitary_matrix_data>*
                                get_serialization_helper() const override
        {
            return serialization_helper<unitary_matrix_data>
                    ::get<unitary_matrix_mult>("unitary_matrix_mult");
        };

        virtual void save(oarchive& os) const override
        {
            matcl::save(os, m_A);
            matcl::save(os, m_B);
            os << (Integer)m_tA;
            os << (Integer)m_tB;
        };
        static unitary_matrix_data* load(iarchive& ar)
        {
            unitary_matrix A, B;
            matcl::load(ar, A);
            matcl::load(ar, B);

            Integer tmp_tA, tmp_tB;
            ar >> tmp_tA;
            ar >> tmp_tB;

            trans_type t_A = (trans_type)tmp_tA;
            trans_type t_B = (trans_type)tmp_tB;
            return new unitary_matrix_mult(A,B,t_A,t_B);
        };
        virtual void save(std::ostream& os) const override
        {
            os << m_A;
            os << " ";
            os << m_B;
            os << " ";
            os << (Integer)m_tA;
            os << " ";
            os << (Integer)m_tB;
            os << " ";
        };

        static unitary_matrix_data* load(std::istream& is)
        {
            unitary_matrix A, B;
            is >> A;
            is >> B;

            Integer tmp_tA, tmp_tB;
            is >> tmp_tA;
            is >> tmp_tB;

            trans_type t_A = (trans_type)tmp_tA;
            trans_type t_B = (trans_type)tmp_tB;
            return new unitary_matrix_mult(A,B,t_A,t_B);
        };
};

class unitary_matrix_nan : public unitary_matrix_data
{
    private:
        Integer     m_M;
        Integer     m_N;
        value_code  m_vc;

    public:
        unitary_matrix_nan(Integer M, Integer N, value_code vc)
            : m_M(M), m_N(N), m_vc(vc)
        {
            m_vc    = matrix_traits::unify_value_types(m_vc, value_code::v_float);

            if (m_vc == value_code::v_object)
                throw error::object_value_type_not_allowed("unitary_matrix_nan");
        };

        virtual ~unitary_matrix_nan() override {};

    public:
        virtual Integer rows() const override
        {
            return m_M;
        };

        virtual Integer cols() const override
        {
            return m_N;
        };
        virtual value_code get_value_code() const override
        {
            return m_vc;
        };
        virtual ti::ti_object get_type() const override
        {
            return ti::ti_object_type(m_vc);
        };

        virtual bool all_finite() const override
        {
            return false;
        };

        virtual bool is_matrix() const override
        { 
            return true;
        };

        virtual bool is_id() const override
        { 
            return false;
        };

        virtual void to_matrix(matcl::Matrix& ret) const override
        {
            ret = make_nan_matrix(m_M, m_N, m_vc);
        };

        virtual data_ptr convert(value_code vc) const override
        {
            data_ptr ret(new unitary_matrix_nan(m_M, m_N, vc));
            return ret;
        };

        virtual void mult_right(matcl::Matrix& ret, const matcl::Matrix& mat, 
                                trans_type t_u) const override
        {
            Integer M       = (t_u == trans_type::no_trans) ? this->rows() : this->cols();
            Integer N       = mat.cols();
            value_code vc   = matrix_traits::unify_value_types(m_vc, mat.get_value_code());
            ret             = make_nan_matrix(M,N,vc);
        };
        
        virtual void mult_left(matcl::Matrix& ret, const matcl::Matrix& mat, 
                                trans_type t_u) const override
        {
            Integer M       = mat.rows();
            Integer N       = (t_u == trans_type::no_trans) ? this->cols() : this->rows();
            value_code vc   = matrix_traits::unify_value_types(m_vc, mat.get_value_code());
            ret             = make_nan_matrix(M,N,vc);
        };

        virtual serialization_helper<unitary_matrix_data>*
                                get_serialization_helper() const override
        {
            return serialization_helper<unitary_matrix_data>
                    ::get<unitary_matrix_nan>("unitary_matrix_nan");
        };

        virtual void save(oarchive& os) const override
        {
            os << m_M;
            os << m_N;
            os << (Integer)m_vc;
        };
        virtual void save(std::ostream& os) const override
        {
            os << " ";
            os << m_M << " ";
            os << m_N << " ";
            os << (Integer)m_vc << " ";
        }

        static unitary_matrix_data* load(std::istream& is)
        {
            Integer M, N, vc_i;

            is >> M;
            is >> N;
            is >> vc_i;

            value_code vc = (value_code)vc_i;
            return new unitary_matrix_nan(M,N,vc);
        };

        static unitary_matrix_data* load(iarchive& ar)
        {
            Integer M, N, vc_i;

            ar >> M;
            ar >> N;
            ar >> vc_i;

            value_code vc = (value_code)vc_i;

            return new unitary_matrix_nan(M,N,vc);
        };
};

static unitary_matrix_data* make_unitary_matrix_from_mat(const Matrix& mat, bool test_finite)
{
    if (test_finite == false)
        return new unitary_matrix_from_mat(mat);

    if (mat.all_finite() == false)
        return new unitary_matrix_nan(mat.rows(), mat.cols(), mat.get_value_code());
    else
        return new unitary_matrix_from_mat(mat);
};

}}

namespace matcl
{

ti::ti_object unitary_matrix_data::get_type() const
{
    return ti::ti_object_type(this->get_value_code());
};

unitary_matrix::unitary_matrix()
    :m_impl(new details::unitary_matrix_from_mat(1.0)), m_trans(trans_type_ext::no_trans)
{}

unitary_matrix::unitary_matrix(const Matrix& mat, bool test_finite)
    : m_impl(details::make_unitary_matrix_from_mat(mat,test_finite))
    , m_trans(trans_type_ext::no_trans)
{}

unitary_matrix::unitary_matrix(Matrix&& mat, bool test_finite)
    : m_impl(details::make_unitary_matrix_from_mat(std::move(mat), test_finite))
    , m_trans(trans_type_ext::no_trans)
{}

unitary_matrix::unitary_matrix(const unitary_matrix_data_ptr& rep)
    :m_impl(rep), m_trans(trans_type_ext::no_trans)
{
    if (!m_impl)
        throw error::uninitialized_object_used("unitary_matrix");
};
unitary_matrix::unitary_matrix(unitary_matrix_data_ptr&& rep)
    :m_impl(std::move(rep)), m_trans(trans_type_ext::no_trans)
{
    if (!m_impl)
        throw error::uninitialized_object_used("unitary_matrix");
};

unitary_matrix::unitary_matrix(const unitary_matrix& mat)
    :m_impl(mat.get_impl()), m_trans(mat.m_trans)
{};
unitary_matrix::unitary_matrix(unitary_matrix&& mat)
    :m_impl(std::move(mat.get_impl())), m_trans(mat.m_trans)
{};

unitary_matrix& unitary_matrix::operator=(const unitary_matrix& mat) &
{
    m_impl  = mat.get_impl();
    m_trans = mat.m_trans;
    return *this;
};
unitary_matrix& unitary_matrix::operator=(unitary_matrix&& mat) &
{
    m_impl  = std::move(mat.get_impl());
    m_trans = mat.m_trans;
    return *this;
};

unitary_matrix::~unitary_matrix()
{};

Integer unitary_matrix::rows() const
{
    return (m_trans == trans_type_ext::no_trans || m_trans == trans_type_ext::conj)? 
                m_impl->rows() : m_impl->cols();
};
Integer unitary_matrix::cols() const
{
    return (m_trans == trans_type_ext::no_trans || m_trans == trans_type_ext::conj)? 
                m_impl->cols() : m_impl->rows();
};
        
Integer unitary_matrix::length() const
{
    Integer r = this->rows();
    Integer c = this->cols();

    if (r == 0 || c == 0)
        return 0;

    return std::max(r,c);
};
                
Real unitary_matrix::numel() const
{
    return Real(this->rows()) * Real(this->cols());
};

bool unitary_matrix::all_finite() const
{
    return m_impl->all_finite();
};

bool unitary_matrix::is_empty() const
{
    Integer r = this->rows();
    Integer c = this->cols();

    if (r == 0 || c == 0)
        return true;
    else
        return false;
};

bool unitary_matrix::is_scalar() const
{
    Integer r = this->rows();
    Integer c = this->cols();

    if (r == 1 && c == 1)
        return true;
    else
        return false;
};

bool unitary_matrix::is_square() const
{
    Integer r = this->rows();
    Integer c = this->cols();

    if (r == c)
        return true;
    else
        return false;
};

bool unitary_matrix::is_vector() const
{
    Integer r = this->rows();
    Integer c = this->cols();

    if (r == 1 || c == 1 || r == 0 || c == 0)
        return true;
    else
        return false;
};

value_code unitary_matrix::get_value_code() const
{
    return m_impl->get_value_code();
};

ti::ti_object unitary_matrix::get_type() const
{
    return m_impl->get_type();
};

matcl::Matrix unitary_matrix::to_matrix() const
{
    matcl::Matrix ret;
    m_impl->to_matrix(ret);
    
    switch (m_trans)
    {
        case trans_type_ext::no_trans:
            return ret;
        case trans_type_ext::conj:
            return conj(ret);
        case trans_type_ext::trans:
            return trans(ret);
        case trans_type_ext::conj_trans:
            return ctrans(ret);
        default:
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
    };
};

unitary_matrix unitary_matrix::from_nan(Integer M, Integer N, value_code vc)
{
    impl_ptr impl(new details::unitary_matrix_nan(M,N,vc));
    return unitary_matrix(impl);
};

unitary_matrix unitary_matrix::convert(value_code new_val_code) const
{
    if (this->get_value_code() == new_val_code)
        return *this;

    impl_ptr ret = m_impl->convert(new_val_code);

    unitary_matrix out = unitary_matrix(ret);
    out.m_trans = this->m_trans;
    return out;
};

unitary_matrix matcl::ctrans(const unitary_matrix& mat)
{
    unitary_matrix ret(mat);
    ret.set_trans(details::trans_manip::link_trans(mat.get_trans(), trans_type_ext::conj_trans));
    return ret;
};

unitary_matrix matcl::trans(const unitary_matrix& mat)
{
    unitary_matrix ret(mat);
    ret.set_trans(details::trans_manip::link_trans(mat.get_trans(), trans_type_ext::trans));
    return ret;
};
unitary_matrix matcl::trans(const unitary_matrix& mat, trans_type t)
{
    unitary_matrix ret(mat);
    ret.set_trans(details::trans_manip::link_trans(mat.get_trans(), details::trans_manip::convert_trans(t)));
    return ret;
};

unitary_matrix matcl::conj(const unitary_matrix& mat)
{
    unitary_matrix ret(mat);
    ret.set_trans(details::trans_manip::link_trans(mat.get_trans(), trans_type_ext::conj));
    return ret;
};

matcl::Matrix matcl::operator*(const unitary_matrix& A, const matcl::Matrix& B)
{
    return mmul(A,B,trans_type::no_trans,trans_type::no_trans);
};
matcl::Matrix matcl::operator*(const unitary_matrix& A, matcl::Matrix&& B)
{
    return mmul(A,std::move(B),trans_type::no_trans,trans_type::no_trans);
};
matcl::Matrix matcl::operator*(const matcl::Matrix& A, const unitary_matrix& B)
{
    return mmul(A,B,trans_type::no_trans,trans_type::no_trans);
};
matcl::Matrix matcl::operator*(matcl::Matrix&& A, const unitary_matrix& B)
{
    return mmul(A,std::move(B),trans_type::no_trans,trans_type::no_trans);
};
unitary_matrix matcl::operator*(const unitary_matrix& A, const unitary_matrix& B)
{
    return mmul(A,B,trans_type::no_trans,trans_type::no_trans);
};

static void mmul_impl(Matrix& ret, const unitary_matrix& A, const matcl::Matrix& B, trans_type t_A0)
{
    trans_type_ext t_A = details::trans_manip::link_trans(A.get_trans(), 
                            details::trans_manip::convert_trans(t_A0));

    if (t_A == trans_type_ext::no_trans)
        return A.get_impl()->mult_right(ret, B, trans_type::no_trans);
    else if (t_A == trans_type_ext::trans)
        return A.get_impl()->mult_right(ret, B, trans_type::trans);
    else if (t_A == trans_type_ext::conj_trans)
        return A.get_impl()->mult_right(ret, B, trans_type::conj_trans);

    value_code vc_A = A.get_impl()->get_value_code();
    if (details::is_real_code(vc_A) == true)
        return A.get_impl()->mult_right(ret, B, trans_type::no_trans);

    // we need to eval conj(A) * op(B) = X => A * op(conj(B)) = conj(X)
    value_code vc_B = B.get_value_code();
    if (details::is_real_code(vc_B) == true)
    {
        A.get_impl()->mult_right(ret, B, trans_type::no_trans);
        ret = conj(ret);
        return;
    }
    else
    {
        A.get_impl()->mult_right(ret, conj(B), trans_type::no_trans);
        ret = conj(ret);
        return;
    }
};

static void mmul_impl(Matrix& ret, const matcl::Matrix& A, const unitary_matrix& B, trans_type t_B0)
{
    trans_type_ext t_B = details::trans_manip::link_trans(B.get_trans(), 
                            details::trans_manip::convert_trans(t_B0));

    if (t_B == trans_type_ext::no_trans)
        return B.get_impl()->mult_left(ret, A, trans_type::no_trans);
    else if (t_B == trans_type_ext::trans)
        return B.get_impl()->mult_left(ret, A, trans_type::trans);
    else if (t_B == trans_type_ext::conj_trans)
        return B.get_impl()->mult_left(ret, A, trans_type::conj_trans);

    value_code vc_B = B.get_impl()->get_value_code();
    if (details::is_real_code(vc_B) == true)
        return B.get_impl()->mult_left(ret, A, trans_type::no_trans);

    // we need to eval A * conj(B) = X => conj(A) * B = conj(X)
    value_code vc_A = A.get_value_code();
    if (details::is_real_code(vc_A) == true)
    {
        B.get_impl()->mult_left(ret, A, trans_type::no_trans);
        ret = conj(ret);
        return;
    }
    else
    {
        B.get_impl()->mult_left(ret, conj(A), trans_type::no_trans);
        ret = conj(ret);
        return;
    }
};

matcl::Matrix matcl::mmul(const unitary_matrix& A, const matcl::Matrix& B0, trans_type t_A, trans_type t_B)
{
    //increase refcont
    Matrix B(B0);

    if (t_B != trans_type::no_trans)
        B = trans(B,t_B);

    Matrix ret;
    mmul_impl(ret, A, B, t_A);
    return ret;
};
matcl::Matrix matcl::mmul(const unitary_matrix& A, matcl::Matrix&& B0, trans_type t_A, trans_type t_B)
{
    //increase refcont
    Matrix B(std::move(B0));

    if (t_B != trans_type::no_trans)
        B = trans(B,t_B);

    Matrix ret;
    mmul_impl(ret, A, B, t_A);
    return ret;
};
matcl::Matrix matcl::mmul(const matcl::Matrix& A0, const unitary_matrix& B, trans_type t_A, trans_type t_B)
{
    //increase refcont
    Matrix A(A0);

    if (t_A != trans_type::no_trans)
        A = trans(A,t_A);

    Matrix ret;
    mmul_impl(ret, A, B, t_B);
    return ret;
};
matcl::Matrix matcl::mmul(matcl::Matrix&& A0, const unitary_matrix& B, trans_type t_A, trans_type t_B)
{
    //increase refcont
    Matrix A(std::move(A0));

    if (t_A != trans_type::no_trans)
        A = trans(A,t_A);

    Matrix ret;
    mmul_impl(ret, A, B, t_B);
    return ret;
};

bool is_matrix(const unitary_matrix& A)
{
    return A.get_impl()->is_matrix();
};

bool is_id(const unitary_matrix& A)
{
    return A.get_impl()->is_id();
};

unitary_matrix matcl::mmul(const unitary_matrix& A, const unitary_matrix& B, trans_type t_A, trans_type t_B)
{
    bool isv    = A.all_finite() && B.all_finite();

    if (isv == false)
    {
        Integer M       = (t_A == trans_type::no_trans) ? A.rows() : A.cols();
        Integer N       = (t_B == trans_type::no_trans) ? B.cols() : B.rows();
        value_code vc   = matrix_traits::unify_value_types(A.get_value_code(), B.get_value_code());

        return unitary_matrix::from_nan(M,N,vc);
    };

    if (is_id(A) == true)
        return trans(B,t_B);

    if (is_id(B) == true)
        return trans(A,t_A);

    if (is_matrix(A) == true)
        return unitary_matrix(mmul(A.to_matrix(), B, t_A, t_B),false);

    if (is_matrix(B) == true)
        return unitary_matrix(mmul(A, B.to_matrix(), t_A, t_B),false);

    unitary_matrix ret;
    using impl_ptr = unitary_matrix::impl_ptr;
    ret.set_impl(impl_ptr(new details::unitary_matrix_mult(A, B, t_A, t_B)));
    return ret;
};

void matcl::save(oarchive& ar,const unitary_matrix& mat)
{
    Integer tr  = (Integer)mat.get_trans();

    ar << tr;
    
    serialization_helper<unitary_matrix_data>* out = mat.get_impl()->get_serialization_helper();

    out->save(ar, mat.get_impl().get());
};
void matcl::load(iarchive& ar,unitary_matrix& mat)
{
    Integer tr;
    ar >> tr;

    using impl_ptr = unitary_matrix::impl_ptr;

    //TODO: is this needed?
    //serialization_helper<unitary_matrix_data>;

    impl_ptr ptr(serialization_helper<unitary_matrix_data>::load(ar));
    
    mat = unitary_matrix(ptr);
    mat.set_trans((trans_type_ext)tr);
};

std::ostream& matcl::operator<<(std::ostream& os, const unitary_matrix& mat)
{
    os << " ";
    Integer tr  = (Integer)mat.get_trans();
    os << tr << " ";

    serialization_helper<unitary_matrix_data>* out = mat.get_impl()->get_serialization_helper();

    out->save(os,mat.get_impl().get());
    return os;
};

std::istream& matcl::operator>>(std::istream& is, unitary_matrix& mat)
{
    Integer tr;
    is >> tr;

    using impl_ptr = unitary_matrix::impl_ptr;
    impl_ptr ptr(serialization_helper<unitary_matrix_data>::load(is));
    
    mat = unitary_matrix(ptr);
    mat.set_trans((trans_type_ext)tr);

    return is;
};

unitary_matrix matcl::rand_unitary(Integer K)
{
    Matrix X = matcl::randn(K,K);

    unitary_matrix ret;
    details::rand_unitary_impl(ret,X);
    return ret;
};
unitary_matrix matcl::frand_unitary(Integer K)
{
    Matrix X = matcl::frandn(K,K);

    unitary_matrix ret;
    details::rand_unitary_impl(ret,X);
    return ret;
};
unitary_matrix matcl::crand_unitary(Integer K)
{
    Matrix X = matcl::crandn(K,K);

    unitary_matrix ret;
    details::rand_unitary_impl(ret,X);
    return ret;
};
unitary_matrix matcl::fcrand_unitary(Integer K)
{
    Matrix X = matcl::fcrandn(K,K);

    unitary_matrix ret;
    details::rand_unitary_impl(ret,X);
    return ret;
};
unitary_matrix matcl::rand_unitary(Integer K, value_code vc)
{
    Matrix X = matcl::randn(K,K,vc);

    unitary_matrix ret;
    details::rand_unitary_impl(ret,X);
    return ret;
};

};
