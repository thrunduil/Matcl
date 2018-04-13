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

#if 0
TODO

#include "matcl-linalg/linear_eq/ksp_solver.h"
#include "matcl-linalg/iterative_linear_eq/ksp_solver_impl.h"
#include "matcl-linalg/petsc_utils/petsc_option.h"
#include "matcl-linalg/linear_eq/linsolve_object.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-linalg/special_matrices/struct_flag_linalg.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/lib_functions/func_matrix.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-linalg/linear_eq/preconditioner.h"

namespace matcl { namespace details
{

struct ksp_solver_linsolve : public linsolve_obj_data
{
    ksp_solver   m_impl;

    ksp_solver_linsolve(const ksp_solver& impl)
        :m_impl(impl)
    {};

    virtual ~ksp_solver_linsolve(){};

    virtual bool is_modified() const override           { return false; };
    virtual bool is_direct() const override             { return false; };
    virtual Integer rows() const override               { return m_impl.rows(); };
    virtual Integer cols() const override               { return m_impl.cols(); };
    virtual value_code get_value_code() const override  { return m_impl.get_value_code(); };
    virtual bool all_finite() const override            { return m_impl.all_finite(); };

    virtual bool is_hermitian() const override
    {
        if (m_impl.is_from_matrix() == false)
            return false;

        bool is_real        = matcl::matrix_traits::is_real(get_value_code());
        const Matrix& tmp   = m_impl.get_matrix();
        return tmp.get_struct().is_hermitian(tmp.is_square(), is_real);
    };
    virtual bool is_posdef() const override
    {
        if (m_impl.is_from_matrix() == false)
            return false;

        return matcl::is_posdef(m_impl.get_matrix().get_struct());
    }

    virtual matcl::Matrix solve(const Matrix& X, trans_type tA) const override
    {
        Matrix Y;
        bool has_sol;
        tie(Y,has_sol) = m_impl.solve(X, tA);

        if (has_sol == false)
            throw error::iterative_solver_not_converged();

        return Y;
    };
    virtual matcl::Matrix solve(Matrix&& X, trans_type tA) const override
    {
        Matrix Y;
        bool has_sol;
        tie(Y,has_sol) = m_impl.solve(std::move(X), tA);

        if (has_sol == false)
            throw error::iterative_solver_not_converged();

        return Y;
    };

    virtual Matrix mmul_right(const Matrix& X, trans_type t) const override
    {
        if (m_impl.is_from_matrix())
            return mmul(m_impl.get_matrix(), X, t);
        else
            return m_impl.get_linear_operator().mmul_right(X, t);
    };

    virtual Matrix mmul_right(Matrix&& X, trans_type t) const override
    {
        if (m_impl.is_from_matrix())
            return mmul(m_impl.get_matrix(), std::move(X), t);
        else
            return m_impl.get_linear_operator().mmul_right(std::move(X), t);
    };

    virtual Matrix mmul_left(const Matrix& X, trans_type t) const override
    {
        if (m_impl.is_from_matrix())
            return mmul(X, m_impl.get_matrix(), trans_type::no_trans, t);
        else
            return mmul_left_op(X, t);
    }
    virtual Matrix mmul_left(Matrix&& X, trans_type t) const override
    {
        if (m_impl.is_from_matrix())
            return mmul(std::move(X), m_impl.get_matrix(), trans_type::no_trans, t);
        else
        {
            Matrix X2(std::move(X));
            return mmul_left_op(X2, t);
        }
    }

    Matrix mmul_left_op(const Matrix& X, trans_type t) const
    {
        if (t == trans_type::no_trans)
        {
            //X * A = Y -> A'*X' = Y'
            Matrix Y = mmul_right(ctrans(X), trans_type::conj_trans);
            return ctrans(std::move(Y));
        }
        else
        {
            //X * op(A) = Y -> A*op(X) = op(Y)
            Matrix Y = mmul_right(trans(X,t), trans_type::no_trans);
            return trans(std::move(Y), t);
        };
    };

    virtual matcl::Matrix base_matrix() const override
    {
        if (m_impl.is_from_matrix())
            return m_impl.get_matrix();

        Matrix I    = matcl::speye(this->cols(), this->cols(), this->get_value_code());
        return mmul_right(std::move(I), trans_type::no_trans);
    };

    virtual data_ptr convert(value_code new_val_code) const override
    {
        //conversion is not allowed;

        (void)new_val_code;

        const linsolve_obj_data* base_ptr = this;

        data_ptr ptr(const_cast<linsolve_obj_data*>(base_ptr)->shared_from_this());
        return ptr;
    };

    virtual Real log_det() const override
    {
        throw error::logdet_not_available();
    };
};

}};

namespace matcl
{

//------------------------------------------------------------------------
//                      ksp_solver
//------------------------------------------------------------------------
ksp_solver::ksp_solver()
    :m_impl(new impl_type(Matrix(1.0),options()))
{};

ksp_solver::ksp_solver(const Matrix& A, const options& opts)
    :m_impl(new impl_type(A, opts))
{}
ksp_solver::ksp_solver(const linear_operator& A, const options& opts)
    :m_impl(new impl_type(A, opts))
{}
ksp_solver::ksp_solver(const linsolve_obj& A)
{
    precond p(A);
    *this = p;
}


ksp_solver::ksp_solver(const Matrix& A, const Matrix& B, const options& opts)
    :m_impl(new impl_type(A, B, opts))
{}
ksp_solver::ksp_solver(const linear_operator& A, const Matrix& B, const options& opts)
    :m_impl(new impl_type(A, B, opts))
{}

ksp_solver::~ksp_solver()
{}

ksp_solver& ksp_solver::operator()(const Matrix& A, const options& opts)
{
    *this = ksp_solver(A, opts);
    return *this;
}
ksp_solver& ksp_solver::operator()(const linear_operator& A, const options& opts)
{
    *this = ksp_solver(A, opts);
    return *this;
}
ksp_solver& ksp_solver::operator()(const Matrix& A, const Matrix& B, const options& opts)
{
    *this = ksp_solver(A, B, opts);
    return *this;
}
ksp_solver& ksp_solver::operator()(const linear_operator& A, const Matrix& B, const options& opts)
{
    *this = ksp_solver(A, B, opts);
    return *this;
}

ksp_solver& ksp_solver::operator()(const Matrix& A)
{
    *this = ksp_solver(A, this->options_set());
    return *this;
}
ksp_solver& ksp_solver::operator()(const linear_operator& A)
{
    *this = ksp_solver(A, this->options_set());
    return *this;
}
ksp_solver& ksp_solver::operator()(const Matrix& A, const Matrix& B)
{
    *this = ksp_solver(A, B, this->options_set());
    return *this;
}
ksp_solver& ksp_solver::operator()(const linear_operator& A, const Matrix& B)
{
    *this = ksp_solver(A, B, this->options_set());
    return *this;
}

void ksp_solver::reset_operator(const Matrix& A)
{
    m_impl->reset_operator(A);
};
void ksp_solver::reset_operator(const linear_operator& A)
{
    m_impl->reset_operator(A);
}

void ksp_solver::rebuild_precond(const Matrix& B)
{
    m_impl->rebuild_preconditioner(B, this->options_set());
}
void ksp_solver::rebuild_precond(const Matrix& B, const options& opts)
{
    m_impl->rebuild_preconditioner(B, opts);
}

void ksp_solver::set_notifier(const notifier& notif)
{
    m_impl->set_notifier(notif);
}

void ksp_solver::set_nullspace_right(const Matrix& Nr, bool with_const)
{
    return m_impl->set_nullspace_right(Nr, with_const);
};
void ksp_solver::set_nullspace_left(const Matrix& Nr, bool with_const)
{
    return m_impl->set_nullspace_left(Nr, with_const);
};
void ksp_solver::remove_nullspace_right()
{
    return m_impl->remove_nullspace_right();
};
void ksp_solver::remove_nullspace_left()
{
    return m_impl->remove_nullspace_left();
};

void ksp_solver::set_options(const options& opts)
{
    m_impl->reset_options(opts);
}

void ksp_solver::set_tolerances(Integer maxit, Real rtol, Real atol, Real dtol)
{
    m_impl->set_tolerances(maxit, rtol, atol, dtol);
};

void ksp_solver::set_precond(const linsolve_obj& precond)
{
    m_impl->set_preconditioner(precond);
};
void ksp_solver::set_precond(const ksp_solver& precond)
{
    m_impl->add_missing(precond.options_missing());
    m_impl->set_preconditioner(precond);
};
void ksp_solver::set_precond(const linear_operator& precond)
{
    m_impl->set_preconditioner(precond);
}
options ksp_solver::options_missing() const
{
    return m_impl->options_missing();
};
options ksp_solver::options_set() const
{
    return m_impl->options_set();
};
KSP ksp_solver::get_KSP() const
{
    return m_impl->get_KSP();
};

std::string ksp_solver::solver_name() const
{
    return m_impl->get_solver();
};

std::string ksp_solver::precond_name() const
{
    return m_impl->get_precond();
};

std::string ksp_solver::solver_name(petsc_solver sol)
{
    return petsc::converter_solver().to_string((Integer)sol);
};

petsc_solver ksp_solver::solver_code(const std::string& sol_name)
{
    return petsc::converter_solver().to_code(sol_name);
};

tuple<Matrix,bool> ksp_solver::solve(const matcl::Matrix& b, trans_type tA) const
{
    return m_impl->solve(*this, b, b, tA, false);
}

tuple<Matrix,bool> ksp_solver::solve(const matcl::Matrix& b, const Matrix& start, trans_type tA) const
{
    return m_impl->solve(*this, b, start, tA, true);
}

linear_operator ksp_solver::get_linear_operator() const
{
    return m_impl->get_linear_operator();
}

Matrix ksp_solver::get_matrix() const
{
    return m_impl->get_matrix();
}

bool ksp_solver::is_from_matrix() const
{
    return m_impl->is_from_matrix() == false;
};

Matrix ksp_solver::get_precond_matrix() const
{
    return m_impl->get_preconditioner_matrix();
};

bool ksp_solver::is_precond_from_matrix() const
{
    return m_impl->is_preconditioner_from_matrix();
};

linsolve_obj ksp_solver::to_linsolve_object2() const
{
    m_impl->finalize_build();
    using rep_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(rep_ptr(new details::ksp_solver_linsolve(*this)));
};

value_code ksp_solver::get_value_code() const
{
    return m_impl->get_value_code();
}

Integer ksp_solver::rows() const
{
    return m_impl->rows();
};

Integer ksp_solver::cols() const
{
    return m_impl->cols();
};

bool ksp_solver::all_finite() const
{
    return m_impl->all_finite();
};

void ksp_solver::set_option(const std::string& opt, const std::string& val)
{
    return m_impl->set_option(opt, val);
};

void ksp_solver::set_id(const std::string& id)
{
    m_impl->set_id(id);
};

const std::string& ksp_solver::get_id() const
{
    return m_impl->get_id();
};

//------------------------------------------------------------------------
//                      ksp_solver_ls
//------------------------------------------------------------------------
ksp_solver_ls::ksp_solver_ls()
    :ksp_solver()
{
    m_impl->set_ls();
};

ksp_solver_ls::ksp_solver_ls(const Matrix& A, const options& opts)
    :ksp_solver(A, opts)
{
    m_impl->set_ls();
};
ksp_solver_ls::ksp_solver_ls(Matrix&& A, const options& opts)
    :ksp_solver(std::move(A), std::move(opts))
{
    m_impl->set_ls();
};
ksp_solver_ls::ksp_solver_ls(const linear_operator& A, const options& opts)
    :ksp_solver(A, opts)
{
    m_impl->set_ls();
};

ksp_solver_ls::~ksp_solver_ls()
{};

Matrix ksp_solver_ls::get_solution_std() const
{
    return m_impl->get_solution_std();
};
}

#endif