/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/matrix/matrix_concat.h"
#include "matcl-matrep/container/matrix_cons_data.h"
#include "matcl-internals/error/error_check_basic.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-matrep/container/matrix2.inl"
#include "matcl-matrep/visitors/insert_visitor.inl"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

#include <list>

namespace matcl
{

value_code details::mat_cons_data::link_value_types(value_code new_vt, value_code this_vt, bool initialized) const
{
    if (initialized == false)
        return new_vt;

    if (new_vt == this_vt)
        return this_vt;

    if (new_vt > this_vt)
        return link_value_types(this_vt,new_vt, true);

    //new_vt < this_vt
    switch (this_vt)
    {
        case value_code::v_integer:
            return value_code::v_integer;
        case value_code::v_real: 
            return value_code::v_real;
        case value_code::v_complex:
            return value_code::v_complex;
        case value_code::v_object:       
            return value_code::v_object;

        case value_code::v_float:
        {
            return (new_vt == value_code::v_integer)? value_code::v_real : value_code::v_float;
        }
        case value_code::v_float_complex:
        {
            return (new_vt == value_code::v_integer || new_vt == value_code::v_real)? 
                    value_code::v_complex : value_code::v_float_complex;
        }
        default:
        {
            matcl_assert(0,"unknown case");
            return value_code::v_integer;
        }
    };
};

details::inplace_type details::data_container::check_for_inplace_build(bool is_row, bool is_sparse, 
                                      matcl::value_code vt, ti::ti_object ti) const
{
    const struct_code sc = is_sparse? struct_code::struct_sparse : struct_code::struct_dense;

    switch(m_type)
    {
        case scal_int:
        case scal_real:
        case scal_float:
        case scal_complex:
        case scal_fcomplex:
            return inplace_type::cannot_continue;

        case type_matrix:
        {
            bool valid = m_matrix.is_unique();
            valid  = valid && m_matrix.get_value_code()== vt;
            valid  = valid && (vt != value_code::v_object || m_matrix.get_type() == ti);
            valid  = valid && m_matrix.get_struct_code() == sc;

            if (valid == true)
                return inplace_type::can_inplace;
            else
            {
                if (is_row == true && m_matrix.cols() == 0)
                    return inplace_type::can_continue;
                else if (m_matrix.rows() == 0)
                    return inplace_type::can_continue;
            };

            return inplace_type::cannot_continue;
        }
        case type_col:
        case type_row:
            return inplace_type::cannot_continue;
    };

    return inplace_type::cannot_continue;
};

mat_row::mat_row()
:m_rows(0),m_cols(0),m_nnz(0)
{};

mat_row::mat_row(const mat_row& mr)
:m_rows(mr.m_rows),m_cols(mr.m_cols),m_nnz(mr.m_nnz),m_data(mr.m_data)
{};

mat_row::mat_row(mat_row&& mr)
:m_rows(mr.m_rows),m_cols(mr.m_cols),m_nnz(mr.m_nnz),m_data(std::move(mr.m_data))
{};

mat_row& mat_row::operator=(const mat_row& mr)
{
    m_rows	= mr.m_rows;
    m_cols	= mr.m_cols;
    m_nnz	= mr.m_nnz;
    m_data	= mr.m_data;

    return *this;
};

mat_row& mat_row::operator=(mat_row&& mr)
{
    m_rows	= mr.m_rows;
    m_cols	= mr.m_cols;
    m_nnz	= mr.m_nnz;
    m_data	= std::move(mr.m_data);

    return *this;
};

mat_row::~mat_row()
{};

bool mat_row::is_initialized() const
{
    return m_data ? m_data->is_initialized() : false;
};

ti::ti_object mat_row::get_type() const
{
    return m_data? m_data->get_type() : ti::predefined::get_ti_int();
};

matcl::value_code mat_row::get_value_type() const
{
    return m_data? m_data->get_value_type() : value_code::v_integer;
};

void mat_row::make_unique()
{
    if (m_data.unique() == false && m_data.get() != nullptr)
        m_data = ref_ptr(new details::mat_cons_data(*m_data.get()));
};

template<class S>
void mat_row::add_scalar_impl(const S& val)
{
    if (!m_data)
        m_data = ref_ptr(new details::mat_cons_data());

    if (m_rows == 0)
    {
        m_rows = 1;
        m_cols += 1;
        m_nnz  += 1;
        m_data->add_scalar(val);
        return;
    };

    error::check_horzcat(m_rows, m_cols, 1, 1);
    m_cols += 1;
    m_nnz  += 1;
    m_data->add_scalar(val);

    return;
};

void mat_row::add_mat(const Matrix& mat)
{
    if (!m_data)
        m_data = ref_ptr(new details::mat_cons_data());

    Integer mat_rows        = mat.rows();
    Integer mat_cols        = mat.cols();
    Integer mat_nnz         = mat.structural_nnz();

    if (m_rows == 0)
    {
        m_rows  = mat_rows;
        m_cols  += mat_cols;
        m_nnz   += mat_nnz;
        m_data->add_matrix(mat);
        return;
    };

    error::check_horzcat(m_rows, m_cols, mat_rows, mat_cols);
    m_cols  += mat_cols;
    m_nnz   += mat_nnz;
    m_data->add_matrix(mat);
};

void mat_row::add_col(const mat_col& mc)
{
    if (!mc.m_data)
        return;

    if (!m_data)
        m_data = ref_ptr(new details::mat_cons_data());

    if (m_rows == 0)
    {
        m_rows = mc.m_rows;
        m_cols += mc.m_cols;
        m_nnz  += mc.m_nnz;
        m_data->add_col(mc);
        return;
    };

    error::check_horzcat(m_rows, m_cols, mc.m_rows, mc.m_cols);
    m_cols += mc.m_cols;
    m_nnz  += mc.m_nnz;
    m_data->add_col(mc);
};

void mat_row::add_row(const mat_row& mr)
{
    const details::mat_cons_data* mr_data = mr.m_data.get();
    if (!mr_data)
        return;

    using matrix_vector = details::mat_cons_data::matrix_vector;

    const matrix_vector& mr_vector = mr_data->m_vector;
    using iterator = matrix_vector::const_iterator;

    iterator pos = mr_vector.begin();
    while (pos != mr_vector.end())
    {
        switch (pos->m_type)
        {
            case details::data_container::type_matrix:
            {
                add_mat(pos->m_matrix);
                break;
            }
            case details::data_container::type_col:
            {
                add_col(pos->m_col);
                break;
            }
            case details::data_container::type_row:
            {
                add_row(pos->m_row);
                break;
            }
            case details::data_container::scal_int:
            {
                add_int(pos->m_int);
                break;
            }
            case details::data_container::scal_real:
            {
                add_real(pos->m_real);
                break;
            }
            case details::data_container::scal_float:
            {
                add_float(pos->m_float);
                break;
            }
            case details::data_container::scal_complex:
            {
                add_complex(pos->m_complex[0],pos->m_complex[1]);
                break;
            }
            case details::data_container::scal_fcomplex:
            {
                add_fcomplex(pos->m_fcomplex[0],pos->m_fcomplex[1]);
                break;
            }
            default:
            {
                matcl_assert(0,"unknown case");
                throw;
            }
        };		
        ++pos;
    };
    return;
};

mat_row& mat_row::add(const Matrix& mat)
{
    make_unique();
    add_mat(mat);
    return *this;
};

mat_row& mat_row::add(Matrix&& mat0)
{
    make_unique();

    Matrix mat(std::move(mat0));
    add_mat(std::move(mat));
    return *this;
};

mat_row& mat_row::operator,(const Matrix& mat)
{
    return add(mat);
};

mat_row& mat_row::operator,(Matrix&& mat0)
{
    return add(std::move(mat0));
};

mat_row& mat_row::add(const mat_row& mr)
{
    make_unique();
    add_row(mr);
    return *this;
};

mat_row& mat_row::add(mat_row&& mr)
{
    make_unique();
    add_row(std::move(mr));
    return *this;
};

mat_row& mat_row::add(const mat_col& mc)
{
    make_unique();
    add_col(mc);
    return *this;
};

mat_row& mat_row::add(mat_col&& mc)
{
    make_unique();
    add_col(std::move(mc));
    return *this;
};

mat_row& mat_row::operator,(const mat_row& mr)
{
    return add(mr);
};

mat_row& mat_row::operator,(mat_row&& mr)
{
    return add(std::move(mr));
};

mat_row& mat_row::operator,(const mat_col& mc)
{
    return add(mc);
};

mat_row& mat_row::operator,(mat_col&& mc)
{
    return add(std::move(mc));
};

void mat_row::add_int(Integer val)
{
    add_scalar_impl(val);
};

void mat_row::add_real(Real val)
{
    add_scalar_impl(val);
};

void mat_row::add_float(Float val)
{
    add_scalar_impl(val);
};

void mat_row::add_complex(Real re, Real im)
{
    add_scalar_impl(Complex(re,im));
};

void mat_row::add_fcomplex(Float re, Float im)
{
    add_scalar_impl(Float_complex(re,im));
};

mat_col::mat_col()
    :m_rows(0),m_cols(0),m_nnz(0)
{};

mat_col::~mat_col()
{};

mat_col::mat_col(const mat_col& mr)
:m_rows(mr.m_rows),m_cols(mr.m_cols),m_nnz(mr.m_nnz),m_data(mr.m_data)
{};

mat_col::mat_col(mat_col&& mr)
:m_rows(mr.m_rows),m_cols(mr.m_cols),m_nnz(mr.m_nnz),m_data(std::move(mr.m_data))
{};

mat_col& mat_col::operator=(const mat_col& mr)
{
    m_rows	= mr.m_rows;
    m_cols	= mr.m_cols;
    m_nnz	= mr.m_nnz;
    m_data	= mr.m_data;

    return *this;
};

mat_col& mat_col::operator=(mat_col&& mr)
{
    m_rows	= mr.m_rows;
    m_cols	= mr.m_cols;
    m_nnz	= mr.m_nnz;
    m_data	= std::move(mr.m_data);

    return *this;
};

ti::ti_object mat_col::get_type() const
{
    return m_data? m_data->get_type() : ti::predefined::get_ti_int();
};

matcl::value_code mat_col::get_value_type() const
{
    return m_data? m_data->get_value_type() : value_code::v_integer;
};

void mat_col::make_unique()
{
    if (m_data.unique() == false && m_data.get() != nullptr)
        m_data = ref_ptr(new details::mat_cons_data(*m_data.get()));
};

bool mat_col::is_initialized() const
{
    return m_data ? m_data->is_initialized() : false;
};

void mat_col::add_mat(const Matrix& mat)
{
    if (!m_data)
        m_data = ref_ptr(new details::mat_cons_data());

    Integer mat_rows        = mat.rows();
    Integer mat_cols        = mat.cols();
    Integer mat_nnz         = mat.structural_nnz();

    if (m_cols == 0)
    {
        m_cols = mat_cols;
        m_rows += mat_rows;
        m_nnz  += mat_nnz;
        m_data->add_matrix(mat);
        return;
    };

    error::check_vertcat(m_rows, m_cols, mat_rows, mat_cols);
    m_rows += mat_rows;
    m_nnz  += mat_nnz;
    m_data->add_matrix(mat);
};

void mat_col::add_col(const mat_col& mc)
{
    const details::mat_cons_data* mc_data = mc.m_data.get();
    if (!mc_data)
        return;

    using matrix_vector = details::mat_cons_data::matrix_vector;

    const matrix_vector& mc_vector = mc_data->m_vector;
    using iterator = matrix_vector::const_iterator;

    iterator pos = mc_vector.begin();
    while (pos != mc_vector.end())
    {
        switch (pos->m_type)
        {
            case details::data_container::type_matrix:
            {
                add_mat(pos->m_matrix);
                break;
            }
            case details::data_container::type_col:
            {
                add_col(pos->m_col);
                break;
            }
            case details::data_container::type_row:
            {
                add_row(pos->m_row);
                break;
            }
            case details::data_container::scal_int:
            {
                add_int(pos->m_int);
                break;
            }
            case details::data_container::scal_real:
            {
                add_real(pos->m_real);
                break;
            }
            case details::data_container::scal_float:
            {
                add_float(pos->m_float);
                break;
            }
            case details::data_container::scal_complex:
            {
                add_complex(pos->m_complex[0],pos->m_complex[1]);
                break;
            }
            case details::data_container::scal_fcomplex:
            {
                add_fcomplex(pos->m_fcomplex[0],pos->m_fcomplex[1]);
                break;
            }
            default:
            {
                matcl_assert(0,"unknown case");
            }
        };		
        ++pos;
    };
    return;
};

void mat_col::add_row(const mat_row& mr)
{
    if (!mr.m_data)
        return;

    if (!m_data)
        m_data = ref_ptr(new details::mat_cons_data());

    if (m_cols == 0)
    {
        m_cols = mr.m_cols;
        m_rows += mr.m_rows;
        m_nnz  += mr.m_nnz;
        m_data->add_row(mr);
        return;
    };

    error::check_vertcat(m_rows, m_cols, mr.m_rows, mr.m_cols);
    m_rows += mr.m_rows;
    m_nnz  += mr.m_nnz;
    m_data->add_row(mr);
};

template<class S>
void mat_col::add_scalar_impl(const S& val)
{
    if (!m_data)
        m_data = ref_ptr(new details::mat_cons_data());

    if (m_cols == 0)
    {
        m_cols = 1;
        m_rows += 1;
        m_nnz  += 1;
        m_data->add_scalar(val);
        return;
    };

    error::check_vertcat(m_rows, m_cols, 1, 1);
    m_rows += 1;
    m_nnz  += 1;
    m_data->add_scalar(val);
    return;
};

mat_col& mat_col::add(const Matrix& mat)
{
    make_unique();
    add_mat(mat);
    return *this;
};

mat_col& mat_col::add(Matrix&& mat0)
{
    make_unique();

    Matrix mat(std::move(mat0));
    add_mat(std::move(mat));
    return *this;
};

mat_col& mat_col::operator,(const Matrix& mat)
{
    return add(mat);
};

mat_col& mat_col::operator,(Matrix&& mat0)
{
    return add(std::move(mat0));
};

mat_col& mat_col::add(const mat_row& mr)
{
    make_unique();
    add_row(mr);
    return *this;
};

mat_col& mat_col::add(mat_row&& mr)
{
    make_unique();
    add_row(std::move(mr));
    return *this;
};

mat_col& mat_col::add(const mat_col& mc)
{
    make_unique();
    add_col(mc);
    return *this;
};

mat_col& mat_col::add(mat_col&& mc)
{
    make_unique();
    add_col(std::move(mc));
    return *this;
};

void mat_col::add_int(Integer val)
{
    add_scalar_impl(val);
};

void mat_col::add_real(Real val)
{
    add_scalar_impl(val);
};

void mat_col::add_float(Float val)
{
    add_scalar_impl(val);
};

void mat_col::add_complex(Real re, Real im)
{
    add_scalar_impl(Complex(re,im));
};

void mat_col::add_fcomplex(Float re, Float im)
{
    add_scalar_impl(Float_complex(re,im));
};

namespace details
{
    template<class out,class in>
    struct get_scalar
    {
        using ti_type = ti::ti_object;
        static out eval(ti_type, const in& val)
        {
            return gr::converter<out,in>::eval(val);
        };
    };

    template<class in>
    struct get_scalar<Object,in>
    {
        static Object eval(ti::ti_object ti, const in& val)
        {
            return gr::converter<Object,in>::eval(val,ti);
        };
    };

    template<class out>
    struct get_scalar<out,out>
    {
        using ti_type = ti::ti_object;
        static out eval(ti_type, const out& val)
        {
            return val;
        };
    };

    template<>
    struct get_scalar<Object,Object>
    {
        using ti_type = ti::ti_object;
        static Object eval(ti_type ti, const Object& val)
        {
            return gr::converter<Object,Object>::eval(val,ti);
        };
    };

    struct insert_functor:public details::extract_type_switch<void,insert_functor,true>
    {
        template<class M1, class mat_type>
        static void eval(const Matrix& , const M1& mat,mat_type& out, Integer r, Integer c)
        {
            using struct_type       = typename M1::struct_type;
            using M1_value          = typename M1::value_type;
            using struct_type_ret   = typename mat_type::struct_type;
            using matrix_value      = typename mat_type::value_type;

            static const bool is_valid = ((int)type_to_code<M1_value>::value 
                                          <= (int)type_to_code<matrix_value>::value);
            return insert_visitor_impl<mat_type,struct_type_ret,M1,struct_type,is_valid>
                        ::eval(out,r,c,mat);
        };

        template<class M1, class mat_type>
        static void eval_scalar(const Matrix& , const M1& mat,mat_type& out, Integer r, Integer c)
        {
            using M1_value          = M1;
            using struct_type_ret   = typename mat_type::struct_type;
            using matrix_value      = typename mat_type::value_type;

            static const bool is_valid = ((int)type_to_code<M1_value>::value 
                                          <= (int)type_to_code<matrix_value>::value);
            return insert_visitor_scal_impl<mat_type,struct_type_ret,M1,is_valid>
                        ::eval(out,r,c,mat);
        };
    };

    template<class matrix_type>
    struct sparse_matrix_constructor_row
    {
        using value_type = typename matrix_type::value_type;

        static void eval(matrix_type& out, const mat_row& mc, Integer col_start, Integer offset)
        {
            using matrix_vector = details::mat_cons_data::matrix_vector;
            using iterator      = matrix_vector::const_iterator;

            const matrix_vector& mc_vector = mc.m_data->m_vector;            

            Integer row_start = 0;

            iterator pos = mc_vector.begin() + offset;
            while (pos != mc_vector.end())
            {
                switch (pos->m_type)
                {
                    case details::data_container::type_matrix:
                    {
                        const Matrix& mat   = pos->m_matrix;
                        Integer c           = mat.cols();
                        add_matrix(out,row_start,col_start,mat);

                        col_start           += c;
                        break;
                    }
                    case details::data_container::scal_int:
                    {
                        Integer tmp         = pos->m_int;
                        add_scal(out,row_start,col_start,get_scalar<value_type,Integer>::eval(mc.get_type(),tmp));
                        col_start           += 1;
                        break;
                    }
                    case details::data_container::scal_real:
                    {
                        const Real& tmp     = pos->m_real;
                        add_scal(out,row_start,col_start,get_scalar<value_type,Real>::eval(mc.get_type(),tmp));
                        col_start           += 1;
                        break;
                    }
                    case details::data_container::scal_float:
                    {
                        const Float& tmp    = pos->m_float;
                        add_scal(out,row_start,col_start,get_scalar<value_type,Float>::eval(mc.get_type(),tmp));
                        col_start           += 1;
                        break;
                    }
                    case details::data_container::scal_complex:
                    {
                        const Complex& tmp  = Complex(pos->m_complex[0],pos->m_complex[1]);
                        add_scal(out,row_start,col_start,get_scalar<value_type,Complex>::eval(mc.get_type(),tmp));
                        col_start           += 1;
                        break;
                    }
                    case details::data_container::scal_fcomplex:
                    {
                        const Float_complex& tmp = Float_complex(pos->m_fcomplex[0],pos->m_fcomplex[1]);
                        add_scal(out,row_start,col_start,get_scalar<value_type,Float_complex>::eval(mc.get_type(),tmp));
                        col_start           += 1;
                        break;
                    }
                    case details::data_container::type_col:
                    {
                        Matrix mat          = pos->m_col.to_matrix();
                        Integer c           = mat.cols();
                        add_matrix(out,row_start,col_start,mat);

                        col_start           += c;
                        break;
                    }
                    case details::data_container::type_row:
                    {
                        matcl_assert(0,"this case should be removed");
                    }
                    default:
                    {
                        matcl_assert(0,"unknown case");
                    }
                };		

                ++pos;

                if (pos != mc_vector.end())
                {
                    raw::details::sparse_ccs<value_type>& rep = out.rep();
                    if (col_start+1 <= out.cols())
                        rep.ptr_c()[col_start+1] = rep.ptr_c()[col_start];
                };
            };

            return;
        };

        private:
            static void add_matrix(matrix_type& out, Integer row_start, Integer col_start, const Matrix& mat)
            {
                return insert_functor
                    ::make<const Matrix&,matrix_type&>(mat,out,row_start,col_start);
            };

            static void add_scal(matrix_type& out, Integer row_start, Integer col_start, const value_type& mat)
            {
                raw::details::sparse_ccs<value_type>& rep = out.rep();

                Integer k       = rep.ptr_c()[col_start+1]++;
                rep.ptr_r()[k]  = row_start;
                mrd::reset_helper(rep.ptr_x()[k],mat);
            };
    };

    template<class matrix_type>
    struct sparse_matrix_constructor_col
    {
        using value_type        = typename matrix_type::value_type;
        using const_matrix_type = mr::const_matrix<matrix_type>;

        static void eval(matrix_type& out, const mat_col& mc)
        {
            if (out.cols() == 1)
            {
                eval_col(out,mc);
                return;
            };

            matcl::mat_code mat_type = matrix_traits::mat_type_info_type<matrix_type>::matrix_code;

            using matrix_vector     = details::mat_cons_data::matrix_vector;            

            const matrix_vector& mc_vector = mc.m_data->m_vector;
            using iterator = matrix_vector::const_iterator;			

            std::list<const_matrix_type> mat_vec;

            iterator pos = mc_vector.begin();
            while(pos != mc_vector.end())
            {
                switch (pos->m_type)
                {
                    case details::data_container::type_matrix:
                    {
                        const Matrix& mat = pos->m_matrix;
                        Matrix new_mat = convert(mat,mat_type);
                        mat_vec.push_back(const_matrix_type(new_mat.get_impl<matrix_type>()));
                        break;
                    }
                    case details::data_container::type_row:
                    {
                        Matrix mat = pos->m_row.to_matrix();
                        Matrix new_mat = convert(mat,mat_type);
                        mat_vec.push_back(new_mat.get_impl<matrix_type>());
                        break;
                    }
                    case details::data_container::scal_int:
                    case details::data_container::scal_real:
                    case details::data_container::scal_float:
                    case details::data_container::scal_complex:
                    case details::data_container::scal_fcomplex:
                    case details::data_container::type_col:
                    {
                        matcl_assert(0,"this case should be removed");
                    }
                    default:
                    {
                        matcl_assert(0,"unknown case");
                    }
                };
                ++pos;
            };

            add_matrices(out, mat_vec);
            return;
        };

        private:
            static void eval_col(matrix_type& out, const mat_col& mc)
            {
                using matrix_vector = details::mat_cons_data::matrix_vector;
                using iterator      = matrix_vector::const_iterator;

                const matrix_vector& mc_vector = mc.m_data->m_vector;                

                iterator pos = mc_vector.begin();
                Integer row_start = 0;                

                while(pos != mc_vector.end())
                {
                    switch (pos->m_type)
                    {
                        case details::data_container::type_matrix:
                        {
                            const Matrix& mat   = pos->m_matrix;
                            Integer r           = mat.rows();
                            add_matrix(out,row_start,0,mat);
                            row_start           +=r;
                            break;
                        }
                        case details::data_container::type_row:
                        {
                            Matrix mat          = pos->m_row.to_matrix();
                            Integer r           = mat.rows();
                            add_matrix(out,row_start,0,mat);
                            row_start           +=r;
                            break;
                        }
                        case details::data_container::scal_int:
                        {
                            Integer tmp         = pos->m_int;
                            add_scal(out,row_start,0,get_scalar<value_type,Integer>::eval(mc.get_type(),tmp));
                            row_start           +=1;
                            break;
                        }
                        case details::data_container::scal_real:
                        {
                            add_scal(out,row_start,0,get_scalar<value_type,Real>::eval(mc.get_type(),pos->m_real));
                            row_start           +=1;
                            break;
                        }
                        case details::data_container::scal_float:
                        {
                            add_scal(out,row_start,0,get_scalar<value_type,Float>::eval(mc.get_type(),pos->m_float));
                            row_start           +=1;
                            break;
                        }
                        case details::data_container::scal_complex:
                        {
                            const Complex& tmp  = Complex(pos->m_complex[0],pos->m_complex[1]);
                            add_scal(out,row_start,0,get_scalar<value_type,Complex>::eval(mc.get_type(),tmp));
                            row_start           +=1;
                            break;
                        }
                        case details::data_container::scal_fcomplex:
                        {
                            const Float_complex& tmp = Float_complex(pos->m_fcomplex[0],pos->m_fcomplex[1]);
                            add_scal(out,row_start,0,get_scalar<value_type,Float_complex>::eval(mc.get_type(),tmp));
                            row_start   +=1;
                            break;
                        }				
                        case details::data_container::type_col:
                        {
                            matcl_assert(0,"this case should be removed");
                        }
                        default:
                        {
                            matcl_assert(0,"unknown case");
                        }
                    };
                    ++pos;
                };

                return;
            };
        
            static void add_matrix(matrix_type& out, Integer row_start, Integer col_start, const Matrix& mat)
            {
                return insert_functor
                    ::make<const Matrix&,matrix_type&>(mat,out,row_start,col_start);
            };
            
            static void add_scal(matrix_type& out, Integer row_start, Integer col_start, const value_type& mat)
            {
                raw::details::sparse_ccs<value_type>& rep = out.rep();

                Integer k       = rep.ptr_c()[col_start+1]++;
                rep.ptr_r()[k]  = row_start;
                mrd::reset_helper(rep.ptr_x()[k],mat);
            }
            
            static void add_matrices(matrix_type& out, const std::list<const_matrix_type>& mat_vec)
            {
                insert_sparse_cols<matrix_type>(out, mat_vec);
                return;
            };
    };

    template<class matrix_type>
    struct dense_matrix_constructor_row
    {
        using value_type = typename matrix_type::value_type;

        static void eval(matrix_type& out, Integer row_start, Integer col_start, const mat_row& mc,
                         Integer offset)
        {
            using matrix_vector = details::mat_cons_data::matrix_vector;
            using iterator      = matrix_vector::const_iterator ;

            const matrix_vector& mc_vector = mc.m_data->m_vector;            

            iterator pos = mc_vector.begin() + offset;

            while (pos != mc_vector.end())
            {
                switch (pos->m_type)
                {
                    case details::data_container::type_matrix:
                    {
                        const Matrix& mat = pos->m_matrix;
                        add_matrix(out,row_start,col_start,mat);

                        col_start += imult(mat.cols(),out.ld());
                        break;
                    }
                    case details::data_container::scal_int:
                    {
                        Integer tmp = pos->m_int;
                        add_scal(out,row_start,col_start,get_scalar<value_type,Integer>::eval(mc.get_type(),tmp));
                        col_start += out.ld();
                        break;
                    }
                    case details::data_container::scal_real:
                    {
                        const Real& tmp = pos->m_real;
                        add_scal(out,row_start,col_start,get_scalar<value_type,Real>::eval(mc.get_type(),tmp));
                        col_start += out.ld();
                        break;
                    }
                    case details::data_container::scal_float:
                    {
                        const Float& tmp = pos->m_float;
                        add_scal(out,row_start,col_start,get_scalar<value_type,Float>::eval(mc.get_type(),tmp));
                        col_start += out.ld();
                        break;
                    }
                    case details::data_container::scal_complex:
                    {
                        const Complex& tmp = Complex(pos->m_complex[0],pos->m_complex[1]);
                        add_scal(out,row_start,col_start,get_scalar<value_type,Complex>::eval(mc.get_type(),tmp));
                        col_start += out.ld();
                        break;
                    }
                    case details::data_container::scal_fcomplex:
                    {
                        const Float_complex& tmp = Float_complex(pos->m_fcomplex[0],pos->m_fcomplex[1]);
                        add_scal(out,row_start,col_start,get_scalar<value_type,Float_complex>::eval(mc.get_type(),tmp));
                        col_start += out.ld();
                        break;
                    }
                    case details::data_container::type_col:
                    {
                        dense_matrix_constructor_col<matrix_type>::eval(out,row_start,col_start,pos->m_col,0);
                        col_start += imult(pos->m_col.m_cols,out.ld());
                        break;
                    }
                    case details::data_container::type_row:
                    {
                        matcl_assert(0,"this case should be removed");
                    }
                    default:
                    {
                        matcl_assert(0,"unknown case");
                    }
                };		
                ++pos;
            };
            return;
        };

        private:
            static void add_matrix(matrix_type& out, Integer row_start, Integer col_start, const Matrix& mat)
            {
                return insert_functor
                    ::make<const Matrix&,matrix_type&>(mat,out,row_start,col_start);
            };
        
            static void add_scal(matrix_type& out, Integer row_start, Integer col_start, const value_type& mat)
            {
                mrd::reset_helper(out.ptr()[row_start+col_start],mat);
            };
    };

    template<class matrix_type>
    struct dense_matrix_constructor_col
    {
        using value_type = typename matrix_type::value_type;

        static void eval(matrix_type& out, Integer row_start, Integer col_start, const mat_col& mc,
                         Integer offset)
        {
            using matrix_vector = details::mat_cons_data::matrix_vector;
            using iterator      = matrix_vector::const_iterator;

            const matrix_vector& mc_vector = mc.m_data->m_vector;            

            iterator pos = mc_vector.begin() + offset;            

            while (pos != mc_vector.end())
            {
                switch (pos->m_type)
                {
                    case details::data_container::type_matrix:
                    {
                        const Matrix& mat   = pos->m_matrix;
                        Integer r           = mat.rows();
                        add_matrix(out,row_start,col_start,mat);

                        row_start += r;
                        break;
                    }
                    case details::data_container::scal_int:
                    {
                        Integer tmp = pos->m_int;
                        add_scal(out,row_start,col_start,get_scalar<value_type,Integer>::eval(mc.get_type(),tmp));
                        row_start += 1;
                        break;
                    }
                    case details::data_container::scal_real:
                    {
                        const Real& tmp = pos->m_real;
                        add_scal(out,row_start,col_start,get_scalar<value_type,Real>::eval(mc.get_type(),tmp));
                        row_start += 1;
                        break;
                    }
                    case details::data_container::scal_float:
                    {
                        const Float& tmp = pos->m_float;
                        add_scal(out,row_start,col_start,get_scalar<value_type,Float>::eval(mc.get_type(),tmp));
                        row_start += 1;
                        break;
                    }
                    case details::data_container::scal_fcomplex:
                    {
                        const Float_complex& tmp = Float_complex(pos->m_fcomplex[0],pos->m_fcomplex[1]);
                        add_scal(out,row_start,col_start,get_scalar<value_type,Float_complex>::eval(mc.get_type(),tmp));
                        row_start += 1;
                        break;
                    }
                    case details::data_container::scal_complex:
                    {
                        const Complex& tmp = Complex(pos->m_complex[0],pos->m_complex[1]);
                        add_scal(out,row_start,col_start,get_scalar<value_type,Complex>::eval(mc.get_type(),tmp));
                        row_start += 1;
                        break;
                    }
                    case details::data_container::type_row:
                    {
                        dense_matrix_constructor_row<matrix_type>::eval(out,row_start,col_start,pos->m_row,0);
                        row_start += pos->m_row.m_rows;
                        break;
                    }
                    case details::data_container::type_col:
                    {
                        matcl_assert(0,"this case should be removed");
                    }
                    default:
                    {
                        matcl_assert(0,"unknown case");
                    }
                };		
                ++pos;
            };
            return;
        };

        private:
            static void add_matrix(matrix_type& out, Integer row_start, Integer col_start,const Matrix& mat)
            {
                return insert_functor
                    ::make<const Matrix&,matrix_type&>(mat,out,row_start,col_start);
            };
        
            static void add_scal(matrix_type& out, Integer row_start, Integer col_start, const value_type& mat)
            {
                mrd::reset_helper(out.ptr()[row_start+col_start],mat);
            };
    };
};

const Matrix mat_row::to_matrix() const
{
    if (!m_data)
        return Matrix(raw::integer_sparse(ti::ti_empty(),0,0),true);

    size_t n_elem = m_data->m_vector.size();
    matcl::value_code m_value_code = this->m_data->get_value_type();

    if (m_rows == 0 || m_cols == 0 || m_nnz == 0 || n_elem == 0)
    {
        if (m_value_code == value_code::v_integer)
        {
            raw::integer_sparse out(ti::ti_empty(),m_rows,m_cols);
            return Matrix(out,false);
        }
        else if (m_value_code == value_code::v_object)
        {
            raw::object_sparse out(m_data->get_type(),m_rows,m_cols);
            return Matrix(out,false);
        }
        else if (m_value_code == value_code::v_float 
                 || m_value_code == value_code::v_float_complex)
        {
            raw::float_sparse out(ti::ti_empty(),m_rows,m_cols);
            return Matrix(out,false);
        }
        else
        {
            raw::real_sparse out(ti::ti_empty(),m_rows,m_cols);
            return Matrix(out,false);
        };
    };

    if (n_elem == 1)
    {
        const details::data_container& dc = m_data->m_vector[0];
        switch (dc.m_type)
        {
            case details::data_container::type_matrix:
            {
                Matrix out(dc.m_matrix);
                return out;
            }
            case details::data_container::type_col:
            {
                return dc.m_col.to_matrix();
            }
            case details::data_container::type_row:
            {
                return dc.m_row.to_matrix();
            }
            case details::data_container::scal_int:
            {
                return dc.m_int;
            }
            case details::data_container::scal_real:
            {
                return dc.m_real;
            }
            case details::data_container::scal_float:
            {
                return dc.m_float;
            }
            case details::data_container::scal_complex:
            {
                return Complex(dc.m_complex[0],dc.m_complex[1]);
            }
            case details::data_container::scal_fcomplex:
            {
                return Float_complex(dc.m_fcomplex[0],dc.m_fcomplex[1]);
            }
            default:
            {
                matcl_assert(0,"unknown case");
                throw;
            }
        };
    };

    Matrix out  = const_cast<mat_row*>(this)->build_matrix_inplace();
    return out;
};

Matrix mat_row::build_matrix_inplace()
{
    Real density    = Real(m_nnz)/(m_rows+1.)/(m_cols+1.);
    bool is_sparse  = (density < optim_params::max_sparse_density_min);

    matcl::value_code m_value_code = this->m_data->get_value_type();

    //only unique mat_row object can be updated inplace. Since mat_row cannot be
    //shared between threads, build_matrix() cannot be called by other threads.
    if (m_data.unique() == false)
        return build_matrix_new();

    auto& vec   = m_data->m_vector;
    size_t n    = vec.size();

    Integer mat_to_update = -1;

    for (size_t i = 0; i < n; ++i)
    {
        details::inplace_type it;
        it = vec[i].check_for_inplace_build(true, is_sparse, m_value_code, m_data->get_type());

        switch (it)
        {
            case details::inplace_type::can_inplace:
                mat_to_update = (Integer)i;
                goto exit_for_label;
            case details::inplace_type::can_continue:
                continue;
            default:
                return build_matrix_new();
        };
    };    

    if (mat_to_update == -1)
        return build_matrix_new();

  exit_for_label:

    Matrix mat = std::move(vec[mat_to_update].m_matrix);
    mat.set_struct(struct_flag());

    switch (m_value_code)
    {
        case value_code::v_integer:
        {
            if (is_sparse)
            {
                using MT    = raw::integer_sparse;
                MT& tmp     = mat.get_impl_unique<MT>();

                Integer col_start = tmp.cols();

                tmp.prepare_for_concat(this->m_rows, this->m_cols, this->m_nnz);
                details::sparse_matrix_constructor_row<MT>::eval(tmp,*this,col_start,mat_to_update+1);

                mat = Matrix(tmp,true);
            }
            else
            {
                using MT        = raw::integer_dense;
                MT& tmp         = mat.get_impl_unique<MT>();

                Integer col_start = imult(tmp.ld(), tmp.cols());

                tmp.prepare_for_concat(this->m_rows, this->m_cols);
                details::dense_matrix_constructor_row<MT>::eval(tmp,0,col_start,*this,mat_to_update+1);

                mat = Matrix(tmp,true);
            }
            break;
        }
        case value_code::v_float:
        {
            if (is_sparse)
            {
                using MT    = raw::float_sparse;
                MT& tmp     = mat.get_impl_unique<MT>();

                Integer col_start = tmp.cols();

                tmp.prepare_for_concat(this->m_rows, this->m_cols, this->m_nnz);
                details::sparse_matrix_constructor_row<MT>::eval(tmp,*this,col_start,mat_to_update+1);

                mat = Matrix(tmp,true);
            }
            else
            {
                using MT    = raw::float_dense;
                MT& tmp     = mat.get_impl_unique<MT>();

                Integer col_start = imult(tmp.ld(), tmp.cols());

                tmp.prepare_for_concat(this->m_rows, this->m_cols);
                details::dense_matrix_constructor_row<MT>::eval(tmp,0,col_start,*this,mat_to_update+1);

                mat = Matrix(tmp,true);
            }
            break;
        }
        case value_code::v_real:
        {
            if (is_sparse)
            {
                using MT    = raw::real_sparse;
                MT& tmp     = mat.get_impl_unique<MT>();

                Integer col_start = tmp.cols();

                tmp.prepare_for_concat(this->m_rows, this->m_cols, this->m_nnz);
                details::sparse_matrix_constructor_row<MT>::eval(tmp,*this,col_start, mat_to_update+1);

                mat = Matrix(tmp,true);
            }
            else
            {
                using MT    = raw::real_dense;
                MT& tmp     = mat.get_impl_unique<MT>();

                Integer col_start = imult(tmp.ld(), tmp.cols());

                tmp.prepare_for_concat(this->m_rows, this->m_cols);
                details::dense_matrix_constructor_row<MT>::eval(tmp,0,col_start,*this,mat_to_update+1);

                mat = Matrix(tmp,true);
            }
            break;
        }
        case value_code::v_float_complex:
        {
            if (is_sparse)
            {
                using MT    = raw::float_complex_sparse;
                MT& tmp     = mat.get_impl_unique<MT>();

                Integer col_start = tmp.cols();

                tmp.prepare_for_concat(this->m_rows, this->m_cols, this->m_nnz);
                details::sparse_matrix_constructor_row<MT>::eval(tmp,*this,col_start,mat_to_update+1);

                mat = Matrix(tmp,true);
            }
            else
            {
                using MT    = raw::float_complex_dense;
                MT& tmp     = mat.get_impl_unique<MT>();

                Integer col_start = imult(tmp.ld(), tmp.cols());

                tmp.prepare_for_concat(this->m_rows, this->m_cols);
                details::dense_matrix_constructor_row<MT>::eval(tmp,0,col_start,*this,mat_to_update+1);

                mat = Matrix(tmp,true);
            }
            break;
        }
        case value_code::v_complex:
        {
            if (is_sparse)
            {
                using MT    = raw::complex_sparse;
                MT& tmp     = mat.get_impl_unique<MT>();

                Integer col_start = tmp.cols();

                tmp.prepare_for_concat(this->m_rows, this->m_cols, this->m_nnz);
                details::sparse_matrix_constructor_row<MT>::eval(tmp,*this,col_start,mat_to_update+1);

                mat = Matrix(tmp,true);
            }
            else
            {
                using MT    = raw::complex_dense;
                MT& tmp     = mat.get_impl_unique<MT>();

                Integer col_start = imult(tmp.ld(), tmp.cols());

                tmp.prepare_for_concat(this->m_rows, this->m_cols);
                details::dense_matrix_constructor_row<MT>::eval(tmp,0,col_start,*this,mat_to_update+1);

                mat = Matrix(tmp,true);
            }
            break;
        }
        case value_code::v_object:
        {
            if (is_sparse)
            {
                using MT    = raw::object_sparse;
                MT& tmp     = mat.get_impl_unique<MT>();

                Integer col_start = tmp.cols();

                tmp.prepare_for_concat(this->m_rows, this->m_cols, this->m_nnz);
                details::sparse_matrix_constructor_row<MT>::eval(tmp,*this,col_start, mat_to_update+1);

                mat = Matrix(tmp,true);
            }
            else
            {
                using MT    = raw::object_dense;
                MT& tmp     = mat.get_impl_unique<MT>();

                Integer col_start = imult(tmp.ld(), tmp.cols());

                tmp.prepare_for_concat(this->m_rows, this->m_cols);
                details::dense_matrix_constructor_row<MT>::eval(tmp,0,col_start,*this,mat_to_update+1);

                mat = Matrix(tmp,true);
            }
            break;
        }
        default:
        {
            matcl_assert(0,"unknown value type code");
            throw;
        }
    };

    //matrix is updated, m_data must be modified in order to forbit other inplace updates    
    this->m_data = ref_ptr(new details::mat_cons_data());
    this->m_data->add_matrix(mat);

    return mat;
};

Matrix mat_row::build_matrix_new() const
{
    Real density    = Real(m_nnz)/(m_rows+1.)/(m_cols+1.);
    bool is_sparse  = (density < optim_params::max_sparse_density_min);

    matcl::value_code m_value_code = this->m_data->get_value_type();

    switch (m_value_code)
    {
        case value_code::v_integer:
        {
            if (is_sparse)
            {
                raw::integer_sparse out(ti::ti_empty(),m_rows,m_cols,m_nnz);
                details::sparse_matrix_constructor_row<raw::integer_sparse>::eval(out,*this,0,0);
                return Matrix(out,true);
            }
            else
            {
                raw::integer_dense out(ti::ti_empty(),m_rows,m_cols);
                details::dense_matrix_constructor_row<raw::integer_dense>::eval(out,0,0,*this,0);
                return Matrix(out,true);
            }
        }
        case value_code::v_real:
        {
            if (is_sparse)
            {
                raw::real_sparse out(ti::ti_empty(),m_rows,m_cols,m_nnz);
                details::sparse_matrix_constructor_row<raw::real_sparse>::eval(out,*this,0,0);
                return Matrix(out,true);
            }
            else
            {
                raw::real_dense out(ti::ti_empty(),m_rows,m_cols);
                details::dense_matrix_constructor_row<raw::real_dense>::eval(out,0,0,*this,0);
                return Matrix(out,true);
            }
        }
        case value_code::v_float:
        {
            if (is_sparse)
            {
                raw::float_sparse out(ti::ti_empty(),m_rows,m_cols,m_nnz);
                details::sparse_matrix_constructor_row<raw::float_sparse>::eval(out,*this,0,0);
                return Matrix(out,true);
            }
            else
            {
                raw::float_dense out(ti::ti_empty(),m_rows,m_cols);
                details::dense_matrix_constructor_row<raw::float_dense>::eval(out,0,0,*this,0);
                return Matrix(out,true);
            }
        }
        case value_code::v_complex:
        {
            if (is_sparse)
            {
                raw::complex_sparse out(ti::ti_empty(),m_rows,m_cols,m_nnz);
                details::sparse_matrix_constructor_row<raw::complex_sparse>::eval(out,*this,0,0);
                return Matrix(out,true);
            }
            else
            {
                raw::complex_dense out(ti::ti_empty(),m_rows,m_cols);
                details::dense_matrix_constructor_row<raw::complex_dense>::eval(out,0,0,*this,0);
                return Matrix(out,true);
            }
        }
        case value_code::v_float_complex:
        {
            if (is_sparse)
            {
                raw::float_complex_sparse out(ti::ti_empty(),m_rows,m_cols,m_nnz);
                details::sparse_matrix_constructor_row<raw::float_complex_sparse>::eval(out,*this,0,0);
                return Matrix(out,true);
            }
            else
            {
                raw::float_complex_dense out(ti::ti_empty(),m_rows,m_cols);
                details::dense_matrix_constructor_row<raw::float_complex_dense>::eval(out,0,0,*this,0);
                return Matrix(out,true);
            }
        };
        case value_code::v_object:
        {
            if (is_sparse)
            {
                raw::object_sparse out(m_data->get_type(),m_rows,m_cols,m_nnz);
                details::sparse_matrix_constructor_row<raw::object_sparse>::eval(out,*this,0,0);
                return Matrix(out,true);
            }
            else
            {
                raw::object_dense out(m_data->get_type(),m_rows,m_cols);
                details::dense_matrix_constructor_row<raw::object_dense>::eval(out,0,0,*this,0);
                return Matrix(out,true);
            }
        }
        default:
        {
            matcl_assert(0,"unknown value type code");
            throw;
        }
    };
};

const Matrix mat_col::to_matrix() const
{
    if (!m_data)
        return Matrix(raw::integer_sparse(ti::ti_empty(),0,0),true);

    size_t n_elem = m_data->m_vector.size();
    matcl::value_code m_value_code = this->m_data->get_value_type();

    if (m_rows == 0 || m_cols == 0 || m_nnz == 0 || n_elem == 0)
    {
        if (m_value_code == value_code::v_integer)
        {
            raw::integer_sparse out(ti::ti_empty(),m_rows,m_cols);
            return Matrix(out,true);
        }
        else if (m_value_code == value_code::v_object)
        {
            raw::object_sparse out(m_data->get_type(),m_rows,m_cols);
            return Matrix(out,true);
        }
        else if (m_value_code == value_code::v_float 
                 || m_value_code == value_code::v_float_complex)
        {
            raw::float_sparse out(ti::ti_empty(),m_rows,m_cols);
            return Matrix(out,false);
        }
        else
        {
            raw::real_sparse out(ti::ti_empty(),m_rows,m_cols);
            return Matrix(out,true);
        };
    };

    if (n_elem == 1)
    {
        const details::data_container& dc = m_data->m_vector[0];
        switch (dc.m_type)
        {
            case details::data_container::type_matrix:
            {
                Matrix out(dc.m_matrix);
                return out;
            }
            case details::data_container::type_col:
            {
                return dc.m_col.to_matrix();
            }
            case details::data_container::type_row:
            {
                return dc.m_row.to_matrix();
            }
            case details::data_container::scal_int:
            {
                return dc.m_int;
            }
            case details::data_container::scal_real:
            {
                return dc.m_real;
            }
            case details::data_container::scal_float:
            {
                return dc.m_float;
            }
            case details::data_container::scal_complex:
            {
                return Complex(dc.m_complex[0],dc.m_complex[1]);
            }
            case details::data_container::scal_fcomplex:
            {
                return Float_complex(dc.m_fcomplex[0],dc.m_fcomplex[1]);
            }
            default:
            {
                matcl_assert(0,"unknown case");
                throw;
            }
        };
    };

    Matrix out = const_cast<mat_col*>(this)->build_matrix_inplace();
    return out;
};

Matrix mat_col::build_matrix_inplace()
{
    Real density    = Real(m_nnz)/(m_rows+1.)/(m_cols+1.);
    bool is_sparse  = (density < optim_params::max_sparse_density_min);

    if (is_sparse == true)
    {
        //sparse matrices cannot be updated inplace in column oriented way
        return build_matrix_new();
    };

    matcl::value_code m_value_code = this->m_data->get_value_type();

    //only unique mat_row object can be updated inplace. Since mat_row cannot be
    //shared between threads, build_matrix() cannot be called by other threads.
    if (m_data.unique() == false)
        return build_matrix_new();

    auto& vec       = m_data->m_vector;
    size_t n        = vec.size();

    Integer mat_to_update = -1;

    for (size_t i = 0; i < n; ++i)
    {
        details::inplace_type it;
        it = vec[i].check_for_inplace_build(true, is_sparse, m_value_code, m_data->get_type());

        switch (it)
        {
            case details::inplace_type::can_inplace:
                mat_to_update = (Integer)i;
                goto exit_for_label;
            case details::inplace_type::can_continue:
                continue;
            default:
                return build_matrix_new();
        };
    };    

    if (mat_to_update == -1)
        return build_matrix_new();

  exit_for_label:

    Matrix mat = std::move(vec[mat_to_update].m_matrix);

    switch (m_value_code)
    {
        case value_code::v_integer:
        {
            using MT    = raw::integer_dense;
            MT& tmp     = mat.get_impl_unique<MT>();

            Integer row_start = tmp.rows();

            tmp.prepare_for_concat(this->m_rows, this->m_cols);            
            details::dense_matrix_constructor_col<MT>::eval(tmp,row_start,0,*this,mat_to_update+1);
            mat = Matrix(tmp,true);

            break;
        }
        case value_code::v_real:
        {
            using MT    = raw::real_dense;
            MT& tmp     = mat.get_impl_unique<MT>();

            Integer row_start = tmp.rows();

            tmp.prepare_for_concat(this->m_rows, this->m_cols);            
            details::dense_matrix_constructor_col<MT>::eval(tmp,row_start,0,*this,mat_to_update+1);
            mat = Matrix(tmp,true);

            break;
        }
        case value_code::v_float:
        {
            using MT    = raw::float_dense;
            MT& tmp     = mat.get_impl_unique<MT>();

            Integer row_start = tmp.rows();

            tmp.prepare_for_concat(this->m_rows, this->m_cols);            
            details::dense_matrix_constructor_col<MT>::eval(tmp,row_start,0,*this,mat_to_update+1);
            mat = Matrix(tmp,true);

            break;
        }
        case value_code::v_complex:
        {
            using MT    = raw::complex_dense;
            MT& tmp     = mat.get_impl_unique<MT>();

            Integer row_start = tmp.rows();

            tmp.prepare_for_concat(this->m_rows, this->m_cols);            
            details::dense_matrix_constructor_col<MT>::eval(tmp,row_start,0,*this,mat_to_update+1);
            mat = Matrix(tmp,true);

            break;
        }
        case value_code::v_float_complex:
        {
            using MT    = raw::float_complex_dense;
            MT& tmp     = mat.get_impl_unique<MT>();

            Integer row_start = tmp.rows();

            tmp.prepare_for_concat(this->m_rows, this->m_cols);            
            details::dense_matrix_constructor_col<MT>::eval(tmp,row_start,0,*this,mat_to_update+1);
            mat = Matrix(tmp,true);

            break;
        }
        case value_code::v_object:
        {
            using MT    = raw::object_dense;
            MT& tmp     = mat.get_impl_unique<MT>();

            Integer row_start = tmp.rows();

            tmp.prepare_for_concat(this->m_rows, this->m_cols);            
            details::dense_matrix_constructor_col<MT>::eval(tmp,row_start,0,*this,mat_to_update+1);
            mat = Matrix(tmp,true);

            break;
        }
        default:
        {
            matcl_assert(0,"unknown value type code");
            throw;
        }
    };

    //matrix is updated, m_data must be modified in order to forbit other inplace updates
    this->m_data = ref_ptr(new details::mat_cons_data());
    this->m_data->add_matrix(mat);

    return mat;
}

Matrix mat_col::build_matrix_new() const
{
    Real density = Real(m_nnz)/(m_rows+1.)/(m_cols+1.);
    bool is_sparse = (density < optim_params::max_sparse_density_min);

    matcl::value_code m_value_code = this->m_data->get_value_type();

    switch (m_value_code)
    {
        case value_code::v_integer:
        {
            if (is_sparse)
            {
                raw::integer_sparse out(ti::ti_empty(),m_rows,m_cols,m_nnz);
                if (m_nnz == 0)
                    return Matrix(out,true);

                details::sparse_matrix_constructor_col<raw::integer_sparse>::eval(out,*this);
                return Matrix(out,true);
            }
            else
            {
                raw::integer_dense out(ti::ti_empty(),m_rows,m_cols);
                details::dense_matrix_constructor_col<raw::integer_dense>::eval(out,0,0,*this,0);
                return Matrix(out,true);
            }
        }
        case value_code::v_float:
        {
            if (is_sparse)
            {
                raw::float_sparse out(ti::ti_empty(),m_rows,m_cols,m_nnz);
                if (m_nnz == 0)
                    return Matrix(out,true);

                details::sparse_matrix_constructor_col<raw::float_sparse>::eval(out,*this);
                return Matrix(out,true);
            }
            else
            {
                raw::float_dense out(ti::ti_empty(),m_rows,m_cols);
                details::dense_matrix_constructor_col<raw::float_dense>::eval(out,0,0,*this,0);
                return Matrix(out,true);
            }
        }
        case value_code::v_real:
        {
            if (is_sparse)
            {
                raw::real_sparse out(ti::ti_empty(),m_rows,m_cols,m_nnz);
                if (m_nnz == 0)
                    return Matrix(out,true);

                details::sparse_matrix_constructor_col<raw::real_sparse>::eval(out,*this);
                return Matrix(out,true);
            }
            else
            {
                raw::real_dense out(ti::ti_empty(),m_rows,m_cols);
                details::dense_matrix_constructor_col<raw::real_dense>::eval(out,0,0,*this,0);
                return Matrix(out,true);
            }
        }
        case value_code::v_float_complex:
        {
            if (is_sparse)
            {
                raw::float_complex_sparse out(ti::ti_empty(),m_rows,m_cols,m_nnz);
                if (m_nnz == 0)
                    return Matrix(out,true);

                details::sparse_matrix_constructor_col<raw::float_complex_sparse>::eval(out,*this);
                return Matrix(out,true);
            }
            else
            {
                raw::float_complex_dense out(ti::ti_empty(),m_rows,m_cols);
                details::dense_matrix_constructor_col<raw::float_complex_dense>::eval(out,0,0,*this,0);
                return Matrix(out,true);
            }
        }
        case value_code::v_complex:
        {
            if (is_sparse)
            {
                raw::complex_sparse out(ti::ti_empty(),m_rows,m_cols,m_nnz);
                if (m_nnz == 0)
                    return Matrix(out,true);

                details::sparse_matrix_constructor_col<raw::complex_sparse>::eval(out,*this);
                return Matrix(out,true);
            }
            else
            {
                raw::complex_dense out(ti::ti_empty(),m_rows,m_cols);
                details::dense_matrix_constructor_col<raw::complex_dense>::eval(out,0,0,*this,0);
                return Matrix(out,true);
            }
        }
        case value_code::v_object:
        {
            if (is_sparse)
            {
                raw::object_sparse out(m_data->get_type(),m_rows,m_cols,m_nnz);
                if (m_nnz == 0)
                    return Matrix(out,true);

                details::sparse_matrix_constructor_col<raw::object_sparse>::eval(out,*this);
                return Matrix(out,true);
            }
            else
            {
                raw::object_dense out(m_data->get_type(),m_rows,m_cols);
                details::dense_matrix_constructor_col<raw::object_dense>::eval(out,0,0,*this,0);
                return Matrix(out,true);
            }
        }
        default:
        {
            matcl_assert(0,"unknown value type code");
            throw;
        }
    };
};

template MATCL_MATREP_EXPORT void mat_row::add_scalar_impl<Integer>(const Integer& val);
template MATCL_MATREP_EXPORT void mat_row::add_scalar_impl<Real>(const Real& val);
template MATCL_MATREP_EXPORT void mat_row::add_scalar_impl<Float>(const Float& val);
template MATCL_MATREP_EXPORT void mat_row::add_scalar_impl<Complex>(const Complex& val);
template MATCL_MATREP_EXPORT void mat_row::add_scalar_impl<Float_complex>(const Float_complex& val);
template MATCL_MATREP_EXPORT void mat_row::add_scalar_impl<Object>(const Object& val);

template MATCL_MATREP_EXPORT void mat_col::add_scalar_impl<Integer>(const Integer& val);
template MATCL_MATREP_EXPORT void mat_col::add_scalar_impl<Real>(const Real& val);
template MATCL_MATREP_EXPORT void mat_col::add_scalar_impl<Float>(const Float& val);
template MATCL_MATREP_EXPORT void mat_col::add_scalar_impl<Complex>(const Complex& val);
template MATCL_MATREP_EXPORT void mat_col::add_scalar_impl<Float_complex>(const Float_complex& val);
template MATCL_MATREP_EXPORT void mat_col::add_scalar_impl<Object>(const Object& val);

};