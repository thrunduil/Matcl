/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-matrep/general/config.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/lib_functions/func_matrix.h"
#include "matcl-matrep/matrix/struct_flag_ext.h"
#include "matcl-internals/base/optim_params.h"

namespace matcl { namespace details
{

//----------------------------------------------------------------------
//                          matrix_info
//----------------------------------------------------------------------
class matrix_info
{
    private:
        Integer				    m_rows;
        Integer				    m_cols;
        matcl::struct_code	    m_str;
        matcl::value_code	    m_val_type;
        struct_flag             m_flag;
        Real				    m_nnz;

        matrix_info();

    public:        
        matrix_info(const Matrix& mat);

        bool                is_real_matrix() const;
        static matrix_info  make();
        static matrix_info	link_info(const matrix_info&, const matrix_info&);
        static Real			estimate_cost(const matrix_info&, const matrix_info&);
        static Real			estimate_cost_nnz(const matrix_info&, const matrix_info&);
        static Real			struct_mult_cost(matcl::struct_code, matcl::struct_code);
        static Real			estim_nnz(const matrix_info&, const matrix_info&);

        static matcl::struct_code   link_struct(matcl::struct_code s1, matcl::struct_code s2);
        static matcl::value_code    link_value(matcl::value_code v1, matcl::value_code v2);
};

matrix_info::matrix_info()
:m_rows(0),m_cols(0),m_str(struct_code::struct_scalar), m_val_type(value_code::v_integer), m_nnz(0)
{};

matrix_info matrix_info::make()
{
    return matrix_info();
};

bool matrix_info::is_real_matrix() const
{
    return m_val_type == value_code::v_integer || matrix_traits::is_float_real(m_val_type);
};
matrix_info::matrix_info(const Matrix& mat)
:m_rows(mat.rows()),m_cols(mat.cols()),m_str(mat.get_struct_code()), 
 m_val_type(mat.get_value_code()), m_nnz(mat.structural_nnz()), m_flag(mat.get_struct())
{};

matrix_info matrix_info::link_info(const matrix_info& mi1, const matrix_info& mi2)
{
    bool is_sq_1    = (mi1.m_rows == mi1.m_cols);
    bool is_sq_2    = (mi2.m_rows == mi2.m_cols);
    bool is_sq_ret  = (mi1.m_rows == mi2.m_cols);

    matrix_info out;
    out.m_rows      = mi1.m_rows;
    out.m_cols      = mi2.m_cols;
    out.m_str       = link_struct(mi1.m_str,mi2.m_str);
    out.m_val_type  = link_value(mi1.m_val_type,mi2.m_val_type);
    out.m_nnz       = estim_nnz(mi1,mi2);
    out.m_flag      = predefined_struct_ext
                        ::mult_struct(mi1.m_flag,mi2.m_flag, trans_type::no_trans, trans_type::no_trans,
                                      mi1.is_real_matrix(), mi2.is_real_matrix(), is_sq_1, is_sq_2, is_sq_ret);

    return out;
};

matcl::value_code matrix_info::link_value(matcl::value_code v1, matcl::value_code v2)
{
    return matrix_traits::unify_value_types(v1,v2);
};

matcl::struct_code matrix_info::link_struct(matcl::struct_code s1, matcl::struct_code s2)
{
    return matcl::struct_code(std::min(s1,s2));
};

Real matrix_info::estim_nnz(const matrix_info& mi1, const matrix_info& mi2)
{
    if (mi1.m_cols == 0 || mi2.m_cols == 0 || mi1.m_rows == 0 || mi2.m_rows == 0)
    {
        return 0;
    };

    if(mi1.m_flag.is_diag())
    {
        if (mi1.m_flag.is_id())
            return mi2.m_nnz;

        if (mi1.m_rows == mi1.m_cols) 
            return mi2.m_nnz;
    };

    if(mi2.m_flag.is_diag())
    {
        if (mi2.m_flag.is_id())
            return mi1.m_nnz;

        if (mi2.m_rows == mi2.m_cols) 
            return mi1.m_nnz;
    };

    Real rb = Real(mi2.m_nnz)/Real(mi2.m_cols);
    Real da = Real(mi1.m_nnz)/Real(mi1.m_rows)/Real(mi1.m_cols);
    Real db = Real(mi2.m_nnz)/Real(mi2.m_rows)/Real(mi2.m_cols);
    Real pb0 = std::pow(1. - db,Real(mi2.m_rows));

    Real nnz  = 0;
    if (pb0 != 1.)
    {
        rb      = std::min(rb/(1-pb0),Real(mi2.m_rows));
        Real ed = (1. - std::pow(1. - da,rb))*(1-pb0);

        nnz     = mi1.m_rows*mi2.m_cols*ed;
    };
    matcl_assert(nnz/mi1.m_rows/mi2.m_cols <= 1,"error in estimation of number of nonzeros");
    return nnz;
};

Real matrix_info::estimate_cost(const matrix_info& mi1, const matrix_info& mi2)
{
    Real out = estimate_cost_nnz(mi1,mi2);
    if (out == 0)
        return 0;

    value_code vc = link_value(mi1.m_val_type, mi1.m_val_type);

    switch (vc)
    {
        case value_code::v_integer:
            out = out * 2.;
            break;
        case value_code::v_float:
            out = out / 2.;
            break;
        case value_code::v_real:
            break;
        case value_code::v_float_complex:
            out = out * 2.;
            break;
        case value_code::v_complex:
            out = out * 4.;
            break;
        case value_code::v_object:
            out = out * 100.;
            break;
    };

    out				= out * struct_mult_cost(mi1.m_str, mi2.m_str);
    return out;
};

Real matrix_info::estimate_cost_nnz(const matrix_info& mi1, const matrix_info& mi2)
{
    if (mi2.m_nnz == 0 || mi1.m_nnz == 0)
        return 0;

    if (mi2.m_flag.is_id() == true || mi1.m_flag.is_id() == true)
        return 0;

    if (mi1.m_flag.is_diag() == true)
        return mi2.m_nnz;

    if (mi2.m_flag.is_diag() == true)
        return mi1.m_nnz;

    Real out = 0;
    Real rb = Real(mi2.m_nnz)/mi2.m_cols;
    Real da = Real(mi1.m_nnz)/Real(mi1.m_rows)/Real(mi1.m_cols);
    Real db = Real(mi2.m_nnz)/Real(mi2.m_rows)/Real(mi2.m_cols);
    Real pb0 = std::pow(1. - db,Real(mi2.m_rows));

    if (pb0 == 1.)
        return out;

    rb				= std::min(rb/(1-pb0),(Real)mi2.m_rows);

    out             = da*rb*mi2.m_cols;
    return out;
};

Real matrix_info::struct_mult_cost(matcl::struct_code s1, matcl::struct_code s2)
{
    switch (s1)
    {
        case struct_code::struct_dense:
            switch (s2)
            {
                case struct_code::struct_dense:
                    return optim_params::mult_cost_dense_dense;
                case struct_code::struct_sparse:
                    return optim_params::mult_cost_dense_sparse;
                case struct_code::struct_banded:
                    return optim_params::mult_cost_dense_banded;
            };
            break;
        case struct_code::struct_sparse:
            switch (s2)
            {
                case struct_code::struct_dense:
                    return optim_params::mult_cost_sparse_dense;
                case struct_code::struct_sparse:
                    return optim_params::mult_cost_sparse_sparse;
                case struct_code::struct_banded:
                    return optim_params::mult_cost_sparse_banded;
            };
            break;
        case struct_code::struct_banded:
            switch (s2)
            {
                case struct_code::struct_dense:
                    return optim_params::mult_cost_banded_dense;
                case struct_code::struct_sparse:
                    return optim_params::mult_cost_banded_sparse;
                case struct_code::struct_banded:
                    return optim_params::mult_cost_banded_banded;
            };
            break;
    };
    matcl_assert(0,"unknown case");
    throw;
};

//----------------------------------------------------------------------
//                          scalar
//----------------------------------------------------------------------
class scalar
{
    private:
        Matrix                  m_scal;
        value_code              m_value_code;
        bool                    m_set;

        scalar(const Matrix& scal, value_code vc_ret);

    public:
        scalar(value_code vc_ret);

        scalar                  mult_scalar(const Matrix& other) const;
        Matrix                  mult_matrix(const Matrix& other) const;
        Matrix                  to_matrix() const;
};

scalar::scalar(value_code vc_ret)
    :m_set(false), m_value_code(vc_ret)
{};

scalar::scalar(const Matrix& scal, value_code vc_ret)
    :m_set(true), m_scal(scal), m_value_code(vc_ret)
{
    value_code vc0  = scal.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vc0, vc_ret);

    if (vc0 != vc_ret)
        m_scal      = convert_value(m_scal, vc);
};

scalar scalar::mult_scalar(const Matrix& other) const
{
    if (m_set == false)
        return scalar(other, m_value_code);
    else
        return scalar(m_scal*other,m_value_code);
};

Matrix scalar::mult_matrix(const Matrix& other) const
{
    if (m_set == false)
        return other;
    else
        return m_scal*other;
};

Matrix scalar::to_matrix() const
{
    if (m_set == false)
        return 1;
    else
        return m_scal;
};

//----------------------------------------------------------------------
//                          matrix_chain
//----------------------------------------------------------------------
class matrix_chain
{
    private:
        std::vector<Matrix>     m_list;
        raw::integer_dense		m_indices;
        scalar                  m_scal;
        bool                    scaled;
        value_code              m_vc_return;

    public:
        matrix_chain(Integer n, value_code vc_ret);
        
        void                add(const Matrix& m);
        Matrix	            make();

        static value_code   make_value_code(value_code v1, value_code v2);
        static value_code   make_value_code(value_code v1, value_code v2, value_code v3);

    private:
        Matrix	            make(Integer i, Integer j);
        Matrix              mult(const Matrix& A, const Matrix& B, value_code vc_ret);
        void	            factor();
};

matrix_chain::matrix_chain(Integer n, value_code vc_ret)
:m_indices(ti::ti_empty()), m_vc_return(vc_ret), m_scal(vc_ret)
{
    m_list.reserve(n);
    scaled      = false;
};

value_code matrix_chain::make_value_code(value_code v1, value_code v2)
{
    value_code vc   = matrix_traits::unify_value_types(v1, v2);
    vc              = matrix_traits::real_value_type(vc);
    return vc;
}

value_code matrix_chain::make_value_code(value_code v1, value_code v2, value_code v3)
{
    value_code vc   = make_value_code(v1, v2);
    vc              = make_value_code(vc, v3);
    return vc;
}

void matrix_chain::add(const Matrix& m)
{
    if (m.is_scalar())
    {
        m_scal = m_scal.mult_scalar(m);
        return;
    }

    if (m_list.size() > 0)
    {
        error::check_mul(m_list.back().rows(),m_list.back().cols(),m.rows(),m.cols(),
                         trans_type::no_trans,trans_type::no_trans);
    };
    
    m_list.push_back(m);
};

Matrix matrix_chain::make()
{
    if (m_list.size() == 0)
        return m_scal.to_matrix();

    if (m_list.size() == 1)
        return m_scal.mult_matrix(m_list[0]);

    if (m_list.size() == 2)
    {
        Integer nz1 = m_list[0].structural_nnz();
        Integer nz2 = m_list[1].structural_nnz();

        if (nz1 < nz2)
            return mult(m_scal.mult_matrix(m_list[0]),m_list[1],m_vc_return);
        else
            return mult(m_list[0],m_scal.mult_matrix(m_list[1]),m_vc_return);
    };

    factor();
    return make(1,static_cast<Integer>(m_list.size()));
};

void matrix_chain::factor()
{
    Integer N       = static_cast<Integer>(m_list.size());

    m_indices.assign_to_fresh(raw::integer_dense(ti::ti_empty(),N,N));
    raw::real_dense m_cost = raw::real_dense(ti::ti_empty(),N,N);

    Real* ptr_cost = m_cost.ptr();
    Integer* ptr_i = m_indices.ptr();

    for (Integer i = 0, pos = 0; i < N; ++i, pos += m_cost.ld() + 1)
    {
        ptr_cost[pos] = 0;
    };

    using vec_info  = std::vector<matrix_info>;
    using mat_info  = std::vector<vec_info>;

    mat_info m_info(N);

    for (Integer i = 0; i < N; ++i)
    {
        const Matrix& A = m_list[i];
        matrix_info inf_A(A);

        vec_info v_info(N,matrix_info::make());
        v_info[i] = inf_A;

        for (Integer k = i+1; k < N; ++k)
        {
            const Matrix& B = m_list[k];
            matrix_info inf_B(B);

            inf_A = matrix_info::link_info(inf_A,inf_B);
            v_info[k] = inf_A;
        };

        m_info[i] = v_info;
    };

    Integer ldc = m_cost.ld();
    Integer ldi = m_indices.ld();

    for (Integer l = 1; l < N; ++l)
    {
        for (Integer i = 0; i < N-l; ++i)
        {            
            Integer j = i + l;
            ptr_cost[i+j*ldc] = constants::inf();

            for (Integer k = i; k < j; ++k)
            {
                Real cost = matrix_info::estimate_cost(m_info[i][k],m_info[k+1][j]);
                Real q = ptr_cost[i+k*ldc] + ptr_cost[k+1+j*ldc] + cost;

                if (q < ptr_cost[i+j*ldc])
                {
                    ptr_cost[i+j*ldc] = q;
                    ptr_i[i+j*ldi] = k+1;
                };
            };
        };
    };
    return;
};

Matrix matrix_chain::make(Integer i, Integer j)
{
    if (i==j)
    {
        if (scaled)
            return Matrix(m_list[i-1]);
        else
        {
            scaled = true;
            return m_scal.mult_matrix(m_list[i-1]);
        };
    };

    Integer k   = m_indices.ptr()[i-1+(j-1)*m_indices.ld()];
    Matrix X    = make(i,k);
    Matrix Y    = make(k+1,j);

    return mult(X,Y,m_vc_return);
};

Matrix matrix_chain::mult(const Matrix& A, const Matrix& B, value_code vc_ret)
{
    if (vc_ret != value_code::v_integer && (A.get_value_code() == value_code::v_integer) 
        && (B.get_value_code() == value_code::v_integer))
    {
        if (A.structural_nnz() < B.structural_nnz())
        {
            matcl::mat_code mt  = matrix_traits::get_matrix_type(value_code::v_real,A.get_struct_code());
            return convert(A,mt)*B;
        }
        else
        {
            matcl::mat_code mt  = matrix_traits::get_matrix_type(value_code::v_real,B.get_struct_code());
            return A*convert(B,mt);
        };
    }
    else
    {
        return A*B;
    };
};

};

//----------------------------------------------------------------------
//                          chain_mult
//----------------------------------------------------------------------
Matrix matcl::chain_mult(const Matrix& A1, const Matrix& A2, const Matrix& A3)
{
    value_code vc   = details::matrix_chain::make_value_code(A1.get_value_code(), 
                        A2.get_value_code(), A3.get_value_code());

    details::matrix_chain mc(3, vc);
    mc.add(A1); mc.add(A2); mc.add(A3);

    return mc.make();
};

Matrix matcl::chain_mult(const Matrix& A1, const Matrix& A2, const Matrix& A3, const Matrix& A4)
{
    value_code vc   = details::matrix_chain::make_value_code(A1.get_value_code(), 
                        A2.get_value_code(), A3.get_value_code());
    vc              = details::matrix_chain::make_value_code(vc, A4.get_value_code()); 

    details::matrix_chain mc(4, vc);
    mc.add(A1); mc.add(A2); mc.add(A3); mc.add(A4);

    return mc.make();
};

Matrix matcl::chain_mult(const Matrix& A1, const Matrix& A2, const Matrix& A3, const Matrix& A4, 
                        const Matrix& A5)
{
    value_code vc   = details::matrix_chain::make_value_code(A1.get_value_code(), 
                        A2.get_value_code(), A3.get_value_code());
    vc              = details::matrix_chain::make_value_code(vc, A4.get_value_code(),
                        A5.get_value_code()); 

    details::matrix_chain mc(5, vc);
    mc.add(A1); mc.add(A2); mc.add(A3); mc.add(A4); mc.add(A5);

    return mc.make();
};

Matrix matcl::chain_mult(const Matrix& A1, const Matrix& A2, const Matrix& A3, const Matrix& A4, 
                        const Matrix& A5, const Matrix& A6)
{
    value_code vc   = details::matrix_chain::make_value_code(A1.get_value_code(), 
                        A2.get_value_code(), A3.get_value_code());
    vc              = details::matrix_chain::make_value_code(vc, A4.get_value_code(),
                        A5.get_value_code()); 
    vc              = details::matrix_chain::make_value_code(vc, A6.get_value_code()); 

    details::matrix_chain mc(6, vc);
    mc.add(A1); mc.add(A2); mc.add(A3); mc.add(A4); mc.add(A5); mc.add(A6);

    return mc.make();
};

};
