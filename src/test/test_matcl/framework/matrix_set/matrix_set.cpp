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

#include "matrix_set.h"
#include "matcl-core/IO/logger.h"
#include "test/test_matcl/framework/matrix_utils.h"
#include "test_options.h"

#pragma warning( push )
#pragma warning(disable:4702)	//  unreachable code
    #include <boost/lexical_cast.hpp>
#pragma warning( pop )

#ifndef _MSC_VER

#include <cstdarg>
namespace
{
    int sprintf_s(char * _DstBuf, size_t _SizeInBytes, const char * _Format, ...)
    {
        va_list args;
        va_start (args, _Format);
        int n = vsnprintf(_DstBuf, _SizeInBytes, _Format, args);
        va_end (args);

        return n;
    }
}
#endif

namespace matcl { namespace test
{

static const Integer max_int = 1000;

static std::string val_type_to_str(matcl::value_code val_type)
{
    switch (val_type)
    {
        case value_code::v_integer:
        {
            return "integer";
        }
        case value_code::v_float:
        {
            return "float";
        }
        case value_code::v_real:
        {
            return "real";
        }
        case value_code::v_float_complex:
        {
            return "float_complex";
        }
        case value_code::v_complex:
        {
            return "complex";
        }
    };
    
    return "";
};

static std::string struct_type_to_str(matcl::struct_code str_type)
{
    switch (str_type)
    {
        case struct_code::struct_dense:
        {
            return "dense";
            break;
        }
        case struct_code::struct_banded:
        {
            return "banded";
            break;
        }
        case struct_code::struct_sparse:
        {
            return "sparse";
        }
    };
    return "";
};

//=================================================================================
//                      MATRIX SET
//=================================================================================
Real matrix_set::make(const Matrix& mat,unary_function* func, bool show_partial_res, int code) const
{
    if (show_partial_res)
    {
        std::string mat_name = make_mat_name(mat,code);
        matcl::out_stream << mat_name + "\n";
    };
    Real val = 0.;
    try
    {
        val = func->eval(mat,show_partial_res, code);
        if (show_partial_res)
        {			
            if (func->is_error())
                matcl::out_stream  << std::string("    ") + func->get_error() + "\n";
            else
                matcl::out_stream  << std::string("    ") + boost::lexical_cast<std::string>(val) + "\n";
        };
    }
    catch(const std::exception& ex)
    {
        val = 1.; // exception here means test fail
        if (show_partial_res)
            matcl::out_stream  << std::string("    ") + ex.what() + "\n";
    };

    return val;
};

Real matrix_set::make(const Scalar& mat,unary_function* func, bool show_partial_res, int code) const
{
    if (show_partial_res)
    {
        std::string mat_name = make_scal_name(mat,code);
        matcl::out_stream << mat_name + "\n";
    };

    Real val = 0.;
    
    try
    {
        val = func->eval_scal(mat,show_partial_res,code);
        if (show_partial_res)
        {			
            if (func->is_error())
                matcl::out_stream  << std::string("    ") + func->get_error() + "\n";
            else
                matcl::out_stream  << std::string("    ") + boost::lexical_cast<std::string>(val) + "\n";
        };
    }
    catch(const std::exception& ex)
    {
        val = 1.; // exception here means test fail
        if (show_partial_res)
            matcl::out_stream  << std::string("    ") + ex.what() + "\n";
    };

    return val;
};

std::string matrix_set::make_mat_name(const Matrix& mat, int code) const
{
    Integer m   = mat.rows();
    Integer n   = mat.cols();

    char buff[100];

    matcl::value_code val_type	= mat.get_value_code();
    matcl::struct_code str_type	= mat.get_struct_code();

    std::string val_str;
    std::string str_str;

    val_str = val_type_to_str(val_type);
    str_str = struct_type_to_str(str_type);

    sprintf_s(buff,100,"%dx%d %s %s matrix, code: %d",m,n,val_str.c_str(),str_str.c_str(),code);

    std::string out(buff);
    return out;
};

std::string matrix_set::make_scal_name(const Scalar& mat, int code) const
{
    char buff[100];

    matcl::value_code val_type	= mat.get_value_code();

    std::string val_str = val_type_to_str(val_type);

    sprintf_s(buff,100,"%s scalar, code: %d",val_str.c_str(),code);

    std::string out(buff);
    return out;
};

Real matrix_set::make(unary_function* func, options opts, const std::vector<Matrix>& matrices,
                      const std::vector<Scalar>& scalars) const
{
    Real out                = 0.;
    int size                = (int)matrices.size();
    bool show_partial_res   = opts.show_partial_res;
    bool show_memleaks      = opts.show_memleaks;

    for (int code = opts.first_matrix_code; code < size; code++)
    {
        long N      = matcl::details::no_existing_objects();
        Real tmp    = matrix_set::make(matrices[code],func,show_partial_res,code);
        N           = matcl::details::no_existing_objects() - N - func->n_new_objects();

        if (show_partial_res == false && tmp != 0.)
            matrix_set::make(matrices[code],func,true,code);

        if (show_memleaks && N != 0)
            matcl::out_stream  << "memleaks\n";

        out += tmp;
    };

    size = (int)scalars.size();

    for (int code = opts.first_scalar_code; code < size; code++)
    {
        Real tmp = matrix_set::make(scalars[code],func,show_partial_res,code);

        if (show_partial_res == false && tmp != 0.)
            matrix_set::make(scalars[code],func,true,code);

        out += tmp;
    };

    if (show_partial_res)
        matcl::out_stream  << std::string("total error: ") + boost::lexical_cast<std::string>(out) + "\n";

    return out/(matrices.size()+scalars.size());
};

void matrix_set::add_matrices(std::vector<Matrix>& mat, Integer m, Integer n)
{
    add_matrices_dense(mat,m,n);

    for (size_t i = 0; i < m_sparse_info.size(); ++i)
        add_matrices_sparse(mat,m,n,m_sparse_info[i]);

    for (size_t i = 0; i < m_band_info.size(); ++i)
        add_matrices_band(mat,m,n,m_band_info[i].first,m_band_info[i].second);
};

void matrix_set::add_matrices_dense(std::vector<Matrix>& mc, Integer m, Integer n)
{
    if (m == 1 && n == 1)
    {
        mc.push_back(m_rand->rand_scalar_int());
        mc.push_back(m_rand->rand_scalar_real());
        mc.push_back(m_rand->rand_scalar_float());
        mc.push_back(m_rand->rand_scalar_compl());
        mc.push_back(m_rand->rand_scalar_fcompl());
    };

    //int
    {
        Matrix tmp = m_rand->rand_dense_int(m,n);
        Matrix m1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_dense_int(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //real
    {
        Matrix tmp = m_rand->rand_dense_real(m,n);
        Matrix m1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_dense_real(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //complex
    {
        Matrix tmp = m_rand->rand_dense_compl(m,n);
        Matrix m1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_dense_compl(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //float
    {
        Matrix tmp = m_rand->rand_dense_float(m,n);
        Matrix m1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_dense_float(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //float complex
    {
        Matrix tmp = m_rand->rand_dense_fcompl(m,n);
        Matrix m1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_dense_fcompl(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };
};

void matrix_set::add_matrices_sparse(std::vector<Matrix>& mc, Integer m, Integer n, Real d)
{
    //int
    {
        Matrix tmp = m_rand->rand_sparse_int(m,n,d);
        Matrix m1 = tmp;

        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_sparse_int(m,n+5,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix m3 = tmp;        

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //real
    {
        Matrix tmp = m_rand->rand_sparse_real(m,n,d);
        Matrix m1 = tmp;

        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_sparse_real(m,n+5,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //complex
    {
        Matrix tmp = m_rand->rand_sparse_compl(m,n,d);
        Matrix m1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_sparse_compl(m,n+5,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //float
    {
        Matrix tmp = m_rand->rand_sparse_float(m,n,d);
        Matrix m1 = tmp;

        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_sparse_float(m,n+5,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //float complex
    {
        Matrix tmp = m_rand->rand_sparse_fcompl(m,n,d);
        Matrix m1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_sparse_fcompl(m,n+5,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

};

void matrix_set::add_matrices_band(std::vector<Matrix>& mc, Integer m, Integer n, Integer ld, Integer ud)
{
    Integer diag_l = (m == 0)? 0 : 1;
    Integer diag_u = (n == 0)? 0 : 1;

    if (ld + diag_l > m)
        return;

    if (ud + diag_u > n)
        return;

    //int
    {
        Matrix tmp = m_rand->rand_band_int(m,n,-ld,ud);
        Matrix m1 = tmp;

        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_band_int(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        tmp = m_rand->rand_band_int(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix m4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
            case 3: mc.push_back(m4); break;
        };
    };

    //real
    {
        Matrix tmp = m_rand->rand_band_real(m,n,-ld,ud);
        Matrix m1 = tmp;

        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_band_real(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        tmp = m_rand->rand_band_real(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix m4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
            case 3: mc.push_back(m4); break;
        };
    };

    //complex
    {
        Matrix tmp = m_rand->rand_band_compl(m,n,-ld,ud);
        Matrix m1 = tmp;

        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_band_compl(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        tmp = m_rand->rand_band_compl(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix m4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
            case 3: mc.push_back(m4); break;
        };
    };

    //float
    {
        Matrix tmp = m_rand->rand_band_float(m,n,-ld,ud);
        Matrix m1 = tmp;

        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_band_float(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        tmp = m_rand->rand_band_float(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix m4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
            case 3: mc.push_back(m4); break;
        };
    };

    //float complex
    {
        Matrix tmp = m_rand->rand_band_fcompl(m,n,-ld,ud);
        Matrix m1 = tmp;

        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_band_fcompl(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        tmp = m_rand->rand_band_fcompl(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix m4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
            case 3: mc.push_back(m4); break;
        };
    };
};

Real matrix_set::make(unary_function* func, options opts) const
{
    return matrix_set::make(func,opts,m_matrices, m_scalars);
};

Integer matrix_set::get_matrix_vector_size() const
{
    return (Integer)m_matrices.size();
};
Matrix matrix_set::get_matrix(int code) const
{
    int s = (int)m_matrices.size();

    if (code >=  s)
        return m_matrices.back();

    if (code < 0)
        return m_matrices.front();

    return m_matrices[code];
};

//=================================================================================
//                      MATRIX SET BIN
//=================================================================================
Real matrix_set_bin::make(const matrix_pair& mat,bin_function* func, bool show_partial_res, int code) const
{
    if (show_partial_res)
    {
        std::string mat_name = make_mat_name(mat,code);
        matcl::out_stream << mat_name + "\n";
    };
    Real val = 0.;
  
    try
    {
        val = func->eval(mat.first,mat.second,code);
        if (show_partial_res)
        {			
            if (func->is_error())
                matcl::out_stream  << std::string("    ") + func->get_error() + "\n";
            else
                matcl::out_stream  << std::string("    ") + boost::lexical_cast<std::string>(val) + "\n";
        };
    }
    catch(const std::exception& ex)
    {
        val = 1.; // exception here means test fail
        if (show_partial_res)       
            matcl::out_stream  << std::string("    ") + ex.what() + "\n";
    };

    return val;
};

Real matrix_set_bin::make(const scalar_pair& mat,bin_function* func, bool show_partial_res, int code) const
{
    if (show_partial_res)
    {
        std::string mat_name = make_scal_name(mat,code);
        matcl::out_stream << mat_name + "\n";
    };

    Real val = 0.;

    try
    {
        val = func->eval_scal(mat.first,mat.second,code);
        if (show_partial_res)
        {			
            if (func->is_error())
                matcl::out_stream  << std::string("    ") + func->get_error() + "\n";
            else
                matcl::out_stream  << std::string("    ") + boost::lexical_cast<std::string>(val) + "\n";
        };
    }
    catch(const std::exception& ex)
    {
        val = 1.; // exception here means test fail
        if (show_partial_res)
            matcl::out_stream  << std::string("    ") + ex.what() + "\n";
    };

    return val;
};

std::string matrix_set_bin::make_mat_name(const matrix_pair& mat, int code) const
{
    Integer m_1 = mat.first.rows();
    Integer n_1 = mat.first.cols();
    Integer m_2 = mat.second.rows();
    Integer n_2 = mat.second.cols();

    char buff[100];

    matcl::value_code val_type_1	= mat.first.get_value_code();
    matcl::struct_code str_type_1	= mat.first.get_struct_code();

    matcl::value_code val_type_2	= mat.second.get_value_code();
    matcl::struct_code str_type_2	= mat.second.get_struct_code();

    std::string val_str_1, val_str_2;
    std::string str_str_1, str_str_2;

    val_str_1 = val_type_to_str(val_type_1);
    val_str_2 = val_type_to_str(val_type_2);
    str_str_1 = struct_type_to_str(str_type_1);
    str_str_2 = struct_type_to_str(str_type_2);

    sprintf_s(buff,100,"%dx%d %s %s matrix and %dx%d %s %s matrix, code: %d",
                m_1,n_1,val_str_1.c_str(),str_str_1.c_str(),
                m_2,n_2,val_str_2.c_str(),str_str_2.c_str(),
                code);

    std::string out(buff);
    return out;
};

std::string matrix_set_bin::make_scal_name(const scalar_pair& mat, int code) const
{
    char buff[100];

    matcl::value_code val_type_1	= mat.first.get_value_code();
    matcl::value_code val_type_2	= mat.second.get_value_code();

    std::string val_str_1, val_str_2;

    val_str_1 = val_type_to_str(val_type_1);
    val_str_2 = val_type_to_str(val_type_2);

    sprintf_s(buff,100,"%s scalar and %s scalar, code: %d",
                val_str_1.c_str(),val_str_2.c_str(),code);

    std::string out(buff);
    return out;
};

Real matrix_set_bin::make(bin_function* func, options opts, const std::vector<matrix_pair>& matrices,
                          const std::vector<scalar_pair>& scalars) const
{
    Real out = 0.;
    bool show_partial_res = opts.show_partial_res;

    int size = (int)matrices.size();	

    for (int code = opts.first_matrix_code; code < size; code++)
    {
        Real tmp = matrix_set_bin::make(matrices[code],func,show_partial_res,code);

        if (show_partial_res == false && tmp != 0.)
            matrix_set_bin::make(matrices[code],func,true,code);

        out += tmp;
    };

    size = (int)scalars.size();	

    for (int code = opts.first_scalar_code; code < size; code++)
    {
        Real tmp = matrix_set_bin::make(scalars[code],func,show_partial_res,code);
    
        if (show_partial_res == false && tmp != 0.)
            matrix_set_bin::make(scalars[code],func,true,code);

        out += tmp;
    };

    if (show_partial_res)
        matcl::out_stream  << std::string("total error: ") + boost::lexical_cast<std::string>(out) + "\n";

    return out/(matrices.size()+scalars.size());
};

Real matrix_set_bin::make(bin_function* func, options opts) const
{
    return matrix_set_bin::make(func,opts,m_matrices,m_scalars);
};

matrix_set_bin::matrix_pair matrix_set_bin::get_matrix(int code) const
{
    if (code > (int)m_matrices.size() )
        return m_matrices.back();

    if (code < 0)
        return m_matrices.front();

    return m_matrices[code];
};
matrix_set_bin::scalar_pair matrix_set_bin::get_scalar(int code) const
{
    if (code > (int)m_scalars.size() || code < 0)
        throw std::runtime_error("invalid scalar code");

    return m_scalars[code];
};

void matrix_set_bin::add_matrices(std::vector<matrix_pair>& mat, Integer m1, Integer n1, 
                                  Integer m2, Integer n2)
{
    add_matrices_dense(mat, m1,n1,m2,n2);

    for (size_t i = 0; i < m_sparse_info.size(); ++i)
        add_matrices_sparse(mat,m1,n1,m2,n2,m_sparse_info[i]);

    for (size_t i = 0; i < m_band_info.size(); ++i)
        add_matrices_band(mat,m1,n1,m2,n2,m_band_info[i].first,m_band_info[i].second);
};

void matrix_set_bin::add_matrices(std::vector<matrix_pair>& mat, Integer m, Integer n, Matrix nm)
{
    add_matrices_dense(mat,nm,m,n);

    for (size_t i = 0; i < m_sparse_info.size(); ++i)
        add_matrices_sparse(mat,nm,m,n,m_sparse_info[i]);

    for (size_t i = 0; i < m_band_info.size(); ++i)
        add_matrices_band(mat,nm,m,n,m_band_info[i].first,m_band_info[i].second);
};

void matrix_set_bin::add_matrices_dense(std::vector<matrix_pair>& mat, Integer m1, Integer n1, 
                                        Integer m2, Integer n2)
{
    if (m1 == 1 && n1 == 1)
    {
        add_matrices(mat, m2, n2, m_rand->rand_scalar_int());
        add_matrices(mat, m2, n2, m_rand->rand_scalar_real());
        add_matrices(mat, m2, n2, m_rand->rand_scalar_float());
        add_matrices(mat, m2, n2, m_rand->rand_scalar_compl());
        add_matrices(mat, m2, n2, m_rand->rand_scalar_fcompl());
    };

    //int
    {
        Matrix tmp = m_rand->rand_dense_int(m1,n1);
        Matrix t1 = tmp;
        
        tmp.reserve(m1+2,n1+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_dense_int(m1+5,n1+5);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
        };
    };

    //real
    {
        Matrix tmp = m_rand->rand_dense_real(m1,n1);
        Matrix t1 = tmp;
        tmp.reserve(m1+2,n1+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_dense_real(m1+5,n1+5);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
        };
    };

    //complex
    {
        Matrix tmp = m_rand->rand_dense_compl(m1,n1);
        Matrix t1 = tmp;
        tmp.reserve(m1+2,n1+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_dense_compl(m1+5,n1+5);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
        };
    };

    //float
    {
        Matrix tmp = m_rand->rand_dense_float(m1,n1);
        Matrix t1 = tmp;
        tmp.reserve(m1+2,n1+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_dense_float(m1+5,n1+5);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
        };
    };

    //float complex
    {
        Matrix tmp = m_rand->rand_dense_fcompl(m1,n1);
        Matrix t1 = tmp;
        tmp.reserve(m1+2,n1+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_dense_fcompl(m1+5,n1+5);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
        };
    };
}

void matrix_set_bin::add_matrices_sparse(std::vector<matrix_pair>& mat, Integer m1, Integer n1, 
                                         Integer m2, Integer n2, Real d)
{
    //int
    {
        Matrix tmp = m_rand->rand_sparse_int(m1,n1,d);
        Matrix t1 = tmp;
        tmp.reserve(m1+2,n1+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_sparse_int(m1,n1+5,d);
        tmp = tmp(colon(),colon(2,2+n1-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
        };
    };

    //real
    {
        Matrix tmp = m_rand->rand_sparse_real(m1,n1,d);
        Matrix t1 = tmp;
        tmp.reserve(m1+2,n1+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_sparse_real(m1,n1+5,d);
        tmp = tmp(colon(),colon(2,2+n1-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
        };
    };

    //complex
    {
        Matrix tmp = m_rand->rand_sparse_compl(m1,n1,d);
        Matrix t1 = tmp;
        tmp.reserve(m1+2,n1+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_sparse_compl(m1,n1+5,d);
        tmp = tmp(colon(),colon(2,2+n1-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
        };
    };

    //float
    {
        Matrix tmp = m_rand->rand_sparse_float(m1,n1,d);
        Matrix t1 = tmp;
        tmp.reserve(m1+2,n1+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_sparse_float(m1,n1+5,d);
        tmp = tmp(colon(),colon(2,2+n1-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
        };
    };

    //float complex
    {
        Matrix tmp = m_rand->rand_sparse_fcompl(m1,n1,d);
        Matrix t1 = tmp;
        tmp.reserve(m1+2,n1+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_sparse_fcompl(m1,n1+5,d);
        tmp = tmp(colon(),colon(2,2+n1-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
        };
    };
}

void matrix_set_bin::add_matrices_band(std::vector<matrix_pair>& mat, Integer m1, Integer n1, 
                                       Integer m2, Integer n2, Integer ld, Integer ud)
{
    Integer diag_l = (m1 == 0)? 0 : 1;
    Integer diag_u = (n1 == 0)? 0 : 1;

    if (ld + diag_l > m1)
        return;

    if (ud + diag_u > n1)
        return;

    //int
    {
        Matrix tmp = m_rand->rand_band_int(m1,n1,-ld,ud);
        Matrix t1 = tmp;

        tmp.reserve_band(m1+2,n1+2,-(ld+2),ud+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_band_int(m1+5,n1+5,-ld,ud);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        Matrix t3 = tmp;

        tmp = m_rand->rand_band_int(m1+5,n1+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix t4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
            case 3: add_matrices(mat, m2, n2, t4); break;
        };
    };

    //real
    {
        Matrix tmp = m_rand->rand_band_real(m1,n1,-ld,ud);
        Matrix t1 = tmp;
        tmp.reserve_band(m1+2,n1+2,-(ld+2),ud+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_band_real(m1+5,n1+5,-ld,ud);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        Matrix t3 = tmp;

        tmp = m_rand->rand_band_real(m1+5,n1+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix t4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
            case 3: add_matrices(mat, m2, n2, t4); break;
        };
    };

    //complex
    {
        Matrix tmp = m_rand->rand_band_compl(m1,n1,-ld,ud);
        Matrix t1 = tmp;
        tmp.reserve_band(m1+2,n1+2,-(ld+2),ud+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_band_compl(m1+5,n1+5,-ld,ud);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        Matrix t3 = tmp;

        tmp = m_rand->rand_band_compl(m1+5,n1+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix t4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
            case 3: add_matrices(mat, m2, n2, t4); break;
        };
    };

    //float
    {
        Matrix tmp = m_rand->rand_band_float(m1,n1,-ld,ud);
        Matrix t1 = tmp;
        tmp.reserve_band(m1+2,n1+2,-(ld+2),ud+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_band_float(m1+5,n1+5,-ld,ud);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        Matrix t3 = tmp;

        tmp = m_rand->rand_band_float(m1+5,n1+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix t4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
            case 3: add_matrices(mat, m2, n2, t4); break;
        };
    };

    //float complex
    {
        Matrix tmp = m_rand->rand_band_fcompl(m1,n1,-ld,ud);
        Matrix t1 = tmp;
        tmp.reserve_band(m1+2,n1+2,-(ld+2),ud+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_band_fcompl(m1+5,n1+5,-ld,ud);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        Matrix t3 = tmp;

        tmp = m_rand->rand_band_fcompl(m1+5,n1+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m1-1),colon(2,2+n1-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix t4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: add_matrices(mat, m2, n2, t1); break;
            case 1: add_matrices(mat, m2, n2, t2); break;
            case 2: add_matrices(mat, m2, n2, t3); break;
            case 3: add_matrices(mat, m2, n2, t4); break;
        };
    };
};

void matrix_set_bin::add_matrices_dense(std::vector<matrix_pair>& mat, Matrix nm, Integer m, Integer n)
{
    if (m == 1 && n == 1)
    {
        mat.push_back(matrix_pair(nm,m_rand->rand_scalar_int()));
        mat.push_back(matrix_pair(nm,m_rand->rand_scalar_real()));
        mat.push_back(matrix_pair(nm,m_rand->rand_scalar_float()));
        mat.push_back(matrix_pair(nm,m_rand->rand_scalar_compl()));
        mat.push_back(matrix_pair(nm,m_rand->rand_scalar_fcompl()));
    };

    //int
    {
        Matrix tmp = m_rand->rand_dense_int(m,n);
        Matrix t1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_dense_int(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
        };
    };

    //real
    {
        Matrix tmp = m_rand->rand_dense_real(m,n);
        Matrix t1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_dense_real(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
        };
    };

    //complex
    {
        Matrix tmp = m_rand->rand_dense_compl(m,n);
        Matrix t1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_dense_compl(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
        };
    };

    //float
    {
        Matrix tmp = m_rand->rand_dense_float(m,n);
        Matrix t1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_dense_float(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
        };
    };

    //float complex
    {
        Matrix tmp = m_rand->rand_dense_fcompl(m,n);
        Matrix t1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_dense_fcompl(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
        };
    };
}

void matrix_set_bin::add_matrices_sparse(std::vector<matrix_pair>& mat, Matrix nm, Integer m, Integer n, 
                                         Real d)
{
    //int
    {
        Matrix tmp = m_rand->rand_sparse_int(m,n,d);
        Matrix t1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_sparse_int(m,n+5,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
        };
    };

    //real
    {
        Matrix tmp = m_rand->rand_sparse_real(m,n,d);
        Matrix t1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_sparse_real(m,n+5,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
        };
    };

    //complex
    {
        Matrix tmp = m_rand->rand_sparse_compl(m,n,d);
        Matrix t1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_sparse_compl(m,n+5,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
        };
    };

    //float
    {
        Matrix tmp = m_rand->rand_sparse_float(m,n,d);
        Matrix t1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_sparse_float(m,n+5,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
        };
    };

    //float complex
    {
        Matrix tmp = m_rand->rand_sparse_fcompl(m,n,d);
        Matrix t1 = tmp;
        tmp.reserve(m+2,n+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_sparse_fcompl(m,n+5,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix t3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
        };
    };
};

void matrix_set_bin::add_matrices_band(std::vector<matrix_pair>& mat, Matrix nm, Integer m, Integer n, 
                                      Integer ld, Integer ud)
{
    Integer diag_l = (m == 0)? 0 : 1;
    Integer diag_u = (n == 0)? 0 : 1;

    if (ld + diag_l > m)
        return;

    if (ud + diag_u > n)
        return;

    //int
    {
        Matrix tmp = m_rand->rand_band_int(m,n,-ld,ud);
        Matrix t1 = tmp;
        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_band_int(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix t3 = tmp;

        tmp = m_rand->rand_band_int(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix t4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
            case 3: add(mat,matrix_pair(nm,t4)); break;
        };
    };

    //real
    {
        Matrix tmp = m_rand->rand_band_real(m,n,-ld,ud);
        Matrix t1 = tmp;
        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_band_real(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix t3 = tmp;

        tmp = m_rand->rand_band_real(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix t4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
            case 3: add(mat,matrix_pair(nm,t4)); break;
        };
    };

    //complex
    {
        Matrix tmp = m_rand->rand_band_compl(m,n,-ld,ud);
        Matrix t1 = tmp;
        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_band_compl(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix t3 = tmp;

        tmp = m_rand->rand_band_compl(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix t4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
            case 3: add(mat,matrix_pair(nm,t4)); break;
        };
    };

    //float
    {
        Matrix tmp = m_rand->rand_band_float(m,n,-ld,ud);
        Matrix t1 = tmp;
        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_band_float(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix t3 = tmp;

        tmp = m_rand->rand_band_float(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix t4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
            case 3: add(mat,matrix_pair(nm,t4)); break;
        };
    };

    //float complex
    {
        Matrix tmp = m_rand->rand_band_fcompl(m,n,-ld,ud);
        Matrix t1 = tmp;
        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix t2 = tmp;

        tmp = m_rand->rand_band_fcompl(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix t3 = tmp;

        tmp = m_rand->rand_band_fcompl(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix t4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: add(mat,matrix_pair(nm,t1)); break;
            case 1: add(mat,matrix_pair(nm,t2)); break;
            case 2: add(mat,matrix_pair(nm,t3)); break;
            case 3: add(mat,matrix_pair(nm,t4)); break;
        };
    };
};

void matrix_set_bin::add(std::vector<matrix_pair>& mat,const matrix_pair& mp)
{
    Real prob       = test_options::get_binmat_prob();
    bool instert    = rand() < prob;

    if (instert)
        mat.push_back(mp);
};

};};

