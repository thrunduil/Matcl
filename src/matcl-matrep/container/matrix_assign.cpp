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

#include "matcl-matrep/container/matrix2.inl"
#include "matcl-matrep/visitors/assign_visitor.inl" 
#include "matcl-matrep/visitors/assign_visitor_get.inl" 
#include "matcl-matrep/visitors/assign_visitor_mat.inl" 
#include "matcl-matrep/visitors/assign_visitor_del.inl" 
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-core/error/exception_classes.h"
#include "matcl-internals/error/error_check_basic.h"
#include "matcl-matrep/container/matrix_container.inl"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-matrep/details/extract_type2_switch.h" 
#include "matcl-matrep/details/struct_flag_predefined.h"

namespace matcl
{

namespace md = matcl::details;

namespace details
{

static bool is_tr_view(const colon& c1,const colon& c2)
{
    bool r1, r2;

    Integer m_s1    = 0;
    Integer m_s2    = 0;

    switch(c1.m_flag)
    {
        case colon::t_all:
        {
            r1  = false;
            break;
        }

        case colon::t_range_simple:
        case colon::t_range_mat:
        case colon::t_end:
        {
            r1      = true;
            m_s1    = c1.m_s;

            if (c1.m_i != 1)
                return false;

            break;        
        }
        default:
            return false;
    };

    switch(c2.m_flag)
    {
        case colon::t_all:
        {
            r2 = false; 
            break;
        }
        case colon::t_range_simple:
        case colon::t_range_mat:
        case colon::t_end:
        {
            r2      = true;
            m_s2    = c2.m_s;

            if (c2.m_i != 1)
                return false;

            break;
        }
        default:
        {
            return false;
        };
    };

    if (r1 == true)
    {
        if (r2 == true)
        {
            if (m_s1 == m_s2)   return true;
            else                return false;
        }
        else
        {
            if (m_s1 == 1)      return true;
            else                return false;
        }
    }
    else
    {
        if (r2 == true)
        {
            if (m_s2 == 1)      return true;
            else                return false;
        }
        else
        {
            return true;
        }
    };
};

struct assign_visitor_1: public extract_type2_switch_nc<Matrix&,assign_visitor_1>
{
    Integer ind;

    assign_visitor_1(Integer ind) : ind(ind){};

    template<class M1, class M2>
    Matrix& eval_mat_mat(Matrix& handle, M1& A, const M2& B)
    {
        Integer r = B.rows();
        Integer c = B.cols();

        if (r == 0 && c == 0)
        {
            handle = delrows(vec(handle),ind);
            return handle;
        };

        error::check_assign_1(1,r,c);
        
        using value_type_2 = typename M2::value_type;

        return assign_functor<M1,value_type_2>::eval(handle,A,B(1,1),ind);
    };

    template<class M1, class M2>
    Matrix& eval_mat_scal(Matrix& handle, M1& A, const M2& B)
    {
        using value_type_2 = M2;
        return assign_functor<M1,value_type_2>::eval(handle,A,B,ind);
    };

    template<class M1, class M2>
    Matrix& eval_scal_mat(Matrix& handle, M1&, const M2& B)
    {
        Integer r = B.rows();
        Integer c = B.cols();

        if (r == 0 && c == 0)
        {
            handle = delrows(vec(handle),ind);
            return handle;
        };

        using FullMatrix = raw::Matrix<M1,struct_dense>;
        error::check_assign_1(1,r,c);

        using value_type_2 = typename M2::value_type;
        return assign_functor_scal<M1,value_type_2>::eval(handle,B(1,1),ind);
    };

    template<class M1, class M2>
    Matrix& eval_scal_scal(Matrix& handle, M1&, const M2& B)
    {
        using value_type_2 = M2;
        return assign_functor_scal<M1,value_type_2>::eval(handle,B,ind);
    };
};

struct assign_visitor_2: public extract_type2_switch_nc<Matrix&,assign_visitor_2>
{
    Integer ind_1, ind_2;

    assign_visitor_2(Integer ind_1, Integer ind_2) : ind_1(ind_1), ind_2(ind_2){};

    template<class M1, class M2>
    Matrix& eval_mat_mat(Matrix& handle, M1& A, const M2& B)
    {
        error::check_assign_2(1,1,B.rows(),B.cols());
        using value_type_2 = typename M2::value_type;

        return assign_functor<M1,value_type_2>::eval(handle,A,B(1,1),ind_1,ind_2);
    };

    template<class M1, class M2>
    Matrix& eval_mat_scal(Matrix& handle, M1& A, const M2& B)
    {
        using value_type_2 = M2;
        return assign_functor<M1,value_type_2>::eval(handle,A,B,ind_1,ind_2);
    };

    template<class M1, class M2>
    Matrix& eval_scal_mat(Matrix& handle, M1&, const M2& B)
    {
        using FullMatrix = raw::Matrix<M1,struct_dense>;
        error::check_assign_2(1,1,B.rows(),B.cols());
        using value_type_2 = typename M2::value_type;

        return assign_functor_scal<M1,value_type_2>::eval(handle,B(1,1),ind_1,ind_2);
    };

    template<class M1, class M2>
    Matrix& eval_scal_scal(Matrix& handle, M1&, const M2& B)
    {
        using value_type_2 = M2;
        return assign_functor_scal<M1,value_type_2>::eval(handle,B,ind_1,ind_2);
    };
};

struct assign_visitor_3: public extract_type2_switch_nc<Matrix&,assign_visitor_3>
{
    const colon* ind;

    assign_visitor_3(const colon* ind) : ind(ind){};

    template<class M1, class M2>
    Matrix& eval_mat_mat(Matrix& handle, M1& A, const M2& B)
    {
        Integer r = B.rows(), c = B.cols();
        if ( r == 1 && c == 1)
        {
            using value_type_2 = typename M2::value_type;
            return assign_functor<M1,value_type_2>::eval(handle,A,B(1,1),*ind);
        };
        return assign_functor_mat<M1,M2>::eval(handle,B,*ind);

    };

    template<class M1, class M2>
    Matrix& eval_mat_scal(Matrix& handle, M1& A, const M2& B)
    {
        using value_type_2 = M2;
        return assign_functor<M1,value_type_2>::eval(handle,A,B,*ind);
    };

    template<class M1, class M2>
    Matrix& eval_scal_mat(Matrix& handle, M1& , const M2& B)
    {
        using FullMatrix = raw::Matrix<M1,struct_dense>;

        Integer r = B.rows(), c = B.cols();
        if ( r == 1 && c == 1)
        {
            using value_type_2 = typename M2::value_type;
            handle = matcl::full(handle);
            return assign_functor<FullMatrix,value_type_2>
                    ::eval(handle,handle.get_impl_unique<FullMatrix>(),B(1,1),*ind);
        };
        handle = matcl::full(handle);
        return assign_functor_mat<FullMatrix,M2>::eval(handle,B,*ind);
    };

    template<class M1, class M2>
    Matrix& eval_scal_scal(Matrix& handle, M1&, const M2& B)
    {
        using value_type_2 = M2;
        return assign_functor_scal<M1,value_type_2>::eval(handle,B,*ind);
    };
};

struct assign_visitor_4: public extract_type2_switch_nc<Matrix&,assign_visitor_4>
{
    const colon* ind_1;
    const colon* ind_2;

    assign_visitor_4(const colon* ind_1, const colon* ind_2)
        : ind_1(ind_1), ind_2(ind_2){};

    template<class M1, class M2>
    Matrix& eval_mat_mat(Matrix& handle, M1& A, const M2& B)
    {
        Integer r = B.rows(), c = B.cols();
        if ( r == 1 && c == 1)
        {
            using value_type_2 = typename M2::value_type;
            return assign_functor<M1,value_type_2>::eval(handle,A,B(1,1),*ind_1,*ind_2);
        };
        return assign_functor_mat<M1,M2>::eval(handle,B,*ind_1,*ind_2);
    };

    template<class M1, class M2>
    Matrix& eval_mat_scal(Matrix& handle, M1& A, const M2& B)
    {
        using value_type_2 = M2;
        return assign_functor<M1,value_type_2>::eval(handle,A,B,*ind_1,*ind_2);
    };

    template<class M1, class M2>
    Matrix& eval_scal_mat(Matrix& handle, M1& , const M2& B)
    {
        using FullMatrix = raw::Matrix<M1,struct_dense>;

        Integer r = B.rows(), c = B.cols();
        if ( r == 1 && c == 1)
        {
            using value_type_2 = typename M2::value_type;
            handle = matcl::full(handle);
            return assign_functor<FullMatrix,value_type_2>
                ::eval(handle,handle.get_impl_unique<FullMatrix>(),B(1,1),*ind_1,*ind_2);
        };

        handle = matcl::full(handle);
        return assign_functor_mat<FullMatrix,M2>::eval(handle,B,*ind_1,*ind_2);
    };

    template<class M1, class M2>
    Matrix& eval_scal_scal(Matrix& handle, M1&, const M2& B)
    {
        using value_type_2 = M2;
        return assign_functor_scal<M1,value_type_2>::eval(handle,B,*ind_1,*ind_2);
    };
};

struct assign_visitor_4_drop: public extract_type_switch<Matrix&,assign_visitor_4_drop,false>
{
    template<class Mat>
    static Matrix& eval(Matrix& handle, Mat& A, const colon* ind_1, const colon* ind_2, Real tol)
    {
        return assign_functor_drop<Mat>::eval(handle,A,*ind_1,*ind_2,tol);
    }

    template<class T>
    static Matrix& eval_scalar(Matrix& handle, T&, const colon*, const colon*, Real)
    {
        return handle;
    };
};

struct assign_visitor_4_add: public extract_type_switch<Matrix&,assign_visitor_4_add,false>
{
    template<class Mat>
    static Matrix& eval(Matrix& handle, Mat& A, const colon* ind_1, const colon* ind_2)
    {
        return assign_functor_add<Mat>::eval(handle,A,*ind_1,*ind_2);
    }

    template<class T>
    static Matrix& eval_scalar(Matrix& handle, T&, const colon*, const colon*)
    {
        return handle;
    };
};

struct assign_visitor_3_drop: public extract_type_switch<Matrix&,assign_visitor_3_drop,false>
{
    template<class Mat>
    static Matrix& eval(Matrix& handle, Mat& A, const colon* ind_1, Real tol)
    {
        return assign_functor_drop<Mat>::eval(handle,A,*ind_1,tol);
    }

    template<class T>
    static Matrix& eval_scalar(Matrix& handle, T&, const colon*, Real)
    {
        return handle;
    };
};

struct assign_visitor_3_add: public extract_type_switch<Matrix&,assign_visitor_3_add,false>
{
    template<class Mat>
    static Matrix& eval(Matrix& handle, Mat& A, const colon* ind_1)
    {
        return assign_functor_add<Mat>::eval(handle,A,*ind_1);
    }

    template<class T>
    static Matrix& eval_scalar(Matrix& handle, T&, const colon*)
    {
        return handle;
    };
};

struct assign_visitor_diag_drop: public extract_type_switch<Matrix&,assign_visitor_diag_drop,false>
{
    template<class Mat>
    static Matrix& eval(Matrix& handle, Mat& A, Integer d, Real tol)
    {
        return assign_functor_drop<Mat>::eval_diag(handle,A,d,tol);
    }

    template<class T>
    static Matrix& eval_scalar(Matrix& handle, T&, Integer, Real)
    {
        return handle;
    };
};

struct assign_visitor_diag_add: public extract_type_switch<Matrix&,assign_visitor_diag_add,false>
{
    template<class Mat>
    static Matrix& eval(Matrix& handle, Mat& A, Integer d)
    {
        return assign_functor_add<Mat>::eval_diag(handle,A,d);
    }

    template<class T>
    static Matrix& eval_scalar(Matrix& handle, T&, Integer)
    {
        return handle;
    };
};

struct assign_visitor_1_drop: public extract_type_switch<Matrix&,assign_visitor_1_drop,false>
{
    template<class Mat>
    static Matrix& eval(Matrix& handle, Mat& A,Integer ind1, Real tol)
    {        
        return assign_functor_drop<Mat>::eval(handle,A,ind1,tol);
    }

    template<class T>
    static Matrix& eval_scalar(Matrix& handle, T&,Integer, Real)
    {
        return handle;
    };
};

struct assign_visitor_1_add: public extract_type_switch<Matrix&,assign_visitor_1_add,false>
{
    template<class Mat>
    static Matrix& eval(Matrix& handle, Mat& A,Integer ind1)
    {
        return assign_functor_add<Mat>::eval(handle,A,ind1);
    }

    template<class T>
    static Matrix& eval_scalar(Matrix& handle, T&,Integer)
    {
        return handle;
    };
};

struct assign_visitor_2_drop: public extract_type_switch<Matrix&,assign_visitor_2_drop,false>
{
    template<class Mat>
    static Matrix& eval(Matrix& handle, Mat& A,Integer ind1, Integer ind2, Real tol)
    {
        return assign_functor_drop<Mat>::eval(handle,A,ind1, ind2, tol);
    }

    template<class T>
    static Matrix& eval_scalar(Matrix& handle, T&,Integer, Integer, Real)
    {
        return handle;
    };
};

struct assign_visitor_2_add: public extract_type_switch<Matrix&,assign_visitor_2_add,false>
{
    template<class Mat>
    static Matrix& eval(Matrix& handle, Mat& A,Integer ind1, Integer ind2)
    {
        return assign_functor_add<Mat>::eval(handle,A,ind1, ind2);
    }

    template<class T>
    static Matrix& eval_scalar(Matrix& handle, T&,Integer, Integer)
    {
        return handle;
    };
};

struct assign_visitor_diag: public extract_type2_switch_nc<Matrix&,assign_visitor_diag>
{
    Integer d;

    assign_visitor_diag(Integer d) : d(d){};

    template<class M1, class M2>
    Matrix& eval_mat_mat(Matrix& handle, M1& A, const M2& B)
    {
        error::check_diag(d,A.rows(),A.cols());

        Integer r = B.rows(), c = B.cols();
        if ( r == 1 && c == 1)
        {
            using value_type_2 = typename M2::value_type;
            return assign_functor<M1,value_type_2>::eval_diag(handle,A,B(1,1),d);
        };

        return assign_functor_mat<M1,M2>::eval_diag(handle,B,d);
    };

    template<class M1, class M2>
    Matrix& eval_mat_scal(Matrix& handle, M1& A, const M2& B)
    {
        error::check_diag(d,A.rows(),A.cols());

        using value_type_2 = M2;
        return assign_functor<M1,value_type_2>::eval_diag(handle,A,B,d);
    };
    
    template<class M1, class M2>
    Matrix& eval_scal_mat(Matrix& handle, M1&, const M2& B)
    {
        using FullMatrix = raw::Matrix<M1,struct_dense>;

        error::check_diag(d,1,1);

        error::check_assign_1(1,B.rows(),B.cols());
        using value_type_2 = typename M2::value_type;
        return assign_functor_scal<M1,value_type_2>::eval(handle,B(1,1),1);
    };
    
    template<class M1, class M2>
    Matrix& eval_scal_scal(Matrix& handle, M1&, const M2& B)
    {
        error::check_diag(d,1,1);

        using value_type_2 = M2;
        return assign_functor_scal<M1,value_type_2>::eval(handle,B,1);
    };
};

};

Matrix& matcl::sub_matrix::operator=(const Matrix& mat) const &&
{
    details::prepare_for_assign(*this->m_matrix, mat);

    if (m_colon_2)
        return details::assign_visitor_4(m_colon_1,m_colon_2).make(*m_matrix,mat);
    else if (m_colon_1)
        return details::assign_visitor_3(m_colon_1).make(*m_matrix,mat);
    else
        return details::assign_visitor_diag(m_d).make(*m_matrix,mat);
};

Matrix& matcl::sub_matrix_1::operator=(const Matrix& mat) const &&
{
    return details::assign_visitor_1(m_ind_1).make(*m_matrix,mat);
};

Matrix& matcl::sub_matrix_2::operator=(const Matrix& mat) const &&
{
    return details::assign_visitor_2(m_ind_1,m_ind_2).make(*m_matrix,mat);
};

Matrix& matcl::sub_matrix::drop_sparse(Real tol) const &&
{
    if (m_colon_2)
        return details::assign_visitor_4_drop::make<Matrix&>(*m_matrix,m_colon_1,m_colon_2,tol);
    else if (m_colon_1)
        return details::assign_visitor_3_drop::make<Matrix&>(*m_matrix, m_colon_1,tol);
    else
        return details::assign_visitor_diag_drop::make<Matrix&>(*m_matrix,m_d,tol);
};

Matrix& matcl::sub_matrix_1::drop_sparse(Real tol) const &&
{
    return details::assign_visitor_1_drop::make<Matrix&>(*m_matrix,m_ind_1,tol);
};

Matrix& matcl::sub_matrix_2::drop_sparse(Real tol) const &&
{
    return details::assign_visitor_2_drop::make<Matrix&>(*m_matrix,m_ind_1,m_ind_2,tol);
};

Matrix& matcl::sub_matrix::add_sparse() const &&
{
    if (m_colon_2)
        return details::assign_visitor_4_add::make<Matrix&>(*m_matrix,m_colon_1,m_colon_2);
    else if (m_colon_1)
        return details::assign_visitor_3_add::make<Matrix&>(*m_matrix,m_colon_1);
    else
        return details::assign_visitor_diag_add::make<Matrix&>(*m_matrix,m_d);
};

Matrix& matcl::sub_matrix_1::add_sparse() const &&
{
    return details::assign_visitor_1_add::make<Matrix&>(*m_matrix,m_ind_1);
};

Matrix& matcl::sub_matrix_2::add_sparse() const &&
{
    return details::assign_visitor_2_add::make<Matrix&>(*m_matrix,m_ind_1,m_ind_2);
};

void details::prepare_for_assign(const Matrix& sub1, const Matrix& rhs1, const Matrix& rhs2)
{
    const matrix_base& sub1_base    = details::matrix_data_accesser::get_base(sub1);       

    // if sub1 is scalar: do nothing
    if (sub1_base.m_type <= mat_code::integer_dense)
        return;

    if (sub1_base.is_effective_unique_mat())
        return;

    long count  = sub1_base.m_value.m_mat.m_refcount->get_count();

    // if sub1 is unique: do nothing
    if (count == 1)
        return;

    // if sub1 count > 3: do nothing, copying rhs1, rhs2 will not help
    if (count > 3)
        return;

    const matrix_base& sub2_base    = details::matrix_data_accesser::get_base(rhs1);
    const matrix_base& sub3_base    = details::matrix_data_accesser::get_base(rhs2);

    bool is_scal_2                  = (sub2_base.m_type <= mat_code::integer_dense);
    bool is_scal_3                  = (sub3_base.m_type <= mat_code::integer_dense);

    bool the_same_1                 = false;
    bool the_same_2                 = false;

    if (is_scal_2 == false && sub1_base.m_value.m_mat.m_refcount == sub2_base.m_value.m_mat.m_refcount)
        the_same_1  = true;

    if (is_scal_3 == false && sub1_base.m_value.m_mat.m_refcount == sub3_base.m_value.m_mat.m_refcount)
        the_same_2  = true;

    if (the_same_1 == true && the_same_2 == true)
    {
        //refcount pointers are the same; after making unique the rhs1, rhs2 matrices, 
        // refcount of sub1 should drop to one
        rhs1.make_unique();
        rhs2.make_unique();

        return;    
    };

    if (count == 2 && (the_same_1 == true || the_same_2 == true))
    {
        //refcount pointers are the same; after making unique the rhs1 or rhs2 matrix, 
        // refcount of sub1 should drop to one
        if (the_same_1)
            rhs1.make_unique();

        if (the_same_2)
            rhs2.make_unique();

        return;    
    };
};

void details::prepare_for_assign(const matcl::Matrix& sub1, const matcl::Matrix& sub2)
{
    const matrix_base& sub1_base    = details::matrix_data_accesser::get_base(sub1);
    
    // if sub1 is scalar: do nothing
    if (sub1_base.m_type <= mat_code::integer_dense)
        return;

    if (sub1_base.is_effective_unique_mat())
        return;

    long count  = sub1_base.m_value.m_mat.m_refcount->get_count();

    // if sub1 is unique: do nothing
    // if sub1 count > 2: do nothing, copying sub2 will not help
    if (count != 2)
        return;

    const matrix_base& sub2_base    = details::matrix_data_accesser::get_base(sub2);

    // if sub2 is scalar: do nothing
    if (sub2_base.m_type <= mat_code::integer_dense)
        return;

    // if refcount ptr of sub1 and sub2 are different, then do nothing
    if (sub1_base.m_value.m_mat.m_refcount != sub2_base.m_value.m_mat.m_refcount)
        return;

    //refcount pointers are the same; after making unique the sub2 matrix, 
    // refcount of sub1 should drop to one
    sub2.make_unique();
    return;    
}

Matrix& matcl::sub_matrix::operator=(const sub_matrix& mat0) const &&
{    
    return std::move(*this).operator=(mat0.to_matrix());
};

Matrix& matcl::sub_matrix::operator=(const sub_matrix_1& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

Matrix& matcl::sub_matrix::operator=(const sub_matrix_2& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

Matrix& matcl::sub_matrix_1::operator=(const sub_matrix& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

Matrix& matcl::sub_matrix_1::operator=(const sub_matrix_1& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

Matrix& matcl::sub_matrix_1::operator=(const sub_matrix_2& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

Matrix& matcl::sub_matrix_2::operator=(const sub_matrix& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

Matrix& matcl::sub_matrix_2::operator=(const sub_matrix_1& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

Matrix& matcl::sub_matrix_2::operator=(const sub_matrix_2& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

Matrix& sub_matrix::assign_unique(const Matrix& mat) const &&
{
    m_matrix->mark_unique(true);
    Matrix& ret = (std::move(*this) = mat);
    m_matrix->mark_unique(false);
    return ret;
}

Matrix& sub_matrix::assign_unique(const sub_matrix& mat) const &&
{
    m_matrix->mark_unique(true);
    Matrix& ret = (std::move(*this) = mat);
    m_matrix->mark_unique(false);
    return ret;
}

Matrix& sub_matrix::assign_unique(const sub_matrix_1& mat) const &&
{
    m_matrix->mark_unique(true);
    Matrix& ret = (std::move(*this) = mat);
    m_matrix->mark_unique(false);
    return ret;
}

Matrix& sub_matrix::assign_unique(const sub_matrix_2& mat) const &&
{
    m_matrix->mark_unique(true);
    Matrix& ret = (std::move(*this) = mat);
    m_matrix->mark_unique(false);
    return ret;
}

Matrix& sub_matrix_1::assign_unique(const Matrix& mat) const &&
{
    m_matrix->mark_unique(true);
    Matrix& ret = (std::move(*this) = mat);
    m_matrix->mark_unique(false);
    return ret;
}

Matrix& sub_matrix_1::assign_unique(const sub_matrix& mat) const &&
{
    m_matrix->mark_unique(true);
    Matrix& ret = (std::move(*this) = mat);
    m_matrix->mark_unique(false);
    return ret;
}

Matrix& sub_matrix_1::assign_unique(const sub_matrix_1& mat) const &&
{
    m_matrix->mark_unique(true);
    Matrix& ret = (std::move(*this) = mat);
    m_matrix->mark_unique(false);
    return ret;
}

Matrix& sub_matrix_1::assign_unique(const sub_matrix_2& mat) const &&
{
    m_matrix->mark_unique(true);
    Matrix& ret = (std::move(*this) = mat);
    m_matrix->mark_unique(false);
    return ret;
}

Matrix& sub_matrix_2::assign_unique(const Matrix& mat) const &&
{
    m_matrix->mark_unique(true);
    Matrix& ret = (std::move(*this) = mat);
    m_matrix->mark_unique(false);
    return ret;
}

Matrix& sub_matrix_2::assign_unique(const sub_matrix& mat) const &&
{
    m_matrix->mark_unique(true);
    Matrix& ret = (std::move(*this) = mat);
    m_matrix->mark_unique(false);
    return ret;
}

Matrix& sub_matrix_2::assign_unique(const sub_matrix_1& mat) const &&
{
    m_matrix->mark_unique(true);
    Matrix& ret = (std::move(*this) = mat);
    m_matrix->mark_unique(false);
    return ret;
}

Matrix& sub_matrix_2::assign_unique(const sub_matrix_2& mat) const &&
{
    m_matrix->mark_unique(true);
    Matrix& ret = (std::move(*this) = mat);
    m_matrix->mark_unique(false);
    return ret;
}

const Matrix Matrix::operator()(Integer i,Integer j) const
{
    switch(m_type)
    {
        case mat_code::integer_scalar:
        {
            error::check_index(i,j,1,1);
            return m_value.val_int;
        }
        case mat_code::real_scalar:
        {
            error::check_index(i,j,1,1);
            return m_value.val_real;
        }
        case mat_code::float_scalar:
        {
            error::check_index(i,j,1,1);
            return m_value.val_float;
        }
        case mat_code::complex_scalar:
        {
            error::check_index(i,j,1,1);
            return Complex(m_value.val_complex[0],m_value.val_complex[1]);
        }
        case mat_code::float_complex_scalar:
        {
            error::check_index(i,j,1,1);
            return Float_complex(m_value.val_fcomplex[0],m_value.val_fcomplex[1]);
        }
        case mat_code::object_scalar:
        {
            error::check_index(i,j,1,1);
            return get_object();
        }
        case mat_code::integer_dense:
        {
            using cont_type = details::matrix_container<Integer,struct_dense>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::real_dense:
        {
            using cont_type = details::matrix_container<Real,struct_dense>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::float_dense:
        {
            using cont_type = details::matrix_container<Float,struct_dense>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::complex_dense:
        {
            using cont_type = details::matrix_container<Complex,struct_dense>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::float_complex_dense:
        {
            using cont_type = details::matrix_container<Float_complex,struct_dense>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::object_dense:
        {
            using cont_type = details::matrix_container<Object,struct_dense>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::integer_sparse:
        {
            using cont_type = details::matrix_container<Integer,struct_sparse>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::real_sparse:	
        {
            using cont_type = details::matrix_container<Real,struct_sparse>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::float_sparse:	
        {
            using cont_type = details::matrix_container<Float,struct_sparse>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::complex_sparse:
        {
            using cont_type = details::matrix_container<Complex,struct_sparse>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::float_complex_sparse:
        {
            using cont_type = details::matrix_container<Float_complex,struct_sparse>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::object_sparse:
        {
            using cont_type = details::matrix_container<Object,struct_sparse>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::integer_band:	
        {
            using cont_type = details::matrix_container<Integer,struct_banded>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::real_band:	
        {
            using cont_type = details::matrix_container<Real,struct_banded>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::float_band:	
        {
            using cont_type = details::matrix_container<Float,struct_banded>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::complex_band:
        {
            using cont_type = details::matrix_container<Complex,struct_banded>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::float_complex_band:
        {
            using cont_type = details::matrix_container<Float_complex,struct_banded>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        case mat_code::object_band:
        {
            using cont_type = details::matrix_container<Object,struct_banded>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i,j);
        }
        default:
        {
            matcl_assert(0,"unknown case");
            throw;
        }
    };	
};

const Matrix Matrix::operator()(Integer i) const
{
    switch(m_type)
    {
        case mat_code::integer_scalar:
        {
            error::check_index(i, 1);
            return m_value.val_int;
        }
        case mat_code::real_scalar:
        {
            error::check_index(i, 1);
            return m_value.val_real;
        }
        case mat_code::float_scalar:
        {
            error::check_index(i, 1);
            return m_value.val_float;
        }
        case mat_code::complex_scalar:
        {
            error::check_index(i, 1);
            return Complex(m_value.val_complex[0],m_value.val_complex[1]);
        }
        case mat_code::float_complex_scalar:
        {
            error::check_index(i, 1);
            return Float_complex(m_value.val_fcomplex[0],m_value.val_fcomplex[1]);
        }
        case mat_code::object_scalar:
        {
            error::check_index(i, 1);
            return get_object();
        }
        case mat_code::integer_dense:
        {
            using cont_type = details::matrix_container<Integer,struct_dense>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::real_dense:
        {
            using cont_type = details::matrix_container<Real,struct_dense>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::float_dense:
        {
            using cont_type = details::matrix_container<Float,struct_dense>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::complex_dense:
        {
            using cont_type = details::matrix_container<Complex,struct_dense>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::float_complex_dense:
        {
            using cont_type = details::matrix_container<Float_complex,struct_dense>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::object_dense:
        {
            using cont_type = details::matrix_container<Object,struct_dense>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::integer_sparse:
        {
            using cont_type = details::matrix_container<Integer,struct_sparse>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::real_sparse:	
        {
            using cont_type = details::matrix_container<Real,struct_sparse>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::float_sparse:	
        {
            using cont_type = details::matrix_container<Float,struct_sparse>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::complex_sparse:
        {
            using cont_type = details::matrix_container<Complex,struct_sparse>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::float_complex_sparse:
        {
            using cont_type = details::matrix_container<Float_complex,struct_sparse>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::object_sparse:
        {
            using cont_type = details::matrix_container<Object,struct_sparse>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::integer_band:	
        {
            using cont_type = details::matrix_container<Integer,struct_banded>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::real_band:	
        {
            using cont_type = details::matrix_container<Real,struct_banded>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::float_band:	
        {
            using cont_type = details::matrix_container<Float,struct_banded>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::complex_band:
        {
            using cont_type = details::matrix_container<Complex,struct_banded>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::float_complex_band:
        {
            using cont_type = details::matrix_container<Float_complex,struct_banded>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        case mat_code::object_band:
        {
            using cont_type = details::matrix_container<Object,struct_banded>;
            return static_cast<const cont_type*>(m_value.m_mat.m_mat_ptr)->get_elem(i);
        }
        default:
        {
            matcl_assert(0,"unknown case");
            throw;
        }
    };	
};

namespace details
{
    struct get_functor_3 : public details::extract_type_switch<void,get_functor_3,true>
    {
        template<class T>
        static void eval(const Matrix& handle, const T& mat, const colon& c1, Matrix& ret)
        {
            details::get_sub_functor<T>::eval(handle,mat,c1, ret);
        };

        template<class T>
        static void eval_scalar(const Matrix& handle, const T& mat, const colon& c1, Matrix& ret)
        {
            details::get_sub_functor_scal<T>::eval(handle,mat,c1, ret);
        };
    };

    struct get_functor_4 : public details::extract_type_switch<void,get_functor_4,true>
    {
        template<class T>
        static void eval(const Matrix& handle, const T& mat, const colon& c1, const colon& c2, Matrix& ret)
        {
            details::get_sub_functor<T>::eval(handle,mat,c1,c2, ret);
        };

        template<class T>
        static void eval_scalar(const Matrix& handle, const T& mat, const colon& c1, const colon& c2, Matrix& ret)
        {
            details::get_sub_functor_scal<T>::eval(handle,mat,c1,c2, ret);
        };
    };
};

const Matrix Matrix::operator()(const colon& c1) const
{
    Matrix ret;
    details::get_functor_3::make<const Matrix&,const colon&>(*this,c1,ret);
    return ret;
}

const Matrix Matrix::operator()(const colon& c1,const colon& c2) const
{
    Matrix ret;
    details::get_functor_4::make<const Matrix&,const colon&,const colon&>(*this,c1,c2,ret);

    if (details::is_tr_view(c1,c2))
    {
        ret.set_struct(md::predefined_struct::get_rectangle_view(this->get_struct(), ret.is_square()));
    }
    return ret;
};

namespace details
{
    struct delrows_functor : public details::extract_type_switch<void,delrows_functor,true>
    {
        template<class T>
        static void eval(const Matrix& handle, const T& mat, Matrix& ret, const colon& c1, bool rvalue)
        {
            return details::del_rows_functor<T>::eval(ret,handle,mat,c1,rvalue);
        };

        template<class T>
        static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret, const colon& c1, bool)
        {
            return details::del_functor_scal<T>::eval(ret,handle,mat,c1,true);
        };
    };

    struct delcols_functor : public details::extract_type_switch<void,delcols_functor,true>
    {
        template<class T>
        static void eval(const Matrix& handle, const T& mat, Matrix& ret, const colon& c1, bool rvalue)
        {
            return details::del_cols_functor<T>::eval(ret,handle,mat,c1,rvalue);
        };

        template<class T>
        static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret, const colon& c1, bool)
        {
            return details::del_functor_scal<T>::eval(ret,handle,mat,c1,false);
        };
    };

    struct delrowscols_functor : public details::extract_type_switch<void,delrowscols_functor,true>
    {
        template<class T>
        static void eval(const Matrix& handle, const T& mat, Matrix& ret, const colon& c1,
                         const colon& c2, bool rvalue)
        {
            return details::del_rowscols_functor<T>::eval(ret,handle,mat,c1,c2,rvalue);
        };

        template<class T>
        static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret, const colon& c1,
                                const colon& c2, bool)
        {
            return details::del_functor_scal<T>::eval2(ret,handle,mat,c1,c2);
        };
    };

};

const Matrix Matrix::delrows(const colon& c) const &
{
    Matrix ret;
    details::delrows_functor::make<const Matrix&>(*this, ret, c, false);
    return ret;
};

const Matrix Matrix::delrows(const colon& c) const &&
{
    Matrix ret;
    details::delrows_functor::make<const Matrix&>(*this, ret, c, true);
    return ret;
};

const Matrix Matrix::delcols(const colon& c) const &
{
    Matrix ret;
    details::delcols_functor::make<const Matrix&>(*this, ret, c, false);
    return ret;
};

const Matrix Matrix::delcols(const colon& c) const &&
{
    Matrix ret;
    details::delcols_functor::make<const Matrix&>(*this, ret, c, true);
    return ret;
};

const Matrix Matrix::delrowscols(const colon& c1, const colon& c2) const &
{
    Matrix ret;
    details::delrowscols_functor::make<const Matrix&>(*this,ret, c1,c2, false);
    return ret;
};

const Matrix Matrix::delrowscols(const colon& c1, const colon& c2) const &&
{
    Matrix ret;
    details::delrowscols_functor::make<const Matrix&>(*this, ret, c1, c2, true);
    return ret;
};

const Matrix matcl::sub_matrix::to_matrix() const
{
    const Matrix& tmp = *m_matrix;
    if (m_colon_2)
        return tmp(*m_colon_1,*m_colon_2);
    else if(m_colon_1)
        return tmp(*m_colon_1);
    else
        return get_diag(tmp,m_d);
};

const Matrix matcl::sub_matrix_1::to_matrix() const
{
    const Matrix& tmp = *m_matrix;
    return tmp(m_ind_1);
};

const Matrix matcl::sub_matrix_2::to_matrix() const
{
    const Matrix& tmp = *m_matrix;
    return tmp(m_ind_1,m_ind_2);
};

Matrix& matcl::sub_matrix::get_view_info(Integer& first_row, Integer& rows, 
                            Integer& first_col, Integer& cols) const
{
    if (m_colon_1 == nullptr || m_colon_2 == nullptr)
        throw error::unable_to_create_view(error::unable_create_view_reason::invalid_colon_type);

    md::colon_info ci;
    md::make_index(m_matrix->rows(),m_matrix->cols(),*m_colon_1,*m_colon_2,ci);

    if (ci.r_flag == 0 || ci.c_flag == 0)
        throw error::unable_to_create_view(error::unable_create_view_reason::step_not_one);

    if (ci.r_step != 1 || ci.c_step != 1)
        throw error::unable_to_create_view(error::unable_create_view_reason::step_not_one);
    
    first_row   = ci.r_start;
    first_col   = ci.c_start;
    rows        = ci.r_size;
    cols        = ci.c_size;

    return *m_matrix;
};

Matrix& matcl::sub_matrix_1::get_view_info(Integer& first_row, Integer& rows, 
                            Integer& first_col, Integer& cols) const
{
    Integer r = m_matrix->rows();
    Integer c = m_matrix->cols();

    Integer i, j;
    details::pos2ind(m_ind_1,r,i,j);

    if (i < 0 || j < 0 || i >= r || j >= c)
        throw error::invalid_index(i+1, j+1, r, c);

    first_row   = i;
    first_col   = j;
    rows        = 1;
    cols        = 1;

    return *m_matrix;
};

Matrix& matcl::sub_matrix_2::get_view_info(Integer& first_row, Integer& rows, 
                            Integer& first_col, Integer& cols) const
{
    Integer r = m_matrix->rows();
    Integer c = m_matrix->cols();

    if (m_ind_1 < 1 || m_ind_2 < 1 || m_ind_1 > r || m_ind_2 > c)
        throw error::invalid_index(m_ind_1, m_ind_2, r, c);

    first_row   = m_ind_1;
    first_col   = m_ind_2;
    rows        = 1;
    cols        = 1;

    return *m_matrix;
};

template Matrix& matcl::sub_matrix::assign_scalar<Integer>(const Integer& val) const;
template Matrix& matcl::sub_matrix::assign_scalar<Real>(const Real& val) const;
template Matrix& matcl::sub_matrix::assign_scalar<Float>(const Float& val) const;
template Matrix& matcl::sub_matrix::assign_scalar<Complex>(const Complex& val) const;
template Matrix& matcl::sub_matrix::assign_scalar<Float_complex>(const Float_complex& val) const;
template Matrix& matcl::sub_matrix::assign_scalar<Object>(const Object& val) const;

template Matrix& matcl::sub_matrix_1::assign_scalar<Integer>(const Integer& val) const;
template Matrix& matcl::sub_matrix_1::assign_scalar<Real>(const Real& val) const;
template Matrix& matcl::sub_matrix_1::assign_scalar<Float>(const Float& val) const;
template Matrix& matcl::sub_matrix_1::assign_scalar<Complex>(const Complex& val) const;
template Matrix& matcl::sub_matrix_1::assign_scalar<Float_complex>(const Float_complex& val) const;
template Matrix& matcl::sub_matrix_1::assign_scalar<Object>(const Object& val) const;

template Matrix& matcl::sub_matrix_2::assign_scalar<Integer>(const Integer& val) const;
template Matrix& matcl::sub_matrix_2::assign_scalar<Real>(const Real& val) const;
template Matrix& matcl::sub_matrix_2::assign_scalar<Float>(const Float& val) const;
template Matrix& matcl::sub_matrix_2::assign_scalar<Complex>(const Complex& val) const;
template Matrix& matcl::sub_matrix_2::assign_scalar<Float_complex>(const Float_complex& val) const;
template Matrix& matcl::sub_matrix_2::assign_scalar<Object>(const Object& val) const;

};