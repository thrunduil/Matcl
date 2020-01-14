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

#pragma once

#include "matcl-matrep/container/matrix_container.h"
#include "matcl-matrep/container/matrix2.inl"
#include "matcl-internals/base/utils.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl { namespace details
{

namespace gr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

inline matrix_container_base::matrix_container_base(matcl::mat_code type)
:m_type(type),refcount_str(1)
{};

template<class val_type, class str_type>
inline matrix_container<val_type,str_type>::matrix_container(const Matrix& mat)
:matrix_container_base(matrix_traits::mat_type_info_type_2<val_type,str_type>::matrix_code)
,m_matrix(mat, Matrix::copy_is_safe())
{};

template<class val_type, class str_type>
inline matrix_container<val_type,str_type>::matrix_container(Matrix&& mat)
:matrix_container_base(matrix_traits::mat_type_info_type_2<val_type,str_type>::matrix_code)
,m_matrix(std::move(mat))
{};

template<class val_type, class str_type>
void * matrix_container<val_type,str_type>::operator new(size_t )
{
    using matrix_type = raw::Matrix<val_type,str_type>;
    void* ptr = container_pool::malloc<(int)type_to_code<matrix_type>::value>();
    return ptr;
};

template<class val_type, class str_type>
void   matrix_container<val_type,str_type>::operator delete(void *p, size_t)
{
    if (p)
        free_container((matrix_container_base*)p);
};

template<class val_type, class str_type>
struct get_diag_info
{
    using Matrix = raw::Matrix<val_type,str_type>;

    static Integer eval_struct_udiags(const Matrix& m, bool use_flag)
    {
        Integer c = m.cols();
        if (c <= 1)
            return 0;

        if (use_flag)
        {
            struct_flag::diag_type dt = m.get_struct().get_udiags();

            if (dt == struct_flag::zero)
                return 0;

            if (c <= 2)
                return 1;

            if (dt != struct_flag::general)
                return 1;
        };

        return c - 1;
    };

    static Integer eval_struct_ldiags(const Matrix& m, bool use_flag)
    {
        Integer r = m.rows();
        if (r <= 1)
            return 0;

        if (use_flag)
        {
            struct_flag::diag_type dt = m.get_struct().get_ldiags();

            if (dt == struct_flag::zero)
                return 0;

            if (r <= 2)
                return 1;

            if (dt != struct_flag::general)
                return 1;
        };

        return r - 1;
    };
};

template<class val_type>
struct get_diag_info<val_type,struct_banded>
{
    using Matrix = raw::Matrix<val_type,struct_banded>;

    static Integer eval_struct_udiags(const Matrix& m, bool use_flag)
    {
        Integer nud = m.number_superdiagonals();

        if (nud == 0)
            return 0;

        if (use_flag)
        {
            struct_flag::diag_type dt = m.get_struct().get_udiags();

            if (dt == struct_flag::zero)
                return 0;

            if (nud == 1)
                return 1;

            if (dt != struct_flag::general)
                return 1;
        };

        return nud;
    };

    static Integer eval_struct_ldiags(const Matrix& m, bool use_flag)
    {
        Integer nld = m.number_subdiagonals();

        if (nld == 0)
            return 0;

        if (use_flag)
        {
            struct_flag::diag_type dt = m.get_struct().get_ldiags();

            if (dt == struct_flag::zero)
                return 0;

            if (nld == 1)
                return 1;

            if (dt != struct_flag::general)
                return 1;
        };

        return nld;
    };
};

template<class val_type, class str_type>
inline Integer matrix_container<val_type,str_type>::struct_ldiags_impl(bool use_flag) const
{
    return get_diag_info<val_type,str_type>::eval_struct_ldiags(m_matrix, use_flag);
};

template<class val_type, class str_type>
inline Integer matrix_container<val_type,str_type>::struct_udiags_impl(bool use_flag) const
{
    return get_diag_info<val_type,str_type>::eval_struct_udiags(m_matrix, use_flag);
};

template<class val_type, class str_type>
bool matrix_container<val_type,str_type>::is_scalar_true_impl() const
{
    Integer r = m_matrix.rows(), c = m_matrix.cols();

    //if (r == 0 && c == 0)
    //    return false;

    error::check_scalar(r,c);
    return mrd::cast_bool_helper<val_type>::eval(m_matrix(1,1));
};

template<class val_type, class str_type>
struct get_elem_helper
{};

template<class val_type>
struct get_elem_helper<val_type,struct_dense>
{
    using Matrix = typename raw::Matrix<val_type,struct_dense>;

    static const val_type& eval(const Matrix& m, Integer i,Integer j)
    {
        error::check_index(i, j, m.rows(), m.cols());
        return m.ptr()[i-1+(j-1)*m.ld()];
    };

    static const val_type& eval(const Matrix& m,Integer i)
    {
        error::check_index(i, m.size());

        Integer j;
        pos2ind(i,m.rows(),i,j);

        return m.ptr()[i+j*m.ld()];
    };
};

template<class val_type>
struct get_elem_helper<val_type,struct_sparse>
{
    using Matrix = typename raw::Matrix<val_type,struct_sparse>;

    static val_type eval(const Matrix& m, Integer i,Integer j)
    {
        error::check_index(i, j, m.rows(), m.cols());

        const raw::details::sparse_ccs<val_type>& rep = m.rep();
        Integer k;
        if (rep.has_element(i-1,j-1,k) == false)
            return default_value<val_type>(m.get_type());

        return val_type(rep.ptr_x()[k]);
    };

    static val_type eval(const Matrix& m,Integer i)
    {
        Integer r = m.rows();
        if (r == 0)
            error::check_index(i, 0);

        Integer j;
        pos2ind(i,r,i,j);

        return eval(m,i+1,j+1);
    };
};

template<class val_type>
struct get_elem_helper<val_type,struct_banded>
{
    using Matrix = typename raw::Matrix<val_type,struct_banded>;

    static val_type eval(const Matrix& mat, Integer i,Integer j)
    {
        Integer r = mat.rows();
        Integer c = mat.cols();
        --i;
        --j;

        if (i < 0 || j < 0 || i >= r || j >= c)
            throw error::invalid_index(i+1, j+1, r, c);

        if ( i < mat.first_row(j) || i > mat.last_row(j) )
            return default_value<val_type>(mat.get_type());

        return val_type(mat.rep_ptr()[mat.element_pos(i,j)]);
    };

    static val_type eval(const Matrix& mat,Integer i)
    {
        Integer r = mat.rows();
        if (r == 0)
            error::check_index(i, 0);

        Integer j;
        pos2ind(i,r,i,j);
        return eval(mat,i+1,j+1);
    };
};

template<class val_type, class str_type>
inline Matrix matrix_container<val_type,str_type>::get_elem(Integer i) const
{
    return get_elem_helper<val_type,str_type>::eval(m_matrix,i);
};

template<class val_type, class str_type>
inline Matrix matrix_container<val_type,str_type>::get_elem(Integer i,Integer j) const
{
    return get_elem_helper<val_type,str_type>::eval(m_matrix,i,j);
};

template<class val_type, class str_type>
matcl::Matrix matrix_container<val_type,str_type>::get_scalar_impl() const
{
    Integer r = rows();
    Integer c = cols();

    if (r != 1 || c != 1)
        throw error::scalar_required(r,c);

    return get_elem(1,1);
};

template<class val_type, class str_type>
inline void matrix_container<val_type,str_type>::set_struct_impl(const struct_flag& fl) const
{
    return m_matrix.set_struct(fl);
};

template<class val_type, class str_type>
inline void matrix_container<val_type,str_type>::add_struct_impl(const struct_flag& fl) const
{
    return m_matrix.add_struct(fl);
};

template<class val_type, class str_type>
inline bool matrix_container<val_type,str_type>::is_same_matrix_impl(const matrix_container_base* other) const
{
    const matrix_container* other_exact = static_cast<const matrix_container*>(other);
    return m_matrix.is_same_matrix(other_exact->m_matrix);
};

template<class val_type, class str_type>
inline typename matrix_container<val_type,str_type>::Matrix
matrix_container<val_type,str_type>::reserve_impl(Integer r, Integer c) const
{
    return m_matrix.reserve(r,c);
};

template<class val_type, class str_type>
struct resize_functor
{
    using Matrix = raw::Matrix<val_type,str_type>;

    static Matrix eval(const Matrix& mat, Integer r, Integer c)
    {
        if (mat.is_unique())
            return const_cast<Matrix&>(mat).resize(r,c);
        else
            return mat.resize(r,c);
    };
};

template<class val_type, class str_type>
inline typename matrix_container<val_type,str_type>::Matrix
matrix_container<val_type,str_type>::resize_impl(Integer r, Integer c) const
{
    return resize_functor<val_type,str_type>::eval(m_matrix,r,c);
};

template<class c_type, class val_ty, class struct_ty>
struct matrix_container_cons
{
    using type = matrix_container<val_ty,struct_ty>;
};

template<class c_type, class val_ty, class struct_ty>
struct matrix_container_cons<const c_type,val_ty,struct_ty>
{
    using type = const matrix_container<val_ty,struct_ty>;
};

template<class ret, class derived>
struct eval_type
{
    template<class c_type>
    static ret make(c_type* ptr)
    {
        switch(ptr->get_type())
        {
            case mat_code::integer_dense:
            {
                using type = typename matrix_container_cons<c_type,Integer,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::real_dense:
            {
                using type = typename matrix_container_cons<c_type,Real,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::float_dense:
            {
                using type = typename matrix_container_cons<c_type,Float,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::complex_dense:
            {
                using type = typename matrix_container_cons<c_type,Complex,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::float_complex_dense:
            {
                using type = typename matrix_container_cons<c_type,Float_complex,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::object_dense:
            {
                using type = typename matrix_container_cons<c_type,Object,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::integer_sparse:
            {
                using type = typename matrix_container_cons<c_type,Integer,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::real_sparse:
            {
                using type = typename matrix_container_cons<c_type,Real,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::float_sparse:
            {
                using type = typename matrix_container_cons<c_type,Float,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::complex_sparse:
            {
                using type = typename matrix_container_cons<c_type,Complex,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::float_complex_sparse:
            {
                using type = typename matrix_container_cons<c_type,Float_complex,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::object_sparse:
            {
                using type = typename matrix_container_cons<c_type,Object,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::integer_band:
            {
                using type = typename matrix_container_cons<c_type,Integer,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::real_band:
            {
                using type = typename matrix_container_cons<c_type,Real,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::float_band:
            {
                using type = typename matrix_container_cons<c_type,Float,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::complex_band:
            {
                using type = typename matrix_container_cons<c_type,Complex,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::float_complex_band:
            {
                using type = typename matrix_container_cons<c_type,Float_complex,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            case mat_code::object_band:
            {
                using type = typename matrix_container_cons<c_type,Object,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr));
            }
            default:
                matcl_assert(0,"invalid case");
                throw;
        };
    };

    template<class Arg,class c_type>
    static ret make(c_type* ptr, typename details::hide_type<Arg>::type arg)
    {
        switch(ptr->get_type())
        {
            case mat_code::integer_dense:
            {
                using type = typename matrix_container_cons<c_type, Integer,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::real_dense:
            {
                using type = typename matrix_container_cons<c_type, Real,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::float_dense:
            {
                using type = typename matrix_container_cons<c_type, Float,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::complex_dense:
            {
                using type = typename matrix_container_cons<c_type, Complex,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::float_complex_dense:
            {
                using type = typename matrix_container_cons<c_type, Float_complex,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::object_dense:
            {
                using type = typename matrix_container_cons<c_type, Object,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::integer_sparse:
            {
                using type = typename matrix_container_cons<c_type, Integer,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::real_sparse:
            {
                using type = typename matrix_container_cons<c_type, Real,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::float_sparse:
            {
                using type = typename matrix_container_cons<c_type, Float,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::complex_sparse:
            {
                using type = typename matrix_container_cons<c_type, Complex,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::float_complex_sparse:
            {
                using type = typename matrix_container_cons<c_type, Float_complex,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::object_sparse:
            {
                using type = typename matrix_container_cons<c_type, Object,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::integer_band:
            {
                using type = typename matrix_container_cons<c_type, Integer,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::real_band:
            {
                using type = typename matrix_container_cons<c_type, Real,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::float_band:
            {
                using type = typename matrix_container_cons<c_type, Float,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::complex_band:
            {
                using type = typename matrix_container_cons<c_type, Complex,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::float_complex_band:
            {
                using type = typename matrix_container_cons<c_type, Float_complex,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            case mat_code::object_band:
            {
                using type = typename matrix_container_cons<c_type, Object,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr),arg);
            }
            default:
                matcl_assert(0,"invalid case");
                throw;
        };
    };		

    template<class Arg1, class Arg2,class c_type>
    static ret make(c_type* ptr, typename details::hide_type<Arg1>::type arg1, 
                                 typename details::hide_type<Arg2>::type arg2)
    {
        switch(ptr->get_type())
        {
            case mat_code::integer_dense:
            {
                using type = typename matrix_container_cons<c_type, Integer,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::real_dense:
            {
                using type = typename matrix_container_cons<c_type, Real,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::float_dense:
            {
                using type = typename matrix_container_cons<c_type, Float,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::complex_dense:
            {
                using type = typename matrix_container_cons<c_type, Complex,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::float_complex_dense:
            {
                using type = typename matrix_container_cons<c_type, Float_complex,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::object_dense:
            {
                using type = typename matrix_container_cons<c_type, Object,struct_dense>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::integer_sparse:
            {
                using type = typename matrix_container_cons<c_type, Integer,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::real_sparse:
            {
                using type = typename matrix_container_cons<c_type, Real,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::float_sparse:
            {
                using type = typename matrix_container_cons<c_type, Float,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::complex_sparse:
            {
                using type = typename matrix_container_cons<c_type, Complex,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::float_complex_sparse:
            {
                using type = typename matrix_container_cons<c_type, Float_complex,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::object_sparse:
            {
                using type = typename matrix_container_cons<c_type, Object,struct_sparse>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::integer_band:
            {
                using type = typename matrix_container_cons<c_type, Integer,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::real_band:
            {
                using type = typename matrix_container_cons<c_type, Real,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::float_band:
            {
                using type = typename matrix_container_cons<c_type, Float,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::complex_band:
            {
                using type = typename matrix_container_cons<c_type, Complex,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::float_complex_band:
            {
                using type = typename matrix_container_cons<c_type, Float_complex,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            case mat_code::object_band:
            {
                using type = typename matrix_container_cons<c_type, Object,struct_banded>::type;
                return derived::eval(static_cast<type*>(ptr),arg1,arg2);
            }
            default:
                matcl_assert(0,"invalid case");
                throw;
        };
    };		
};

struct eval_type_rows : public eval_type<Integer,eval_type_rows>
{
    template<class T>
    static Integer eval(T* mat)                     { return mat->rows_impl(); };
};

struct eval_type_cols : public eval_type<Integer,eval_type_cols>
{
    template<class T>
    static Integer eval(T* mat)                     { return mat->cols_impl(); };
};

struct eval_type_struct_nnz : public eval_type<Integer,eval_type_struct_nnz>
{
    template<class T>
    static Integer eval(T* mat)                     { return mat->struct_nnz_impl(); };
};

struct eval_type_struct_ldiags : public eval_type<Integer,eval_type_struct_ldiags>
{
    template<class T>
    static Integer eval(T* mat, bool use_flags)     { return mat->struct_ldiags_impl(use_flags); };
};

struct eval_type_struct_udiags : public eval_type<Integer,eval_type_struct_udiags>
{
    template<class T>
    static Integer eval(T* mat, bool use_flags)     { return mat->struct_udiags_impl(use_flags); };
};

struct eval_type_size : public eval_type<int_tup_2,eval_type_size>
{
    template<class T>
    static int_tup_2 eval(T* mat)                   { return mat->size_impl(); };
};

struct eval_type_is_scalar_true : public eval_type<bool,eval_type_is_scalar_true>
{
    template<class T>
    static bool eval(T* mat)                        { return mat->is_scalar_true_impl(); };
};

struct eval_type_get_ti : public eval_type<ti::ti_object,eval_type_get_ti>
{
    template<class T>
    static ti::ti_object eval(T* mat)               
    { 
        using val_type = typename T::Matrix::value_type;
        return ti::convert_ti_object<val_type>(mat->get_ti_impl()); 
    };
};

struct eval_type_disp : public eval_type<void,eval_type_disp>
{
    template<class T>
    static void eval(T* mat, const disp_stream_ptr& os) { return mat->disp_impl(os); };
};

struct eval_type_get_struct : public eval_type<struct_flag&,eval_type_get_struct>
{
    template<class T>
    static struct_flag& eval(T* mat)                { return mat->get_struct_impl(); };
};

struct eval_type_get_struct_const : public eval_type<const struct_flag&,eval_type_get_struct_const>
{
    template<class T>
    static const struct_flag& eval(const T* mat)    { return mat->get_struct_impl(); };
};

struct eval_type_set_struct : public eval_type<void,eval_type_set_struct>
{
    template<class T>
    static void eval(const T* mat, struct_flag fl)  { mat->set_struct_impl(fl); };
};

struct eval_type_add_struct : public eval_type<void,eval_type_add_struct>
{
    template<class T>
    static void eval(const T* mat, struct_flag fl)  { mat->add_struct_impl(fl); };
};

struct eval_type_get_scalar : public eval_type<matcl::Matrix,eval_type_get_scalar>
{
    template<class T>
    static matcl::Matrix eval(T* mat)				{ return mat->get_scalar_impl(); };
};

struct eval_type_reserve : public eval_type<matcl::Matrix,eval_type_reserve>
{
    template<class T>
    static matcl::Matrix eval(T* mat, Integer r, Integer c)  
    { 
        return matcl::Matrix(mat->reserve_impl(r,c),false); 
    };
};

struct eval_type_reserve_band : public eval_type<matcl::Matrix,eval_type_reserve_band>
{
    template<class V, class S>
    static matcl::Matrix eval(const matrix_container<V,struct_dense>* mat, const S& arg)  
    { 
        return matcl::Matrix(mat->reserve_impl(arg.r,arg.c), false); 
    };

    template<class V, class S>
    static matcl::Matrix eval(const matrix_container<V,struct_sparse>* mat, const S& arg)  
    { 
        return matcl::Matrix(mat->reserve_impl(arg.r,arg.c), false); 
    };

    template<class V, class S>
    static matcl::Matrix eval(const matrix_container<V,struct_banded>* mat, const S& arg)  
    { 
        using BM = raw::Matrix<V,struct_banded>;
        return matcl::Matrix(mat->get().reserve(arg.r,arg.c,arg.fd,arg.ld),false); 
    };
};

struct eval_type_resize : public eval_type<matcl::Matrix,eval_type_resize>
{
    template<class T>
    static matcl::Matrix eval(T* mat, Integer r, Integer c)  
    { 
        return matcl::Matrix(mat->resize_impl(r,c), false); 
    };
};

struct eval_type_resize_band : public eval_type<matcl::Matrix,eval_type_resize_band>
{
    template<class V, class S>
    static matcl::Matrix eval(const matrix_container<V,struct_dense>* mat, const S& arg)  
    { 
        return matcl::Matrix(mat->resize_impl(arg.r,arg.c), false); 
    };

    template<class V, class S>
    static matcl::Matrix eval(const matrix_container<V,struct_sparse>* mat, const S& arg)  
    { 
        return matcl::Matrix(mat->resize_impl(arg.r,arg.c), false); 
    };

    template<class V, class S>
    static matcl::Matrix eval(const matrix_container<V,struct_banded>* mat, const S& arg)  
    { 
        using BM = raw::Matrix<V,struct_banded>;
        return matcl::Matrix(mat->get().resize(arg.r,arg.c,arg.fd,arg.ld),false); 
    };
};

};};