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

#include "matcl-matrep/matrix/matrix_rep_dense.h"
#include "matcl-matrep/matrix/matrix_rep_sparse.h"
#include "matcl-matrep/matrix/matrix_rep_band.h"
#include "matcl-matrep/matrix/matrix_rep_functions.h"
#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/details/matrix.inl"
#include "matcl-matrep/details/matrix_details_subs.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-core/details/integer.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"

#define MATCL_MATREP_EXPORT_REP MATCL_MATREP_EXPORT

namespace matcl
{

template<class T>
dense_matrix<T,true>::dense_matrix()
    :m_matrix(T(0))
{
    const mat_type& mat= m_matrix.impl_unique<mat_type>();
    init_from_rep(mat);
};

template<class T>
dense_matrix<T,false>::dense_matrix()
    :base_type()
{};

template<class T>
dense_matrix<T,true>::dense_matrix(bool val)
    :m_matrix(matcl::Matrix(T(val)))
{
    const mat_type& mat= m_matrix.impl_unique<mat_type>();
    init_from_rep(mat);
};

template<class T>
dense_matrix<T,true>::dense_matrix(const T& val)
    :m_matrix(matcl::Matrix(val))
{
    const mat_type& mat= m_matrix.impl_unique<mat_type>();
    init_from_rep(mat);
};

template<class T>
void dense_matrix<T,true>::update_rep()
{
    const mat_type& mat = m_matrix.impl<mat_type>();
    m_matrix            = Matrix(mat,false);
    init_from_rep(mat);
};

template<class T>
void dense_matrix<T,true>::init_from_rep(const mat_type& mat)
{
    m_rows      = mat.rows();
    m_cols      = mat.cols();
    m_max_rows  = mat.max_rows();
    m_max_cols  = mat.max_cols();
    m_ld        = mat.ld();
    m_size      = mat.size();
    m_ptr       = const_cast<T*>(mat.ptr());
    m_flag      = const_cast<struct_flag*>(&mat.get_struct());
};

template<class T>
dense_matrix<T,true>::dense_matrix(const matcl::Matrix& m)
{
    const mat_type& mat = m.impl<mat_type>();
    m_matrix            = Matrix(mat,false);
    init_from_rep(mat);
};

template<class T>
dense_matrix<T,true>::dense_matrix(matcl::Matrix&& m0)
{
    Matrix m            = std::move(m0);
    const mat_type& mat = m.impl<mat_type>();
    m_matrix            = Matrix(mat,false);
    init_from_rep(mat);
};

template<class T>
dense_matrix<T,true>::dense_matrix(matcl::Matrix& m, str_make_unique)
{
    mat_type mat= m.impl_unique<mat_type>();

    //assignment to m_matrix must be done after getting impl type
    //otherwise m cannot be unique, and memory will be copied always
    m_matrix    = m;

    init_from_rep(mat);
};

template<class T>
dense_matrix<T,true>::dense_matrix(matcl::Matrix&& m, str_make_unique)
{
    mat_type mat= m.impl_unique<mat_type>();

    //assignment to m_matrix must be done after getting impl type
    //otherwise m cannot be unique, and memory will be copied always
    m_matrix    = std::move(m);

    init_from_rep(mat);
};

template<class T>
dense_matrix<T,false>::dense_matrix(matcl::Matrix& m)
    :base_type(m, str_make_unique())
{};

template<class T>
dense_matrix<T,false>::dense_matrix(matcl::Matrix&& m)
    :base_type(m, str_make_unique())
{};

template<class T>
dense_matrix<T,true>::dense_matrix(const sub_dense_matrix& sub)
    :dense_matrix(sub.to_matrix())
{};

template<class T>
dense_matrix<T,true>::dense_matrix(const dense_row<T>& con)
    :dense_matrix(con.to_matrix())
{};

template<class T>
dense_matrix<T,true>::dense_matrix(const dense_col<T>& con)
    :dense_matrix(con.to_matrix())
{};

template<class T>
dense_matrix<T,true>::dense_matrix(const dense_matrix& mat)
    :m_matrix(mat.m_matrix)
{
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
dense_matrix<T,true>::dense_matrix(dense_matrix&& mat)
    :m_matrix(std::move(mat.m_matrix))
{
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
dense_matrix<T,true>::dense_matrix(const sparse_matrix<T>& mat)
    :dense_matrix(mat.to_matrix())
{};

template<class T>
dense_matrix<T,true>::dense_matrix(sparse_matrix<T>&& mat)
    :dense_matrix(mat.to_matrix())
{};

template<class T>
dense_matrix<T,true>::dense_matrix(const band_matrix<T>& mat)
    :dense_matrix(mat.to_matrix())
{};

template<class T>
dense_matrix<T,true>::dense_matrix(band_matrix<T>&& mat)
    :dense_matrix(mat.to_matrix())
{};

template<class T>
dense_matrix<T,false>::dense_matrix(dense_matrix<T,false>&& mat)
    :dense_matrix(std::move(mat.m_matrix))
{
    init_from_rep(m_matrix.get_impl<mat_type>());
}

template<class T>
dense_matrix<T,false>::dense_matrix(dense_matrix<T,true>&& mat)
    :dense_matrix(mat.m_matrix)
{}

template<class T>
dense_matrix<T,false>::dense_matrix(dense_matrix<T,true>& mat)
    :dense_matrix(mat.m_matrix)
{}

template<class T>
dense_matrix<T,true>& dense_matrix<T,true>::operator=(const dense_matrix& mat) &
{
    m_matrix = mat.m_matrix;
    init_from_rep(m_matrix.get_impl<mat_type>());
    return *this;
};

template<class T>
dense_matrix<T,true>& dense_matrix<T,true>::operator=(dense_matrix&& mat) &
{
    m_matrix = std::move(mat.m_matrix);
    init_from_rep(m_matrix.get_impl<mat_type>());
    return *this;
};

template<class T>
dense_matrix<T,false>& dense_matrix<T,false>::operator=(dense_matrix&& mat) &
{
    m_matrix = std::move(mat.m_matrix);
    init_from_rep(m_matrix.get_impl<mat_type>());
    return *this;
};

template<class T>
dense_matrix<T,true>::~dense_matrix()
{};

template<class T>
dense_matrix<T,false>::~dense_matrix()
{};

template<class T>
void dense_matrix<T,true>::set_struct(const struct_flag& sc) const
{
    m_flag->set(sc);
};

template<class T>
void dense_matrix<T,true>::add_struct(const struct_flag& sc) const
{
    m_flag->add(sc);
};

template<class T>
ti::ti_object dense_matrix<T,true>::get_type() const
{
    return m_matrix.get_type();
};

template<class T>
const typename dense_matrix<T,true>::dense_matrix 
dense_matrix<T,true>::delrows(const colon& c) const &
{
    return dense_matrix(m_matrix.delrows(c));
}

template<class T>
const typename dense_matrix<T,true>::dense_matrix 
dense_matrix<T,true>::delrows(const colon& c) const &&
{
    return dense_matrix(std::move(m_matrix).delrows(c));
}

template<class T>
const typename dense_matrix<T,true>::dense_matrix
dense_matrix<T,true>::delcols(const colon& c) const &
{
    return dense_matrix(m_matrix.delcols(c));
};

template<class T>
const typename dense_matrix<T,true>::dense_matrix
dense_matrix<T,true>::delcols(const colon& c) const &&
{
    return dense_matrix(std::move(m_matrix).delcols(c));
};

template<class T>
const typename dense_matrix<T,true>::dense_matrix 
dense_matrix<T,true>::delrowscols(const colon& c1, const colon& c2) const &
{
    return dense_matrix(m_matrix.delrowscols(c1, c2));
};

template<class T>
const typename dense_matrix<T,true>::dense_matrix 
dense_matrix<T,true>::delrowscols(const colon& c1, const colon& c2) const &&
{
    return dense_matrix(std::move(m_matrix).delrowscols(c1, c2));
};

template<class T>
const typename dense_matrix<T,true>::dense_matrix
dense_matrix<T,true>::operator()(const colon& r) const
{
    return dense_matrix(m_matrix(r));
};

template<class T>
const typename dense_matrix<T,true>::dense_matrix
dense_matrix<T,true>::operator()(const colon& r, const colon& c) const
{
    return dense_matrix(m_matrix(r,c));
};

template<class T>
const typename dense_matrix<T,true>::dense_matrix
dense_matrix<T,true>::diag(Integer d) const
{
    return dense_matrix(m_matrix.diag(d));
};

template<class T>
Integer dense_matrix<T,true>::length() const
{
    Integer r = rows();
    Integer c = cols();

    if (r == 0 || c == 0)
        return 0;

    return (r > c) ? r : c;
};

template<class T>
Integer dense_matrix<T,true>::structural_nnz() const
{
    return m_matrix.structural_nnz();
};

template<class T>
Integer dense_matrix<T,true>::structural_ldiags(bool use_flags) const
{
    return m_matrix.structural_ldiags(use_flags);
};

template<class T>
Integer dense_matrix<T,true>::structural_udiags(bool use_flags) const
{
    return m_matrix.structural_udiags(use_flags);
};

template<class T>
Real dense_matrix<T,true>::numel() const
{
    return Real(rows()) * Real(cols());
};

template<class T>
bool dense_matrix<T,true>::all_finite() const
{
    return m_matrix.all_finite();
};

template<class T>
void dense_matrix<T,true>::throw_error_single_index(Integer r, Integer size) const
{
    throw error::invalid_single_index(r,size);
};

template<class T>
void dense_matrix<T,true>::throw_error_double_index(Integer r, Integer c, Integer rows, Integer cols) const
{
    throw error::invalid_double_index(r, c, rows, cols);
};

template<class T>
const typename dense_matrix<T,true>::dense_matrix
dense_matrix<T,true>::clone() const
{
    return dense_matrix(m_matrix.clone());
};

template<class T>
const typename dense_matrix<T,true>::dense_matrix&
dense_matrix<T,true>::make_unique() const
{
    m_matrix.make_unique();
    return *this;
};

template<class T>
void dense_matrix<T,true>::resize(Integer r, Integer c)
{
    m_matrix.resize(r,c);
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
void dense_matrix<T,true>::reserve(Integer r, Integer c)
{
    m_matrix.reserve(r,c);
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
void dense_matrix<T,true>::resize_band(Integer r, Integer c, Integer ld, Integer ud)
{
    (void)ld;
    (void)ud;
    m_matrix.resize(r,c);
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
void dense_matrix<T,true>::reserve_band(Integer r, Integer c, Integer ld, Integer ud)
{
    (void)ld;
    (void)ud;
    m_matrix.reserve(r,c);
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
sub_dense_matrix<T> dense_matrix<T,true>::operator()(const colon& r)
{
    sub_dense_matrix ret(this,r);
    return ret;
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::zeros(Integer r, Integer c)
{
    return dense_matrix<T>(matcl::zeros(r,c, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::zeros(ti::ti_object ti, Integer r, Integer c)
{
    return dense_matrix<T>(matcl::zeros(ti, r,c));
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::ones(Integer r, Integer c)
{
    return dense_matrix<T>(matcl::ones(r,c, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::ones(ti::ti_object ti, Integer r, Integer c)
{
    return dense_matrix<T>(matcl::ones(ti, r,c));
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::eye(Integer r, Integer c)
{
    return dense_matrix<T>(matcl::eye(r,c, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::eye(Integer r)
{
    return dense_matrix<T>(matcl::eye(r,r, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::eye(ti::ti_object ti, Integer r, Integer c)
{
    return dense_matrix<T>(matcl::eye(ti, r,c));
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::eye(ti::ti_object ti, Integer r)
{
    return dense_matrix<T>(matcl::eye(ti, r,r));
};

template<class T>
struct make_helper{};

template<>
struct make_helper<Integer>
{
    using T         = Integer;
    using mat_type  = dense_matrix<T>;
    
    static mat_type eval(Integer rows,Integer cols, const T *arr)
    {
        return dense_matrix<T>(matcl::make_integer_dense(rows, cols, arr));
    };

    static mat_type eval(Integer rows,Integer cols, const T *arr, Integer ld)
    {
        return dense_matrix<T>(matcl::make_integer_dense(rows, cols, arr, ld));
    };

    static mat_type eval_foreign(Integer rows,Integer cols, T *arr, Integer ld)
    {
        return dense_matrix<T>(matcl::make_dense_foreign(rows, cols, arr, ld));
    };

    static mat_type eval_noinit(Integer rows,Integer cols, T*& arr)
    {
        return dense_matrix<T>(matcl::make_integer_dense_noinit(rows, cols, arr));
    };

    static mat_type eval(const T& val, Integer rows,Integer cols)
    {
        return dense_matrix<T>(matcl::make_integer_dense(val, rows, cols));
    };

    static mat_type eval_range(T s, T e)
    {
        return dense_matrix<T>(matcl::irange(s,e));
    };

    static mat_type eval_range(T s, T i, T e)
    {
        return dense_matrix<T>(matcl::irange(s,i, e));
    };
};

template<>
struct make_helper<Real>
{
    using T         = Real;
    using mat_type  = dense_matrix<T>;
    
    static mat_type eval(Integer rows,Integer cols, const T *arr)
    {
        return dense_matrix<T>(matcl::make_real_dense(rows, cols, arr));
    };

    static mat_type eval(Integer rows,Integer cols, const T *arr, Integer ld)
    {
        return dense_matrix<T>(matcl::make_real_dense(rows, cols, arr, ld));
    };

    static mat_type eval_foreign(Integer rows,Integer cols, T *arr, Integer ld)
    {
        return dense_matrix<T>(matcl::make_dense_foreign(rows, cols, arr, ld));
    };

    static mat_type eval_noinit(Integer rows,Integer cols, T*& arr)
    {
        return dense_matrix<T>(matcl::make_real_dense_noinit(rows, cols, arr));
    };

    static mat_type eval(const T& val, Integer rows,Integer cols)
    {
        return dense_matrix<T>(matcl::make_real_dense(val, rows, cols));
    };

    static mat_type eval_randn(Integer r, Integer c, const rand_state& rand_ptr)
    {
        return dense_matrix<T>(matcl::randn(r,c,rand_ptr));
    };

    static mat_type eval_range(T s, T e)
    {
        return dense_matrix<T>(matcl::range(s,e));
    };

    static mat_type eval_range(T s, T i, T e)
    {
        return dense_matrix<T>(matcl::range(s,i, e));
    };

    static mat_type eval_linspace(T s, T e, Integer n)
    {
        return dense_matrix<T>(matcl::linspace(s,e,n));
    };

    static mat_type eval_logspace(T s, T e, Integer n)
    {
        return dense_matrix<T>(matcl::logspace(s,e,n));
    };
};

template<>
struct make_helper<Float>
{
    using T         = Float;
    using mat_type  = dense_matrix<T>;
    
    static mat_type eval(Integer rows,Integer cols, const T *arr)
    {
        return dense_matrix<T>(matcl::make_float_dense(rows, cols, arr));
    };

    static mat_type eval(Integer rows,Integer cols, const T *arr, Integer ld)
    {
        return dense_matrix<T>(matcl::make_float_dense(rows, cols, arr, ld));
    };

    static mat_type eval_foreign(Integer rows,Integer cols, T *arr, Integer ld)
    {
        return dense_matrix<T>(matcl::make_dense_foreign(rows, cols, arr, ld));
    };

    static mat_type eval_noinit(Integer rows,Integer cols, T*& arr)
    {
        return dense_matrix<T>(matcl::make_float_dense_noinit(rows, cols, arr));
    };

    static mat_type eval(const T& val, Integer rows,Integer cols)
    {
        return dense_matrix<T>(matcl::make_float_dense(val, rows, cols));
    };

    static mat_type eval_randn(Integer r, Integer c, const rand_state& rand_ptr)
    {
        return dense_matrix<T>(matcl::frandn(r,c,rand_ptr));
    };

    static mat_type eval_range(T s, T e)
    {
        return dense_matrix<T>(matcl::frange(s,e));
    };

    static mat_type eval_range(T s, T i, T e)
    {
        return dense_matrix<T>(matcl::frange(s,i, e));
    };

    static mat_type eval_linspace(T s, T e, Integer n)
    {
        return dense_matrix<T>(matcl::flinspace(s,e,n));
    };

    static mat_type eval_logspace(T s, T e, Integer n)
    {
        return dense_matrix<T>(matcl::flogspace(s,e,n));
    };
};

template<>
struct make_helper<Complex>
{
    using T         = Complex;
    using mat_type  = dense_matrix<T>;
    
    static mat_type eval(Integer rows,Integer cols, const T *arr)
    {
        return dense_matrix<T>(matcl::make_complex_dense(rows, cols, arr));
    };

    static mat_type eval(Integer rows,Integer cols, const T *arr, Integer ld)
    {
        return dense_matrix<T>(matcl::make_complex_dense(rows, cols, arr, ld));
    };

    static mat_type eval_foreign(Integer rows,Integer cols, T *arr, Integer ld)
    {
        return dense_matrix<T>(matcl::make_dense_foreign(rows, cols, arr, ld));
    };

    static mat_type eval_noinit(Integer rows,Integer cols, T*& arr)
    {
        return dense_matrix<T>(matcl::make_complex_dense_noinit(rows, cols, arr));
    };

    static mat_type eval(const T& val, Integer rows,Integer cols)
    {
        return dense_matrix<T>(matcl::make_complex_dense(val, rows, cols));
    };

    static mat_type eval_randn(Integer r, Integer c, const rand_state& rand_ptr)
    {
        return dense_matrix<T>(matcl::crandn(r,c,rand_ptr));
    };

    static mat_type eval_complex(Integer r, Integer c, const Real* re, const Real* im)
    {
        return dense_matrix<T>(matcl::make_complex_dense(r,c,re,im));
    };
};

template<>
struct make_helper<Float_complex>
{
    using T         = Float_complex;
    using mat_type  = dense_matrix<T>;
    
    static mat_type eval(Integer rows,Integer cols, const T *arr)
    {
        return dense_matrix<T>(matcl::make_float_complex_dense(rows, cols, arr));
    };

    static mat_type eval(Integer rows,Integer cols, const T *arr, Integer ld)
    {
        return dense_matrix<T>(matcl::make_float_complex_dense(rows, cols, arr, ld));
    };

    static mat_type eval_foreign(Integer rows,Integer cols, T *arr, Integer ld)
    {
        return dense_matrix<T>(matcl::make_dense_foreign(rows, cols, arr, ld));
    };

    static mat_type eval(Integer rows,Integer cols, T *arr, Integer ld)
    {
        return dense_matrix<T>(matcl::make_dense_foreign(rows, cols, arr, ld));
    };

    static mat_type eval_noinit(Integer rows,Integer cols, T*& arr)
    {
        return dense_matrix<T>(matcl::make_float_complex_dense_noinit(rows, cols, arr));
    };

    static mat_type eval(const T& val, Integer rows,Integer cols)
    {
        return dense_matrix<T>(matcl::make_float_complex_dense(val, rows, cols));
    };

    static mat_type eval_randn(Integer r, Integer c, const rand_state& rand_ptr)
    {
        return dense_matrix<T>(matcl::fcrandn(r,c,rand_ptr));
    };

    static mat_type eval_complex(Integer r, Integer c, const Float* re, const Float* im)
    {
        return dense_matrix<T>(matcl::make_float_complex_dense(r,c,re,im));
    };
};

template<>
struct make_helper<Object>
{
    using T         = Object;
    using mat_type  = dense_matrix<T>;
    
    static mat_type eval(ti::ti_object ti,Integer rows,Integer cols, const T *arr)
    {
        return dense_matrix<T>(matcl::make_object_dense(ti, rows, cols, arr));
    };

    static mat_type eval(ti::ti_object ti,Integer rows,Integer cols, const T *arr, Integer ld)
    {
        return dense_matrix<T>(matcl::make_object_dense(ti, rows, cols, arr, ld));
    };

    static mat_type eval_foreign(ti::ti_object ti,Integer rows,Integer cols, T *arr, Integer ld)
    {
        return dense_matrix<T>(matcl::make_dense_foreign(ti, rows, cols, arr, ld));
    };

    static mat_type eval_noinit(ti::ti_object ti,Integer rows,Integer cols, T*& arr)
    {
        return dense_matrix<T>(matcl::make_object_dense_noinit(ti, rows, cols, arr));
    };

    static mat_type eval(const T& val, Integer rows,Integer cols)
    {
        return dense_matrix<T>(matcl::make_object_dense(val, rows, cols));
    };
};

template<class T>
dense_matrix<T> dense_matrix<T>::rand(Integer r, Integer c, const rand_state& rand_ptr)
{
    return dense_matrix<T>(matcl::rand(r,c, details::value_to_code<T>::value, rand_ptr));
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::randn(Integer r, Integer c, const rand_state& rand_ptr)
{
    return make_helper<T>::eval_randn(r,c,rand_ptr);
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::range(T s, T e)
{
    return make_helper<T>::eval_range(s,e);
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::range(T s, T i, T e)
{
    return make_helper<T>::eval_range(s,i,e);
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::linspace(typename details::real_type<T>::type s, 
                                          typename details::real_type<T>::type e, Integer n)
{
    return make_helper<T>::eval_linspace(s,e,n);
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::logspace(typename details::real_type<T>::type s, 
                                          typename details::real_type<T>::type e, Integer n)
{
    return make_helper<T>::eval_logspace(s,e,n);
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::make_complex(Integer rows,Integer cols, 
                        const typename details::real_type<T>::type *ar_r,
                        const typename details::real_type<T>::type* ar_i)
{
    return make_helper<T>::eval_complex(rows,cols,ar_r,ar_i);
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::make(Integer rows,Integer cols)
{
    return dense_matrix<T>(matcl::make_dense_matrix(rows, cols, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::make(ti::ti_object ti, Integer rows,Integer cols)
{
    return dense_matrix<T>(matcl::make_object_dense(ti, rows, cols));
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::make(Integer rows,Integer cols, const T *arr)
{
    return make_helper<T>::eval(rows, cols, arr);    
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::make(Integer rows,Integer cols, const T *arr, Integer ld)
{
    return make_helper<T>::eval(rows, cols, arr, ld);    
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::make(ti::ti_object ti, Integer rows,Integer cols, const T *arr)
{
    return make_helper<T>::eval(ti, rows, cols, arr);    
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::make(ti::ti_object ti, Integer rows,Integer cols, const T *arr, Integer ld)
{
    return make_helper<T>::eval(ti, rows, cols, arr, ld);    
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::make_foreign(Integer rows,Integer cols, T* arr, Integer ld)
{
    return make_helper<T>::eval_foreign(rows, cols, arr, ld);    
}

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::make_foreign(ti::ti_object ti, Integer rows,Integer cols, Object* arr, 
                                Integer ld)
{
    return make_helper<T>::eval_foreign(ti, rows, cols, arr, ld);    
};

template<class T>
dense_matrix<T> dense_matrix<T>::make(const T& val, Integer rows, Integer cols)
{
    return make_helper<T>::eval(val, rows, cols);    
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::make_noinit(Integer rows,Integer cols, T*& ptr_data)
{
    return make_helper<T>::eval_noinit(rows, cols, ptr_data);    
};

template<class T>
template<class Enable>
dense_matrix<T> dense_matrix<T>::make_noinit(ti::ti_object ti, Integer rows,Integer cols, T*& ptr_data)
{
    return make_helper<T>::eval_noinit(ti, rows, cols, ptr_data);    
};

template<class T>
dense_matrix<T> dense_matrix<T>::diag(const dense_matrix& v, Integer d)
{
    return dense_matrix<T>(matcl::diag(v.to_matrix(),d));
};

template<class T> 
dense_matrix<T> dense_matrix<T>::diags(const dense_matrix& A, const Matrix &d, Integer r, Integer c)
{
    return dense_matrix<T>(matcl::diags(A.to_matrix(),d,r,c));
};

template<class T>
sub_dense_matrix<T> dense_matrix<T,true>::operator()(const colon& r, const colon& c)
{
    sub_dense_matrix ret(this,r,c);
    return ret;
};

template<class T>
sub_dense_matrix<T> dense_matrix<T,true>::diag(Integer d)
{
    sub_dense_matrix ret(d,this);
    return ret;
};

template<class T>
const typename sub_dense_matrix<T>::matrix_type 
sub_dense_matrix<T>::to_matrix() const
{
    const Matrix& tmp = Matrix(*m_matrix);
    if (m_colon_2)
        return matrix_type(tmp(*m_colon_1,*m_colon_2));
    else if(m_colon_1)
        return matrix_type(tmp(*m_colon_1));
    else
        return matrix_type(get_diag(tmp,m_d));
};

template<class T>
typename sub_dense_matrix<T>::matrix_type& 
sub_dense_matrix<T>::operator=(const matrix_type& mat) const &&
{
    if (m_colon_2)
    {
        m_matrix->m_matrix(*m_colon_1,*m_colon_2) = mat.to_matrix();

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    }
    else if (m_colon_1)
    {
        m_matrix->m_matrix(*m_colon_1) = mat.to_matrix();

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    }
    else
    {
        m_matrix->m_matrix.diag(m_d) = mat.to_matrix();

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    };
}

template<class T>
typename sub_dense_matrix<T>::matrix_type& 
sub_dense_matrix<T>::operator=(const T& mat) const &&
{
    if (m_colon_2)
    {
        m_matrix->m_matrix(*m_colon_1,*m_colon_2) = mat;

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    }
    else if (m_colon_1)
    {
        m_matrix->m_matrix(*m_colon_1) = mat;

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    }
    else
    {
        m_matrix->m_matrix.diag(m_d) = mat;

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    };
}

template<class T>
typename sub_dense_matrix<T>::matrix_type& 
sub_dense_matrix<T>::operator=(const sub_dense_matrix<T>& mat0) const &&
{
    if (m_colon_2)
    {
        m_matrix->m_matrix(*m_colon_1,*m_colon_2) = mat0.m_matrix->m_matrix;

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    }
    else if (m_colon_1)
    {
        m_matrix->m_matrix(*m_colon_1) = mat0.m_matrix->m_matrix;

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    }
    else
    {
        m_matrix->m_matrix.diag(m_d) = mat0.m_matrix->m_matrix;

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    };
};

template<class T>
template<class V> 
inline V sub_dense_matrix<T>::get_scalar() const
{
    return Matrix(*this).get_scalar<V>();
};

template<class T>
inline Integer sub_dense_matrix<T>::rows() const
{
    return Matrix(*this).rows();
};

template<class T>
inline Integer sub_dense_matrix<T>::cols() const
{
    return Matrix(*this).cols();
};

template<class T>
inline Integer sub_dense_matrix<T>::length() const
{
    return Matrix(*this).length();
};

template<class T>
inline Real sub_dense_matrix<T>::numel() const
{
    return Matrix(*this).numel();
};

template<class T>
inline bool sub_dense_matrix<T>::all_finite() const
{
    return Matrix(*this).all_finite();
};

template<class T>
inline sub_dense_matrix<T>::operator bool() const
{
    return Matrix(*this).operator bool();
};

template<class T>
void test_compile(const dense_matrix<T>& v)
{
    disp(v);

    {
        std::ostringstream os;
        std::istringstream is;

        os << v;

        dense_matrix<T> tmp;
        is >> tmp;

        oarchive oa(os);
        save(oa, v);

        iarchive ia(is);
        load(ia, tmp);
    };

    using ret   = dense_matrix<T>;
    using ret2  = dense_matrix<typename md::unify_types<T,Float>::type>;
    using ret3  = dense_matrix<typename md::real_type_int_real<T>::type>;
    using ret4  = dense_matrix<typename md::real_type<T>::type>;
    using vect  = std::vector<dense_matrix<T>>;

    ret tmp;
    ret2 tmp2;
    ret3 tmp3;
    ret4 tmp4;

    tmp = delrows(v, colon());
    tmp = delcols(v,colon());
    tmp = delrowscols(v,colon(), colon());
    tmp = horzcat(v,v);
    tmp = vertcat(v,v);
    tmp = blkdiag(v,v);
    tmp = horzcat({v,v});
    tmp = vertcat({v,v});
    tmp = blkdiag({v,v});

    tmp = horzcat(vect{v,v});
    tmp = vertcat(vect{v,v});
    tmp = blkdiag(vect{v,v});

    tmp = repmat(v,1,1);

    tmp = full(v);

    {
        sparse_matrix<T> tmp_mat;
        tmp = full(tmp_mat);
    };
    {
        band_matrix<T> tmp_mat;
        tmp = full(tmp_mat);
    };

    tmp = clone(v);
    tmp = trans(v);
    tmp = ctrans(v);
    tmp = trans(v,trans_type::trans);
    tmp = trans(v,trans_type_ext::trans);
    tmp = vec(v);
    tmp = tril(v);
    tmp = triu(v);
    tmp = flipud(v);
    tmp = fliplr(v);
    tmp = reshape(v,1,1);
    tmp = get_diag(v);
    tmp = rot90(v);
    find(v);
    find(v, *(test_function*)nullptr);
    find(v, *(test_type_function<T>*)nullptr);
    find2(v);
    find2(v,*(test_function*)nullptr);
    find2(v, *(test_type_function<T>*)nullptr);
    find3(v);
    find3(v,*(test_function*)nullptr);
    find3(v, *(test_type_function<T>*)nullptr);
    tmp = sort(v);
    sort2(v);
    tmp = sortrows(v);
    sortrows2(v);
    tmp = sortrows(v,Matrix());
    sortrows2(v,Matrix());
    tmp = sortcols(v);
    sortcols2(v);
    tmp = sortcols(v,Matrix());
    sortcols2(v,Matrix());
    issorted(v);
    tmp = drop_sparse(v,0);
    sparse_matrix<T> ret_sp = convert<sparse_matrix<T>>(v);
    convert_value<Real>(v);
    tmp = convert_object(v,v.get_type());

    nnz(v, 1);
    all(v);
    all(v,*(test_function*)nullptr);
    all(v, *(test_type_function<T>*)nullptr);
    any(v);
    any(v,*(test_function*)nullptr);
    any(v, *(test_type_function<T>*)nullptr);
    tmp = sum(v);
    tmp = prod(v);
    tmp = cumsum(v);
    tmp = cumprod(v);
    tmp = min_d(v);
    tmp = max_d(v);
    min2(v);
    max2(v);
    tmp = min(v);
    tmp = max(v);

    tmp2 = mean(v);
    tmp3 = std(v);
    tmp4 = min_abs_d(v);
    tmp4 = max_abs_d(v);
    tmp4 = min_abs2(v).get<1>();
    tmp4 = max_abs2(v).get<1>();
    tmp4 = min_abs(v);
    tmp4 = max_abs(v);

    nnz_vec(v);
    all_vec(v);
    all_vec(v,*(test_function*)nullptr);
    all_vec(v,*(test_type_function<T>*)nullptr);
    any_vec(v);
    any_vec(v,*(test_function*)nullptr);
    any_vec(v, *(test_type_function<T>*)nullptr);
    sum_vec(v);
    prod_vec(v);
    min_vec(v);
    max_vec(v);
    min2_vec(v);
    max2_vec(v);

    using scal2  = typename md::unify_types<T,Float>::type;
    using scal3  = typename md::real_type_int_real<T>::type;
    using scal4  = typename md::real_type<T>::type;

    scal2 sc2 = mean_vec(v);
    scal3 sc3 = std_vec(v);
    scal4 sc4 = min_abs_vec(v);
    sc4 = max_abs_vec(v);
    sc4 = min_abs2_vec(v).get<1>();
    sc4 = max_abs2_vec(v).get<1>();

    (void)sc2;
    (void)sc3;
};

template class MATCL_MATREP_EXPORT_REP matcl::dense_matrix<Integer,true>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_matrix<Real,true>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_matrix<Float,true>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_matrix<Complex,true>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_matrix<Float_complex,true>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_matrix<Object,true>;

template class MATCL_MATREP_EXPORT_REP matcl::dense_matrix<Integer,false>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_matrix<Real,false>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_matrix<Float,false>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_matrix<Complex,false>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_matrix<Float_complex,false>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_matrix<Object,false>;

template class MATCL_MATREP_EXPORT_REP matcl::sub_dense_matrix<Integer>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_dense_matrix<Real>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_dense_matrix<Float>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_dense_matrix<Complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_dense_matrix<Float_complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_dense_matrix<Object>;

template class MATCL_MATREP_EXPORT_REP matcl::dense_row<Integer>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_row<Real>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_row<Float>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_row<Complex>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_row<Float_complex>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_row<Object>;

template class MATCL_MATREP_EXPORT_REP matcl::dense_col<Integer>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_col<Real>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_col<Float>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_col<Complex>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_col<Float_complex>;
template class MATCL_MATREP_EXPORT_REP matcl::dense_col<Object>;

template MATCL_MATREP_EXPORT_REP
dense_matrix<Integer> dense_matrix<Integer>::zeros(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::zeros(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::zeros(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Complex> dense_matrix<Complex>::zeros(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float_complex> dense_matrix<Float_complex>::zeros(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Object> dense_matrix<Object>::zeros(ti::ti_object, Integer r, Integer c);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Integer> dense_matrix<Integer>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Complex> dense_matrix<Complex>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float_complex> dense_matrix<Float_complex>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Object> dense_matrix<Object>::ones(ti::ti_object, Integer r, Integer c);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Integer> dense_matrix<Integer>::eye(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::eye(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::eye(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Complex> dense_matrix<Complex>::eye(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float_complex> dense_matrix<Float_complex>::eye(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Object> dense_matrix<Object>::eye(ti::ti_object, Integer r, Integer c);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Integer> dense_matrix<Integer>::eye(Integer r);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::eye(Integer r);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::eye(Integer r);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Complex> dense_matrix<Complex>::eye(Integer r);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float_complex> dense_matrix<Float_complex>::eye(Integer r);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Object> dense_matrix<Object>::eye(ti::ti_object, Integer r);

#pragma warning(push)
#pragma warning(disable:5037) //an out-of-line definition of a member of a class template
                              // cannot have default arguments

// is it a VS bug?

template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::randn(Integer r, Integer c, const rand_state& rand_ptr);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::randn(Integer r, Integer c, const rand_state& rand_ptr);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Float_complex> dense_matrix<Float_complex>::randn(Integer r, Integer c, const rand_state& rand_ptr);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Complex> dense_matrix<Complex>::randn(Integer r, Integer c, const rand_state& rand_ptr);

#pragma warning(pop)

template MATCL_MATREP_EXPORT_REP
dense_matrix<Integer> dense_matrix<Integer>::range(Integer s, Integer e);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Integer> dense_matrix<Integer>::range(Integer s, Integer i, Integer e);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::range(Real s, Real e);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::range(Real s, Real i, Real e);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::range(Float s, Float e);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::range(Float s, Float i, Float e);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::linspace(Real s, Real e, Integer n);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::linspace(Float s, Float e, Integer n);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::logspace(Real s, Real e, Integer n);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::logspace(Float s, Float e, Integer n);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Integer> dense_matrix<Integer>::make(Integer rows, Integer cols, const Integer *arr);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Integer> dense_matrix<Integer>::make(Integer rows, Integer cols, const Integer *arr, Integer ld);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Integer> dense_matrix<Integer>::make_foreign(Integer rows, Integer cols, Integer *arr, Integer ld);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::make(Integer rows, Integer cols, const Real *arr);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::make(Integer rows, Integer cols, const Real *arr, Integer ld);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::make_foreign(Integer rows, Integer cols, Real *arr, Integer ld);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::make(Integer rows, Integer cols, const Float *arr);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::make(Integer rows, Integer cols, const Float *arr, Integer ld);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::make_foreign(Integer rows, Integer cols, Float *arr, Integer ld);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Complex> dense_matrix<Complex>::make(Integer rows, Integer cols, const Complex *arr);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Complex> dense_matrix<Complex>::make(Integer rows, Integer cols, const Complex *arr, Integer ld);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Complex> dense_matrix<Complex>::make_foreign(Integer rows, Integer cols, Complex *arr, Integer ld);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Float_complex> dense_matrix<Float_complex>::make(Integer rows, Integer cols, const Float_complex *arr);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float_complex> dense_matrix<Float_complex>::make(Integer rows, Integer cols, const Float_complex *arr,
                                                              Integer ld);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float_complex> dense_matrix<Float_complex>::make_foreign(Integer rows, Integer cols, Float_complex *arr,
                                                              Integer ld);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Object> dense_matrix<Object>::make(ti::ti_object ti, Integer rows, Integer cols, const Object *arr);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Object> dense_matrix<Object>::make(ti::ti_object ti, Integer rows, Integer cols, const Object *arr,
                                                Integer ld);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Object> dense_matrix<Object>::make_foreign(ti::ti_object ti, Integer rows, Integer cols, Object *arr,
                                                Integer ld);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Integer> dense_matrix<Integer>::make(Integer rows, Integer cols);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::make(Integer rows, Integer cols);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::make(Integer rows, Integer cols);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Complex> dense_matrix<Complex>::make(Integer rows, Integer cols);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float_complex> dense_matrix<Float_complex>::make(Integer rows, Integer cols);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Object> dense_matrix<Object>::make(ti::ti_object ti, Integer rows, Integer cols);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Integer> dense_matrix<Integer>::make_noinit(Integer rows, Integer cols, Integer *& arr);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Real> dense_matrix<Real>::make_noinit(Integer rows, Integer cols, Real *& arr);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float> dense_matrix<Float>::make_noinit(Integer rows, Integer cols, Float *& arr);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Complex> dense_matrix<Complex>::make_noinit(Integer rows, Integer cols, Complex *& arr);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Float_complex> dense_matrix<Float_complex>::make_noinit(Integer rows, Integer cols,
                                                                     Float_complex*& arr);
template MATCL_MATREP_EXPORT_REP
dense_matrix<Object> dense_matrix<Object>::make_noinit(ti::ti_object ti, Integer rows, Integer cols,
                                                       Object*& arr);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Complex> dense_matrix<Complex>::make_complex(Integer rows, Integer cols, const Real* re, 
                                                          const Real* im);

template MATCL_MATREP_EXPORT_REP
dense_matrix<Float_complex> dense_matrix<Float_complex>::make_complex(Integer rows, Integer cols, 
                                                                      const Float* re, const Float* im);

template void test_compile(const dense_matrix<Integer>& v);
template void test_compile(const dense_matrix<Real>& v);
template void test_compile(const dense_matrix<Float>& v);
template void test_compile(const dense_matrix<Complex>& v);
template void test_compile(const dense_matrix<Float_complex>& v);
template void test_compile(const dense_matrix<Object>& v);

};
