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

#include "matcl-matrep/matrix/matrix_rep_band.h"
#include "matcl-matrep/matrix/matrix_rep_dense.h"
#include "matcl-matrep/matrix/matrix_rep_sparse.h"
#include "matcl-matrep/matrix/matrix_rep_functions.h"
#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/details/matrix.inl"
#include "matcl-matrep/details/matrix_details_subs.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-core/details/integer.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-internals/base/utils.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"

#define MATCL_MATREP_EXPORT_REP MATCL_MATREP_EXPORT

namespace matcl
{

template<class Mat>
struct convert_to_band
{
    static Matrix eval(const Matrix& m)
    {
        return convert(m, Mat::matrix_code);
    };

    static Matrix eval(Matrix&& m)
    {
        Matrix loc(std::move(m));
        return convert(loc, Mat::matrix_code);
    };
};

//--------------------------------------------------------------------
//              band_matrix<T, true>
//--------------------------------------------------------------------
template<class T>
band_matrix<T,true>::band_matrix()
    :m_matrix(convert_to_band<mat_type>::eval(T(0)))
{
    const mat_type& mat= m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};

template<class T>
band_matrix<T,true>::band_matrix(bool val)
    :m_matrix(convert_to_band<mat_type>::eval(T(val)))
{
    const mat_type& mat= m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};

template<class T>
band_matrix<T,true>::band_matrix(const T& val)
    :m_matrix(convert_to_band<mat_type>::eval(val))
{
    const mat_type& mat= m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};

template<class T>
band_matrix<T,true>::band_matrix(const matcl::Matrix& m)    
    : m_matrix(convert_to_band<mat_type>::eval(m))
{
    const mat_type& mat = m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};

template<class T>
band_matrix<T,true>::band_matrix(matcl::Matrix&& m0)
    : m_matrix(convert_to_band<mat_type>::eval(std::move(m0)))
{
    const mat_type& mat = m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};

template<class T>
band_matrix<T,true>::band_matrix(const sub_band_matrix& sub)
    :band_matrix(sub.to_matrix())
{};

template<class T>
band_matrix<T,true>::band_matrix(const sub_band_matrix_1& sub)
    :band_matrix(sub.to_matrix())
{};

template<class T>
band_matrix<T,true>::band_matrix(const sub_band_matrix_2& sub)
    :band_matrix(sub.to_matrix())
{};


template<class T>
band_matrix<T,true>::band_matrix(const band_matrix& mat)
    :m_matrix(mat.m_matrix)
{
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
band_matrix<T,true>::band_matrix(band_matrix&& mat)
    :m_matrix(std::move(mat.m_matrix))
{
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
band_matrix<T,true>::band_matrix(const dense_matrix& mat)
    :band_matrix(mat.to_matrix())
{};

template<class T>
band_matrix<T,true>::band_matrix(dense_matrix&& mat)
    :band_matrix(std::move(mat.to_matrix()))
{};

template<class T>
band_matrix<T,true>::band_matrix(const sparse_matrix& mat)
    :band_matrix(mat.to_matrix())
{};

template<class T>
band_matrix<T,true>::band_matrix(sparse_matrix&& mat)
    :band_matrix(std::move(mat.to_matrix()))
{};

template<class T>
band_matrix<T,true>& band_matrix<T,true>::operator=(const band_matrix& mat) &
{
    m_matrix = mat.m_matrix;
    init_from_rep(m_matrix.get_impl<mat_type>());
    return *this;
};

template<class T>
band_matrix<T,true>& band_matrix<T,true>::operator=(band_matrix&& mat) &
{
    m_matrix = std::move(mat.m_matrix);
    init_from_rep(m_matrix.get_impl<mat_type>());
    return *this;
};

template<class T>
band_matrix<T,true>::~band_matrix()
{};

template<class T>
Integer band_matrix<T,true>::size() const
{ 
    return imult_c(m_rows, m_cols); 
};

template<class T>
void band_matrix<T,true>::set_struct(const struct_flag& sc) const
{
    m_flag->set(sc);
};

template<class T>
void band_matrix<T,true>::add_struct(const struct_flag& sc) const
{
    m_flag->add(sc);
};

template<class T>
ti::ti_object band_matrix<T,true>::get_type() const
{
    return m_matrix.get_type();
};

template<class T>
const sparse_matrix<T>
band_matrix<T,true>::delrows(const colon& c) const &
{
    return sparse_matrix(m_matrix.delrows(c));
}

template<class T>
const sparse_matrix<T>
band_matrix<T,true>::delrows(const colon& c) const &&
{
    return sparse_matrix(std::move(m_matrix).delrows(c));
}

template<class T>
const sparse_matrix<T>
band_matrix<T,true>::delcols(const colon& c) const &
{
    return sparse_matrix(m_matrix.delcols(c));
};

template<class T>
const sparse_matrix<T>
band_matrix<T,true>::delcols(const colon& c) const &&
{
    return sparse_matrix(std::move(m_matrix).delcols(c));
};

template<class T>
const sparse_matrix<T>
band_matrix<T,true>::delrowscols(const colon& c1, const colon& c2) const &
{
    return sparse_matrix(m_matrix.delrowscols(c1, c2));
};

template<class T>
const sparse_matrix<T>
band_matrix<T,true>::delrowscols(const colon& c1, const colon& c2) const &&
{
    return sparse_matrix(std::move(m_matrix).delrowscols(c1, c2));
};

template<class T>
const typename band_matrix<T,true>::band_matrix
band_matrix<T,true>::operator()(const colon& r) const
{
    return band_matrix(m_matrix(r));
};

template<class T>
const typename band_matrix<T,true>::band_matrix
band_matrix<T,true>::operator()(const colon& r, const colon& c) const
{
    return band_matrix(m_matrix(r,c));
};

template<class T>
T band_matrix<T, true>::operator()(Integer p) const
{
    Integer rows= this->rows();
    Integer r   = (p-1) % rows + 1;
    Integer c   = (p-1) / rows + 1;

    return operator()(r,c);
};

template<class T>
typename band_matrix<T, true>::sub_band_matrix_1
band_matrix<T, true>::operator()(Integer p)
{
    sub_band_matrix_1 ret(this,p);
    return ret;
};

template<class T>
struct get_ti_type
{
    static ti::ti_empty eval(const matcl::Matrix&)      { return ti::ti_empty(); };
};

template<>
struct get_ti_type<Object>
{
    static ti::ti_object eval(const matcl::Matrix& m)   { return m.get_type(); };
};

template<class T>
T band_matrix<T, true>::operator()(Integer r, Integer c) const
{
    if (r < 1 || c < 1 || r > m_rows || c > m_cols)
        throw error::invalid_index_band(r, c, m_rows, m_cols, m_ldiags, m_udiags);

    if (std::max(1, c - m_udiags) > r || r > std::min(m_rows, c + m_ldiags))
       return matcl::details::default_value<T>(get_ti_type<T>::eval(m_matrix));
    else
        return m_base_ptr[m_udiags + r - c + (c-1) * m_base_ld];
};

template<class T>
typename band_matrix<T, true>::sub_band_matrix_2
band_matrix<T, true>::operator()(Integer r, Integer c)
{
    sub_band_matrix_2 ret(this, r, c);
    return ret;
};

template<class T>
const typename band_matrix<T,true>::dense_matrix
band_matrix<T,true>::diag(Integer d) const
{
    return dense_matrix(m_matrix.diag(d));
};

template<class T>
Integer band_matrix<T,true>::length() const
{
    Integer r = rows();
    Integer c = cols();

    if (r == 0 || c == 0)
        return 0;

    return (r > c) ? r : c;
};

template<class T>
Integer band_matrix<T,true>::structural_nnz() const
{
    return m_matrix.structural_nnz();
};

template<class T>
Integer band_matrix<T,true>::structural_ldiags(bool use_flags) const
{
    return m_matrix.structural_ldiags(use_flags);
};

template<class T>
Integer band_matrix<T,true>::structural_udiags(bool use_flags) const
{
    return m_matrix.structural_udiags(use_flags);
};

template<class T>
Real band_matrix<T,true>::numel() const
{
    return Real(rows()) * Real(cols());
};

template<class T>
bool band_matrix<T,true>::all_finite() const
{
    return m_matrix.all_finite();
};

template<class T>
const typename band_matrix<T,true>::band_matrix
band_matrix<T,true>::clone() const
{
    return band_matrix(m_matrix.clone());
};

template<class T>
typename band_matrix<T,true>::band_matrix&
band_matrix<T,true>::make_unique()
{
    m_matrix.make_unique();
    const mat_type& mat = m_matrix.get_impl<mat_type>();
    init_from_rep(mat);

    return *this;
};

template<class T>
void band_matrix<T,true>::resize(Integer r, Integer c)
{
    m_matrix.resize(r,c);
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
void band_matrix<T,true>::reserve(Integer r, Integer c)
{
    m_matrix.reserve(r,c);
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
void band_matrix<T,true>::resize_band(Integer r, Integer c, Integer fd, Integer ld)
{
    m_matrix.resize_band(r,c,fd,ld);
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
void band_matrix<T,true>::reserve_band(Integer r, Integer c, Integer fd, Integer ld)
{
    m_matrix.reserve_band(r,c,fd,ld);
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
sub_band_matrix<T> band_matrix<T,true>::operator()(const colon& r)
{
    sub_band_matrix ret(this, r);
    return ret;
};

template<class T>
sub_band_matrix<T> band_matrix<T,true>::operator()(const colon& r, const colon& c)
{
    sub_band_matrix ret(this, r, c);
    return ret;
};

template<class T>
sub_band_matrix<T> band_matrix<T,true>::diag(Integer d)
{
    sub_band_matrix ret(d, this);
    return ret;
};

template<class T>
void band_matrix<T,true>::update_rep()
{
    m_matrix   = convert(m_matrix, mat_type::matrix_code);
    const mat_type& mat = m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};

template<class T>
void band_matrix<T,true>::init_from_rep(const mat_type& mat)
{
    m_ldiags        = -mat.first_diag();
    m_udiags        = mat.last_diag();
    m_rows          = mat.rows();
    m_cols          = mat.cols();

    m_base_ld       = mat.ld();
    m_base_size     = mat.impl_size();
    m_base_ptr      = const_cast<value_type*>(mat.rep_ptr());
    m_flag          = const_cast<struct_flag*>(&mat.get_struct());
};

template<class T>
band_matrix<T,true>::band_matrix(matcl::Matrix& m, str_make_unique)
{
    if (m.is_unique() == true)
    {
        m_matrix           = convert_to_band<mat_type>::eval(m);
        const mat_type& mat = m_matrix.get_impl<mat_type>();
        init_from_rep(mat);
    }
    else 
    {
        //TODO: convert m to band?
        m_matrix           = convert_to_band<mat_type>::eval(m).make_unique();
        const mat_type& mat = m_matrix.get_impl<mat_type>();

        init_from_rep(mat);
    }    
};

template<class T>
band_matrix<T,true>::band_matrix(matcl::Matrix&& m, str_make_unique)
{
    m_matrix           = convert_to_band<mat_type>::eval(std::move(m)).make_unique();
    const mat_type& mat = m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};

//--------------------------------------------------------------------
//              band_matrix<T, true> static functions
//--------------------------------------------------------------------
template<class T>
template<class Enable>
band_matrix<T,true> band_matrix<T,true>::zeros(Integer r, Integer c,Integer fd, Integer ld)
{
    return band_matrix(bzeros(r,c,fd,ld, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
band_matrix<T,true> band_matrix<T,true>::zeros(ti::ti_object ti, Integer r, Integer c,Integer fd, Integer ld)
{
    return band_matrix(bzeros(ti, r,c,fd,ld));
};

template<class T>
template<class Enable>
band_matrix<T,true> band_matrix<T,true>::ones(Integer r, Integer c)
{
    return band_matrix(bones(r,c,details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
band_matrix<T,true> band_matrix<T,true>::ones(ti::ti_object ti, Integer r, Integer c)
{
    return band_matrix(bones(ti, r,c));
};

template<class T>
template<class Enable>
band_matrix<T,true> band_matrix<T,true>::eye(Integer r, Integer c, Integer fd, Integer ld)
{
    return band_matrix(beye(r,c,fd,ld, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
band_matrix<T,true> band_matrix<T,true>::eye(Integer r, Integer fd, Integer ld)
{
    return band_matrix(beye(r,r,fd,ld, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
band_matrix<T,true> band_matrix<T,true>::eye(ti::ti_object ti, Integer r, Integer c, Integer fd, Integer ld)
{
    return band_matrix(beye(ti, r,c,fd,ld));
};

template<class T>
template<class Enable>
band_matrix<T,true> band_matrix<T,true>::eye(ti::ti_object ti, Integer r, Integer fd, Integer ld)
{
    return band_matrix(beye(ti, r,r,fd,ld));
};

template<class T>
band_matrix<T,true> band_matrix<T,true>::diag(const dense_matrix& v, Integer d)
{
    return band_matrix(bdiag(Matrix(v),d));
};

template<class T>
band_matrix<T,true> band_matrix<T,true>::diags(const dense_matrix& A, const Matrix &d, Integer r, Integer c)
{
    return band_matrix(bdiags(Matrix(A),d,r,c));
};

template<class T>
band_matrix<T> band_matrix<T>::rand(Integer r, Integer c, Integer fd, Integer ld, const rand_state& rand_ptr)
{
    return band_matrix<T>(matcl::rand_band(r,c, fd, ld, details::value_to_code<T>::value, rand_ptr));
};

template<class T>
template<class Enable>
band_matrix<T> band_matrix<T>::randn(Integer r, Integer c, Integer fd, Integer ld, const rand_state& rand_ptr)
{
    return band_matrix<T>(matcl::randn_band(r,c, fd, ld, details::value_to_code<T>::value, rand_ptr));
};

template<class T>
struct make_helper{};

template<>
struct make_helper<Integer>
{
    using T         = Integer;
    using mat_type  = band_matrix<T>;

    static mat_type eval(const T& val,Integer rows,Integer cols, Integer fd, Integer ld)
    {
        return band_matrix<T>(matcl::make_integer_band(val, rows, cols, fd, ld));
    };
};

template<>
struct make_helper<Real>
{
    using T         = Real;
    using mat_type  = band_matrix<T>;

    static mat_type eval(const T& val,Integer rows,Integer cols, Integer fd, Integer ld)
    {
        return band_matrix<T>(matcl::make_real_band(val, rows, cols, fd, ld));
    };
};

template<>
struct make_helper<Float>
{
    using T         = Float;
    using mat_type  = band_matrix<T>;

    static mat_type eval(const T& val,Integer rows,Integer cols, Integer fd, Integer ld)
    {
        return band_matrix<T>(matcl::make_float_band(val, rows, cols, fd, ld));
    };
};

template<>
struct make_helper<Complex>
{
    using T         = Complex;
    using mat_type  = band_matrix<T>;

    static mat_type eval(const T& val,Integer rows,Integer cols, Integer fd, Integer ld)
    {
        return band_matrix<T>(matcl::make_complex_band(val, rows, cols, fd, ld));
    };
};

template<>
struct make_helper<Float_complex>
{
    using T         = Float_complex;
    using mat_type  = band_matrix<T>;
    
    static mat_type eval(const T& val,Integer rows,Integer cols, Integer fd, Integer ld)
    {
        return band_matrix<T>(matcl::make_float_complex_band(val, rows, cols, fd, ld));
    };
};

template<>
struct make_helper<Object>
{
    using T         = Object;
    using mat_type  = band_matrix<T>;
    
    static mat_type eval(const T& val,Integer rows,Integer cols, Integer fd, Integer ld)
    {
        return band_matrix<T>(matcl::make_object_band(val, rows, cols, fd, ld));
    };
};

template<class T>
template<class Enable>
band_matrix<T> band_matrix<T>::make(Integer rows,Integer cols, Integer fd, Integer ld)
{
    return band_matrix<T>(matcl::make_band_matrix(rows, cols, fd, ld, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
band_matrix<T> band_matrix<T>::make(ti::ti_object ti, Integer rows,Integer cols, Integer fd, Integer ld)
{
    return band_matrix<T>(matcl::make_object_band(ti, rows, cols, fd, ld));
};

template<class T>
band_matrix<T> band_matrix<T>::make(const T& val,Integer rows,Integer cols, Integer fd, Integer ld)
{
    return make_helper<T>::eval(val, rows, cols, fd, ld);
};

template<class T>
template<class Enable>
band_matrix<T> band_matrix<T>::make_noinit(Integer rows,Integer cols, Integer fd, Integer ld)
{
    return band_matrix<T>(matcl::make_band_noinit(rows, cols, fd, ld, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
band_matrix<T> band_matrix<T>::make_noinit(ti::ti_object ti, Integer rows,Integer cols, Integer fd, Integer ld)
{
    return band_matrix<T>(matcl::make_object_band_noinit(ti, rows, cols, fd, ld));
};

//--------------------------------------------------------------------
//              band_matrix<T, false>
//--------------------------------------------------------------------
template<class T>
band_matrix<T,false>::band_matrix()
    :base_type()
{};

template<class T>
band_matrix<T,false>::band_matrix(matcl::Matrix& m)
    :base_type(m, str_make_unique())
{};

template<class T>
band_matrix<T,false>::band_matrix(matcl::Matrix&& m)
    :base_type(std::move(m), str_make_unique())
{};

template<class T>
band_matrix<T,false>::band_matrix(band_matrix<T,false>&& mat)
    :band_matrix(std::move(mat.m_matrix))
{}

template<class T>
band_matrix<T,false>::band_matrix(band_matrix<T,true>&& mat)
    :band_matrix(std::move(mat.m_matrix))
{}

template<class T>
band_matrix<T,false>::band_matrix(band_matrix<T,true>& mat)
    :band_matrix(mat.make_unique().m_matrix)
{}

template<class T>
band_matrix<T,false>& band_matrix<T,false>::operator=(band_matrix&& mat) &
{
    m_matrix = std::move(mat.m_matrix);
    init_from_rep(m_matrix.get_impl<mat_type>());
    return *this;
};

template<class T>
band_matrix<T,false>::~band_matrix()
{};

//--------------------------------------------------------------------
//              sub_band_matrix<T>
//--------------------------------------------------------------------
template<class T>
const Matrix
sub_band_matrix<T>::to_raw_matrix() const
{
    const Matrix& tmp = Matrix(*m_matrix);
    if (m_colon_2)
        return tmp(*m_colon_1,*m_colon_2);
    else if(m_colon_1)
        return tmp(*m_colon_1);
    else
        return get_diag(tmp,m_d);
};

template<class T>
const typename sub_band_matrix<T>::matrix_type 
sub_band_matrix<T>::to_matrix() const
{
    return matrix_type(this->to_raw_matrix());
};

template<class T>
typename sub_band_matrix<T>::matrix_type& 
sub_band_matrix<T>::operator=(const matrix_type& mat) const &&
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
typename sub_band_matrix<T>::matrix_type& 
sub_band_matrix<T>::operator=(const sub_band_matrix<T>& mat0) const &&
{
    if (m_colon_2)
    {
        m_matrix->m_matrix(*m_colon_1,*m_colon_2) = mat0.to_raw_matrix();

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    }
    else if (m_colon_1)
    {
        m_matrix->m_matrix(*m_colon_1) = mat0.to_raw_matrix();

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    }
    else
    {
        m_matrix->m_matrix.diag(m_d) = mat0.to_raw_matrix();

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    };
};

template<class T>
typename sub_band_matrix<T>::matrix_type& 
sub_band_matrix<T>::operator=(const T& val) const &&
{
    if (m_colon_2)
    {
        m_matrix->m_matrix(*m_colon_1,*m_colon_2) = val;

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    }
    else if (m_colon_1)
    {
        m_matrix->m_matrix(*m_colon_1) = val;

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    }
    else
    {
        m_matrix->m_matrix.diag(m_d) = val;

        //m_matrix->m_matrix could be made unique, so we need to update rep pointers
        m_matrix->update_rep();
        return *m_matrix;
    };
};

template<class T>
typename sub_band_matrix<T>::matrix_type& 
sub_band_matrix<T>::operator=(const sub_band_matrix_1<T>& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

template<class T>
typename sub_band_matrix<T>::matrix_type& 
sub_band_matrix<T>::operator=(const sub_band_matrix_2<T>& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

template<class T>
T sub_band_matrix_1<T>::to_matrix() const
{
    const matrix_type& mat = *m_matrix;
    return mat(m_ind_1);
};

template<class T>
typename sub_band_matrix_1<T>::matrix_type& 
sub_band_matrix_1<T>::operator=(const matrix_type& mat) const &&
{
    m_matrix->m_matrix(m_ind_1) = mat.to_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_band_matrix_1<T>::matrix_type& 
sub_band_matrix_1<T>::operator=(const T& val) const &&
{
    m_matrix->m_matrix(m_ind_1) = val;

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_band_matrix_1<T>::matrix_type& 
sub_band_matrix_1<T>::operator=(const sub_band_matrix<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1) = Matrix(mat.to_matrix());

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_band_matrix_1<T>::matrix_type& 
sub_band_matrix_1<T>::operator=(const sub_band_matrix_1<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1) = mat.to_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_band_matrix_1<T>::matrix_type& 
sub_band_matrix_1<T>::operator=(const sub_band_matrix_2<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1) = mat.to_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
T sub_band_matrix_2<T>::to_matrix() const
{
    const matrix_type& mat = *m_matrix;
    return mat(m_ind_1,m_ind_2);
};

template<class T>
typename sub_band_matrix_2<T>::matrix_type& 
sub_band_matrix_2<T>::operator=(const matrix_type& mat) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = mat.to_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_band_matrix_2<T>::matrix_type& 
sub_band_matrix_2<T>::operator=(const T& val) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = val;

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_band_matrix_2<T>::matrix_type& 
sub_band_matrix_2<T>::operator=(const sub_band_matrix<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = Matrix(mat.to_matrix());

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_band_matrix_2<T>::matrix_type& 
sub_band_matrix_2<T>::operator=(const sub_band_matrix_1<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = mat.to_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_band_matrix_2<T>::matrix_type& 
sub_band_matrix_2<T>::operator=(const sub_band_matrix_2<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = mat.to_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
template<class V> 
inline V sub_band_matrix<T>::get_scalar() const
{
    return Matrix(*this).get_scalar<V>();
};

template<class T>
inline Integer sub_band_matrix<T>::rows() const
{
    return Matrix(*this).rows();
};

template<class T>
inline Integer sub_band_matrix<T>::cols() const
{
    return Matrix(*this).cols();
};

template<class T>
inline Integer sub_band_matrix<T>::length() const
{
    return Matrix(*this).length();
};

template<class T>
inline Real sub_band_matrix<T>::numel() const
{
    return Matrix(*this).numel();
};

template<class T>
inline bool sub_band_matrix<T>::all_finite() const
{
    return Matrix(*this).all_finite();
};

template<class T>
inline sub_band_matrix<T>::operator bool() const
{
    return Matrix(*this).operator bool();
};

template<class T>
template<class V> 
inline V sub_band_matrix_1<T>::get_scalar() const
{
    return Matrix(*this).get_scalar<V>();
};

template<class T>
inline Integer sub_band_matrix_1<T>::rows() const
{
    return Matrix(*this).rows();
};

template<class T>
inline Integer sub_band_matrix_1<T>::cols() const
{
    return Matrix(*this).cols();
};

template<class T>
inline Integer sub_band_matrix_1<T>::length() const
{
    return Matrix(*this).length();
};

template<class T>
inline Real sub_band_matrix_1<T>::numel() const
{
    return Matrix(*this).numel();
};

template<class T>
inline bool sub_band_matrix_1<T>::all_finite() const
{
    return Matrix(*this).all_finite();
};

template<class T>
inline sub_band_matrix_1<T>::operator bool() const
{
    return Matrix(*this).operator bool();
};

template<class T>
template<class V> 
inline V sub_band_matrix_2<T>::get_scalar() const
{
    return Matrix(*this).get_scalar<V>();
};

template<class T>
inline Integer sub_band_matrix_2<T>::rows() const
{
    return Matrix(*this).rows();
};

template<class T>
inline Integer sub_band_matrix_2<T>::cols() const
{
    return Matrix(*this).cols();
};

template<class T>
inline Integer sub_band_matrix_2<T>::length() const
{
    return Matrix(*this).length();
};

template<class T>
inline Real sub_band_matrix_2<T>::numel() const
{
    return Matrix(*this).numel();
};

template<class T>
inline bool sub_band_matrix_2<T>::all_finite() const
{
    return Matrix(*this).all_finite();
};

template<class T>
inline sub_band_matrix_2<T>::operator bool() const
{
    return Matrix(*this).operator bool();
};

//--------------------------------------------------------------------
//              TESTS
//--------------------------------------------------------------------

template<class T>
void test_compile(const band_matrix<T>& v)
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

    delrows(v, colon());
    delcols(v,colon());
    delrowscols(v,colon(), colon());
    horzcat(v,v);
    vertcat(v,v);
    blkdiag(v,v);
    horzcat({v,v});
    vertcat({v,v});
    blkdiag({v,v});
    repmat(v,1,1);

    band(v);

    {
        sparse_matrix<T> tmp;
        band(tmp);
    };
    {
        dense_matrix<T> tmp;
        band(tmp);
    };

    clone(v);
    trans(v);
    ctrans(v);
    trans(v,trans_type::trans);
    trans(v,trans_type_ext::trans);
    vec(v);
    tril(v);
    triu(v);
    flipud(v);
    fliplr(v);
    reshape(v,1,1);
    get_diag(v);
    rot90(v);
    find(v);
    find(v, *(test_function*)nullptr);
    find(v, *(test_type_function<T>*)nullptr);
    find2(v);
    find2(v,*(test_function*)nullptr);
    find2(v, *(test_type_function<T>*)nullptr);
    find3(v);
    find3(v,*(test_function*)nullptr);
    find3(v, *(test_type_function<T>*)nullptr);
    sort(v);
    sort2(v);
    sortrows(v);
    sortrows2(v);
    sortrows(v,Matrix());
    sortrows2(v,Matrix());
    sortcols(v);
    sortcols2(v);
    sortcols(v,Matrix());
    sortcols2(v,Matrix());
    issorted(v);
    drop_sparse(v,0);
    convert<dense_matrix<T>>(v);
    convert_value<Real>(v);
    convert_object(v,v.get_type());

    nnz(v, 1);
    all(v);
    all(v,*(test_function*)nullptr);
    all(v, *(test_type_function<T>*)nullptr);
    any(v);
    any(v,*(test_function*)nullptr);
    any(v, *(test_type_function<T>*)nullptr);
    sum(v);
    prod(v);
    cumsum(v);
    cumprod(v);
    min_d(v);
    max_d(v);
    min2(v);
    max2(v);
    min(v);
    max(v);

    using ret   = dense_matrix<T>;
    using ret2  = dense_matrix<typename md::unify_types<T,Float>::type>;
    using ret3  = dense_matrix<typename md::real_type_int_real<T>::type>;
    using ret4  = dense_matrix<typename md::real_type<T>::type>;

    ret tmp;
    ret2 tmp2;
    ret3 tmp3;
    ret4 tmp4;

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

};

namespace matcl
{

template class MATCL_MATREP_EXPORT_REP matcl::band_matrix<Integer,true>;
template class MATCL_MATREP_EXPORT_REP matcl::band_matrix<Real,true>;
template class MATCL_MATREP_EXPORT_REP matcl::band_matrix<Float,true>;
template class MATCL_MATREP_EXPORT_REP matcl::band_matrix<Complex,true>;
template class MATCL_MATREP_EXPORT_REP matcl::band_matrix<Float_complex,true>;
template class MATCL_MATREP_EXPORT_REP matcl::band_matrix<Object,true>;

template class MATCL_MATREP_EXPORT_REP matcl::band_matrix<Integer,false>;
template class MATCL_MATREP_EXPORT_REP matcl::band_matrix<Real,false>;
template class MATCL_MATREP_EXPORT_REP matcl::band_matrix<Float,false>;
template class MATCL_MATREP_EXPORT_REP matcl::band_matrix<Complex,false>;
template class MATCL_MATREP_EXPORT_REP matcl::band_matrix<Float_complex,false>;
template class MATCL_MATREP_EXPORT_REP matcl::band_matrix<Object,false>;

template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix<Integer>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix<Real>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix<Float>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix<Complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix<Float_complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix<Object>;

template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix_1<Integer>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix_1<Real>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix_1<Float>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix_1<Complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix_1<Float_complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix_1<Object>;

template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix_2<Integer>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix_2<Real>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix_2<Float>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix_2<Complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix_2<Float_complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_band_matrix_2<Object>;

#pragma warning(push)
#pragma warning(disable:5037) //an out-of-line definition of a member of a class template
                              // cannot have default arguments

// is it a VS bug?

template MATCL_MATREP_EXPORT_REP 
band_matrix<Real> band_matrix<Real>::randn(Integer r, Integer c, Integer fd, Integer ld,
                                        const rand_state& rand_ptr);

template MATCL_MATREP_EXPORT_REP 
band_matrix<Float> band_matrix<Float>::randn(Integer r, Integer c, Integer fd, Integer ld,
                                        const rand_state& rand_ptr);

template MATCL_MATREP_EXPORT_REP 
band_matrix<Complex> band_matrix<Complex>::randn(Integer r, Integer c, Integer fd, Integer ld,
                                        const rand_state& rand_ptr);
template MATCL_MATREP_EXPORT_REP 
band_matrix<Float_complex> band_matrix<Float_complex>::randn(Integer r, Integer c, Integer fd, Integer ld,
                                        const rand_state& rand_ptr);

#pragma warning(pop)

template MATCL_MATREP_EXPORT_REP
band_matrix<Integer> band_matrix<Integer>::make(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Real> band_matrix<Real>::make(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Float> band_matrix<Float>::make(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Complex> band_matrix<Complex>::make(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Float_complex> band_matrix<Float_complex>::make(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Object> band_matrix<Object>::make(ti::ti_object, Integer r, Integer c, Integer fd, Integer ld);

template MATCL_MATREP_EXPORT_REP
band_matrix<Integer> band_matrix<Integer>::make_noinit(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Real> band_matrix<Real>::make_noinit(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Float> band_matrix<Float>::make_noinit(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Complex> band_matrix<Complex>::make_noinit(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Float_complex> band_matrix<Float_complex>::make_noinit(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Object> band_matrix<Object>::make_noinit(ti::ti_object, Integer r, Integer c, Integer fd, Integer ld);

template MATCL_MATREP_EXPORT_REP
band_matrix<Integer> band_matrix<Integer>::zeros(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Real> band_matrix<Real>::zeros(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Float> band_matrix<Float>::zeros(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Complex> band_matrix<Complex>::zeros(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Float_complex> band_matrix<Float_complex>::zeros(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Object> band_matrix<Object>::zeros(ti::ti_object, Integer r, Integer c, Integer fd, Integer ld);

template MATCL_MATREP_EXPORT_REP
band_matrix<Integer> band_matrix<Integer>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
band_matrix<Real> band_matrix<Real>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
band_matrix<Float> band_matrix<Float>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
band_matrix<Complex> band_matrix<Complex>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
band_matrix<Float_complex> band_matrix<Float_complex>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
band_matrix<Object> band_matrix<Object>::ones(ti::ti_object, Integer r, Integer c);

template MATCL_MATREP_EXPORT_REP
band_matrix<Integer> band_matrix<Integer>::eye(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Real> band_matrix<Real>::eye(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Float> band_matrix<Float>::eye(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Complex> band_matrix<Complex>::eye(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Float_complex> band_matrix<Float_complex>::eye(Integer r, Integer c, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Object> band_matrix<Object>::eye(ti::ti_object, Integer r, Integer c, Integer fd, Integer ld);

template MATCL_MATREP_EXPORT_REP
band_matrix<Integer> band_matrix<Integer>::eye(Integer r, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Real> band_matrix<Real>::eye(Integer r, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Float> band_matrix<Float>::eye(Integer r, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Complex> band_matrix<Complex>::eye(Integer r, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Float_complex> band_matrix<Float_complex>::eye(Integer r, Integer fd, Integer ld);
template MATCL_MATREP_EXPORT_REP
band_matrix<Object> band_matrix<Object>::eye(ti::ti_object, Integer r, Integer fd, Integer ld);

template void test_compile(const band_matrix<Integer>& v);
template void test_compile(const band_matrix<Real>& v);
template void test_compile(const band_matrix<Float>& v);
template void test_compile(const band_matrix<Complex>& v);
template void test_compile(const band_matrix<Float_complex>& v);
template void test_compile(const band_matrix<Object>& v);

};
