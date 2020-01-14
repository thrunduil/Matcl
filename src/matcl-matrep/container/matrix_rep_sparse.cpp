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

#include "matcl-matrep/matrix/matrix_rep_sparse.h"
#include "matcl-matrep/matrix/matrix_rep_dense.h"
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

template<class Mat>
struct convert_to_sparse
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
//              sparse_matrix<T, true>
//--------------------------------------------------------------------
template<class T>
sparse_matrix<T,true>::sparse_matrix()
    :m_matrix(convert_to_sparse<mat_type>::eval(T(0)))
{
    const mat_type& mat= m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};

template<class T>
sparse_matrix<T,true>::sparse_matrix(bool val)
    :m_matrix(convert_to_sparse<mat_type>::eval(T(val)))
{
    const mat_type& mat= m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};

template<class T>
sparse_matrix<T,true>::sparse_matrix(const T& val)
    :m_matrix(convert_to_sparse<mat_type>::eval(val))
{
    const mat_type& mat= m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};

template<class T>
sparse_matrix<T,true>::sparse_matrix(const matcl::Matrix& m)
    :m_matrix(convert_to_sparse<mat_type>::eval(m))
{
    const mat_type& mat = m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};

template<class T>
sparse_matrix<T,true>::sparse_matrix(matcl::Matrix&& m0)
    : m_matrix(convert_to_sparse<mat_type>::eval(std::move(m0)))
{
    const mat_type& mat = m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};

template<class T>
sparse_matrix<T,true>::sparse_matrix(const sub_sparse_matrix& sub)
    :sparse_matrix(sub.to_matrix())
{};

template<class T>
sparse_matrix<T,true>::sparse_matrix(const sub_sparse_matrix_1& sub)
    :sparse_matrix(sub.to_matrix())
{};

template<class T>
sparse_matrix<T,true>::sparse_matrix(const sub_sparse_matrix_2& sub)
    :sparse_matrix(sub.to_matrix())
{};

template<class T>
sparse_matrix<T,true>::sparse_matrix(const sparse_row<T>& con)
    :sparse_matrix(con.to_matrix())
{};

template<class T>
sparse_matrix<T,true>::sparse_matrix(const sparse_col<T>& con)
    :sparse_matrix(con.to_matrix())
{};

template<class T>
sparse_matrix<T,true>::sparse_matrix(const sparse_matrix& mat)
    :m_matrix(mat.m_matrix)
{
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
sparse_matrix<T,true>::sparse_matrix(sparse_matrix&& mat)
    :m_matrix(std::move(mat.m_matrix))
{
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
sparse_matrix<T,true>::sparse_matrix(const dense_matrix& mat)
    :sparse_matrix(mat.to_matrix())
{};

template<class T>
sparse_matrix<T,true>::sparse_matrix(dense_matrix&& mat)
    :sparse_matrix(std::move(mat.to_matrix()))
{};

template<class T>
sparse_matrix<T,true>::sparse_matrix(const band_matrix<T>& mat)
    :sparse_matrix(mat.to_matrix())
{};

template<class T>
sparse_matrix<T,true>::sparse_matrix(band_matrix<T>&& mat)
    :sparse_matrix(std::move(mat.to_matrix()))
{};

template<class T>
sparse_matrix<T,true>& sparse_matrix<T,true>::operator=(const sparse_matrix& mat) &
{
    m_matrix = mat.m_matrix;
    init_from_rep(m_matrix.get_impl<mat_type>());
    return *this;
};

template<class T>
sparse_matrix<T,true>& sparse_matrix<T,true>::operator=(sparse_matrix&& mat) &
{
    m_matrix = std::move(mat.m_matrix);
    init_from_rep(m_matrix.get_impl<mat_type>());
    return *this;
};

template<class T>
sparse_matrix<T,true>::~sparse_matrix()
{};

template<class T>
void sparse_matrix<T,true>::set_struct(const struct_flag& sc) const
{
    m_flag->set(sc);
};

template<class T>
void sparse_matrix<T,true>::add_struct(const struct_flag& sc) const
{
    m_flag->add(sc);
};

template<class T>
ti::ti_object sparse_matrix<T,true>::get_type() const
{
    return m_matrix.get_type();
};

template<class T>
const typename sparse_matrix<T,true>::sparse_matrix 
sparse_matrix<T,true>::delrows(const colon& c) const &
{
    return sparse_matrix(m_matrix.delrows(c));
}

template<class T>
const typename sparse_matrix<T,true>::sparse_matrix 
sparse_matrix<T,true>::delrows(const colon& c) const &&
{
    return sparse_matrix(std::move(m_matrix).delrows(c));
}

template<class T>
const typename sparse_matrix<T,true>::sparse_matrix
sparse_matrix<T,true>::delcols(const colon& c) const &
{
    return sparse_matrix(m_matrix.delcols(c));
};

template<class T>
const typename sparse_matrix<T,true>::sparse_matrix
sparse_matrix<T,true>::delcols(const colon& c) const &&
{
    return sparse_matrix(std::move(m_matrix).delcols(c));
};

template<class T>
const typename sparse_matrix<T,true>::sparse_matrix 
sparse_matrix<T,true>::delrowscols(const colon& c1, const colon& c2) const &
{
    return sparse_matrix(m_matrix.delrowscols(c1, c2));
};

template<class T>
const typename sparse_matrix<T,true>::sparse_matrix 
sparse_matrix<T,true>::delrowscols(const colon& c1, const colon& c2) const &&
{
    return sparse_matrix(std::move(m_matrix).delrowscols(c1, c2));
};

template<class T>
const typename sparse_matrix<T,true>::sparse_matrix
sparse_matrix<T,true>::operator()(const colon& r) const
{
    return sparse_matrix(m_matrix(r));
};

template<class T>
const typename sparse_matrix<T,true>::sparse_matrix
sparse_matrix<T,true>::operator()(const colon& r, const colon& c) const
{
    return sparse_matrix(m_matrix(r,c));
};

template<class T>
T sparse_matrix<T, true>::operator()(Integer p) const
{
    Integer rows= this->rows();
    Integer r   = (p-1) % rows + 1;
    Integer c   = (p-1) / rows + 1;

    return operator()(r,c);
};

template<class T>
typename sparse_matrix<T, true>::sub_sparse_matrix_1
sparse_matrix<T, true>::operator()(Integer p)
{
    sub_sparse_matrix_1 ret(this,p);
    return ret;
};

template<class T>
T sparse_matrix<T, true>::operator()(Integer r, Integer c) const
{
    return m_matrix.get_impl<mat_type>()(r,c);
};

template<class T>
typename sparse_matrix<T, true>::sub_sparse_matrix_2
sparse_matrix<T, true>::operator()(Integer r, Integer c)
{
    sub_sparse_matrix_2 ret(this,r,c);
    return ret;
};

template<class T>
const typename sparse_matrix<T,true>::dense_matrix
sparse_matrix<T,true>::diag(Integer d) const
{
    return dense_matrix(m_matrix.diag(d));
};

template<class T>
Integer sparse_matrix<T,true>::length() const
{
    Integer r = rows();
    Integer c = cols();

    if (r == 0 || c == 0)
        return 0;

    return (r > c) ? r : c;
};

template<class T>
Integer sparse_matrix<T,true>::structural_nnz() const
{
    return m_matrix.structural_nnz();
};

template<class T>
Integer sparse_matrix<T,true>::structural_ldiags(bool use_flags) const
{
    return m_matrix.structural_ldiags(use_flags);
};

template<class T>
Integer sparse_matrix<T,true>::structural_udiags(bool use_flags) const
{
    return m_matrix.structural_udiags(use_flags);
};

template<class T>
Real sparse_matrix<T,true>::numel() const
{
    return Real(rows()) * Real(cols());
};

template<class T>
bool sparse_matrix<T,true>::all_finite() const
{
    return m_matrix.all_finite();
};

template<class T>
const typename sparse_matrix<T,true>::sparse_matrix
sparse_matrix<T,true>::clone() const
{
    return sparse_matrix(m_matrix.clone());
};

template<class T>
typename sparse_matrix<T,true>::sparse_matrix&
sparse_matrix<T,true>::make_unique()
{
    m_matrix.make_unique();
    const mat_type& mat = m_matrix.get_impl<mat_type>();
    init_from_rep(mat);

    return *this;
};

template<class T>
void sparse_matrix<T,true>::resize(Integer r, Integer c)
{
    m_matrix.resize(r,c);
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
void sparse_matrix<T,true>::reserve(Integer r, Integer c)
{
    m_matrix.reserve(r,c);
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
void sparse_matrix<T,true>::resize_band(Integer r, Integer c, Integer ld, Integer ud)
{
    (void)ld;
    (void)ud;
    m_matrix.resize(r,c);
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
void sparse_matrix<T,true>::reserve_band(Integer r, Integer c, Integer ld, Integer ud)
{
    (void)ld;
    (void)ud;
    m_matrix.reserve(r,c);
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
sub_sparse_matrix<T> sparse_matrix<T,true>::operator()(const colon& r)
{
    sub_sparse_matrix ret(this,r);
    return ret;
};

template<class T>
sub_sparse_matrix<T> sparse_matrix<T,true>::operator()(const colon& r, const colon& c)
{
    sub_sparse_matrix ret(this,r,c);
    return ret;
};

template<class T>
sub_sparse_matrix<T> sparse_matrix<T,true>::diag(Integer d)
{
    sub_sparse_matrix ret(d,this);
    return ret;
};

template<class T>
Integer sparse_matrix<T,true>::nzmax() const
{
    return m_matrix.get_impl<mat_type>().nzmax();
};

template<class T>
bool sparse_matrix<T,true>::is_sorted() const
{
    return m_matrix.get_impl<mat_type>().rep().is_sorted();
};

template<class T>
bool sparse_matrix<T,true>::is_entry(Integer r, Integer c, Integer &k) const
{
    return m_matrix.get_impl<mat_type>().rep().has_element(r,c,k);
};

template<class T>
void sparse_matrix<T,true>::add_memory(Integer s)
{
    m_matrix.get_impl_unique<mat_type>().rep().add_memory(s);
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
void sparse_matrix<T,true>::sort()
{
    m_matrix.get_impl_unique<mat_type>().rep().sort();
    init_from_rep(m_matrix.get_impl<mat_type>());
};

template<class T>
void sparse_matrix<T,true>::update_rep()
{
    m_matrix   = convert(m_matrix, mat_type::matrix_code);
    const mat_type& mat = m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};

template<class T>
void sparse_matrix<T,true>::init_from_rep(const mat_type& mat)
{
    m_rows      = mat.rows();
    m_cols      = mat.cols();
    m_max_cols  = mat.max_cols();
    m_offset    = mat.rep().offset();
    m_c         = const_cast<Integer*>(mat.rep().ptr_c());
    m_r         = const_cast<Integer**>(mat.rep().ptr_root_r());
    m_x         = const_cast<value_type**>(mat.rep().ptr_root_x());
    m_flag      = const_cast<struct_flag*>(&mat.rep().get_struct());
};

template<class T>
sparse_matrix<T,true>::sparse_matrix(matcl::Matrix& m, str_make_unique)
{
    if (m.is_unique() == true)
    {
        m_matrix           = convert_to_sparse<mat_type>::eval(m);
        const mat_type& mat = m_matrix.get_impl<mat_type>();
        init_from_rep(mat);
    }
    else 
    {
        //TODO: convert m to dense?
        m_matrix           = convert_to_sparse<mat_type>::eval(m).make_unique();
        const mat_type& mat = m_matrix.get_impl<mat_type>();

        init_from_rep(mat);
    }    
};

template<class T>
sparse_matrix<T,true>::sparse_matrix(matcl::Matrix&& m, str_make_unique)
{
    m_matrix           = convert_to_sparse<mat_type>::eval(std::move(m)).make_unique();
    const mat_type& mat = m_matrix.get_impl<mat_type>();
    init_from_rep(mat);
};
//--------------------------------------------------------------------
//              sparse_matrix<T, true> static functions
//--------------------------------------------------------------------
template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::zeros(Integer r, Integer c, Integer nz)
{
    return sparse_matrix<T>(matcl::spzeros(r,c, nz, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::zeros(ti::ti_object ti, Integer r, Integer c, Integer nz)
{
    return sparse_matrix<T>(matcl::spzeros(ti, r,c, nz));
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::ones(Integer r, Integer c)
{
    return sparse_matrix<T>(matcl::spones(r,c, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::ones(ti::ti_object ti, Integer r, Integer c)
{
    return sparse_matrix<T>(matcl::spones(ti, r,c));
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::eye(Integer r, Integer c)
{
    return sparse_matrix<T>(matcl::speye(r,c, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::eye(Integer r)
{
    return sparse_matrix<T>(matcl::speye(r,r, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::eye(ti::ti_object ti, Integer r, Integer c)
{
    return sparse_matrix<T>(matcl::speye(ti, r,c));
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::eye(ti::ti_object ti, Integer r)
{
    return sparse_matrix<T>(matcl::speye(ti, r,r));
};

template<class T>
sparse_matrix<T> sparse_matrix<T>::diag(const dense_matrix& v, Integer d)
{
    return sparse_matrix<T>(spdiag(Matrix(v),d));
};

template<class T> 
sparse_matrix<T> sparse_matrix<T>::diags(const dense_matrix& A, const Matrix &d, Integer r, Integer c)
{
    return sparse_matrix<T>(spdiags(Matrix(A),d,r,c));
};

template<class T>
struct make_helper{};

template<>
struct make_helper<Integer>
{
    using T         = Integer;
    using mat_type  = sparse_matrix<T>;

    static mat_type eval(const Integer *trip_r, const Integer *trip_c, 
                        const T* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax)
    {
        return sparse_matrix<T>(matcl::make_integer_sparse(trip_r, trip_c, trip_x, r, c, nnz, nzmax));
    };
};

template<>
struct make_helper<Real>
{
    using T         = Real;
    using mat_type  = sparse_matrix<T>;

    static mat_type eval(const Integer *trip_r, const Integer *trip_c, 
                        const T* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax)
    {
        return sparse_matrix<T>(matcl::make_real_sparse(trip_r, trip_c, trip_x, r, c, nnz, nzmax));
    };
};

template<>
struct make_helper<Float>
{
    using T         = Float;
    using mat_type  = sparse_matrix<T>;

    static mat_type eval(const Integer *trip_r, const Integer *trip_c, 
                        const T* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax)
    {
        return sparse_matrix<T>(matcl::make_float_sparse(trip_r, trip_c, trip_x, r, c, nnz, nzmax));
    };
};

template<>
struct make_helper<Complex>
{
    using T         = Complex;
    using mat_type  = sparse_matrix<T>;

    static mat_type eval(const Integer *trip_r, const Integer *trip_c, 
                        const T* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax)
    {
        return sparse_matrix<T>(matcl::make_complex_sparse(trip_r, trip_c, trip_x, r, c, nnz, nzmax));
    };

    static mat_type eval(const Integer *trip_r, const Integer *trip_c, 
                        const Real* re, const Real* im, Integer r, Integer c, Integer nnz, Integer nzmax)
    {
        return sparse_matrix<T>(matcl::make_complex_sparse(trip_r, trip_c, re, im, r, c, nnz, nzmax));
    };
};

template<>
struct make_helper<Float_complex>
{
    using T         = Float_complex;
    using mat_type  = sparse_matrix<T>;
    
    static mat_type eval(const Integer *trip_r, const Integer *trip_c, 
                        const T* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax)
    {
        return sparse_matrix<T>(matcl::make_float_complex_sparse(trip_r, trip_c, trip_x, r, c, nnz, nzmax));
    };

    static mat_type eval(const Integer *trip_r, const Integer *trip_c, 
                        const Float* re, const Float* im, Integer r, Integer c, Integer nnz, Integer nzmax)
    {
        return sparse_matrix<T>(matcl::make_float_complex_sparse(trip_r, trip_c, re, im, r, c, nnz, nzmax));
    };
};

template<>
struct make_helper<Object>
{
    using T         = Object;
    using mat_type  = sparse_matrix<T>;
    
    static mat_type eval(ti::ti_object ti, const Integer *trip_r, const Integer *trip_c, 
                        const T* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax)
    {
        return sparse_matrix<T>(matcl::make_object_sparse(ti, trip_r, trip_c, trip_x, r, c, nnz, nzmax));
    };
};

template<class T>
sparse_matrix<T> sparse_matrix<T>::rand(Integer r, Integer c, Real d, const rand_state& rand_ptr)
{
    return sparse_matrix<T>(matcl::sprand(r,c, d, details::value_to_code<T>::value, rand_ptr));
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::randn(Integer r, Integer c, Real d, const rand_state& rand_ptr)
{
    return sparse_matrix<T>(matcl::sprandn(r,c, d, details::value_to_code<T>::value, rand_ptr));
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::make(Integer rows,Integer cols, Integer nz)
{
    return sparse_matrix<T>(matcl::make_sparse_matrix(rows, cols, nz, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::make(ti::ti_object ti, Integer rows,Integer cols, Integer nz)
{
    return sparse_matrix<T>(matcl::make_object_sparse(ti, rows, cols, nz));
};

template<class T>
sparse_matrix<T> sparse_matrix<T>::make(const Matrix& trip_r, const Matrix& trip_c, 
                        const Matrix& trip_x, Integer r, Integer c, Integer nzmax)
{
    return sparse_matrix<T>(matcl::make_sparse_matrix(trip_r, trip_c, trip_x, r,c,nzmax));
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::make_noinit(Integer rows,Integer cols, Integer nzmax)
{
    return sparse_matrix<T>(make_sparse_noinit(rows,cols,nzmax, details::value_to_code<T>::value));
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::make_noinit(ti::ti_object ti, Integer rows,Integer cols, Integer nzmax)
{
    return sparse_matrix<T>(make_object_sparse_noinit(ti, rows,cols,nzmax));
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::make(const Integer *trip_r, const Integer *trip_c, 
                        const T* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax)
{
     return make_helper<T>::eval(trip_r, trip_c, trip_x, r, c, nnz, nzmax);
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::make(ti::ti_object ti, const Integer *trip_r, const Integer *trip_c, 
                        const T* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax)
{
     return make_helper<T>::eval(ti, trip_r, trip_c, trip_x, r, c, nnz, nzmax);
};

template<class T>
template<class Enable>
sparse_matrix<T> sparse_matrix<T>::make_complex(const Integer *trip_r, const Integer *trip_c, 
                            const typename details::real_type<T>::type * trip_re, 
                            const typename details::real_type<T>::type * trip_im, 
                            Integer r, Integer c, Integer nnz, Integer nzmax )
{
    return make_helper<T>::eval(trip_r, trip_c, trip_re, trip_im, r, c, nnz, nzmax);
};

//--------------------------------------------------------------------
//              sparse_matrix<T, false>
//--------------------------------------------------------------------
template<class T>
sparse_matrix<T,false>::sparse_matrix()
    :base_type()
{};

template<class T>
sparse_matrix<T,false>::sparse_matrix(matcl::Matrix& m)
    :base_type(m, str_make_unique())
{};

template<class T>
sparse_matrix<T,false>::sparse_matrix(matcl::Matrix&& m)
    :base_type(std::move(m), str_make_unique())
{};

template<class T>
sparse_matrix<T,false>::sparse_matrix(sparse_matrix<T,false>&& mat)
    :sparse_matrix(std::move(mat.m_matrix))
{}

template<class T>
sparse_matrix<T,false>::sparse_matrix(sparse_matrix<T,true>&& mat)
    :sparse_matrix(std::move(mat.m_matrix))
{}

template<class T>
sparse_matrix<T,false>::sparse_matrix(sparse_matrix<T,true>& mat)
    :sparse_matrix(mat.make_unique().m_matrix)
{}

template<class T>
sparse_matrix<T,false>& sparse_matrix<T,false>::operator=(sparse_matrix&& mat) &
{
    m_matrix = std::move(mat.m_matrix);
    init_from_rep(m_matrix.get_impl<mat_type>());
    return *this;
};

template<class T>
sparse_matrix<T,false>::~sparse_matrix()
{};

//--------------------------------------------------------------------
//              sub_sparse_matrix<T>
//--------------------------------------------------------------------
template<class T>
const Matrix 
sub_sparse_matrix<T>::to_raw_matrix() const
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
const typename sub_sparse_matrix<T>::matrix_type 
sub_sparse_matrix<T>::to_matrix() const
{
    return matrix_type(this->to_raw_matrix());
};

template<class T>
typename sub_sparse_matrix<T>::matrix_type& 
sub_sparse_matrix<T>::operator=(const matrix_type& mat) const &&
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
typename sub_sparse_matrix<T>::matrix_type& 
sub_sparse_matrix<T>::drop_sparse(Real tol) const &&
{
    Matrix& tmp = m_matrix->m_matrix;
    if (m_colon_2)
    {
        tmp(*m_colon_1,*m_colon_2).drop_sparse(tol);
        return *m_matrix;
    }
    else if(m_colon_1)
    {
        tmp(*m_colon_1).drop_sparse(tol);
        return *m_matrix;
    }
    else
    {
        tmp.diag(m_d).drop_sparse(tol);
        return *m_matrix;
    };
};

template<class T>
typename sub_sparse_matrix<T>::matrix_type& 
sub_sparse_matrix<T>::add_sparse() const &&
{
    if (m_colon_2)
    {
        m_matrix->m_matrix(*m_colon_1,*m_colon_2).add_sparse();
        return *m_matrix;
    }
    else if(m_colon_1)
    {
        m_matrix->m_matrix(*m_colon_1).add_sparse();
        return *m_matrix;
    }
    else
    {
        m_matrix->m_matrix.diag(m_d).add_sparse();
        return *m_matrix;
    };
}

template<class T>
typename sub_sparse_matrix<T>::matrix_type& 
sub_sparse_matrix<T>::operator=(const sub_sparse_matrix<T>& mat0) const &&
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
typename sub_sparse_matrix<T>::matrix_type& 
sub_sparse_matrix<T>::operator=(const T& val) const &&
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
typename sub_sparse_matrix<T>::matrix_type& 
sub_sparse_matrix<T>::operator=(const sub_sparse_matrix_1<T>& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

template<class T>
typename sub_sparse_matrix<T>::matrix_type& 
sub_sparse_matrix<T>::operator=(const sub_sparse_matrix_2<T>& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

template<class T>
T sub_sparse_matrix_1<T>::to_matrix() const
{
    const matrix_type& mat = *m_matrix;
    return mat(m_ind_1);
};

template<class T>
typename sub_sparse_matrix_1<T>::matrix_type& 
sub_sparse_matrix_1<T>::drop_sparse(Real tol) const &&
{
    m_matrix->m_matrix(m_ind_1).drop_sparse(tol);
    return *m_matrix;
};

template<class T>
typename sub_sparse_matrix_1<T>::matrix_type& 
sub_sparse_matrix_1<T>::add_sparse() const &&
{
    m_matrix->m_matrix(m_ind_1).add_sparse();
    return *m_matrix;
};

template<class T>
typename sub_sparse_matrix_1<T>::matrix_type& 
sub_sparse_matrix_1<T>::operator=(const matrix_type& mat) const &&
{
    m_matrix->m_matrix(m_ind_1) = mat.to_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_sparse_matrix_1<T>::matrix_type& 
sub_sparse_matrix_1<T>::operator=(const T& val) const &&
{
    m_matrix->m_matrix(m_ind_1) = val;

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_sparse_matrix_1<T>::matrix_type& 
sub_sparse_matrix_1<T>::operator=(const sub_sparse_matrix<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1) = mat.to_raw_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_sparse_matrix_1<T>::matrix_type& 
sub_sparse_matrix_1<T>::operator=(const sub_sparse_matrix_1<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1) = mat.to_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_sparse_matrix_1<T>::matrix_type& 
sub_sparse_matrix_1<T>::operator=(const sub_sparse_matrix_2<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1) = mat.to_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
T sub_sparse_matrix_2<T>::to_matrix() const
{
    const matrix_type& mat = *m_matrix;
    return mat(m_ind_1,m_ind_2);
};

template<class T>
typename sub_sparse_matrix_2<T>::matrix_type& 
sub_sparse_matrix_2<T>::drop_sparse(Real tol) const &&
{
    m_matrix->m_matrix(m_ind_1, m_ind_2).drop_sparse(tol);
    return *m_matrix;
};

template<class T>
typename sub_sparse_matrix_2<T>::matrix_type& 
sub_sparse_matrix_2<T>::add_sparse() const &&
{
    m_matrix->m_matrix(m_ind_1, m_ind_2).add_sparse();
    return *m_matrix;
};

template<class T>
typename sub_sparse_matrix_2<T>::matrix_type& 
sub_sparse_matrix_2<T>::operator=(const matrix_type& mat) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = mat.to_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_sparse_matrix_2<T>::matrix_type& 
sub_sparse_matrix_2<T>::operator=(const T& val) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = val;

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_sparse_matrix_2<T>::matrix_type& 
sub_sparse_matrix_2<T>::operator=(const sub_sparse_matrix<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = mat.to_raw_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_sparse_matrix_2<T>::matrix_type& 
sub_sparse_matrix_2<T>::operator=(const sub_sparse_matrix_1<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = mat.to_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
typename sub_sparse_matrix_2<T>::matrix_type& 
sub_sparse_matrix_2<T>::operator=(const sub_sparse_matrix_2<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = mat.to_matrix();

    //m_matrix->m_matrix could be made unique, so we need to update rep pointers
    m_matrix->update_rep();
    return *m_matrix;
}

template<class T>
template<class V> 
inline V sub_sparse_matrix<T>::get_scalar() const
{
    return Matrix(*this).get_scalar<V>();
};

template<class T>
inline Integer sub_sparse_matrix<T>::rows() const
{
    return Matrix(*this).rows();
};

template<class T>
inline Integer sub_sparse_matrix<T>::cols() const
{
    return Matrix(*this).cols();
};

template<class T>
inline Integer sub_sparse_matrix<T>::length() const
{
    return Matrix(*this).length();
};

template<class T>
inline Real sub_sparse_matrix<T>::numel() const
{
    return Matrix(*this).numel();
};

template<class T>
inline bool sub_sparse_matrix<T>::all_finite() const
{
    return Matrix(*this).all_finite();
};

template<class T>
inline sub_sparse_matrix<T>::operator bool() const
{
    return Matrix(*this).operator bool();
};

template<class T>
template<class V> 
inline V sub_sparse_matrix_1<T>::get_scalar() const
{
    return Matrix(*this).get_scalar<V>();
};

template<class T>
inline Integer sub_sparse_matrix_1<T>::rows() const
{
    return Matrix(*this).rows();
};

template<class T>
inline Integer sub_sparse_matrix_1<T>::cols() const
{
    return Matrix(*this).cols();
};

template<class T>
inline Integer sub_sparse_matrix_1<T>::length() const
{
    return Matrix(*this).length();
};

template<class T>
inline Real sub_sparse_matrix_1<T>::numel() const
{
    return Matrix(*this).numel();
};

template<class T>
inline bool sub_sparse_matrix_1<T>::all_finite() const
{
    return Matrix(*this).all_finite();
};

template<class T>
inline sub_sparse_matrix_1<T>::operator bool() const
{
    return Matrix(*this).operator bool();
};

template<class T>
template<class V> 
inline V sub_sparse_matrix_2<T>::get_scalar() const
{
    return Matrix(*this).get_scalar<V>();
};

template<class T>
inline Integer sub_sparse_matrix_2<T>::rows() const
{
    return Matrix(*this).rows();
};

template<class T>
inline Integer sub_sparse_matrix_2<T>::cols() const
{
    return Matrix(*this).cols();
};

template<class T>
inline Integer sub_sparse_matrix_2<T>::length() const
{
    return Matrix(*this).length();
};

template<class T>
inline Real sub_sparse_matrix_2<T>::numel() const
{
    return Matrix(*this).numel();
};

template<class T>
inline bool sub_sparse_matrix_2<T>::all_finite() const
{
    return Matrix(*this).all_finite();
};

template<class T>
inline sub_sparse_matrix_2<T>::operator bool() const
{
    return Matrix(*this).operator bool();
};

//--------------------------------------------------------------------
//              TESTS
//--------------------------------------------------------------------
template<class T>
void test_compile(const sparse_matrix<T>& v)
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

    using vect  = std::vector<sparse_matrix<T>>;

    delrows(v, colon());
    delcols(v,colon());
    delrowscols(v,colon(), colon());
    horzcat(v,v);
    vertcat(v,v);
    blkdiag(v,v);
    horzcat({v,v});
    vertcat({v,v});
    blkdiag({v,v});
    horzcat(vect{v,v});
    vertcat(vect{v,v});
    blkdiag(vect{v,v});

    repmat(v,1,1);

    sparse(v);

    {
        dense_matrix<T> tmp;
        sparse(tmp);
    };
    {
        band_matrix<T> tmp;
        sparse(tmp);
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

template class MATCL_MATREP_EXPORT_REP matcl::sparse_matrix<Integer,true>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_matrix<Real,true>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_matrix<Float,true>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_matrix<Complex,true>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_matrix<Float_complex,true>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_matrix<Object,true>;

template class MATCL_MATREP_EXPORT_REP matcl::sparse_matrix<Integer,false>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_matrix<Real,false>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_matrix<Float,false>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_matrix<Complex,false>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_matrix<Float_complex,false>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_matrix<Object,false>;

template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix<Integer>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix<Real>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix<Float>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix<Complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix<Float_complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix<Object>;

template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix_1<Integer>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix_1<Real>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix_1<Float>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix_1<Complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix_1<Float_complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix_1<Object>;

template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix_2<Integer>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix_2<Real>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix_2<Float>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix_2<Complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix_2<Float_complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sub_sparse_matrix_2<Object>;

template class MATCL_MATREP_EXPORT_REP matcl::sparse_row<Integer>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_row<Real>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_row<Float>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_row<Complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_row<Float_complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_row<Object>;

template class MATCL_MATREP_EXPORT_REP matcl::sparse_col<Integer>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_col<Real>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_col<Float>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_col<Complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_col<Float_complex>;
template class MATCL_MATREP_EXPORT_REP matcl::sparse_col<Object>;

#pragma warning(push)
#pragma warning(disable:5037) //an out-of-line definition of a member of a class template
                              // cannot have default arguments

// is it a VS bug?

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float> sparse_matrix<Float>::randn(Integer r, Integer c, Real d, const rand_state& rand_ptr);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Real> sparse_matrix<Real>::randn(Integer r, Integer c, Real d, const rand_state& rand_ptr);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float_complex> sparse_matrix<Float_complex>::randn(Integer r, Integer c, Real d, 
                                                                 const rand_state& rand_ptr);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Complex> sparse_matrix<Complex>::randn(Integer r, Integer c, Real d, const rand_state& rand_ptr);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Integer> sparse_matrix<Integer>::make(Integer r, Integer c, Integer nnz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Real> sparse_matrix<Real>::make(Integer r, Integer c, Integer nnz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float> sparse_matrix<Float>::make(Integer r, Integer c, Integer nnz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Complex> sparse_matrix<Complex>::make(Integer r, Integer c, Integer nnz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float_complex> sparse_matrix<Float_complex>::make(Integer r, Integer c, Integer nnz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Object> sparse_matrix<Object>::make(ti::ti_object, Integer r, Integer c, Integer nnz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Integer> sparse_matrix<Integer>::make_noinit(Integer r, Integer c, Integer nnz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Real> sparse_matrix<Real>::make_noinit(Integer r, Integer c, Integer nnz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float> sparse_matrix<Float>::make_noinit(Integer r, Integer c, Integer nnz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Complex> sparse_matrix<Complex>::make_noinit(Integer r, Integer c, Integer nnz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float_complex> sparse_matrix<Float_complex>::make_noinit(Integer r, Integer c, Integer nnz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Object> sparse_matrix<Object>::make_noinit(ti::ti_object, Integer r, Integer c, Integer nnz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Integer> sparse_matrix<Integer>::make(const Integer *trip_r, const Integer *trip_c, 
                            const Integer* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Real> sparse_matrix<Real>::make(const Integer *trip_r, const Integer *trip_c, 
                            const Real* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float> sparse_matrix<Float>::make(const Integer *trip_r, const Integer *trip_c, 
                            const Float* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Complex> sparse_matrix<Complex>::make(const Integer *trip_r, const Integer *trip_c, 
                            const Complex* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float_complex> sparse_matrix<Float_complex>::make(const Integer *trip_r, const Integer *trip_c, 
                            const Float_complex* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Object> sparse_matrix<Object>::make(ti::ti_object ti, const Integer *trip_r, const Integer *trip_c, 
                            const Object* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Complex> sparse_matrix<Complex>::make_complex(const Integer *trip_r, const Integer *trip_c, 
                            const Real * trip_re, const Real* trip_im, 
                            Integer r, Integer c, Integer nnz, Integer nzmax);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float_complex> sparse_matrix<Float_complex>::make_complex(const Integer *trip_r, const Integer *trip_c, 
                            const Float * trip_re, const Float* trip_im, 
                            Integer r, Integer c, Integer nnz, Integer nzmax);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Integer> sparse_matrix<Integer>::zeros(Integer r, Integer c, Integer nz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Real> sparse_matrix<Real>::zeros(Integer r, Integer c, Integer nz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float> sparse_matrix<Float>::zeros(Integer r, Integer c, Integer nz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Complex> sparse_matrix<Complex>::zeros(Integer r, Integer c, Integer nz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float_complex> sparse_matrix<Float_complex>::zeros(Integer r, Integer c, Integer nz);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Object> sparse_matrix<Object>::zeros(ti::ti_object, Integer r, Integer c, Integer nz);

#pragma warning(pop)

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Integer> sparse_matrix<Integer>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Real> sparse_matrix<Real>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float> sparse_matrix<Float>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Complex> sparse_matrix<Complex>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float_complex> sparse_matrix<Float_complex>::ones(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Object> sparse_matrix<Object>::ones(ti::ti_object, Integer r, Integer c);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Integer> sparse_matrix<Integer>::eye(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Real> sparse_matrix<Real>::eye(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float> sparse_matrix<Float>::eye(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Complex> sparse_matrix<Complex>::eye(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float_complex> sparse_matrix<Float_complex>::eye(Integer r, Integer c);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Object> sparse_matrix<Object>::eye(ti::ti_object, Integer r, Integer c);

template MATCL_MATREP_EXPORT_REP
sparse_matrix<Integer> sparse_matrix<Integer>::eye(Integer r);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Real> sparse_matrix<Real>::eye(Integer r);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float> sparse_matrix<Float>::eye(Integer r);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Complex> sparse_matrix<Complex>::eye(Integer r);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Float_complex> sparse_matrix<Float_complex>::eye(Integer r);
template MATCL_MATREP_EXPORT_REP
sparse_matrix<Object> sparse_matrix<Object>::eye(ti::ti_object, Integer r);

template void test_compile(const sparse_matrix<Integer>& v);
template void test_compile(const sparse_matrix<Real>& v);
template void test_compile(const sparse_matrix<Float>& v);
template void test_compile(const sparse_matrix<Complex>& v);
template void test_compile(const sparse_matrix<Float_complex>& v);
template void test_compile(const sparse_matrix<Object>& v);

};
