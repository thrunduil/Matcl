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

#include "matcl-matrep/objects/object_matrix.h"
#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/details/matrix.inl"
#include "matcl-matrep/details/matrix_details_subs.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-matrep/lib_functions/vecfunc.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"

namespace matcl
{

//-----------------------------------------------------------------
//                  object_row
//-----------------------------------------------------------------
template<class T>
object_row<T>::object_row()
{};

template<class T>
object_row<T>::~object_row()
{};

template<class T>
object_row<T>::object_row(const object_row& other)
    :base_type(other)
{};

template<class T>
object_row<T>::object_row(object_row&& other)
    :base_type(std::move(other))
{};

template<class T>
object_row<T>& object_row<T>::operator=(const object_row& rhs)
{
    base_type::operator=(rhs);
    return *this;
};

template<class T>
object_row<T>& object_row<T>::operator=(object_row&& rhs)
{
    base_type::operator=(std::move(rhs));
    return *this;
};

template<class T>
object_row<T>& object_row<T>::add(const T& rhs)
{
    add(object_type<T>(rhs));
    return *this;
};

template<class T>
object_row<T>& object_row<T>::add(T&& rhs)
{
    add(object_type<T>(std::move(rhs)));
    return *this;
};

template<class T>
object_row<T>& object_row<T>::add(const object_type<T>& rhs)
{
    base_type::add(Object(rhs));
    return *this;
};

template<class T>
object_row<T>& object_row<T>::add(object_type<T>&& rhs)
{
    base_type::add(Object(std::move(rhs)));
    return *this;
};

template<class T>
object_row<T>& object_row<T>::add(const object_matrix<T>& rhs)
{
    base_type::add(Matrix(rhs));
    return *this;
};

template<class T>
object_row<T>& object_row<T>::add(object_matrix<T>&& rhs)
{
    base_type::add(Matrix(std::move(rhs)));
    return *this;
};

template<class T>
object_row<T>& object_row<T>::add(const object_row& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
object_row<T>& object_row<T>::add(object_row&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
object_row<T>& object_row<T>::add(const object_col<T>& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
object_row<T>& object_row<T>::add(object_col<T>&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
const object_matrix<T> object_row<T>::to_matrix() const
{
    return object_matrix<T>(base_type::to_matrix());
}

//-----------------------------------------------------------------
//                  object_col
//-----------------------------------------------------------------
template<class T>
object_col<T>::object_col()
{};

template<class T>
object_col<T>::~object_col()
{};

template<class T>
object_col<T>::object_col(const object_col& other)
    :base_type(other)
{};

template<class T>
object_col<T>::object_col(object_col&& other)
    :base_type(std::move(other))
{};

template<class T>
object_col<T>& object_col<T>::operator=(const object_col& rhs)
{
    base_type::operator=(rhs);
    return *this;
};

template<class T>
object_col<T>& object_col<T>::operator=(object_col&& rhs)
{
    base_type::operator=(std::move(rhs));
    return *this;
};

template<class T>
object_col<T>& object_col<T>::add(const T& rhs)
{
    add(object_type<T>(rhs));
    return *this;
};

template<class T>
object_col<T>& object_col<T>::add(T&& rhs)
{
    add(object_type<T>(std::move(rhs)));
    return *this;
};

template<class T>
object_col<T>& object_col<T>::add(const object_type<T>& rhs)
{
    base_type::add(Object(rhs));
    return *this;
};

template<class T>
object_col<T>& object_col<T>::add(object_type<T>&& rhs)
{
    base_type::add(Object(std::move(rhs)));
    return *this;
};

template<class T>
object_col<T>& object_col<T>::add(const object_matrix<T>& rhs)
{
    base_type::add(Matrix(rhs));
    return *this;
};

template<class T>
object_col<T>& object_col<T>::add(object_matrix<T>&& rhs)
{
    base_type::add(Matrix(std::move(rhs)));
    return *this;
};

template<class T>
object_col<T>& object_col<T>::add(const object_row<T>& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
object_col<T>& object_col<T>::add(object_row<T>&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
object_col<T>& object_col<T>::add(const object_col<T>& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
object_col<T>& object_col<T>::add(object_col<T>&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
const object_matrix<T> object_col<T>::to_matrix() const
{
    return object_matrix<T>(base_type::to_matrix());
}

//-----------------------------------------------------------------
//                  object_matrix
//-----------------------------------------------------------------
template<class T>
object_matrix<T>::object_matrix()
    :m_matrix(object_type())
{};

template<class T>
object_matrix<T>::object_matrix(const T& val)
    :m_matrix(object_type(val))
{};

template<class T>
object_matrix<T>::object_matrix(const object_type& v)
    :m_matrix(v)
{};

template<class T>
object_matrix<T>::object_matrix(object_type&& v)
    :m_matrix(std::move(v))
{};

template<class T>
object_matrix<T>::object_matrix(const matcl::Matrix& m)
    :m_matrix(convert_object(m,object_type::get_static_type()))
{};

template<class T>
object_matrix<T>::object_matrix(matcl::Matrix&& m)
    :m_matrix(std::move(m))
{
    m_matrix = convert_object(m_matrix,object_type::get_static_type());
};

template<class T>
object_matrix<T>::object_matrix(const sub_object_matrix& sub)
    :object_matrix(sub.to_matrix())
{};

template<class T>
object_matrix<T>::object_matrix(const sub_object_matrix_1& sub)
    :object_matrix(sub.to_matrix())
{};

template<class T>
object_matrix<T>::object_matrix(const sub_object_matrix_2& sub)
    :object_matrix(sub.to_matrix())
{};

template<class T>
object_matrix<T>::object_matrix(const object_row<T>& sub)
    :object_matrix(sub.to_matrix())
{};

template<class T>
object_matrix<T>::object_matrix(const object_col<T>& sub)
    :object_matrix(sub.to_matrix())
{};

template<class T>
object_matrix<T>::object_matrix(const object_matrix& mat)
    :m_matrix(mat.m_matrix)
{};

template<class T>
object_matrix<T>::object_matrix(object_matrix&& mat)
    :m_matrix(std::move(mat.m_matrix))
{};

template<class T>
object_matrix<T>& object_matrix<T>::operator=(const object_matrix& mat) &
{
    m_matrix = mat.m_matrix;
    return *this;
};

template<class T>
object_matrix<T>& object_matrix<T>::operator=(object_matrix&& mat) &
{
    m_matrix = std::move(mat.m_matrix);
    return *this;
};

template<class T>
object_matrix<T>::~object_matrix()
{};

template<class T>
const struct_flag object_matrix<T>::get_struct() const
{
    return m_matrix.get_struct();
};

template<class T>
void object_matrix<T>::set_struct(const struct_flag& sc) const
{
    m_matrix.set_struct(sc);
};

template<class T>
void object_matrix<T>::add_struct(const struct_flag& sc) const
{
    m_matrix.add_struct(sc);
};

template<class T>
ti::ti_object object_matrix<T>::get_type() const
{
    return m_matrix.get_type();
};

template<class T>
const typename object_matrix<T>::object_matrix 
object_matrix<T>::delrows(const colon& c) const &
{
    return object_matrix(m_matrix.delrows(c));
}

template<class T>
const typename object_matrix<T>::object_matrix 
object_matrix<T>::delrows(const colon& c) const &&
{
    return object_matrix(std::move(m_matrix).delrows(c));
}

template<class T>
const typename object_matrix<T>::object_matrix
object_matrix<T>::delcols(const colon& c) const &
{
    return object_matrix(m_matrix.delcols(c));
};

template<class T>
const typename object_matrix<T>::object_matrix
object_matrix<T>::delcols(const colon& c) const &&
{
    return object_matrix(std::move(m_matrix).delcols(c));
};

template<class T>
const typename object_matrix<T>::object_matrix 
object_matrix<T>::delrowscols(const colon& c1, const colon& c2) const &
{
    return object_matrix(m_matrix.delrowscols(c1, c2));
};

template<class T>
const typename object_matrix<T>::object_matrix 
object_matrix<T>::delrowscols(const colon& c1, const colon& c2) const &&
{
    return object_matrix(std::move(m_matrix).delrowscols(c1, c2));
};

template<class T>
const typename object_matrix<T>::object_matrix
object_matrix<T>::operator()(const colon& r) const
{
    return object_matrix(m_matrix(r));
};

template<class T>
const typename object_matrix<T>::object_matrix
object_matrix<T>::operator()(const colon& r, const colon& c) const
{
    return object_matrix(m_matrix(r,c));
};

template<class T>
object_type<T> object_matrix<T>::operator()(Integer p) const
{
    return object_type(m_matrix(p).get_scalar<Object>(), dynamic::from_object());
};

template<class T>
typename object_matrix<T>::sub_object_matrix_1
object_matrix<T>::operator()(Integer p)
{
    sub_object_matrix_1 ret(this,p);
    return ret;
};

template<class T>
object_type<T> object_matrix<T>::operator()(Integer r, Integer c) const
{
    return object_type(m_matrix(r,c).get_scalar<Object>(), dynamic::from_object());
};

template<class T>
typename object_matrix<T>::sub_object_matrix_2
object_matrix<T>::operator()(Integer r, Integer c)
{
    sub_object_matrix_2 ret(this,r,c);
    return ret;
};

template<class T>
const object_matrix<T> object_matrix<T>::diag(Integer d) const
{
    return object_matrix(m_matrix.diag(d));
};

template<class T>
Integer object_matrix<T>::rows() const
{
    return m_matrix.rows();
};

template<class T>
Integer object_matrix<T>::cols() const
{
    return m_matrix.cols();
};

template<class T>
Integer object_matrix<T>::length() const
{
    Integer r = rows();
    Integer c = cols();

    if (r == 0 || c == 0)
        return 0;

    return (r > c) ? r : c;
};

template<class T>
Integer object_matrix<T>::structural_nnz() const
{
    return m_matrix.structural_nnz();
};

template<class T>
Integer object_matrix<T>::structural_ldiags(bool use_flags) const
{
    return m_matrix.structural_ldiags(use_flags);
};

template<class T>
Integer object_matrix<T>::structural_udiags(bool use_flags) const
{
    return m_matrix.structural_udiags(use_flags);
};

template<class T>
Real object_matrix<T>::numel() const
{
    return Real(rows()) * Real(cols());
};

template<class T>
bool object_matrix<T>::all_finite() const
{
    return m_matrix.all_finite();
};

template<class T>
const typename object_matrix<T>::object_matrix
object_matrix<T>::clone() const
{
    return object_matrix(m_matrix.clone());
};

template<class T>
const typename object_matrix<T>::object_matrix&
object_matrix<T>::make_unique() const
{
    m_matrix.make_unique();
    return *this;
};

template<class T>
void object_matrix<T>::resize(Integer r, Integer c)
{
    m_matrix.resize(r,c);
};

template<class T>
void object_matrix<T>::reserve(Integer r, Integer c)
{
    m_matrix.reserve(r,c);
};

template<class T>
void object_matrix<T>::resize_band(Integer r, Integer c, Integer fd, Integer ld)
{
    m_matrix.resize_band(r,c, fd, ld);
};

template<class T>
void object_matrix<T>::reserve_band(Integer r, Integer c, Integer fd, Integer ld)
{
    m_matrix.reserve_band(r,c,fd,ld);
};

template<class T>
object_type<T> object_matrix<T>::get_scalar() const
{
    matcl::dynamic::object obj = m_matrix.get_scalar<Object>();
    return *reinterpret_cast<const object_type*>(&obj);
};

template<class T>
object_type<T>& object_matrix<T>::get_scalar_unique()
{
    return reinterpret_cast<object_type&>(static_cast<matcl::dynamic::object&>
                (m_matrix.get_scalar_unique<Object>()));
};

template<class T>
const object_type<T>* object_matrix<T>::get_array() const
{
    return reinterpret_cast<const object_type*>(m_matrix.get_array<Object>());
};

template<class T>
object_type<T>* object_matrix<T>::get_array_unique()
{
    return reinterpret_cast<object_type*>(m_matrix.get_array_unique<Object>());
};

template<class T>
sub_object_matrix<T> object_matrix<T>::operator()(const colon& r)
{
    sub_object_matrix ret(this,r);
    return ret;
};

template<class T>
sub_object_matrix<T> object_matrix<T>::operator()(const colon& r, const colon& c)
{
    sub_object_matrix ret(this,r,c);
    return ret;
};

template<class T>
sub_object_matrix<T> object_matrix<T>::diag(Integer d)
{
    sub_object_matrix ret(d,this);
    return ret;
};

template<class T>
const typename sub_object_matrix<T>::matrix_type 
sub_object_matrix<T>::to_matrix() const
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
typename sub_object_matrix<T>::matrix_type& 
sub_object_matrix<T>::drop_sparse(Real tol) const &&
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
typename sub_object_matrix<T>::matrix_type& 
sub_object_matrix<T>::add_sparse() const &&
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
typename sub_object_matrix<T>::matrix_type& 
sub_object_matrix<T>::operator=(const matrix_type& mat) const &&
{
    if (m_colon_2)
    {
        m_matrix->m_matrix(*m_colon_1,*m_colon_2) = mat.to_matrix();
        return *m_matrix;
    }
    else if (m_colon_1)
    {
        m_matrix->m_matrix(*m_colon_1) = mat.to_matrix();
        return *m_matrix;
    }
    else
    {
        m_matrix->m_matrix.diag(m_d) = mat.to_matrix();
        return *m_matrix;
    };
}

template<class T>
typename sub_object_matrix<T>::matrix_type& 
sub_object_matrix<T>::operator=(const sub_object_matrix<T>& mat0) const &&
{
    if (m_colon_2)
    {
        m_matrix->m_matrix(*m_colon_1,*m_colon_2) = mat0.m_matrix->m_matrix;
        return *m_matrix;
    }
    else if (m_colon_1)
    {
        m_matrix->m_matrix(*m_colon_1) = mat0.m_matrix->m_matrix;
        return *m_matrix;
    }
    else
    {
        m_matrix->m_matrix.diag(m_d) = mat0.m_matrix->m_matrix;
        return *m_matrix;
    };
};

template<class T>
typename sub_object_matrix<T>::matrix_type& 
sub_object_matrix<T>::operator=(const T& val) const &&
{
    return std::move(*this).operator=(object_type<T>(val));
};

template<class T>
typename sub_object_matrix<T>::matrix_type& 
sub_object_matrix<T>::operator=(T&& val) const &&
{
    return std::move(*this).operator=(object_type<T>(std::move(val)));
};

template<class T>
typename sub_object_matrix<T>::matrix_type& 
sub_object_matrix<T>::operator=(const object_type<T>& val) const &&
{
    if (m_colon_2)
    {
        m_matrix->m_matrix(*m_colon_1,*m_colon_2) = Object(val);
        return *m_matrix;
    }
    else if (m_colon_1)
    {
        m_matrix->m_matrix(*m_colon_1) = Object(val);
        return *m_matrix;
    }
    else
    {
        m_matrix->m_matrix.diag(m_d) = Object(val);
        return *m_matrix;
    };
};

template<class T>
typename sub_object_matrix<T>::matrix_type& 
sub_object_matrix<T>::operator=(object_type<T>&& val) const &&
{
    if (m_colon_2)
    {
        m_matrix->m_matrix(*m_colon_1,*m_colon_2) = Object(std::move(val));
        return *m_matrix;
    }
    else if (m_colon_1)
    {
        m_matrix->m_matrix(*m_colon_1) = Object(std::move(val));
        return *m_matrix;
    }
    else
    {
        m_matrix->m_matrix.diag(m_d) = Object(std::move(val));
        return *m_matrix;
    };
};

template<class T>
typename sub_object_matrix<T>::matrix_type& 
sub_object_matrix<T>::operator=(const sub_object_matrix_1<T>& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

template<class T>
typename sub_object_matrix<T>::matrix_type& 
sub_object_matrix<T>::operator=(const sub_object_matrix_2<T>& mat0) const &&
{
    return std::move(*this).operator=(mat0.to_matrix());
};

template<class T>
object_type<T> sub_object_matrix_1<T>::to_matrix() const
{
    const matrix_type& mat = *m_matrix;
    return mat(m_ind_1);
};

template<class T>
typename sub_object_matrix_1<T>::matrix_type& 
sub_object_matrix_1<T>::drop_sparse(Real tol) const &&
{
    m_matrix->m_matrix(m_ind_1).drop_sparse(tol);
    return *m_matrix;
};

template<class T>
typename sub_object_matrix_1<T>::matrix_type& 
sub_object_matrix_1<T>::add_sparse() const &&
{
    m_matrix->m_matrix(m_ind_1).add_sparse();
    return *m_matrix;
};

template<class T>
typename sub_object_matrix_1<T>::matrix_type& 
sub_object_matrix_1<T>::operator=(const matrix_type& mat) const &&
{
    m_matrix->m_matrix(m_ind_1) = mat.to_matrix();
    return *m_matrix;
}

template<class T>
typename sub_object_matrix_1<T>::matrix_type& 
sub_object_matrix_1<T>::operator=(const T& val) const &&
{
    return std::move(*this).operator=(object_type<T>(val));
};

template<class T>
typename sub_object_matrix_1<T>::matrix_type& 
sub_object_matrix_1<T>::operator=(T&& val) const &&
{
    return std::move(*this).operator=(object_type<T>(std::move(val)));
};

template<class T>
typename sub_object_matrix_1<T>::matrix_type& 
sub_object_matrix_1<T>::operator=(const object_type<T>& val) const &&
{
    m_matrix->m_matrix(m_ind_1) = Object(val);
    return *m_matrix;
}

template<class T>
typename sub_object_matrix_1<T>::matrix_type& 
sub_object_matrix_1<T>::operator=(object_type<T>&& val) const &&
{
    m_matrix->m_matrix(m_ind_1) = Object(std::move(val));
    return *m_matrix;
}

template<class T>
typename sub_object_matrix_1<T>::matrix_type& 
sub_object_matrix_1<T>::operator=(const sub_object_matrix<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1) = Matrix(mat.to_matrix());
    return *m_matrix;
}

template<class T>
typename sub_object_matrix_1<T>::matrix_type& 
sub_object_matrix_1<T>::operator=(const sub_object_matrix_1<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1) = Object(mat.to_matrix());
    return *m_matrix;
}

template<class T>
typename sub_object_matrix_1<T>::matrix_type& 
sub_object_matrix_1<T>::operator=(const sub_object_matrix_2<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1) = Object(mat.to_matrix());
    return *m_matrix;
}

template<class T>
object_type<T> sub_object_matrix_2<T>::to_matrix() const
{
    const matrix_type& mat = *m_matrix;
    return mat(m_ind_1,m_ind_2);
};

template<class T>
typename sub_object_matrix_2<T>::matrix_type& 
sub_object_matrix_2<T>::drop_sparse(Real tol) const &&
{
    m_matrix->m_matrix(m_ind_1, m_ind_2).drop_sparse(tol);
    return *m_matrix;
};

template<class T>
typename sub_object_matrix_2<T>::matrix_type& 
sub_object_matrix_2<T>::add_sparse() const &&
{
    m_matrix->m_matrix(m_ind_1, m_ind_2).add_sparse();
    return *m_matrix;
};

template<class T>
typename sub_object_matrix_2<T>::matrix_type& 
sub_object_matrix_2<T>::operator=(const matrix_type& mat) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = mat.to_matrix();
    return *m_matrix;
}

template<class T>
typename sub_object_matrix_2<T>::matrix_type& 
sub_object_matrix_2<T>::operator=(const T& val) const &&
{
    return std::move(*this).operator=(object_type<T>(val));
}

template<class T>
typename sub_object_matrix_2<T>::matrix_type& 
sub_object_matrix_2<T>::operator=(T&& val) const &&
{
    return std::move(*this).operator=(object_type<T>(std::move(val)));
}

template<class T>
typename sub_object_matrix_2<T>::matrix_type& 
sub_object_matrix_2<T>::operator=(const object_type<T>& val) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = Object(val);
    return *m_matrix;
}

template<class T>
typename sub_object_matrix_2<T>::matrix_type& 
sub_object_matrix_2<T>::operator=(object_type<T>&& val) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = Object(std::move(val));
    return *m_matrix;
}

template<class T>
typename sub_object_matrix_2<T>::matrix_type& 
sub_object_matrix_2<T>::operator=(const sub_object_matrix<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = Matrix(mat.to_matrix());
    return *m_matrix;
}

template<class T>
typename sub_object_matrix_2<T>::matrix_type& 
sub_object_matrix_2<T>::operator=(const sub_object_matrix_1<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = Object(mat.to_matrix());
    return *m_matrix;
}

template<class T>
typename sub_object_matrix_2<T>::matrix_type& 
sub_object_matrix_2<T>::operator=(const sub_object_matrix_2<T>& mat) const &&
{
    m_matrix->m_matrix(m_ind_1,m_ind_2) = Object(mat.to_matrix());
    return *m_matrix;
}

template<class T>
inline sub_object_matrix<T>::sub_object_matrix(matrix_type* m, const colon& c1)
    :m_matrix(m), m_colon_1(&c1),m_colon_2(nullptr),m_d(0)
{}

template<class T>
inline sub_object_matrix<T>::sub_object_matrix(matrix_type* m, const colon& c1, const colon& c2)
    :m_matrix(m), m_colon_1(&c1),m_colon_2(&c2),m_d(0)
{};

template<class T>
inline sub_object_matrix<T>::sub_object_matrix(Integer diag, matrix_type* m)
    :m_matrix(m), m_colon_1(nullptr), m_colon_2(nullptr), m_d(diag)
{};

template<class T>
inline sub_object_matrix<T>::sub_object_matrix(const sub_object_matrix& m)
    :m_matrix(m.m_matrix),m_colon_1(m.m_colon_1),m_colon_2(m.m_colon_2),m_d(m.m_d)
{};

template<class T>
inline sub_object_matrix<T>::sub_object_matrix(sub_object_matrix&& m)
    :m_matrix(m.m_matrix),m_colon_1(m.m_colon_1),m_colon_2(m.m_colon_2),m_d(m.m_d)
{};

template<class T>
inline sub_object_matrix_1<T>::sub_object_matrix_1(matrix_type* m, Integer r)
    :m_matrix(m), m_ind_1(r)
{}

template<class T>
inline sub_object_matrix_1<T>::sub_object_matrix_1(const sub_object_matrix_1& m)
    :m_matrix(m.m_matrix),m_ind_1(m.m_ind_1)
{};

template<class T>
inline sub_object_matrix_1<T>::sub_object_matrix_1(sub_object_matrix_1&& m)
    :m_matrix(m.m_matrix),m_ind_1(m.m_ind_1)
{};

template<class T>
inline sub_object_matrix_2<T>::sub_object_matrix_2(matrix_type* m, Integer r, Integer c)
    :m_matrix(m), m_ind_1(r), m_ind_2(c)
{}

template<class T>
inline sub_object_matrix_2<T>::sub_object_matrix_2(const sub_object_matrix_2& m)
    :m_matrix(m.m_matrix),m_ind_1(m.m_ind_1),m_ind_2(m.m_ind_2)
{};

template<class T>
inline sub_object_matrix_2<T>::sub_object_matrix_2(sub_object_matrix_2&& m)
    :m_matrix(m.m_matrix),m_ind_1(m.m_ind_1),m_ind_2(m.m_ind_2)
{};

template<class T>
template<class V> 
inline V sub_object_matrix<T>::get_scalar() const
{
    return Matrix(*this).get_scalar<V>();
};

template<class T>
inline Integer sub_object_matrix<T>::rows() const
{
    return Matrix(*this).rows();
};

template<class T>
inline Integer sub_object_matrix<T>::cols() const
{
    return Matrix(*this).cols();
};

template<class T>
inline Integer sub_object_matrix<T>::length() const
{
    return Matrix(*this).length();
};

template<class T>
inline Real sub_object_matrix<T>::numel() const
{
    return Matrix(*this).numel();
};

template<class T>
inline bool sub_object_matrix<T>::all_finite() const
{
    return Matrix(*this).all_finite();
};

template<class T>
inline sub_object_matrix<T>::operator bool() const
{
    return Matrix(*this).operator bool();
};

template<class T>
template<class V> 
inline V sub_object_matrix_1<T>::get_scalar() const
{
    return Matrix(*this).get_scalar<V>();
};

template<class T>
inline Integer sub_object_matrix_1<T>::rows() const
{
    return Matrix(*this).rows();
};

template<class T>
inline Integer sub_object_matrix_1<T>::cols() const
{
    return Matrix(*this).cols();
};

template<class T>
inline Integer sub_object_matrix_1<T>::length() const
{
    return Matrix(*this).length();
};

template<class T>
inline Real sub_object_matrix_1<T>::numel() const
{
    return Matrix(*this).numel();
};

template<class T>
inline bool sub_object_matrix_1<T>::all_finite() const
{
    return Matrix(*this).all_finite();
};

template<class T>
inline sub_object_matrix_1<T>::operator bool() const
{
    return Matrix(*this).operator bool();
};

template<class T>
template<class V> 
inline V sub_object_matrix_2<T>::get_scalar() const
{
    return Matrix(*this).get_scalar<V>();
};

template<class T>
inline Integer sub_object_matrix_2<T>::rows() const
{
    return Matrix(*this).rows();
};

template<class T>
inline Integer sub_object_matrix_2<T>::cols() const
{
    return Matrix(*this).cols();
};

template<class T>
inline Integer sub_object_matrix_2<T>::length() const
{
    return Matrix(*this).length();
};

template<class T>
inline Real sub_object_matrix_2<T>::numel() const
{
    return Matrix(*this).numel();
};

template<class T>
inline bool sub_object_matrix_2<T>::all_finite() const
{
    return Matrix(*this).all_finite();
};

template<class T>
inline sub_object_matrix_2<T>::operator bool() const
{
    return Matrix(*this).operator bool();
};

//-----------------------------------------------------------------
//                  object_matrix
//-----------------------------------------------------------------

template<class T>
inline const Matrix& object_matrix<T>::to_matrix() const &
{ 
    return m_matrix; 
};

template<class T>
inline Matrix&& object_matrix<T>::to_matrix() &&
{ 
    return std::move(m_matrix); 
};

template<class T>
inline object_matrix<T>::operator bool() const
{ 
    return (bool)m_matrix; 
};

template<class T>
inline bool object_matrix<T>::is_empty() const
{ 
    return m_matrix.is_empty(); 
};

template<class T>
inline bool object_matrix<T>::is_scalar() const
{ 
    return m_matrix.is_scalar(); 
};

template<class T>
inline bool object_matrix<T>::is_square() const
{ 
    return m_matrix.is_square(); 
};

template<class T>
inline bool object_matrix<T>::is_vector() const 
{ 
    return m_matrix.is_vector(); 
};

template<class T>
inline bool object_matrix<T>::is_matrix_type() const
{ 
    return m_matrix.is_matrix_type(); 
};

template<class T>
inline bool object_matrix<T>::is_scalar_type() const
{ 
    return m_matrix.is_scalar_type(); 
};

template<class T>
inline bool object_matrix<T>::is_unique() const 
{ 
    return m_matrix.is_unique(); 
}

template<class T>
inline value_code object_matrix<T>::get_value_code() const 
{ 
    return value_code::v_object; 
}

template<class T>
inline struct_code object_matrix<T>::get_struct_code() const 
{ 
    return m_matrix.get_struct_code(); 
}

template<class T>
inline mat_code object_matrix<T>::get_matrix_code() const
{ 
    return m_matrix.get_matrix_code(); 
}

template<class T>
object_matrix<T> object_matrix<T>::zeros(Integer r, Integer c)
{
    return object_matrix<T>(matcl::zeros(object_type::get_static_type(), r, c));
};

template<class T>
object_matrix<T> object_matrix<T>::spzeros(Integer r, Integer c,Integer nnz0)
{
    return object_matrix<T>(matcl::spzeros(object_type::get_static_type(), r, c, nnz0));
};

template<class T>
object_matrix<T> object_matrix<T>::bzeros(Integer r, Integer c,Integer fd, Integer ld)
{
    return object_matrix<T>(matcl::bzeros(object_type::get_static_type(), r, c, fd, ld));
};

template<class T>
object_matrix<T> object_matrix<T>::ones(Integer r, Integer c)
{
    return object_matrix<T>(matcl::ones(object_type::get_static_type(), r, c));
};

template<class T>
object_matrix<T> object_matrix<T>::spones(Integer r, Integer c)
{
    return object_matrix<T>(matcl::spones(object_type::get_static_type(), r, c));
};

template<class T>
object_matrix<T> object_matrix<T>::bones(Integer r, Integer c)
{
    return object_matrix<T>(matcl::bones(object_type::get_static_type(), r, c));
};

template<class T>
object_matrix<T> object_matrix<T>::eye(Integer r, Integer c)
{
    return object_matrix<T>(matcl::eye(object_type::get_static_type(), r, c));
};

template<class T>
object_matrix<T> object_matrix<T>::eye(Integer r)
{
    return object_matrix<T>(matcl::eye(object_type::get_static_type(), r, r));
};

template<class T>
object_matrix<T> object_matrix<T>::speye(Integer r, Integer c)
{
    return object_matrix<T>(matcl::speye(object_type::get_static_type(), r, c));
};

template<class T>
object_matrix<T> object_matrix<T>::speye(Integer r)
{
    return object_matrix<T>(matcl::speye(object_type::get_static_type(), r, r));
};

template<class T>
object_matrix<T> object_matrix<T>::beye(Integer r, Integer c, Integer fd, Integer ld)
{
    return object_matrix<T>(matcl::beye(object_type::get_static_type(), r, c, fd, ld));
};

template<class T>
object_matrix<T> object_matrix<T>::beye(Integer r, Integer fd, Integer ld)
{
    return object_matrix<T>(matcl::beye(object_type::get_static_type(), r, r, fd, ld));
};

template<class T>
object_matrix<T> object_matrix<T>::diag(const object_matrix& v, Integer d)
{
    return object_matrix<T>(matcl::diag(Matrix(v), d));
};

template<class T>
object_matrix<T> object_matrix<T>::spdiag(const object_matrix &v, Integer d)
{
    return object_matrix<T>(matcl::spdiag(Matrix(v), d));
};

template<class T>
object_matrix<T> object_matrix<T>::bdiag(const object_matrix &v, Integer d)
{
    return object_matrix<T>(matcl::bdiag(Matrix(v), d));
};

template<class T>
object_matrix<T> object_matrix<T>::diags(const object_matrix& A, const Matrix &d, Integer r, Integer c)
{
    return object_matrix<T>(matcl::diags(Matrix(A),d,r,c));
};

template<class T>
object_matrix<T> object_matrix<T>::spdiags(const object_matrix &A, const Matrix &d, Integer r, Integer c)
{
    return object_matrix<T>(matcl::spdiags(Matrix(A),d,r,c));
};

template<class T>
object_matrix<T> object_matrix<T>::bdiags(const object_matrix &A, const Matrix &d, Integer r, Integer c)
{
    return object_matrix<T>(matcl::bdiags(Matrix(A),d,r,c));
};

template<class T>
object_matrix<T> object_matrix<T>::make_dense(Integer rows,Integer cols)
{
    return object_matrix<T>(matcl::make_object_dense(object_type::get_static_type(), rows, cols));
};

template<class T>
object_matrix<T> object_matrix<T>::make_dense(Integer rows,Integer cols, const Object *arr)
{
    return object_matrix<T>(matcl::make_object_dense(object_type::get_static_type(), rows, cols, arr));
};

template<class T>
object_matrix<T> object_matrix<T>::make_dense(Integer rows,Integer cols, const Object *arr, Integer ld)
{
    return object_matrix<T>(matcl::make_object_dense(object_type::get_static_type(), rows, cols, arr, ld));
};

template<class T>
object_matrix<T> object_matrix<T>::make_dense_foreign(Integer rows,Integer cols, Object *arr, Integer ld)
{
    return object_matrix<T>(matcl::make_dense_foreign(object_type::get_static_type(), rows, cols, arr, ld));
};

template<class T>
object_matrix<T> object_matrix<T>::make_dense(const T& val, Integer rows, Integer cols)
{
    return object_matrix<T>(matcl::make_object_dense((Object)object_type(val), rows, cols));
};

template<class T>
object_matrix<T> object_matrix<T>::make_dense(const object_type& val, Integer rows, Integer cols)
{
    return object_matrix<T>(matcl::make_object_dense((Object)val, rows, cols));
};

template<class T>
object_matrix<T> object_matrix<T>::make_band(Integer rows,Integer cols, Integer fd, Integer ld)
{
    return object_matrix<T>(matcl::make_object_band((Object)object_type::get_static_type(), 
                                                    rows, cols, fd, ld));
};

template<class T>
object_matrix<T> object_matrix<T>::make_band(const T& val,Integer rows,Integer cols, Integer fd, Integer ld)
{
    return object_matrix<T>(matcl::make_object_band((Object)object_type(val), rows, cols, fd, ld));
};

template<class T>
object_matrix<T> object_matrix<T>::make_band(const object_type& val,Integer rows,Integer cols, 
                     Integer fd, Integer ld)
{
    return object_matrix<T>(matcl::make_object_band((Object)val, rows, cols, fd, ld));
};

template<class T>
object_matrix<T> object_matrix<T>::make_sparse(Integer rows,Integer cols, Integer nzmax)
{
    return object_matrix<T>(matcl::make_object_sparse(object_type::get_static_type(), rows, cols, nzmax));
};

template<class T>
object_matrix<T> object_matrix<T>::make_sparse(const Integer *trip_r, const Integer *trip_c, 
                     const Object* trip_x, Integer r, Integer c, Integer nnz, 
                     Integer nzmax)
{
    return object_matrix<T>(matcl::make_object_sparse(object_type::get_static_type(), trip_r, trip_c,
                                                      trip_x, r, c, nnz, nzmax));
};

template<class T>
object_matrix<T> object_matrix<T>::make_sparse(const Matrix& trip_r, const Matrix& trip_c, 
                     const Matrix& trip_x, Integer r, Integer c, Integer nzmax)
{
    return object_matrix<T>(matcl::make_sparse_matrix(trip_r, trip_c, 
                           Matrix(object_matrix(trip_x)), r, c, nzmax));
};

template<class T>
object_matrix<T> object_matrix<T>::make_dense_noinit(Integer rows,Integer cols, Object*& ptr_data)
{
    return object_matrix<T>(matcl::make_object_dense_noinit(object_type::get_static_type(), rows, cols, ptr_data));
};

template<class T>
object_matrix<T> object_matrix<T>::make_band_noinit(Integer rows,Integer cols, Integer fd, Integer ld)
{
    return object_matrix<T>(matcl::make_object_band_noinit(object_type::get_static_type(), rows, cols, fd, ld));
};

template<class T>
object_matrix<T> object_matrix<T>::make_sparse_noinit(Integer rows,Integer cols, Integer nzmax)
{
    return object_matrix<T>(matcl::make_object_sparse_noinit(object_type::get_static_type(), rows, cols, nzmax));
};


//-----------------------------------------------------------------
//                  IO functions
//-----------------------------------------------------------------
template<class T>
inline void matcl::disp(const object_matrix<T>& m, const disp_stream_ptr& os, const options& opts)
{
    return disp(matcl::Matrix(m), os, opts);
};

template<class T>
inline std::ostream& matcl::operator<<(std::ostream& os, const object_matrix<T>& m)
{
    os << Matrix(m);
    return os;
};

template<class T>
inline std::istream& matcl::operator>>(std::istream& is, object_matrix<T>& m)
{
    Matrix tmp;
    is >> tmp;
    m = object_matrix<T>(tmp);
    return is;
};

template<class T>
inline void matcl::save(oarchive& ar,const object_matrix<T>& mat)
{
    save(ar, matcl::Matrix(mat));
};

template<class T>
inline void matcl::load(iarchive& ar,object_matrix<T>& mat)
{
    Matrix tmp;
    load(ar,tmp);
    mat = object_matrix<T>(tmp);
};

// -----------------------------------------------------------------------
//                 manip functions
// -----------------------------------------------------------------------
template<class T>
inline object_matrix<T> matcl::delrows(const object_matrix<T>& A, const colon& c)
{
    return A.delrows(c);
};

template<class T>
inline object_matrix<T> matcl::delrows(object_matrix<T>&& A, const colon& c)
{
    return std::move(A).delrows(c);
};

template<class T>
inline object_matrix<T> matcl::delcols(const object_matrix<T>& A, const colon& c)
{
    return A.delcols(c);
};

template<class T>
inline object_matrix<T> matcl::delcols(object_matrix<T>&& A, const colon& c)
{
    return std::move(A).delcols(c);
};

template<class T>
inline object_matrix<T> matcl::delrowscols(const object_matrix<T>& A, const colon& c1, const colon& c2)
{
    return A.delrowscols(c1,c2);
};

template<class T>
inline object_matrix<T> matcl::delrowscols(object_matrix<T>&& A, const colon& c1, const colon& c2)
{
    return std::move(A).delrowscols(c1,c2);
};

template<class T>
inline object_matrix<T> matcl::horzcat(const object_matrix<T>& A, const object_matrix<T>& B)
{
    return object_matrix<T>(horzcat(Matrix(A), Matrix(B)));
};

template<class T>
inline object_matrix<T> matcl::horzcat(const std::vector<object_matrix<T>>& mat_list)
{
    mat_row ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return object_matrix<T>(Matrix(ret));
};

template<class T>
inline object_matrix<T> matcl::horzcat(std::initializer_list<object_matrix<T>> mat_list)
{
    mat_row ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return object_matrix<T>(Matrix(ret));
};

template<class T>
inline object_matrix<T> matcl::vertcat(const object_matrix<T>& A, const object_matrix<T>& B)
{
    return object_matrix<T>(vertcat(Matrix(A), Matrix(B)));
};

template<class T>
inline object_matrix<T> matcl::vertcat(const std::vector<object_matrix<T>>& mat_list)
{
    mat_col ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return object_matrix<T>(Matrix(ret));
};

template<class T>
inline object_matrix<T> matcl::vertcat(std::initializer_list<object_matrix<T>> mat_list)
{
    mat_col ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return object_matrix<T>(Matrix(ret));
};

template<class T>
inline object_matrix<T> matcl::blkdiag(const object_matrix<T>& A, const object_matrix<T>& B)
{
    return object_matrix<T>(blkdiag(Matrix(A), Matrix(B)));
};

template<class T>
inline object_matrix<T> matcl::blkdiag(const std::vector<object_matrix<T>>& mat_list)
{
    std::vector<Matrix> vec_mat;
    for( auto it: mat_list)
        vec_mat.push_back(Matrix(it));

    Matrix ret = blkdiag(vec_mat);
    return object_matrix<T>(Matrix(ret));
};

template<class T>
inline object_matrix<T> matcl::blkdiag(std::initializer_list<object_matrix<T>> mat_list)
{
    std::vector<Matrix> vec_mat;
    for( auto it: mat_list)
        vec_mat.push_back(Matrix(it));

    Matrix ret = blkdiag(vec_mat);
    return object_matrix<T>(Matrix(ret));
};

template<class T>
inline object_matrix<T> matcl::repmat(const object_matrix<T>& A, Integer m, Integer n)
{
    return object_matrix<T>(repmat(Matrix(A), m, n));
};

template<class T>
inline object_matrix<T> matcl::sparse(const object_matrix<T>& m)
{
    return object_matrix<T>(sparse(Matrix(m)));
};

template<class T>
inline object_matrix<T> matcl::full(const object_matrix<T>& m)
{
    return object_matrix<T>(full(Matrix(m)));
};

template<class T>
inline object_matrix<T> matcl::band(const object_matrix<T>& m)
{
    return object_matrix<T>(band(Matrix(m)));
};

template<class T>
inline object_matrix<T> matcl::clone(const object_matrix<T>& m)
{
    return object_matrix<T>(clone(Matrix(m)));
};

template<class T>
inline object_matrix<T> matcl::trans(const object_matrix<T>& m)
{
    return object_matrix<T>(trans(Matrix(m)));
};

template<class T>
inline object_matrix<T> matcl::ctrans(const object_matrix<T>& m)
{
    return object_matrix<T>(ctrans(Matrix(m)));
};

template<class T>
inline object_matrix<T> matcl::trans(const object_matrix<T>& m, trans_type t)
{
    return object_matrix<T>(trans(Matrix(m), t));
};

template<class T>
inline object_matrix<T> matcl::trans(const object_matrix<T>& m, trans_type_ext t)
{
    return object_matrix<T>(trans(Matrix(m), t));
};

template<class T>
inline object_matrix<T> matcl::vec(const object_matrix<T>& m)
{
    return object_matrix<T>(vec(Matrix(m)));
};

template<class T>
inline object_matrix<T> matcl::tril(const object_matrix<T>& m, Integer d)
{
    return object_matrix<T>(tril(Matrix(m),d));
};

template<class T>
inline object_matrix<T> matcl::tril(object_matrix<T>&& m, Integer d)
{
    return object_matrix<T>(tril(Matrix(std::move(m)),d));
};

template<class T>
inline object_matrix<T> matcl::triu(const object_matrix<T>& m, Integer d)
{
    return object_matrix<T>(triu(Matrix(m),d));
};

template<class T>
inline object_matrix<T> matcl::triu(object_matrix<T>&& m, Integer d)
{
    return object_matrix<T>(triu(Matrix(std::move(m)),d));
};

template<class T>
inline object_matrix<T> matcl::select_band(const object_matrix<T>& m, Integer ld, Integer ud)
{
    return object_matrix<T>(select_band(Matrix(m), ld, ud));
};

template<class T>
inline object_matrix<T> matcl::flipud(const object_matrix<T>& m)
{
    return object_matrix<T>(flipud(Matrix(m)));
};

template<class T>
inline object_matrix<T> matcl::fliplr(const object_matrix<T>& m)
{
    return object_matrix<T>(fliplr(Matrix(m)));
};

template<class T>
inline object_matrix<T> matcl::reshape(const object_matrix<T>& A, Integer m, Integer n)
{
    return object_matrix<T>(reshape(Matrix(A),m,n));
};

template<class T>
inline object_matrix<T> matcl::get_diag(const object_matrix<T>& m, Integer d)
{
    return object_matrix<T>(get_diag(Matrix(m),d));
};

template<class T>
inline object_matrix<T> matcl::rot90(const object_matrix<T>& m, Integer n)
{
    return object_matrix<T>(rot90(Matrix(m),n));
};

template<class T>
inline dense_matrix<Integer> matcl::find(const object_matrix<T>& v)
{
    return dense_matrix<Integer>(find(Matrix(v)));
};

template<class T>
inline dense_matrix<Integer> matcl::find(const object_matrix<T>& v, const test_function& t)
{
    return dense_matrix<Integer>(find(Matrix(v), t));
};

template<class T>
inline dense_matrix<Integer> matcl::find(const object_matrix<T>& v, const test_object_function<T>& t)
{
    return dense_matrix<Integer>(find(Matrix(v), t.to_test_function()));
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>>
matcl::find2(const object_matrix<T>& v)
{
    Matrix R, C;
    tie(R,C) = find2(Matrix(v));

    dense_matrix<Integer> ret1(R);
    dense_matrix<Integer> ret2(C);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>>
matcl::find2(const object_matrix<T>& v,const test_function& t)
{
    Matrix R, C;
    tie(R,C) = find2(Matrix(v),t);

    dense_matrix<Integer> ret1(R);
    dense_matrix<Integer> ret2(C);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>>
matcl::find2(const object_matrix<T>& v,const test_object_function<T>& t)
{
    Matrix R, C;
    tie(R,C) = find2(Matrix(v),t.to_test_function());

    dense_matrix<Integer> ret1(R);
    dense_matrix<Integer> ret2(C);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>, object_matrix<T>>
matcl::find3(const object_matrix<T>& v)
{
    Matrix R, C, X;
    tie(R,C, X) = find3(Matrix(v));

    dense_matrix<Integer>   ret1(R);
    dense_matrix<Integer>   ret2(C);
    object_matrix<T>        ret3(X);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>, object_matrix<T>>(ret1, ret2, ret3);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>, object_matrix<T>>
matcl::find3(const object_matrix<T>& v, const test_function& t)
{
    Matrix R, C, X;
    tie(R,C, X) = find3(Matrix(v), t);

    dense_matrix<Integer>       ret1(R);
    dense_matrix<Integer>       ret2(C);
    object_matrix<T>            ret3(X);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>, object_matrix<T>>(ret1, ret2, ret3);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>, object_matrix<T>>
matcl::find3(const object_matrix<T>& v, const test_object_function<T>& t)
{
    Matrix R, C, X;
    tie(R,C, X) = find3(Matrix(v), t.to_test_function());

    dense_matrix<Integer>       ret1(R);
    dense_matrix<Integer>       ret2(C);
    object_matrix<T>            ret3(X);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>, object_matrix<T>>(ret1, ret2, ret3);
};

template<class T>
inline object_matrix<T> matcl::sort(const object_matrix<T>& v, int dim, bool asceding)
{
    return object_matrix<T>(sort(Matrix(v), dim, asceding));
};

template<class T>
inline tuple<object_matrix<T>, dense_matrix<Integer>>
matcl::sort2(const object_matrix<T>& v, int dim, bool asceding)
{
    Matrix S, I;
    tie(S,I) = sort2(Matrix(v),dim,asceding);

    object_matrix<T>         ret1(S);
    dense_matrix<Integer>       ret2(I);

    return tuple<object_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline object_matrix<T> matcl::sortrows(const object_matrix<T>& v)
{
    return object_matrix<T>(sortrows(Matrix(v)));
};

template<class T>
inline tuple<object_matrix<T>, dense_matrix<Integer>>
matcl::sortrows2(const object_matrix<T>& v)
{
    Matrix S, I;
    tie(S,I) = sort2(Matrix(v));

    object_matrix<T>        ret1(S);
    dense_matrix<Integer>   ret2(I);

    return tuple<object_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline object_matrix<T> matcl::sortrows(const object_matrix<T>& v, const Matrix& cols)
{
    return object_matrix<T>(sortrows(Matrix(v), cols));
};

template<class T>
inline tuple<object_matrix<T>, dense_matrix<Integer>>
matcl::sortrows2(const object_matrix<T>& v, const Matrix& cols)
{
    Matrix S, I;
    tie(S,I) = sortrows2(Matrix(v),cols);

    object_matrix<T>        ret1(S);
    dense_matrix<Integer>   ret2(I);

    return tuple<object_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
}

template<class T>
inline object_matrix<T> matcl::sortcols(const object_matrix<T>& v)
{
    return object_matrix<T>(sortcols(Matrix(v)));
};

template<class T>
inline tuple<object_matrix<T>, dense_matrix<Integer>>
matcl::sortcols2(const object_matrix<T>& v)
{
    Matrix S, I;
    tie(S,I) = sortcols2(Matrix(v));

    object_matrix<T>         ret1(S);
    dense_matrix<Integer>    ret2(I);

    return tuple<object_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline object_matrix<T> matcl::sortcols(const object_matrix<T>& v, const Matrix& dims)
{
    return object_matrix<T>(sortcols(Matrix(v), dims));
};

template<class T>
inline tuple<object_matrix<T>, dense_matrix<Integer>>
matcl::sortcols2(const object_matrix<T>& v, const Matrix& dims)
{
    Matrix S, I;
    tie(S,I) = sortcols2(Matrix(v),dims);

    object_matrix<T>         ret1(S);
    dense_matrix<Integer>    ret2(I);

    return tuple<object_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline dense_matrix<Integer> matcl::issorted(const object_matrix<T>& v, int dim, bool asceding)
{
    return dense_matrix<Integer>(issorted(Matrix(v),dim,asceding));
};

template<class T>
inline object_matrix<T> matcl::drop_sparse(const object_matrix<T>& A, Real tol)
{
    return object_matrix<T>(drop_sparse(A.to_matrix(), tol));
};

template<class S, class T>
inline object_matrix<S> matcl::convert_object(const object_matrix<T>& A)
{
    return object_matrix<S>(Matrix(A));
};

//--------------------------------------------------------------
//      functions operating on rows or columns
//--------------------------------------------------------------

template<class T>
inline dense_matrix<Integer> matcl::nnz(const object_matrix<T>& A, Integer dim)
{
    return dense_matrix<Integer>(matcl::nnz(matcl::Matrix(A),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::all(const object_matrix<T>& v, int dim)
{
    return dense_matrix<Integer>(matcl::all(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::all(const object_matrix<T>& v,const test_function& t, int dim)
{
    return dense_matrix<Integer>(matcl::all(matcl::Matrix(v),t,dim));
};

template<class T>
inline dense_matrix<Integer> matcl::all(const object_matrix<T>& v, const test_object_function<T>& t, int dim)
{
    return dense_matrix<Integer>(matcl::all(matcl::Matrix(v),t.to_test_function(),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::any(const object_matrix<T>& v, int dim)
{
    return dense_matrix<Integer>(matcl::any(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::any(const object_matrix<T>& v,const test_function& t, int dim)
{
    return dense_matrix<Integer>(matcl::any(matcl::Matrix(v),t,dim));
};

template<class T>
inline dense_matrix<Integer> matcl::any(const object_matrix<T>& v, const test_object_function<T>& t, int dim)
{
    return dense_matrix<Integer>(matcl::any(matcl::Matrix(v),t.to_test_function(),dim));
};

template<class T>
inline object_matrix<T> matcl::sum(const object_matrix<T>& v, int dim)
{
    return object_matrix<T>(matcl::sum(matcl::Matrix(v),dim));
};

template<class T>
inline object_matrix<T> matcl::prod(const object_matrix<T>& v, int dim)
{
    return object_matrix<T>(matcl::prod(matcl::Matrix(v),dim));
};

template<class T>
inline object_matrix<T> matcl::cumsum(const object_matrix<T>& v, int dim)
{
    return object_matrix<T>(matcl::cumsum(matcl::Matrix(v),dim));
};

template<class T>
inline object_matrix<T> matcl::cumprod(const object_matrix<T>& v, int dim)
{
    return object_matrix<T>(matcl::cumprod(matcl::Matrix(v),dim));
};

template<class T>
inline object_matrix<decltype(T()/Integer())>
matcl::mean(const object_matrix<T>& v, int dim)
{
    using ret_type = object_matrix<decltype(T()/Integer())>;
    return ret_type(matcl::mean(matcl::Matrix(v),dim));
};

template<class T>
object_matrix<decltype(decltype(abs(T()))()/Integer())>
inline matcl::std(const object_matrix<T>& v, int dim, bool unbiased)
{
    using ret_type = object_matrix<decltype(decltype(abs(T()))()/Integer())>;
    return ret_type(matcl::std(matcl::Matrix(v),dim,unbiased));
};

template<class T>
inline object_matrix<T> matcl::min_d(const object_matrix<T>& v, int dim)
{
    return object_matrix<T>(matcl::min_d(matcl::Matrix(v),dim));
};

template<class T>
inline object_matrix<T> matcl::max_d(const object_matrix<T>& v, int dim)
{
    return object_matrix<T>(matcl::max_d(matcl::Matrix(v),dim));
};

template<class T>
object_matrix<decltype(abs(T()))>
inline matcl::min_abs_d(const object_matrix<T>& v, int dim)
{
    using ret_type = object_matrix<decltype(abs(T()))>;
    return ret_type(matcl::min_abs_d(matcl::Matrix(v),dim));
};

template<class T>
inline object_matrix<decltype(abs(T()))>
matcl::max_abs_d(const object_matrix<T>& v, int dim)
{
    using ret_type = object_matrix<decltype(abs(T()))>;
    return ret_type(matcl::max_abs_d(matcl::Matrix(v),dim));
};

template<class T>
inline tuple<object_matrix<T>, dense_matrix<Integer>>
matcl::min2(const object_matrix<T>& v, int dim)
{
    using first_type    = object_matrix<T>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::min2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline tuple<object_matrix<T>, dense_matrix<Integer>>
matcl::max2(const object_matrix<T>& v, int dim)
{
    using first_type    = object_matrix<T>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::max2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline tuple<object_matrix<decltype(abs(T()))>, dense_matrix<Integer>>
matcl::min_abs2(const object_matrix<T>& v, int dim)
{
    using first_type    = object_matrix<decltype(abs(T()))>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::min_abs2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline tuple<object_matrix<decltype(abs(T()))>, dense_matrix<Integer>>
matcl::max_abs2(const object_matrix<T>& v, int dim)
{
    using first_type    = object_matrix<decltype(abs(T()))>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::max_abs2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

//--------------------------------------------------------------
//      functions operating on all elements
//--------------------------------------------------------------
template<class T>
inline Integer matcl::nnz_vec(const object_matrix<T>& A)
{
    return matcl::nnz_vec(matcl::Matrix(A));
};

template<class T>
inline bool matcl::all_vec(const object_matrix<T>& v)
{
    return matcl::all_vec(matcl::Matrix(v));
};

template<class T>
inline bool matcl::all_vec(const object_matrix<T>& v,const test_function& t)
{
    return matcl::all_vec(matcl::Matrix(v),t);
};

template<class T>
inline bool matcl::all_vec(const object_matrix<T>& v, const test_object_function<T>& t)
{
    return matcl::all_vec(matcl::Matrix(v),t.to_test_function());
};

template<class T>
inline bool matcl::any_vec(const object_matrix<T>& v)
{
    return matcl::any_vec(matcl::Matrix(v));
};

template<class T>
inline bool matcl::any_vec(const object_matrix<T>& v,const test_function& t)
{
    return matcl::any_vec(matcl::Matrix(v),t);
};

template<class T>
inline bool matcl::any_vec(const object_matrix<T>& v, const test_object_function<T>& t)
{
    return matcl::any_vec(matcl::Matrix(v),t.to_test_function());
};

template<class T>
inline object_type<T> matcl::sum_vec(const object_matrix<T>& v)
{
    return object_type<T>(matcl::sum_vec(matcl::Matrix(v)).get_scalar<Object>(), 
                          dynamic::from_object());
};

template<class T>
inline object_type<T> matcl::prod_vec(const object_matrix<T>& v)
{
    return object_type<T>(matcl::prod_vec(matcl::Matrix(v)).get_scalar<Object>(), 
                          dynamic::from_object());
};

template<class T>
inline object_type<decltype(T()/Integer())>
matcl::mean_vec(const object_matrix<T>& v)
{
    using ret_type = object_type<decltype(T()/Integer())>;
    return ret_type(matcl::mean_vec(matcl::Matrix(v)).get_scalar<Object>(), 
                    dynamic::from_object());
};

template<class T>
inline object_type<decltype(decltype(abs(T()))()/Integer())>
matcl::std_vec(const object_matrix<T>& v, bool unbiased)
{
    using ret_type = object_type<decltype(decltype(abs(T()))()/Integer())>;
    return ret_type(matcl::std_vec(matcl::Matrix(v),unbiased).get_scalar<Object>(), 
                    dynamic::from_object());
};

template<class T>
inline object_type<T> matcl::min_vec(const object_matrix<T>& v)
{
    using ret_type = object_type<T>;
    return ret_type(matcl::min_vec(matcl::Matrix(v)).get_scalar<Object>(), 
                    dynamic::from_object());
};

template<class T>
inline object_type<T> matcl::max_vec(const object_matrix<T>& v)
{
    using ret_type = object_type<T>;
    return ret_type(matcl::max_vec(matcl::Matrix(v)).get_scalar<Object>(), 
                    dynamic::from_object());
};

template<class T>
inline object_type<decltype(abs(T()))>
matcl::min_abs_vec(const object_matrix<T>& v)
{
    using ret_type = object_type<decltype(abs(T()))>;
    return ret_type(matcl::min_abs_vec(matcl::Matrix(v)).get_scalar<Object>(), 
                    dynamic::from_object());
};

template<class T>
inline object_type<decltype(abs(T()))>
matcl::max_abs_vec(const object_matrix<T>& v)
{
    using ret_type = object_type<decltype(abs(T()))>;
    return ret_type(matcl::max_abs_vec(matcl::Matrix(v)).get_scalar<Object>(), 
                    dynamic::from_object());
};

template<class T>
inline tuple<object_type<T>, Integer, Integer >
matcl::min2_vec(const object_matrix<T>& v)
{
    using first_type    = object_type<T>;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::min2_vec(matcl::Matrix(v));

    return ret_type(first_type(ret.get<1>().get_scalar<Object>(), dynamic::from_object()),
                    ret.get<2>(), ret.get<3>());
};

template<class T>
inline tuple<object_type<T>, Integer, Integer >
matcl::max2_vec(const object_matrix<T>& v)
{
    using first_type    = object_type<T>;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::max2_vec(matcl::Matrix(v));

    return ret_type(first_type(ret.get<1>().get_scalar<Object>(), dynamic::from_object()),
                    ret.get<2>(), ret.get<3>());
};

template<class T>
inline tuple<object_type<decltype(abs(T()))>, Integer, Integer >
matcl::min_abs2_vec(const object_matrix<T>& v)
{
    using first_type    = object_type<decltype(abs(T()))>;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::min_abs2_vec(matcl::Matrix(v));

    return ret_type(first_type(ret.get<1>().get_scalar<Object>(), dynamic::from_object()),
                    ret.get<2>(), ret.get<3>());
};

template<class T>
inline tuple<object_type<decltype(abs(T()))>, Integer, Integer >
matcl::max_abs2_vec(const object_matrix<T>& v)
{
    using first_type    = object_type<decltype(abs(T()))>;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::max_abs2_vec(matcl::Matrix(v));

    return ret_type(first_type(ret.get<1>().get_scalar<Object>(), dynamic::from_object()),
                    ret.get<2>(), ret.get<3>());
};

};
