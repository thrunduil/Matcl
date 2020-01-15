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

#include "matcl-internals/container/mat_d.inl"

#include "matcl-internals/container/mat_b.h"
#include "matcl-matrep/details/matrix.inl"
#include "matcl-internals/func/converter.h"
#include "matcl-matrep/func/raw/raw_manip.h"
#include "matcl-internals/base/utils.h"
#include "matcl-internals/base/optim_params.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-core/details/integer.h"
#include "matcl-blas-lapack/level1/level1.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/details/struct_flag_predefined.h"

#pragma warning(push)
#pragma warning(disable:4127) //conditional expression is constant

namespace matcl { namespace raw 
{

namespace gr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

template <class value_type>
Matrix<value_type,struct_dense>::Matrix(tinfo ti)
    : base_type(ti) 
{}

template <class value_type>
Matrix<value_type,struct_dense>::Matrix(tinfo ti,Integer r, Integer c)
    : base_type(ti, r, c) 
{}

template <class value_type>
Matrix<value_type,struct_dense>::Matrix(const Matrix &mat)
    : base_type(mat) 
{}

template <class value_type>
Matrix<value_type,struct_dense>::Matrix(Matrix &&mat)
    : base_type(std::move(mat)) 
{}

template <class value_type>
Matrix<value_type,struct_dense>::Matrix(const base_type &mat)
    : base_type(mat) 
{}

template <class value_type>
Matrix<value_type,struct_dense>::Matrix(base_type &&mat)
    : base_type(std::move(mat)) 
{}

template <class value_type>
Matrix<value_type,struct_dense>::Matrix(tinfo ti,const value_type *arr, Integer r, Integer c)
    : base_type(ti, r, c)
{
    Matrix::read_from(arr, r);	
}

template <class value_type>
Matrix<value_type,struct_dense>::Matrix(tinfo ti,const value_type *arr, Integer r, Integer c, Integer ld)
    : base_type(ti, r, c)
{
    Matrix::read_from(arr, ld);	
}

template <class value_type>
Matrix<value_type,struct_dense>::Matrix(tinfo ti, value_type *arr, Integer r, Integer c, Integer ld, foreign)
    : base_type(ti, r, c, arr, ld)
{}

template <class value_type>
Matrix<value_type,struct_dense>::Matrix(tinfo ti, const value_type &val, Integer r, Integer c)
    : base_type(ti, r, c)
{
    set_to_all(val);
}

#pragma warning(push)
#pragma warning(disable:5037) //an out-of-line definition of a member of a class template
                              // cannot have default arguments

// is it a VS bug?

template<>
template<>
Matrix<Complex,struct_dense>::Matrix(tinfo ti, const Real* xr, const Real* xi, Integer r, Integer c,
        const void*)
: base_type(ti, r, c)
{
    Complex* this_ptr   = m_data.ptr();
    Integer size        = m_data.m_size;

    for(Integer i = 0; i < size; ++i)
        this_ptr[i] = Complex(xr[i],xi[i]);
};

template<>
template<>
Matrix<Float_complex,struct_dense>::Matrix(tinfo ti, const Float* xr, const Float* xi, Integer r, Integer c,
        const void*)
: base_type(ti, r, c)
{
    Float_complex* this_ptr = m_data.ptr();
    Integer size            = m_data.m_size;

    for(Integer i = 0; i < size; ++i)
        this_ptr[i]         = Float_complex(xr[i],xi[i]);
};

#pragma warning(pop)

template <class value_type>
bool Matrix<value_type,struct_dense>::is_same_matrix(const Matrix &m) const
{
    if (this->ptr() == m.ptr() && this->rows() == m.rows() && this->cols() == m.cols())
        return true;
    else
        return false;
}

template<class value_type>
struct set_to_all_impl
{
    static void eval(Matrix<value_type,struct_dense>& mat, const value_type& val)
    {
        value_type* ptr_this    = mat.ptr();
        Integer size            = mat.size();

        level1::set_val<value_type,0>::eval(ptr_this, size, val);
    };
};

template<>
struct set_to_all_impl<Object>
{
    static void eval(Matrix<Object,struct_dense>& mat, const Object& val)
    {
        Object* ptr_this    = mat.ptr();
        Integer size        = mat.size();

        for (Integer i = 0; i < size; ++i)
            details::reset_helper(ptr_this[i], val);
    };
};

template <class value_type>
void Matrix<value_type,struct_dense>::set_to_all(const value_type& val)
{
    return set_to_all_impl<value_type>::eval(*this, val);
};

template <class value_type>
Matrix<value_type,struct_dense> Matrix<value_type,struct_dense>::make_unique(bool keep_bufor) const
{
    if (this->is_unique())
        return *this;
    else
        return this->copy(keep_bufor);
};

template<class value_type>
Matrix<value_type,struct_dense> Matrix<value_type,struct_dense>::copy(bool keep_bufor) const
{
    Integer c   = (keep_bufor == true)? Matrix::max_cols() : Matrix::cols();
    Integer r   = (keep_bufor == true)? Matrix::max_rows() : Matrix::rows();

    Matrix out(Matrix::get_type(), r, c);
    out.m_data.m_rows = Matrix::rows();
    out.m_data.m_cols = Matrix::cols();
    out.m_data.m_size = Matrix::m_data.m_size;
    out.m_data.m_flag = Matrix::m_data.m_flag;

    const value_type* ptr_this  = ptr();
    value_type* ptr_out         = out.ptr();
    Integer this_c              = Matrix::cols();
    Integer this_r              = Matrix::rows();
    Integer out_ld              = out.ld();
    Integer this_ld             = Matrix::ld();

    static const bool Select    = true;
    static const Integer Rows   = 0;
    static const Integer Cols   = 0;
    static const Integer Cont   = 0;

    level1::copy_mat<Select, value_type, value_type, Rows, Cols, Cont>
        ::eval(ptr_out, out_ld, ptr_this, this_ld, this_r, this_c);

    return out;
};

template<class value_type>
Matrix<value_type,struct_dense> Matrix<value_type,struct_dense>::clone(bool keep_bufor) const
{
    Integer c   = (keep_bufor == true)?   Matrix::max_cols() 
                                        : Matrix::cols();
    Integer r   = (keep_bufor == true)?   Matrix::max_rows() 
                                        : Matrix::rows();

    Matrix out(Matrix::get_type(), r, c);
    out.m_data.m_rows = Matrix::rows();
    out.m_data.m_cols = Matrix::cols();
    out.m_data.m_size = Matrix::m_data.m_size;
    out.m_data.m_flag = Matrix::m_data.m_flag;

    const value_type* ptr_this  = ptr();
    value_type* ptr_out         = out.ptr();
    Integer this_c              = Matrix::cols();
    Integer this_r              = Matrix::rows();
    Integer out_ld              = out.ld();
    Integer this_ld             = Matrix::ld();

    if (std::is_same<value_type, Object>::value == true)
    {
        for(Integer j = 0; j < this_c; ++j)
        {
            for(Integer i = 0; i < this_r; ++i)
                mrd::reset_helper(ptr_out[i],mrd::clone_helper<value_type>::eval(ptr_this[i]));

            ptr_out     += out_ld;
            ptr_this    += this_ld;
        };
    }
    else
    {
        static const bool Select  = true;
        static const Integer Rows = 0;
        static const Integer Cols = 0;
        static const Integer Cont = 0;

        level1::copy_mat<Select, value_type,value_type,Rows,Cols,Cont>
            ::eval(ptr_out, out_ld, ptr_this, this_ld, this_r, this_c);
    };

    return out;
};

template<class value_type>
matcl::Matrix Matrix<value_type,struct_dense>::fast_optim() const
{
    if (Matrix::rows() == 1 && Matrix::cols() == 1)
        return *ptr();

    if (Matrix::get_struct().is_diag() == true)
    {
        using BM    = Matrix<value_type,struct_banded>;
        BM tmp      = converter<BM,Matrix>::eval(*this);
        return matcl::Matrix(tmp,false);
    };

    return Matrix::check_change_struct();
};

template<class value_type>
matcl::Matrix Matrix<value_type,struct_dense>::check_change_struct() const
{
    if (Matrix::size() > matcl::optim_params::min_size_to_test_sparsity)
    {
        Real density = Matrix::estim_density();
        if (density < matcl::optim_params::max_sparse_density_min)
        {
            using SparseMatrix  = raw::Matrix<value_type,struct_sparse>;
            return matcl::Matrix(raw::converter<SparseMatrix,Matrix>::eval(*this),false);
        };
    };
    return matcl::Matrix(*this,false);
};

template<class value_type>
Real Matrix<value_type,struct_dense>::estim_density() const
{
    static const Integer n_sample = matcl::optim_params::n_sample_to_test_sparsity;

    Integer step        = Matrix::size()/n_sample;
    step                = (step < 3 ? 3:step);
    const value_type* x = ptr();

    Real nz = 0.5, z = 0.5;

    if (Matrix::ld() == Matrix::rows() || Matrix::rows() == 0)
    {
        Integer size    = Matrix::size();

        for (Integer i = 0; i < size; i += step)
        {
            if (mrd::is_zero(*x) == false)
                ++nz;
            else
                ++z;

            x += step;
        };
    }
    else
    {
        Integer size    = Matrix::size();
        Integer rows    = Matrix::rows();
        Integer ld      = Matrix::ld();

        for (Integer i = 0, pos = 0; i < size; i += step, pos += step)
        {
            while(pos >= rows)
            {
                pos -= rows;
                x   += ld;
            };

            if (mrd::is_zero(x[pos]) == false)
                ++nz;
            else
                ++z;
        };
    };

    return nz/(nz+z);
};

template<class value_type>
Matrix<value_type,struct_banded>
Matrix<value_type,struct_dense>::get_diag_band() const
{
    using BandMatrix = Matrix<value_type,struct_banded>;
    BandMatrix out(Matrix::get_type(), Matrix::m_data.m_rows, Matrix::m_data.m_cols, 0, 0);

    Integer s       = std::min(Matrix::m_data.m_cols, Matrix::m_data.m_rows);
    Integer out_ld  = out.ld();
    Integer this_ld = Matrix::ld();

    value_type * ptr            = out.rep_ptr();
    const value_type * ptr_this = Matrix::ptr();

    for (Integer i = 0; i < s; ++i)
    {
        mrd::reset_helper(*ptr,*ptr_this);

        ptr         += out_ld;
        ptr_this    += this_ld + 1;
    };

    out.get_struct() = this->get_struct();
    return out;
};

template<class value_type>
Matrix<value_type,struct_dense> Matrix<value_type,struct_dense>::make_explicit() const
{
    if (Matrix::ld() == Matrix::rows() || Matrix::rows() == 0 || Matrix::cols() <= 1)
        return *this;

    Matrix out(Matrix::get_type(),Matrix::rows(), Matrix::cols());

    value_type* ptr_out         = out.ptr();
    const value_type* ptr_this  = ptr();

    static const bool Select    = true;
    static const Integer Rows   = 0;
    static const Integer Cols   = 0;
    static const Integer Cont   = 0;

    level1::copy_mat<Select, value_type,value_type,Rows,Cols,Cont>
        ::eval(ptr_out, out.ld(), ptr_this, Matrix::ld(), Matrix::rows(), Matrix::cols());

    return out;
};

template<class value_type>
bool Matrix<value_type,struct_dense>::is_explicit() const
{
    if (Matrix::ld() == Matrix::rows() || Matrix::rows() == 0 || Matrix::cols() <= 1)
        return true;

    return false;
};

template <class value_type>
Matrix<value_type,struct_dense>
Matrix<value_type,struct_dense>::get_diag(Integer d) const
{
    error::check_diag(d, base_type::m_data.m_rows, base_type::m_data.m_cols);

    Integer r   = base_type::m_data.m_rows;
    Integer c   = base_type::m_data.m_cols;
    Integer ld  = Matrix::ld();

    Integer st, s;
    if (d >= 0)
    {
        st  = imult(d,ld);
        s   = (r + d >= c) ? c - d : r;
    }
    else
    {
        st  = - d;
        s   = (r + d >= c) ? c : r + d;
    }

    Matrix res(Matrix::get_type(), s, 1);

    if (s == 0)
        return res;

    value_type * ptr            = res.ptr();
    const value_type * this_ptr = base_type::ptr() + st;

    for (Integer i = 0; i < s; ++i)
    {
        mrd::reset_helper(ptr[i],*this_ptr);
        this_ptr    += ld + 1;
    }

    return res;
};

template<class value_type>
Matrix<value_type,struct_dense> Matrix<value_type,struct_dense>::reshape(Integer r, Integer c) const
{
    error::check_reshape(Matrix::rows(), Matrix::cols(), r, c);

    if (r == Matrix::rows() && c == Matrix::cols())
        return *this;

    if (Matrix::ld() == Matrix::rows() || Matrix::rows() == 0)
    {
        Matrix out(*this);    
        out.m_data.m_rows   = r;
        out.m_data.m_cols   = c;
        out.m_data.m_ld     = std::max(r,1);

        if (Matrix::rows() == r && Matrix::cols() == c)
            out.get_struct() = Matrix::get_struct();
        else
            out.get_struct().reset();

        return out;
    }
    else
    {
        return this->copy().reshape(r,c);
    };    
};

template<class value_type>
Matrix<value_type,struct_dense> 
Matrix<value_type,struct_dense>::reserve(Integer r, Integer c) const
{
    if (r <= Matrix::m_data.m_max_rows && c <= Matrix::m_data.m_max_cols)
        return *this;

    Integer m_c = Matrix::cols();
    Integer m_r = Matrix::rows();

    Matrix out(Matrix::get_type(), std::max(r,Matrix::max_rows()), std::max(c,Matrix::max_cols()));

    out.m_data.m_rows   = m_r;
    out.m_data.m_cols   = m_c;
    out.m_data.m_size   = imult(m_r,m_c);

    out.set_struct(Matrix::get_struct());

    const value_type* this_ptr  = ptr();
    value_type* ptr             = out.ptr();

    static const bool Select    = true;
    static const Integer Rows   = 0;
    static const Integer Cols   = 0;
    static const Integer Cont   = 0;

    level1::copy_mat<Select, value_type,value_type,Rows,Cols,Cont>
        ::eval(ptr, out.ld(), this_ptr, this->ld(), m_r, m_c);

    return out;
};

template<class value_type>
void Matrix<value_type,struct_dense>::prepare_for_concat(Integer r, Integer c)
{
    if (r <= Matrix::m_data.m_max_rows && c <= Matrix::m_data.m_max_cols)
    {
        m_data.m_rows   = r;
        m_data.m_cols   = c;
        m_data.m_size   = imult(r,c);
        m_data.get_struct().reset();
        return;
    };

    Integer new_r = (r == this->rows()) ? r : std::max(this->rows() * 2, r);
    Integer new_c = (c == this->cols()) ? c : std::max(this->cols() * 2, c);

    Matrix out = this->reserve(new_r, new_c);
    
    out.m_data.m_rows = r;
    out.m_data.m_cols = c;
    out.m_data.m_size = imult(r,c);
    out.m_data.get_struct().reset();

    this->assign_to_fresh(out);
};

template<class value_type>
Matrix<value_type,struct_dense> Matrix<value_type,struct_dense>::resize(Integer r, Integer c) const
{
    error::check_resize(r,c);

    if (r <= Matrix::rows() && c <= Matrix::cols())
        return make_view(1,r,1,c);

    if (r <= Matrix::m_data.m_max_rows && c <= Matrix::m_data.m_max_cols)
    {
        Matrix out = copy(true);
        return out.resize(r,c);
    };

    Matrix out = reserve(r,c);
    return out.resize(r,c);
};

template<class value_type>
Matrix<value_type,struct_dense> Matrix<value_type,struct_dense>::resize(Integer r, Integer c)
{
    error::check_resize(r,c);

    if (r <= Matrix::rows() && c <= Matrix::cols())
        return make_view(1, r, 1, c);

    Matrix out = reserve(r,c);

    value_type* ptr = out.ptr();

    Integer m_c     = Matrix::cols();
    Integer m_r     = Matrix::rows();
    Integer out_ld  = out.ld();

    value_type Z    = matcl::details::default_value<value_type>(Matrix::get_type());

    if (m_r < r)
    {
        for (Integer j = 0; j < m_c; ++j)
        {
            level1::set_val<value_type,0>::eval(ptr + m_r, r-m_r, Z);

            ptr += out_ld;
        }; 
    }
    else
    {
        ptr += imult(m_c,out.ld());
    };  

    for (Integer j = m_c; j < c; ++j)
    {
        level1::set_val<value_type,0>::eval(ptr, r, Z);
        ptr += out_ld;
    }; 

    out.m_data.m_rows = r;
    out.m_data.m_cols = c;
    out.m_data.m_size = imult(r,c);

    bool is_sym = (r == c);
    out.m_data.get_struct() = md::predefined_struct::get_resize(out.m_data.get_struct(), is_sym,
                                                                is_real_matrix(out));

    return out;
};

template<class value_type>
Matrix<value_type,struct_dense> 
Matrix<value_type,struct_dense>::make_view(Integer r_start, Integer r_end) const
{
    Integer rows        = Matrix::rows();
    Integer cols        = Matrix::cols();
    Integer ld          = Matrix::ld();

    bool can_view       = ld == rows || rows == 0 || cols <= 1;

    if (can_view == false)
    {
        //view can be created if first and last index is in the same column
        Integer c1      = (r_start - 1) / rows;
        Integer c2      = (r_end - 1) / rows;

        if (c1 == c2)
            can_view    = true;
    };

    if (can_view == true)
    {
        Matrix out(*this);    

        out.m_data.m_rows = r_end-r_start+1;
        out.m_data.m_cols = 1;
        out.m_data.m_size = out.m_data.m_rows;
        
        out.get_struct().reset();

        if (Matrix::rows() != 0)
        {
            Integer pos = r_start - 1;
            Integer c = pos / rows;
            Integer r = pos % rows;
            out.m_data.m_ptr = const_cast<value_type*>(ptr() + r + c*ld);
        };

        return out;
    };

    Integer r = r_end-r_start+1;
    Matrix out(Matrix::get_type(),r,1);

    value_type* ptr_out         = out.ptr();
    const value_type* ptr_this  = ptr();

    for (Integer i = 0, pos = r_start - 1; i < r; ++i, ++pos)
    {
        while (pos >= rows)
        {
            pos         -= rows;
            ptr_this    += ld;
        };

        mrd::reset_helper(ptr_out[i],ptr_this[pos]);
    };

    return out;
};

template<class value_type>
Matrix<value_type,struct_dense>
Matrix<value_type,struct_dense>::make_view(Integer r_start, Integer r_end, 
                                           Integer c_start, Integer c_end) const
{
    Matrix out(*this);    
    out.m_data.m_rows   = r_end-r_start+1;
    out.m_data.m_cols   = c_end-c_start+1;
    out.m_data.m_size   = imult(out.m_data.m_rows,out.m_data.m_cols);
    out.m_data.m_ptr    = const_cast<value_type*>(ptr() + r_start - 1 + (c_start - 1)*Matrix::ld());

    if (r_start == c_start)
    {
        bool is_sym = (r_end == c_end);
        out.m_data.get_struct() = md::predefined_struct::get_rectangle_view(out.m_data.get_struct(), is_sym);
    }
    else
    {
        out.m_data.get_struct().reset();
    }

    return out;
};

template<class value_type>
bool Matrix<value_type,struct_dense>::all_finite() const
{
    return mr::all_finite_helper<Matrix>::eval(*this);
};

template<class val_type>
void Matrix<val_type,struct_dense>::serialize(oarchive_impl & ar, const unsigned int ) const
{
    ar << boost::serialization::base_object<base_type>(const_cast<Matrix&>(*this));
};

template<class val_type>
void Matrix<val_type,struct_dense>::serialize(iarchive_impl & ar, const unsigned int )
{
    ar >> boost::serialization::base_object<base_type>(*this);
};

template<class value_type>
void Matrix<value_type,struct_dense>::destroy_data()
{
    Matrix<value_type,struct_dense>::m_data.destroy_data();
};

template<class value_type>
void Matrix<value_type,struct_dense>::assign_to_fresh(const Matrix& mat)
{
    base_type::assign_to_fresh(mat);
};

template<class value_type>
void Matrix<value_type,struct_dense>::assign_to_fresh(Matrix&& mat)
{
    base_type::assign_to_fresh(std::move(mat));
};

};};

template class matcl::raw::Matrix<matcl::Integer,matcl::struct_dense>;
template class matcl::raw::Matrix<matcl::Float,matcl::struct_dense>;
template class matcl::raw::Matrix<matcl::Real,matcl::struct_dense>;
template class matcl::raw::Matrix<matcl::Float_complex,matcl::struct_dense>;
template class matcl::raw::Matrix<matcl::Complex,matcl::struct_dense>;
template class matcl::raw::Matrix<matcl::Object,matcl::struct_dense>;

#pragma warning(pop)
