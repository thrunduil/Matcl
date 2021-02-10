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

#include "matcl-internals/container/mat_d.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/error/error_check_basic.h"
#include "matcl-matrep/details/matrix.inl"
#include "matcl-matrep/func/raw/raw_manip.h"
#include "matcl-internals/base/utils.h"
#include "matcl-internals/base/optim_params.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-matrep/details/struct_flag_predefined.h"

namespace matcl { namespace raw 
{

namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

template <class value_type>
void sparse_matrix_base<value_type>::destroy_data()
{
    m_data.destroy_data();
};

template <class value_type>
sparse_matrix_base<value_type>::sparse_matrix_base(type_info ti, Integer r, Integer c)
    :m_data(ti,r, c, 0)
{}

template <class value_type>
sparse_matrix_base<value_type>::sparse_matrix_base(type_info ti, Integer r, Integer c, Integer n)
    :m_data(ti,r, c, n)
{}

template <class value_type>
sparse_matrix_base<value_type>::sparse_matrix_base(type_info ti, const Integer *ri, const Integer *ci,
            const value_type *xv, Integer r, Integer c, Integer nnz)
    :m_data(ti)
{
    construct(ri, ci, xv, r, c, nnz, nnz);
}

template <class value_type>
sparse_matrix_base<value_type>::sparse_matrix_base(type_info ti, const Integer *ri, const Integer *ci,
            const value_type *xv, Integer r, Integer c, Integer nnz, Integer nzmax)
    :m_data(ti)
{
    construct(ri, ci, xv, r, c, nnz, nzmax);
}

template <class value_type>
sparse_matrix_base<value_type>::sparse_matrix_base(const sparse_matrix_base &m)
    :m_data(m.m_data)
{}

template <class value_type>
sparse_matrix_base<value_type>::sparse_matrix_base(sparse_matrix_base &&m)
    :m_data(std::move(m.m_data))
{}

template <class value_type>
void sparse_matrix_base<value_type>::construct(const Integer * ri, const Integer *ci,
        const value_type *xv, Integer r, Integer c, Integer nnz, Integer nzmax)
{
    sparse_ccs d(get_type(), r, c, std::max(nnz,nzmax));

    if ((nzmax == 0) || (nnz == 0))
    {
        m_data.assign_to_fresh(d);
        return;
    }

    matcl::pod_workspace<Integer> p(c + 1, 0);

    const Integer * ptr_ri  = ri;
    const Integer * ptr_ci  = ci;
    Integer * ptr_p         = &p[0];

    for (Integer i = 0; i < nnz; ++i)
    {
        error::check_index(*ptr_ri, *ptr_ci, r, c);
        ++ptr_p[*ptr_ci];

        ++ptr_ci;
        ++ptr_ri;
    }

    ptr_p = &p[1];
    Integer * d_c   = d.ptr_c() + 1;
    Integer * d_r   = d.ptr_r();
    value_type* d_x = d.ptr_x();

    for (Integer k = 0, i = 1; i <= c; ++i)
    {
        *ptr_p      += k;
        k	        = *ptr_p;
        *(d_c++)    = *(ptr_p++);
    }

    ptr_ci = ci;
    ptr_ri = ri;
    ptr_p = &p[0];
    const value_type * ptr_xv = xv;
    for (Integer i = 0; i < nnz; ++i)
    {
        Integer k   = ptr_p[*ptr_ci - 1]++;
        d_r[k]      = *ptr_ri - 1;
        mrd::reset_helper(d_x[k],*ptr_xv);

        ++ptr_ci;
        ++ptr_ri;
        ++ptr_xv;
    }

    d.sort();
    d.add_duplications();
    m_data.assign_to_fresh(d);
}

template <class value_type>
template<class RealType>
void sparse_matrix_base<value_type>::construct2(const Integer *, const Integer *,
        const RealType *, const RealType* , Integer , Integer , Integer , Integer )
{};

template <>
template<class RealType>
void sparse_matrix_base<Complex>::construct2(const Integer *ri, const Integer *ci,
        const RealType *xr, const RealType* xi, Integer r, Integer c, Integer nnz, Integer nzmax)
{
    sparse_ccs d(get_type(), r, c, std::max(nnz,nzmax));

    if ((nzmax == 0) || (nnz == 0))
    {
        m_data.assign_to_fresh(d);
        return;
    }
    const Integer * ptr_ri = ri;
    const Integer * ptr_ci = ci;
    matcl::pod_workspace<Integer> p(c + 1, 0);

    Integer * ptr_p = &p[0];
    for (Integer i = 0; i < nnz; ++i)
    {
        error::check_index(*ptr_ri, *ptr_ci, r, c);
        ++ptr_p[*ptr_ci];

        ++ptr_ci;
        ++ptr_ri;
    }

    ptr_p = &p[1];
    Integer * d_c = d.ptr_c() + 1;
    Integer * d_r = d.ptr_r();
    Complex * d_x = d.ptr_x();

    for (Integer k = 0, i = 1; i <= c; ++i)
    {
        *ptr_p      += k;
        k	        = *ptr_p;
        *(d_c++)    = *(ptr_p++);
    }

    ptr_ci = ci;
    ptr_ri = ri;
    const Real * ptr_xr = xr;
    const Real * ptr_xi = xi;
    ptr_p = &p[0];

    for (Integer i = 0; i < nnz; ++i)
    {
        Integer k = ptr_p[*ptr_ci - 1]++;
        d_r[k] = *ptr_ri - 1;
        d_x[k] = Complex(*ptr_xr,*ptr_xi);

        ++ptr_ci;
        ++ptr_ri;
        ++ptr_xr;
        ++ptr_xi;
    }

    d.sort();
    d.add_duplications();
    m_data.assign_to_fresh(d);
}

template <>
template<class RealType>
void sparse_matrix_base<Float_complex>::construct2(const Integer *ri, const Integer *ci,
        const RealType *xr, const RealType* xi, Integer r, Integer c, Integer nnz, Integer nzmax)
{
    sparse_ccs d(get_type(), r, c, std::max(nnz,nzmax));

    if ((nzmax == 0) || (nnz == 0))
    {
        m_data.assign_to_fresh(d);
        return;
    }

    const Integer * ptr_ri = ri;
    const Integer * ptr_ci = ci;
    matcl::pod_workspace<Integer> p(c + 1, 0);

    Integer * ptr_p = &p[0];
    for (Integer i = 0; i < nnz; ++i)
    {
        error::check_index(*ptr_ri, *ptr_ci, r, c);
        ++ptr_p[*ptr_ci];

        ++ptr_ci;
        ++ptr_ri;
    }

    ptr_p = &p[1];
    Integer * d_c = d.ptr_c() + 1;
    Integer * d_r = d.ptr_r();
    Float_complex * d_x = d.ptr_x();

    for (Integer k = 0, i = 1; i <= c; ++i)
    {
        *ptr_p      += k;
        k	        = *ptr_p;
        *(d_c++)    = *(ptr_p++);
    }

    ptr_ci = ci;
    ptr_ri = ri;
    const Float * ptr_xr = xr;
    const Float * ptr_xi = xi;
    ptr_p = &p[0];

    for (Integer i = 0; i < nnz; ++i)
    {
        Integer k = ptr_p[*ptr_ci - 1]++;
        d_r[k] = *ptr_ri - 1;
        d_x[k] = Float_complex(*ptr_xr,*ptr_xi);

        ++ptr_ci;
        ++ptr_ri;
        ++ptr_xr;
        ++ptr_xi;
    }

    d.sort();
    d.add_duplications();
    m_data.assign_to_fresh(d);
}

template <class value_type>
void sparse_matrix_base<value_type>::assign_to_fresh(const sparse_matrix_base &m)
{
    m_data.assign_to_fresh(m.m_data);
}

template <class value_type>
void sparse_matrix_base<value_type>::assign_to_fresh(sparse_matrix_base &&m)
{
    m_data.assign_to_fresh(std::move(m.m_data));
}

template <class value_type>
sparse_matrix_base<value_type> 
sparse_matrix_base<value_type>::copy(bool keep_maxcol) const
{
    return sparse_matrix_base(m_data.copy(keep_maxcol));
};

template <class value_type>
typename sparse_matrix_base<value_type>::DenseMatrix
sparse_matrix_base<value_type>::full() const
{
    Integer r = rows(), c = cols();
    value_type Z = matcl::details::default_value<value_type>(get_type());
    DenseMatrix tmp(get_type(), Z, r, c);

    if (nnz() == 0)
        return tmp;

    value_type* ptr_res     = tmp.ptr();
    const Integer* ptr_c    = m_data.ptr_c();	
    const Integer* ptr_r    = m_data.ptr_r() + m_data.offset();	
    const value_type* ptr_x = m_data.ptr_x() + m_data.offset();

    for (Integer j = 0; j < c; ++j)
    {		
        Integer dr = ptr_c[1] - ptr_c[0];
        ++ptr_c;

        for (Integer k = 0; k < dr; ++k)
            mrd::reset_helper(ptr_res[*(ptr_r++)],*(ptr_x++));			

        ptr_res += r;
    }

    tmp.set_struct(get_struct());
    return tmp;
}

template <class value_type>
typename sparse_matrix_base<value_type>::SparseMatrix
sparse_matrix_base<value_type>::get_diag(Integer d) const
{
    Integer r = rows(), c = cols();
    error::check_diag(d, r, c);

    Integer s, c_first, c_last, r_first;
    if (d > 0)
    {
        s = std::min(c - d, r);
        c_first = d;
        c_last = s+d;
        r_first = 0;
    }
    else
    {
        s = std::min(r + d, c);
        c_first = 0;
        c_last = s;
        r_first = -d;
    };

    SparseMatrix res(get_type(),s,1,s);

    if (s == 0)
        return res;

    const sparse_ccs& Ad = rep();
    sparse_ccs& rd = res.rep();

    const value_type* Ad_x	= Ad.ptr_x();

    Integer* d_c			= rd.ptr_c();
    Integer* d_r			= rd.ptr_r();
    value_type* d_x			= rd.ptr_x();

    d_c[0]					= 0;

    Integer nz = 0, p = 0;

    for (Integer j = c_first, pos = 0; j < c_last; ++j, ++pos, ++r_first)
    {
        if (Ad.has_element(r_first,j,p))
        {
            *(d_r++)			= pos;
            mrd::reset_helper(*(d_x++),Ad_x[p]);

            ++nz;
        };
    };

    d_c[1] = nz;

    rd.add_memory(-1);
    return res;
};

template<class val_type>
matcl::Matrix fast_optim_impl(const Matrix<val_type,struct_sparse>& A)
{
    if (A.rows() == 1 && A.cols() == 1)
        return matcl::Matrix(A(1,1));

    if (A.get_struct().is_diag() == true)
        return matcl::Matrix(A.get_diag_band(),false);

    if (Real(A.nnz())/(A.rows()+1)/(A.cols()+1) > optim_params::max_sparse_density_max)
        return matcl::Matrix(A.full(),false);

    return matcl::Matrix(A,false);
};

template<class val_type>
typename sparse_matrix_base<val_type>::BandMatrix 
sparse_matrix_base<val_type>::get_diag_band() const
{
    Integer m_rows = rows(), m_cols = cols();
    BandMatrix out(get_type(), m_rows, m_cols, 0, 0);

    Integer s = std::min(m_cols, m_rows);

    value_type * ptr = out.rep_ptr();

    const sparse_ccs& Ad = rep();
    const value_type* Ad_x	= Ad.ptr_x();


    Integer p;
    value_type Z = matcl::details::default_value<val_type>(get_type());

    for (Integer i = 0; i < s; ++i)
    {
        if (Ad.has_element(i,i,p))
            mrd::reset_helper(*ptr,Ad_x[p]);
        else
            mrd::reset_helper(*ptr,Z);

        //out.ld() == 1
        ++ptr;
    };

    out.get_struct() = this->get_struct();
    return out;
};

template<class val_type>
const sparse_matrix_base<val_type>
sparse_matrix_base<val_type>::reserve(Integer, Integer c) const
{
    if (c <= max_cols())
        return *this;

    if (rep().ptr_c() == 0)
    {
        sparse_ccs out      = rep();
        out.m_max_cols      = c;
        return sparse_matrix_base(out);
    };

    Integer m_c = cols(), m_nz = nnz();

    sparse_matrix_base out(get_type(), rows(), c, m_nz);

    out.m_data.m_cols       = m_c;
    out.set_struct(get_struct());

    sparse_ccs& rep              = out.rep();
    const sparse_ccs& this_rep   = this->rep();
    
    Integer* d_c            = rep.ptr_c();
    Integer* d_r            = rep.ptr_r();
    val_type* d_x           = rep.ptr_x();

    Integer off             = this_rep.offset();
    const Integer* Ad_c     = this_rep.ptr_c();
    const Integer* Ad_r     = this_rep.ptr_r() + off;
    const val_type* Ad_x    = this_rep.ptr_x() + off;

    for (Integer i = 0; i <= m_c; ++i)
        d_c[i]  = Ad_c[i] - off;

    for (Integer i = 0; i < m_nz; ++i)
    {
        d_r[i]  = Ad_r[i];
        mrd::reset_helper(d_x[i],Ad_x[i]);
    };

    return out;
};

template<class val_type>
void sparse_matrix_base<val_type>::prepare_for_concat(Integer r, Integer c, Integer nnz)
{
    if (this->m_data.m_c == nullptr)
    {
        this->m_data.assign_to_fresh(sparse_ccs(this->get_type(), r, c, nnz));
        return;
    };

    if (nnz <= this->m_data.nzmax() && c <= this->m_data.max_cols())
    {
        //enough place, only increase size
        Integer old_c = this->cols();

        this->m_data.m_rows = std::max(r, this->m_data.m_rows);
        this->m_data.m_cols = c;

        Integer* ptr_c = this->m_data.ptr_c();
        for (Integer i = old_c + 1; i <= c; ++i)
            ptr_c[i] = ptr_c[i-1];

        return;
    };

    Integer old_c = this->cols();

    if (nnz > this->m_data.nzmax())
    {
        //m_r, m_x must be expanded
        this->m_data.add_memory(this->m_data.nzmax() + nnz);
    };

    if (c > this->m_data.max_cols())
    {
        //m_c must be expanded
        this->m_data.add_columns(this->m_data.cols() + c);
    };

    //update info and column pointers
    this->m_data.m_rows = std::max(r, this->m_data.m_rows);
    this->m_data.m_cols = c;

    Integer* ptr_c = this->m_data.ptr_c();
    for (Integer i = old_c + 1; i <= c; ++i)
        ptr_c[i] = ptr_c[i-1];

    return;
};

template<class val_type>
sparse_matrix_base<val_type> sparse_matrix_base<val_type>::resize(Integer r, Integer c)
{
    error::check_resize(r,c);

    if (r < rows())
        return resize_remrows(r,c);

    if (r == rows() && c == cols())
        return *this;

    if (c <= cols())
    {
        sparse_ccs out  = rep();
        out.m_rows      = r;
        out.m_cols      = c;
        out.get_struct()= md::predefined_struct::get_rectangle_view(out.get_struct(), false);
        return sparse_matrix_base(out);
    };

    sparse_matrix_base out = reserve(r,c);

    Integer m_c = cols();

    out.m_data.m_rows = r;
    out.m_data.m_cols = c;

    Integer* d_c = out.rep().ptr_c();

    if (d_c)
    {
        Integer nz = d_c[m_c];

        for (Integer i = m_c + 1; i <= c; ++i)
            d_c[i] = nz;
    };

    bool is_sym = (r == c);
    out.m_data.get_struct() = md::predefined_struct::get_resize(out.m_data.get_struct(), is_sym,
                                                                is_real_matrix(out));

    return out;
}

template<class val_type>
const sparse_matrix_base<val_type>
sparse_matrix_base<val_type>::resize(Integer r, Integer c) const
{
    error::check_resize(r,c);

    if (r == rows() && c == cols())
        return *this;

    if (r < rows() && nnz() > 0)
        return resize_remrows(r,c);

    if (c <= cols())
    {
        sparse_ccs out  = rep();
        out.m_rows      = r;
        out.m_cols      = c;
        out.get_struct()= md::predefined_struct::get_rectangle_view(out.get_struct(), false);
        return sparse_matrix_base(out);
    };

    if (c <= max_cols())
    {
        sparse_matrix_base out = copy(true);
        return out.resize(r,c);
    };

    sparse_matrix_base out = reserve(r,c);
    return out.resize(r,c);
};

template<class val_type>
void sparse_matrix_base<val_type>::change_number_rows(Integer r)
{
    m_data.m_rows      = r;
    m_data.get_struct()= md::predefined_struct::get_rectangle_view(m_data.get_struct(), false);
}

template<class val_type>
void sparse_matrix_base<val_type>::change_number_cols(Integer c)
{
    m_data.m_cols      = c;
    m_data.get_struct()= md::predefined_struct::get_rectangle_view(m_data.get_struct(), false);
}

template<class val_type>
sparse_matrix_base<val_type> 
sparse_matrix_base<val_type>::resize_remrows(Integer r, Integer c) const
{
    sparse_matrix_base out(get_type(), r, c, nnz());

    sparse_ccs& rep              = out.rep();
    const sparse_ccs& this_rep   = this->rep();

    Integer* d_c            = rep.ptr_c();
    Integer* d_r            = rep.ptr_r();
    val_type* d_x           = rep.ptr_x();

    const Integer* Ad_c     = this_rep.ptr_c();
    const Integer* Ad_r     = this_rep.ptr_r();
    const val_type* Ad_x    = this_rep.ptr_x();

    Integer mc              = std::min(c,cols());
    Integer nz              = 0;
    for (Integer j = 0; j < mc; ++j)
    {
        d_c[j]              = nz;
        for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
        {
            Integer pos     = Ad_r[k];

            if (pos >= r)
                continue;

            d_r[nz]         = pos;
            mrd::reset_helper(d_x[nz],Ad_x[k]);
            ++nz;
        };
    };

    for (Integer i = mc; i <= c; ++i)
        d_c[i]              = nz;

    rep.add_memory(-1);
    bool is_sym = (r == c);
    rep.get_struct() = md::predefined_struct::get_resize(this->get_struct(), is_sym, is_real_matrix(rep));

    return out;
}

template<class val_type>
const sparse_matrix_base<val_type>
sparse_matrix_base<val_type>::make_view(Integer c_start, Integer c_end) const
{
    if (c_start == 1 && c_end == cols())
        return *this;

    sparse_matrix_base out(rep());
    out.rep().m_cols        = c_end - c_start + 1;
    out.rep().m_max_cols    = rep().m_cols;   
    if (rep().m_c)
    {
        out.rep().m_c       = rep().m_c + c_start - 1;
        out.rep().m_offset  = out.rep().m_c[0];
    };

    if (c_start == 1)
        out.get_struct()    = md::predefined_struct::get_rectangle_view(get_struct(), false);
    else
        out.get_struct().reset();

    return out;
};

template <class value_type>
sparse_matrix_base<value_type> sparse_matrix_base<value_type>::clone(bool keep_maxcol) const
{
    return sparse_matrix_base(m_data.clone(keep_maxcol));
};

template <class value_type>
void sparse_matrix_base<value_type>::serialize(oarchive_impl & ar, const unsigned int ) const
{
    ar << m_data;
};

template <class value_type>
void sparse_matrix_base<value_type>::serialize(iarchive_impl & ar, const unsigned int )
{
    ar >> m_data;
};

template<class value_type_>
matcl::Matrix Matrix<value_type_,struct_sparse>::fast_optim() const
{
    return fast_optim_impl(*this);
};

template<class V, class W>
struct test_drop
{
    V       m_val;

    test_drop(V v)              : m_val(v) {};

    bool eval(const W& elem)
    { 
        auto aelem  = mrd::abs_helper<W>::eval(elem);
        using AW    = decltype(aelem);
        return (bool)mrd::leq_helper<AW,V>::eval(aelem, m_val); 
    };
};

template<class V>
struct test_zero
{
    bool eval(const V& elem) { return mrd::is_zero(elem); };
};


template<class val_type, class test>
static const_matrix<Matrix<val_type,struct_sparse>>
drop_impl(const Matrix<val_type,struct_sparse>& A, test t)
{
    using sparse_ccs = details::sparse_ccs<val_type>;

    const sparse_ccs& d = A.rep();
    Integer off         = d.offset();
    const Integer* d_c	= d.ptr_c();
    const Integer* dro	= d.ptr_r() + off;
    const Integer* dr	= d.ptr_r();
    const val_type* dxo	= d.ptr_x() + off;        
    const val_type* dx	= d.ptr_x();

    Integer n           = d.cols();
    Integer nz          = d.nnz();

    Integer pos;
    for (pos = 0; pos < nz; ++pos)
    {
        if (t.eval(dxo[pos]))
            break;
    };

    if (pos == nz)
        return A;

    Matrix<val_type,struct_sparse> ret(A.get_type(),A.rows(),A.cols(),nz);

    Integer* ret_c      = ret.rep().ptr_c();
    Integer* ret_r      = ret.rep().ptr_r();
    val_type* ret_x	    = ret.rep().ptr_x();

    for (Integer j = 0; j < pos; ++j)
    {
        ret_r[j]        = dro[j];
        mrd::reset_helper(ret_x[j],dxo[j]);
    };

    nz = 0;
    Integer i;
    for (i = 0; i < n; ++i)
    {
        ret_c[i]        = nz;

        Integer last    = d_c[i+1] - off;
        if (last <= pos)
        {
            nz          = last;
            continue;
        };

        //current column contains zero
        nz              = pos;
        for (Integer k = pos+1; k < last; ++k)
        {
            const val_type& V  = dxo[k];
            if (t.eval(V) == false)
            {
                ret_r[nz]   = dro[k];
                mrd::reset_helper(ret_x[nz],V);
                ++nz;
            };
        };

        ++i;
        break;
    };

    for (Integer j = i; j < n; ++j)
    {
        ret_c[j]        = nz;

        for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
        {
            const val_type& V  = dx[k];
            if (t.eval(V) == false)
            {
                ret_r[nz]   = dr[k];
                mrd::reset_helper(ret_x[nz],V);
                ++nz;
            };
        };
    };
    
    ret_c[n] = nz;
    return ret;
};

template<class value_type_>
const_matrix<Matrix<value_type_,struct_sparse>>
Matrix<value_type_,struct_sparse>::drop() const
{
    return drop_impl(*this,test_zero<value_type_>());
};

template<class value_type_>
const_matrix<Matrix<value_type_,struct_sparse>>
Matrix<value_type_,struct_sparse>::drop(Real tol) const
{
    if (tol == 0.)
        return drop_impl(*this,test_zero<value_type_>());

    test_drop<Real,value_type_> td(mrd::abs_helper<Real>::eval(tol));
    return drop_impl(*this,td);
};

template <class value_type>
bool Matrix<value_type,struct_sparse>::is_same_matrix(const Matrix &m) const
{
    //it is enough to check column pointers
    if (this->rep().ptr_c() == m.rep().ptr_c() && this->rows() == m.rows() && this->cols() == m.cols())
        return true;
    else
        return false;
}

template<class value_type>
bool Matrix<value_type,struct_sparse>::all_finite() const
{
    return mr::all_finite_helper<Matrix>::eval(*this);
};

};};

template class matcl::raw::sparse_matrix_base<matcl::Integer>;
template class matcl::raw::sparse_matrix_base<matcl::Float>;
template class matcl::raw::sparse_matrix_base<matcl::Real>;
template class matcl::raw::sparse_matrix_base<matcl::Float_complex>;
template class matcl::raw::sparse_matrix_base<matcl::Complex>;
template class matcl::raw::sparse_matrix_base<matcl::Object>;

template void matcl::raw::sparse_matrix_base<matcl::Complex>::construct2<matcl::Real>(const Integer* r_ind, 
                    const Integer* c_ind, const Real* xr, const Real* xi, Integer r, Integer c, 
                    Integer nz, Integer nzmax);        

template void matcl::raw::sparse_matrix_base<matcl::Float_complex>
                ::construct2<matcl::Float>(const Integer* r_ind, const Integer* c_ind, const Float* xr, 
                                           const Float* xi, Integer r, Integer c, Integer nz, Integer nzmax);  

template class matcl::raw::Matrix<matcl::Integer, matcl::struct_sparse>;
template class matcl::raw::Matrix<matcl::Float, matcl::struct_sparse>;
template class matcl::raw::Matrix<matcl::Real, matcl::struct_sparse>;
template class matcl::raw::Matrix<matcl::Float_complex, matcl::struct_sparse>;
template class matcl::raw::Matrix<matcl::Complex, matcl::struct_sparse>;
template class matcl::raw::Matrix<matcl::Object, matcl::struct_sparse>;
