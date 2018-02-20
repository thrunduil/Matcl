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

#include "matcl-internals/container/sparse_ccs.h"
#include "matcl-internals/container/mat_base.h"
#include "matcl-scalar/details/enablers.h"

#pragma warning(push)
#pragma warning(disable:4251) //class needs to have dll-interface

namespace matcl { namespace raw 
{

template<class value_type_>
class MATCL_MATREP_EXPORT sparse_matrix_base
{
    public:
        using value_type    = value_type_;
        using struct_type   = struct_sparse;

        using SparseMatrix  = Matrix<value_type,struct_sparse>;
        using DenseMatrix   = Matrix<value_type,struct_dense>;
        using BandMatrix    = Matrix<value_type,struct_banded>;
        using sparse_ccs    = details::sparse_ccs<value_type>;
        using type_info     = ti::ti_type<value_type>;
        using refcount_str  = matcl::details::refcount_str_default;

    public:
        sparse_matrix_base(type_info ti) : m_data(ti)    {}
        sparse_matrix_base(type_info ti, Integer r, Integer c);
        sparse_matrix_base(type_info ti, Integer r, Integer c, Integer nzmax);
        sparse_matrix_base(type_info ti, const Integer* r_ind, const Integer* c_ind, 
                            const value_type* x_ind, Integer r, Integer c, Integer nnz);

        template<class Real_T>
        sparse_matrix_base(type_info ti, const Integer* ri, const Integer* ci, 
                            const Real_T* xr, const Real_T* xi, Integer r, Integer c, Integer nnz,
                            typename matcl::details::enable_if_val_complex<value_type,Real_T>::type = 0);

        sparse_matrix_base(type_info ti, const Integer* r_ind, const Integer* c_ind, const value_type* x_ind, 
                            Integer r,Integer c, Integer nz, Integer nzmax);

        template<class Real_T>
        sparse_matrix_base(type_info ti, const Integer* ri, const Integer* ci, const Real_T* xr, const Real_T* xi, 
                            Integer r, Integer c, Integer nnz, Integer nzmax,
                            typename matcl::details::enable_if_val_complex<value_type,Real_T>::type = 0);

        sparse_matrix_base(const sparse_matrix_base&&) = delete;
        sparse_matrix_base(sparse_matrix_base&&);
        sparse_matrix_base(const sparse_matrix_base&);

        explicit sparse_matrix_base(const sparse_ccs &d) : m_data(d) {}

        void                    mark_unique(bool unique)    { m_data.mark_unique(unique); }        
        bool                    is_unique() const           { return m_data.is_unique(); };
        refcount_str*           get_refstr() const          { return m_data.get_refstr(); };

        void                    destroy_data();

        inline value_type       operator()(Integer i, Integer j) const;
        
        SparseMatrix            get_diag(Integer = 0) const;    

        inline Integer          rows() const;
        inline Integer          cols() const;
        inline Integer          length() const;

        //total length of column pointer (including unused and reserved elements)
        inline Integer          total_cols() const;

        //number of available elements in column pointer
        inline Integer          max_cols() const;
        inline Integer          nnz() const;
        inline Integer          nzmax() const;

        sparse_ccs&             rep()                           { return m_data; }
        const sparse_ccs&       rep() const                     { return m_data; }

        type_info               get_type() const                  { return m_data.get_type(); };
        const struct_flag&      get_struct() const              { return m_data.get_struct(); };
        struct_flag&            get_struct()                    { return m_data.get_struct(); };
        void                    set_struct(struct_flag f) const { m_data.get_struct().set(f); };
        void                    add_struct(struct_flag f) const { m_data.get_struct().add(f); };

        //only for raw matrix not stored in a container, for example allocated on stack
        void                    assign_to_fresh(const sparse_matrix_base&);
        void                    assign_to_fresh(sparse_matrix_base&&);

        sparse_matrix_base      copy(bool keep_maxcol = false) const;
        sparse_matrix_base      clone(bool keep_maxcol = false) const;
        sparse_matrix_base      make_unique(bool keep_bufor = false) const;

        sparse_matrix_base      reserve(Integer r, Integer c) const;
        sparse_matrix_base      resize(Integer r, Integer c) const;
        sparse_matrix_base      resize(Integer r, Integer c);
        sparse_matrix_base      make_view(Integer c_start, Integer c_end) const;        

                                //not thread safe, matrix should be unique
        void                    prepare_for_concat(Integer r, Integer c, Integer nnz);

        void                    serialize(oarchive_impl & ar, const unsigned int version) const;
        void                    serialize(iarchive_impl & ar, const unsigned int version);

    //internal use
    public:
        //change number rows or columns without modifying data; new dimesions must be valid
        //new number of rows and columns cannot be higher than current
        void                    change_number_rows(Integer r);
        void                    change_number_cols(Integer c);

     protected:
        sparse_ccs              m_data;
        void                    construct(const Integer* r_ind, const Integer* c_ind, const value_type* x_ind, 
                                                Integer r, Integer c, Integer nz, Integer nzmax);

        DenseMatrix             full() const;
        BandMatrix              get_diag_band() const;

    private:
        template<class RealType>
        void                    construct2(const Integer* r_ind, const Integer* c_ind, const RealType* xr, 
                                            const RealType* xi, Integer r, Integer c, Integer nz, Integer nzmax);        
        sparse_matrix_base      resize_remrows(Integer r, Integer c) const;

        template<class val_type>
        friend matcl::Matrix    fast_optim_impl(const Matrix<val_type,struct_sparse>& A);

        sparse_matrix_base&     operator=(const sparse_matrix_base&) = delete;
};

template<class value_type_>
class MATCL_MATREP_EXPORT Matrix<value_type_,struct_sparse> : public sparse_matrix_base<value_type_>
{
    public:
        using base_type     = sparse_matrix_base<value_type_>;
        using struct_type   = struct_sparse;
        using value_type    = value_type_;

    public:
        Matrix(type_info ti)            
                : base_type(ti){ ; }

        Matrix(type_info ti,Integer r, Integer c)                
                : base_type(ti,r, c) {}

        Matrix(type_info ti,Integer r, Integer c, Integer nzmax)    
                : base_type(ti, r, c, nzmax) {}

        Matrix(type_info ti, const Integer *ri, const Integer *ci, const value_type *x,
                Integer r, Integer c, Integer nnz)
                : base_type(ti, ri, ci, x, r, c, nnz) {}

        Matrix(type_info ti,const Integer *ri, const Integer *ci, const value_type *x,
                Integer r, Integer c, Integer nnz, Integer nzmax)
                : base_type(ti, ri, ci, x, r, c, nnz, nzmax) {}

        template<class Real_T>
        Matrix(type_info ti,const Integer *ri, const Integer *ci, const Real_T *xr, const Real_T *xi,
                Integer r, Integer c, Integer nnz, 
                typename matcl::details::enable_if_val_complex<value_type, Real_T>::type = 0)
                : base_type(ti,ri, ci, xr, xi, r, c, nnz)  {}

        template<class Real_T>
        Matrix(type_info ti, const Integer *ri, const Integer *ci, const Real_T *xr, const Real_T *xi,
                Integer r, Integer c, Integer nnz, Integer nzmax, 
                typename matcl::details::enable_if_val_complex<value_type, Real_T>::type = 0)
                : base_type(ti,ri, ci, xr, xi, r, c, nnz, nzmax) {}

        Matrix(const base_type &&m) = delete;
        Matrix(base_type &&m)        
                : base_type(std::move(m)) {}
        Matrix(const base_type &m)        
                : base_type(m) {}

        Matrix(const Matrix &&m) = delete;
        Matrix(Matrix &&m)        
                : base_type(std::move(m)) {}
        Matrix(const Matrix &m)        
                : base_type(m) {}

        matcl::Matrix       fast_optim() const;
        Matrix              drop() const;
        Matrix              drop(Real tol) const;

        bool                is_same_matrix(const Matrix& other) const;
        bool                all_finite() const;

    private:
        Matrix&             operator=(const Matrix&) = delete;
};

template<class V>
inline bool is_real_matrix(const sparse_matrix_base<V>&)
{
    return md::is_float_real_scalar<V>::value || std::is_same<V,Integer>::value;
};

};};

#pragma warning(pop)

#include "matcl-internals/container/mat_s.inl"
