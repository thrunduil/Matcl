/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2018
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

#include "matcl-matrep/general/config.h"
#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-internals/container/mat_base.h"
#include "matcl-matrep/details/enablers.h"

namespace matcl { namespace raw 
{
    
template<class value_type_>
class MATCL_MATREP_EXPORT Matrix<value_type_,struct_dense> 
    : public dense_matrix_base<value_type_>
{
    public:
        using struct_type   = struct_dense;
        using value_type    = value_type_;
        using base_type     = dense_matrix_base<value_type_>;
        using refcount_str  = matcl::details::refcount_str_default;
        using BandMatrix    = Matrix<value_type,struct_banded>;
        using tinfo         = ti::ti_type<value_type>;

        using this_type     = Matrix<value_type_, struct_dense>;

        static const matcl::mat_code matrix_code
                            = md::type_to_code<this_type>::value;

        struct foreign{};

    private:
        //not available
        Matrix(const BandMatrix& mat);

    public:

        Matrix(tinfo ti);
        Matrix(tinfo ti,Integer r, Integer c);
                
        Matrix(Matrix &&mat);        
        Matrix(base_type &&mat);

        // mark copying as safe; copy is safe if returned object is not modified
        // unless is unique
        struct copy_is_safe{};

        //TODO: remove this
        struct copy_is_safe_TODO : copy_is_safe{};

        Matrix(const base_type &m, copy_is_safe)
                : Matrix(m) {}; 

        Matrix(const Matrix &m, copy_is_safe)
                : Matrix(m) {}; 

        template<class Real_T>
        Matrix(tinfo ti,const Real_T* xr, const Real_T* xi, Integer r, Integer c,
            typename matcl::details::enable_if_val_complex<value_type,Real_T>::type = 0);

        Matrix(tinfo ti,const value_type *arr, Integer r, Integer c);
        Matrix(tinfo ti,const value_type *arr, Integer r, Integer c, Integer ld);
        Matrix(tinfo ti, value_type *arr, Integer r, Integer c, Integer ld, foreign);
        Matrix(tinfo ti,const value_type &val, Integer r, Integer c);
        
        refcount_str*       get_refstr() const      { return base_type::get_refstr(); };

        inline
        const value_type&	operator()(Integer i, Integer j) const;
        Matrix			    get_diag(Integer d = 0) const;

        //only for raw matrix not stored in a container, for example allocated on stack
        void                assign_to_fresh(const Matrix&);
        void                assign_to_fresh(Matrix&&);
        
        Matrix              reshape(Integer r, Integer c) const;   
        Matrix              make_view(Integer r_start, Integer r_end) const;
        Matrix              make_view(Integer r_start, Integer r_end, Integer c_start, Integer c_end) const;
        matcl::Matrix       fast_optim() const;
        Matrix              copy(bool keep_bufor = false) const;
        Matrix              clone(bool keep_bufor = false) const;
        Matrix              make_unique(bool keep_bufor = false) const;
        bool                is_same_matrix(const Matrix& other) const;
        bool                all_finite() const;

        bool                is_explicit() const;
        Matrix              make_explicit() const;
        void                destroy_data();   

        Matrix              reserve(Integer r, Integer c) const;
        Matrix              resize(Integer r, Integer c) const;
        Matrix              resize(Integer r, Integer c);
        void                prepare_for_concat(Integer r, Integer c);

        value_type*         ptr()                   { return base_type::ptr(); };
        const value_type*   ptr() const             { return base_type::ptr(); };

        void			    serialize(oarchive_impl & ar, const unsigned int version) const;
        void				serialize(iarchive_impl & ar, const unsigned int version);

    private:
        Matrix&             operator=(const Matrix&) = delete;

        BandMatrix          get_diag_band() const;
        matcl::Matrix       check_change_struct() const;
        Real                estim_density() const;
        
        void                set_to_all(const value_type& val);               

        using base_type::ptr;

    public:
        //TODO
        Matrix(const Matrix &mat);
        Matrix(const base_type &mat);
};

};};

#include "matcl-internals/container/mat_d.inl"
