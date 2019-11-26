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

#include "matcl-matrep/general/config.h"
#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-matrep/details/refcount.h"
#include "matcl-matrep/matrix/struct_flag.h"
#include "matcl-core/IO/archive.h"
#include "matcl-matrep/objects/details/type_info_object.h"

namespace matcl { namespace raw
{
    template<class V> class sparse_matrix_base;
};};

namespace matcl { namespace raw { namespace details
{

template<class value_type>
class MATCL_MATREP_EXPORT sparse_ccs
{
    private:
        using tinfo         = ti::ti_type<value_type>;
        using refcount_str  = matcl::details::refcount_str_default;

    public:
        sparse_ccs(tinfo ti);
        sparse_ccs(tinfo ti,Integer r, Integer c, Integer nzmax);
        sparse_ccs(const sparse_ccs&);
        sparse_ccs(sparse_ccs&&);
        sparse_ccs(const sparse_ccs&&) = delete;
        ~sparse_ccs();                

        void                mark_unique(bool is_unique);
        bool                is_unique() const;
        //Integer           refcount() const;        
        sparse_ccs          copy(bool keep_maxcol) const;        
        sparse_ccs          clone(bool keep_maxcol) const;        
        sparse_ccs          make_unique(bool keep_maxcol) const;        
        void                destroy_data();

        void                add_memory(Integer s = 0);
        void                add_columns(Integer nc);
        void                add_duplications();

        bool                has_element(Integer i, Integer j, Integer &k) const;
        value_type          get_elem(Integer i, Integer j) const;
        void                add_element(Integer i, Integer j, Integer k);
        void                remove_element(Integer i, Integer j, Integer k);

        void                sort();
        bool                is_sorted() const;

        refcount_str*       get_refstr() const      { return m_refcount; };
        const Integer*      ptr_c() const           { return m_c; }
        const Integer*      ptr_r() const           { return *m_r; }
        const value_type*   ptr_x() const           { return *m_x; }
        Integer*            ptr_c()                 { return m_c; }
        Integer*            ptr_r()                 { return *m_r; }
        value_type*         ptr_x()                 { return *m_x; }

        Integer**           ptr_root_r()            { return m_r; }
        value_type**        ptr_root_x()            { return m_x; }

        const Integer* const *      ptr_root_r() const      { return m_r; }
        const value_type* const*    ptr_root_x() const      { return m_x; }

        //position of first element in rows and value pointer
        Integer             offset() const          { return m_offset; };
        //position of first element in column pointer
        Integer             offset_cols() const     { return Integer(m_c - m_c_root[0]); };

        Integer             rows() const            { return m_rows; }
        Integer             cols() const            { return m_cols; }
        //number of available elements in column pointer
        Integer             max_cols() const        { return m_max_cols - offset_cols(); }
        //total length of column pointer (including unused and reserved elements)
        Integer             total_cols() const      { return m_max_cols; }
        Integer             nzmax() const           { return *m_nzmax-m_offset; }
        Integer             nnz() const             { return m_c? m_c[m_cols]-m_offset : 0; }
        tinfo               get_type() const        { return m_ti; };
        struct_flag&        get_struct() const      { return m_flag; };

        void                assign_to_fresh(const sparse_ccs& o);
        void                assign_to_fresh(sparse_ccs&& o);

        void                serialize(oarchive_impl & ar, const unsigned int version) const;
        void                serialize(iarchive_impl & ar, const unsigned int version);

    private:
        void                alloc(Integer, Integer, Integer);
        void                resize(Integer, Integer, Integer);

        refcount_str*       m_refcount;
        Integer             m_rows;
        Integer             m_cols;
        Integer             m_max_cols;
        Integer             m_offset;
        Integer *           m_nzmax;
        Integer **          m_c_root;
        Integer *           m_c;
        Integer **          m_r;
        value_type **       m_x;
        tinfo               m_ti;  
        mutable struct_flag m_flag;
        bool                m_effective_unique;

        sparse_ccs&         operator=(const sparse_ccs&) = delete;
        sparse_ccs&         operator=(sparse_ccs&&) = delete;

        friend class matcl::raw::sparse_matrix_base<value_type>;
};

template<class V>
bool is_real_matrix(const sparse_ccs<V>&)
{
    return md::is_float_real_scalar<V>::value || std::is_same<V,Integer>::value;
};

};};};
