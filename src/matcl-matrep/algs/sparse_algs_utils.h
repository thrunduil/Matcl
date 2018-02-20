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

#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-internals/base/utils.h"
#include "matcl-matrep/utils/workspace.h"
#include <list>
#include <map>

namespace matcl { namespace algorithm { namespace details
{

enum sort_type
{
    not_sorted, sorted_increasing, sorted_decreasing
};

sort_type   is_sorted(const raw::integer_dense& vec);
Integer     number_dupl(const raw::integer_dense& vec);
Integer     number_dupl(const raw::integer_dense& vec_r, const raw::integer_dense& vec_c);
void        sort_rows_cols(Integer* ptr_r, Integer* ptr_c, Integer size);
void        sort_rows_cols(Integer* ptr_r, Integer* ptr_c, Integer* indices, Integer size);
template<class T>
void        sort_rows_cols_stable(Integer* ptr_r, Integer* ptr_c, T* ptr_x, Integer size);

void init_row_selector(const raw::integer_dense& ri,Integer val, matcl::pod_workspace<Integer>& v_work_row,
                       Integer& r_min, Integer& r_max);

class row_map
{
    private:
        using list_type         = std::list<Integer>;
        using map_type          = std::map<Integer,list_type>;
        using map_value_type    = map_type::value_type;

        map_type m_map;

    public:
        row_map(){};

        void insert(Integer row, Integer pos)
        {
            using iterator  = map_type::iterator;
            iterator elem   = m_map.find(row);

            if (elem == m_map.end())
            {
                list_type m_list;
                m_list.push_back(pos);
                m_map.insert(elem,map_value_type(row,m_list));
            }
            else
            {
                elem->second.push_back(pos);
            };
        };

        template<class T>
        void add_rows(Integer row,const T& value,T* d_x,Integer* d_r,Integer &nz) const
        {
            using iterator  = map_type::const_iterator;
            iterator elem   = m_map.find(row);
            const list_type& m_list = elem->second;

            using list_iterator = list_type::const_iterator;
            list_iterator pos   = m_list.begin();

            while(pos != m_list.end())
            {
                mrd::reset_helper(d_x[nz],value);
                d_r[nz] = *pos;
                ++nz;
                ++pos;
            };
        };
        
        void clear()
        {
            m_map.clear();
        };
};

};};};