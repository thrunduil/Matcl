/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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
#include "matcl-core/details/integer.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/const_matrix.h"

namespace matcl { namespace details
{

struct colon_info
{
    private:
        using Mat_I     = raw::integer_dense;
        using Mat_Ic    = raw::const_matrix<Mat_I>;

    private:
        //always vectors
        Mat_Ic*     m_ri;
        Mat_Ic*     m_ci;

    public:
        Integer r_start, r_step, r_end, r_size, r_rep_size, r_flag;
        Integer c_start, c_step, c_end, c_size, c_rep_size, c_flag;

        colon_info()
            :m_ri(nullptr), m_ci(nullptr)
        {
            r_start = r_step = r_end = r_size = r_rep_size = r_flag = 0;
            c_start = c_step = c_end = c_size = c_rep_size = c_flag = 0;
        };

        ~colon_info();

        Integer         rows() const;
        Integer         cols() const;
        Integer         rep_rows() const                { return r_rep_size; };
        Integer         rep_cols() const                { return c_rep_size; };
        bool            is_double_mat_colon() const     { return m_ci != nullptr; };

        const Mat_I&    get_rim_2() const               { return m_ri->get(); };
        const Mat_I&    get_cim_2() const               { return m_ci->get(); };

        const Integer*  get_ci_2_ptr() const            { return m_ci ? m_ci->get().ptr() : nullptr; };
        const Integer*  get_ri_2_ptr() const            { return m_ri ? m_ri->get().ptr() : nullptr; };

        Integer         row_index_2(Integer row) const;
        Integer         col_index_2(Integer col) const;

        const Mat_I&    get_rim_1() const               { return m_ri->get(); };
        const Mat_I&    get_rim_r() const               { return m_ri->get(); };
        const Mat_I&    get_rim_c() const               { return m_ci->get(); };
        const Integer*  get_ri_1_ptr() const            { return m_ri ? m_ri->get().ptr() : nullptr; };
        const Integer*  get_ri_r_ptr() const            { return m_ri ? m_ri->get().ptr() : nullptr; };
        const Integer*  get_ri_c_ptr() const            { return m_ci ? m_ci->get().ptr() : nullptr; };

        void            set_ci(const raw::integer_dense& m);
        void            set_ri(const raw::integer_dense& m);
        void            set_ri_2(const Mat_I& r, const Mat_I& c);

    private:
        colon_info(const colon_info&) = delete;
        colon_info& operator=(const colon_info&) = delete;
};

inline Integer colon_info::rows() const
{
    if (r_flag == 0)
        return m_ri->get().size();

    return r_size;
};

inline Integer colon_info::cols() const
{
    if (c_flag == 0)
        return m_ci->get().size();

    return c_size;
};

inline Integer colon_info::row_index_2(Integer row) const
{
    if (r_flag == 0)
        return m_ri->get().ptr()[row-1];

    return r_start + imult(r_step,row-1);
};

inline Integer colon_info::col_index_2(Integer col) const
{
    if (c_flag == 0)
        return m_ci->get().ptr()[col-1];

    return c_start + imult(c_step,col-1);
};

};};
