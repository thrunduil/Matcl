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

#include "matcl-matrep/matrix/matrix.h"

#include <vector>
#include <algorithm>

namespace matcl { namespace details
{

enum class inplace_type : int
{ 
    can_inplace, can_continue, cannot_continue 
};

class data_container
{
    public:
        enum type { type_matrix, type_col, type_row, scal_int, scal_float, scal_real, 
                    scal_complex, scal_fcomplex };

    public:
        Matrix      m_matrix;
        mat_col     m_col;
        mat_row     m_row;

        union
        {
            Integer m_int;
            Real    m_real;
            Float   m_float;
            Real    m_complex[2];
            Float   m_fcomplex[2];
        };

        type        m_type;

        data_container(const Matrix& m)		:m_matrix(m), m_type(type_matrix){};
        data_container(const mat_col& mc)	:m_col(mc), m_type(type_col){};
        data_container(const mat_row& mr)	:m_row(mr), m_type(type_row){};
        data_container(Integer val)			:m_int(val), m_type(scal_int){};
        data_container(Float val)			:m_float(val), m_type(scal_float){};
        data_container(Real val)			:m_real(val), m_type(scal_real){};

        data_container(const Complex& val)	:m_type(scal_complex)
        {	
            m_complex[0] = matcl::real(val);
            m_complex[1] = matcl::imag(val);
        };

        data_container(const Float_complex& val)
            :m_type(scal_fcomplex)
        {	
            m_fcomplex[0] = matcl::real(val);
            m_fcomplex[1] = matcl::imag(val);
        };

        inplace_type    check_for_inplace_build(bool is_row, bool is_sparse, matcl::value_code vt,
                                            ti::ti_object ti) const;
};

class mat_cons_data 
{
    public:
        using matrix_vector = std::vector<data_container>;

    public:
        mat_cons_data();

        mat_cons_data(const mat_cons_data&);

        void	            add_col(const mat_col& mc);
        void	            add_row(const mat_row& mr);
        void	            add_matrix(const Matrix& mc);

        void	            add_scalar(Integer val);
        void	            add_scalar(Real val);
        void	            add_scalar(Float val);
        void	            add_scalar(Complex val);
        void	            add_scalar(Float_complex val);
        void	            add_scalar(const Object& val);
        void	            destroy();

        bool                is_initialized() const;
        ti::ti_object       get_type() const;
        matcl::value_code   get_value_type() const;

    public:
        matrix_vector       m_vector;

    private:
        ti::ti_object       m_ti;
        matcl::value_code   m_vt;
        bool                m_initialized;

    private:		
        mat_cons_data& operator=(const mat_cons_data&)  = delete;

        value_code          link_value_types(value_code new_vt, value_code this_vt, 
                                             bool initialized) const;
};

inline mat_cons_data::mat_cons_data()
:m_ti(ti::predefined::get_ti_int()),m_vt(value_code::v_integer), m_initialized(false)
{};

inline mat_cons_data::mat_cons_data(const mat_cons_data& other)
:m_vector(other.m_vector),m_ti(other.m_ti),m_vt(other.m_vt), m_initialized(other.m_initialized)
{};

inline void mat_cons_data::add_col(const mat_col& mc)
{
    if (mc.is_initialized() == true)
        m_vt = link_value_types(mc.get_value_type(), this->get_value_type(), m_initialized);

    m_initialized = m_initialized || mc.is_initialized();

    if (m_vt == value_code::v_object)
        m_ti = ti::unify_ti<ti::ti_object>(m_ti,mc.get_type());

    m_vector.push_back(data_container(mc));
};

inline void mat_cons_data::add_row(const mat_row& mr)
{
    if (mr.is_initialized() == true)
        m_vt = link_value_types(mr.get_value_type(), this->get_value_type(), m_initialized);

    m_initialized = m_initialized || mr.is_initialized();

    if (m_vt == value_code::v_object)
        m_ti = ti::unify_ti<ti::ti_object>(m_ti,mr.get_type());
    
    m_vector.push_back(data_container(mr));
};

inline void mat_cons_data::add_matrix(const Matrix& m)
{
    m_vt = link_value_types(m.get_value_code(), this->get_value_type(), m_initialized);
    m_initialized = true;

    if (m_vt == value_code::v_object)
        m_ti = ti::unify_ti<ti::ti_object>(m_ti,m.get_type());
    
    m_vector.push_back(data_container(m));
};

inline void mat_cons_data::add_scalar(Integer val)
{
    m_vt = link_value_types(value_code::v_integer, this->get_value_type(), m_initialized);
    m_initialized = true;

    if (m_vt == value_code::v_object)
        m_ti = ti::unify_ti<ti::ti_object>(m_ti,ti::ti_int());
    
    m_vector.push_back(data_container(val));
};

inline void mat_cons_data::add_scalar(Real val)
{
    m_vt = link_value_types(value_code::v_real, this->get_value_type(), m_initialized);
    m_initialized = true;
    if (m_vt == value_code::v_object)
        m_ti = ti::unify_ti<ti::ti_object>(m_ti,ti::ti_real());

    m_vector.push_back(data_container(val));
};

inline void mat_cons_data::add_scalar(Float val)
{
    m_vt = link_value_types(value_code::v_float, this->get_value_type(), m_initialized);
    m_initialized = true;
    if (m_vt == value_code::v_object)
        m_ti = ti::unify_ti<ti::ti_object>(m_ti,ti::ti_float());

    m_vector.push_back(data_container(val));
};

inline void mat_cons_data::add_scalar(Complex val)
{
    m_vt = link_value_types(value_code::v_complex, this->get_value_type(),m_initialized);
    m_initialized = true;
    if (m_vt == value_code::v_object)
        m_ti = ti::unify_ti<ti::ti_object>(m_ti,ti::ti_compl());

    m_vector.push_back(data_container(val));
};

inline void mat_cons_data::add_scalar(Float_complex val)
{
    m_vt = link_value_types(value_code::v_float_complex, this->get_value_type(), m_initialized);
    m_initialized = true;
    if (m_vt == value_code::v_object)
        m_ti = ti::unify_ti<ti::ti_object>(m_ti,ti::ti_float_compl());

    m_vector.push_back(data_container(val));
};

inline void mat_cons_data::add_scalar(const Object& val)
{
    m_vt = value_code::v_object;
    m_initialized = true;
    m_ti = ti::unify_ti<ti::ti_object>(m_ti,val.get_type());
    m_vector.push_back(data_container(val));
};

inline void mat_cons_data::destroy()
{};

inline ti::ti_object mat_cons_data::get_type() const
{
    return m_ti;
};

inline matcl::value_code mat_cons_data::get_value_type() const
{
    return m_vt;
};

inline bool mat_cons_data::is_initialized() const
{
    return m_initialized;
};

};};