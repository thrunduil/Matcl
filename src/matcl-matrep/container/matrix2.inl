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

#include "matcl-matrep/details/matrix.inl"
#include "matcl-matrep/container/matrix_container.inl"
#include "matcl-matrep/visitors/assign_visitor.inl" 


namespace matcl
{

template<class M>
inline const M& details::get_functor<M>::eval(const matrix_base& mat)
{
	using value_type    = typename M::value_type;
	using struct_type   = typename M::struct_type;
	using type          = details::matrix_container<value_type,struct_type>;

	return static_cast<type*>(mat.m_value.m_mat.m_mat_ptr)->get();
};

template<class M>
inline M& details::get_functor<M>::eval(matrix_base& mat)
{
	using value_type    = typename M::value_type;
	using struct_type   = typename M::struct_type;
	using type          = details::matrix_container<value_type,struct_type>;

	return static_cast<type*>(mat.m_value.m_mat.m_mat_ptr)->get();
};

template<class M>
inline details::rvalue_holder<M> details::get_functor<M>::eval_move(matrix_base& mat)
{
	using value_type    = typename M::value_type;
	using struct_type   = typename M::struct_type;
	using type          = details::matrix_container<value_type,struct_type>;

    rvalue_holder<M> tmp(static_cast<type*>(mat.m_value.m_mat.m_mat_ptr));
    mat.m_type = mat_code::integer_scalar;
    
    return std::move(tmp);
};

template<class M>
inline M& details::get_functor<M>::get_ref(matrix_container_base* cont)
{
	using value_type    = typename M::value_type;
	using struct_type   = typename M::struct_type;
	using type          = details::matrix_container<value_type,struct_type>;

    return static_cast<type*>(cont)->get();
};

template<class M>
inline void details::get_functor<M>::destroy(matrix_container_base* cont)
{
	using value_type    = typename M::value_type;
	using struct_type   = typename M::struct_type;
	using type          = details::matrix_container<value_type,struct_type>;

    if (cont)
    {
        static_cast<type*>(cont)->call_destructor();
        free_container(cont);
    };
};

};